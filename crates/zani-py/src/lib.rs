#![allow(non_local_definitions)]
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::Path;

// Import the pure Rust engine from your workspace core crate
use zani::{io, CompressionStrategy, Database, ZaniEngine};

fn parse_strategy(s: &str) -> PyResult<CompressionStrategy> {
    match s.to_lowercase().as_str() {
        "auto" => Ok(CompressionStrategy::Auto),
        "fast" => Ok(CompressionStrategy::Fast),
        "dfast" => Ok(CompressionStrategy::Dfast),
        "greedy" => Ok(CompressionStrategy::Greedy),
        "lazy" => Ok(CompressionStrategy::Lazy),
        "lazy2" => Ok(CompressionStrategy::Lazy2),
        "btlazy2" => Ok(CompressionStrategy::BtLazy2),
        "btopt" => Ok(CompressionStrategy::BtOpt),
        "btultra" => Ok(CompressionStrategy::BtUltra),
        "btultra2" => Ok(CompressionStrategy::BtUltra2),
        _ => Err(PyValueError::new_err(format!("Invalid strategy: '{}'", s))),
    }
}

// ==========================================
// THE DATABASE WRAPPER (Reader & Writer)
// ==========================================

/// PyDatabase is the Python-facing wrapper for the compiled Zstandard database.
/// It acts exactly like a native Python object, but holds pure Rust memory.
///
/// **Architectural Note: Separation of Concerns**
/// `Database` is intentionally separate from `Engine`. This class is strictly a **storage and I/O** 
/// construct (parsing FASTAs, building Zstd dictionaries, and serializing/deserializing bytes). 
/// By not holding transient execution state (like thread counts or batch sizes), multiple different 
/// `Engine` instances can safely read from the exact same `Database` simultaneously across different 
/// threads without locking or mutation.
#[pyclass(name = "Database", module = "zani._zani_rs")]
pub struct PyDatabase {
    pub inner: Database,
}

#[pymethods]
impl PyDatabase {
    /// Initialize an empty database from Python.
    ///
    /// Returns:
    ///     Database: A new, empty database.
    ///
    /// Example:
    ///     >>> import zani
    ///     >>> db = zani.Database()
    #[new]
    pub fn new() -> Self {
        Self {
            inner: Database::new(),
        }
    }

    /// Number of genomes currently loaded in the database.
    ///
    /// Returns:
    ///     int: The number of compiled dictionaries in the database.
    ///
    /// Example:
    ///     >>> len(db)
    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Reads a FASTA/FASTQ file, compiles the sequences, and adds them to the database.
    ///
    /// **Note:** All DNA sequences are automatically **reverse‑complemented** before being added to the
    /// internal `Database`. This ensures that the algorithm always works on the 5'→3' strand, regardless
    /// of the orientation of the input FASTA/FASTQ file.
    ///
    /// Args:
    ///     filepath (str): Path to the `.fna` or `.fastq` file.
    ///     level (int): Zstandard compression level (typically 1-19).
    ///     strategy (str): LZ77 match finding strategy.
    ///     concat (bool): If True, concatenates all sequences in the file into a single genome.
    ///
    /// Raises:
    ///     PyFileNotFoundError: If the FASTA file is missing.
    #[pyo3(signature = (filepath, level=3, strategy="lazy2", concat=true))]
    pub fn add_fasta(
        &mut self,
        filepath: &str,
        level: i32,
        strategy: &str,
        concat: bool,
    ) -> PyResult<()> {
        let path = Path::new(filepath);
        if !path.exists() {
            return Err(PyFileNotFoundError::new_err(format!(
                "File not found: {}",
                filepath
            )));
        }

        let strat = parse_strategy(strategy)?;
        self.inner.add_fasta(path, level, strat, concat);
        Ok(())
    }

    /// Compiles a raw in-memory sequence and adds it to the database.
    ///
    /// Args:
    ///     identifier (bytes): The raw bytes of the genome name.
    ///     sequence (bytes): The raw bytes of the sequence.
    ///     level (int): Zstandard compression level (typically 1-19).
    ///     strategy (str): LZ77 match finding strategy.
    #[pyo3(signature = (identifier, sequence, level=3, strategy="lazy2"))]
    pub fn add_sequence(
        &mut self,
        identifier: &[u8],
        sequence: &[u8],
        level: i32,
        strategy: &str,
    ) -> PyResult<()> {
        let strat = parse_strategy(strategy)?;
        self.inner.add_sequence(identifier, sequence, level, strat);
        Ok(())
    }

    /// Saves the compiled database to disk as a binary file.
    ///
    /// Args:
    ///     filepath (str): Path to save the `.zani` database to.
    ///
    /// Raises:
    ///     IOError: If writing to the disk fails.
    pub fn save(&self, filepath: &str) -> PyResult<()> {
        self.inner
            .save_to_disk(filepath)
            .map_err(|e| PyIOError::new_err(format!("Failed to write database to disk: {}", e)))
    }

    /// Load a previously compiled database from disk.
    ///
    /// Args:
    ///     filepath (str): Path to the compiled `.zani` file.
    ///     level (int): The original compression level used during creation.
    ///     strategy (str): The original compression strategy.
    ///
    /// Returns:
    ///     Database: The loaded database object.
    ///
    /// Raises:
    ///     IOError: If reading from the disk fails.
    #[staticmethod]
    #[pyo3(signature = (filepath, level=3, strategy="lazy2"))]
    pub fn load(filepath: &str, level: i32, strategy: &str) -> PyResult<Self> {
        let strat = parse_strategy(strategy)?;
        let db = Database::load_from_disk(filepath, level, strat)
            .map_err(|e| PyIOError::new_err(format!("Failed to load database from disk: {}", e)))?;

        Ok(Self { inner: db })
    }
}

// Custom error type for missing files
pyo3::create_exception!(
    zani._zani_rs,
    PyFileNotFoundError,
    pyo3::exceptions::PyFileNotFoundError
);

// ==========================================
// THE EXECUTION ENGINE (Releasing the GIL)
// ==========================================

/// Executes the multithreaded Rayon matrix and writes the TSV to disk.
/// This completely bypasses the Python interpreter and drops the GIL.
///
/// **Architectural Note: Separation of Concerns**
/// `Engine` is intentionally separate from `Database`. It is purely a **compute and scheduling** 
/// construct. Keeping it separate prevents `Database` from becoming a monolithic "god object" and avoids
/// ambiguity when comparing two databases (e.g., target vs query). It also allows you to configure an 
/// execution context once and apply it immutably to multiple isolated databases.
#[pyclass(name = "Engine", module = "zani._zani_rs")]
pub struct PyEngine {
    pub compression_level: i32,
    pub batch_size: usize,
    pub threads: usize,
    pub strategy: zani::CompressionStrategy,
}

#[pymethods]
impl PyEngine {
    /// Initializes the computation Engine.
    ///
    /// Args:
    ///     compression_level (int): Zstandard compression level to use during execution (1-19).
    ///     batch_size (int): Number of rows to process before flushing to disk.
    ///     threads (int): Number of threads. 0 = auto-detect all cores.
    ///     strategy (str): Compression strategy.
    ///
    /// Returns:
    ///     Engine: A highly optimized execution engine.
    #[new]
    #[pyo3(signature = (compression_level=3, batch_size=10_000, threads=0, strategy="lazy2"))]
    pub fn new(
        compression_level: i32,
        batch_size: usize,
        threads: usize,
        strategy: &str,
    ) -> PyResult<Self> {
        let strat = parse_strategy(strategy)?;
        Ok(Self {
            compression_level,
            batch_size,
            threads,
            strategy: strat,
        })
    }

    /// Computes the All-vs-All pairwise distance matrix for a database.
    ///
    /// Releases the Python GIL and executes purely in Rust across all CPU cores.
    ///
    /// Args:
    ///     db (Database): The target database to compute against itself.
    ///     output_filepath (str): Path to output the resulting TSV matrix.
    ///
    /// Raises:
    ///     IOError: If writing the TSV fails.
    pub fn all_vs_all(&self, py: Python, db: &PyDatabase, output_filepath: &str) -> PyResult<()> {
        self.search(py, db, db, output_filepath)
    }

    /// Computes the structural distances between a query database and a target database.
    ///
    /// Releases the Python GIL and executes purely in Rust across all CPU cores.
    ///
    /// Args:
    ///     db (Database): The reference target database.
    ///     queries (Database): The query database.
    ///     output_filepath (str): Path to output the resulting TSV matrix.
    ///
    /// Raises:
    ///     IOError: If writing the TSV fails.
    ///     ValueError: If either database is empty.
    pub fn search(
        &self,
        py: Python,
        db: &PyDatabase,
        queries: &PyDatabase,
        output_filepath: &str,
    ) -> PyResult<()> {
        if db.inner.is_empty() || queries.inner.is_empty() {
            return Err(PyValueError::new_err(
                "Cannot run matrix: Database or Queries are empty.",
            ));
        }

        let out_path = Path::new(output_filepath);
        let engine = ZaniEngine::new()
            .with_level(self.compression_level)
            .with_batch_size(self.batch_size)
            .with_threads(self.threads)
            .with_strategy(self.strategy);

        // RELEASE THE PYTHON GIL!
        let write_result = py.allow_threads(|| {
            let (tx, rx) = std::sync::mpsc::sync_channel(100);

            std::thread::scope(|s| {
                s.spawn(|| {
                    engine.query_matrix_batched(&db.inner, &queries.inner, tx);
                });

                io::write_tsv(rx, &db.inner.names, &queries.inner.names, Some(out_path))
            })
        });

        write_result
            .map_err(|e| PyIOError::new_err(format!("Failed to write TSV output to disk: {}", e)))
    }
}

// ==========================================
// MODULE EXPORT
// ==========================================

#[pymodule]
fn _zani_rs(py: Python, m: &PyModule) -> PyResult<()> {
    // Add our custom exceptions
    m.add("PyFileNotFoundError", py.get_type::<PyFileNotFoundError>())?;

    // Add the Database Class
    m.add_class::<PyDatabase>()?;

    // Add the Engine Class
    m.add_class::<PyEngine>()?;

    Ok(())
}
impl Default for PyDatabase {
    fn default() -> Self {
        Self::new()
    }
}
