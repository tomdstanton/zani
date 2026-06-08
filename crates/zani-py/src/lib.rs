#![allow(non_local_definitions)]
use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use std::path::Path;

// Import the pure Rust engine from your workspace core crate
use zani::{io, Database, ZaniEngine};

// ==========================================
// THE DATABASE WRAPPER (Reader & Writer)
// ==========================================

/// PyDatabase is the Python-facing wrapper for the compiled Zstandard database.
/// It acts exactly like a native Python object, but holds pure Rust memory.
#[pyclass(name = "Database", module = "zani._zani_rs")]
pub struct PyDatabase {
    pub inner: Database,
}

#[pymethods]
impl PyDatabase {
    /// Initialize an empty database from Python.
    /// Example: `db = zani.Database()`
    #[new]
    pub fn new() -> Self {
        Self {
            inner: Database::new(),
        }
    }

    /// Number of genomes currently loaded in the database.
    /// Example: `len(db)`
    pub fn __len__(&self) -> usize {
        self.inner.len()
    }

    /// Read a FASTA/FASTQ file, compile the Zstandard dictionaries, and add to the database.
    pub fn add_fasta(&mut self, filepath: &str, level: i32, concat: bool) -> PyResult<()> {
        let path = Path::new(filepath);
        if !path.exists() {
            return Err(PyFileNotFoundError::new_err(format!(
                "File not found: {}",
                filepath
            )));
        }

        self.inner.add_fasta(path, level, concat);
        Ok(())
    }

    /// Writer Class Method: Save the compiled C-struct database to disk as a binary `.zani` file.
    pub fn save(&self, filepath: &str) -> PyResult<()> {
        self.inner.save_to_disk(filepath).map_err(|e| {
            PyIOError::new_err(format!("Failed to write database to disk: {}", e))
        })
    }

    /// Reader Class Method: Load a previously compiled `.zani` database from disk.
    /// Because it relies on the C-FFI, we must re-instantiate the ZSTD dictionaries.
    #[staticmethod]
    pub fn load(filepath: &str, level: i32) -> PyResult<Self> {
        let db = Database::load_from_disk(filepath, level).map_err(|e| {
            PyIOError::new_err(format!("Failed to load database from disk: {}", e))
        })?;

        Ok(Self { inner: db })
    }
}


// Custom error type for missing files
pyo3::create_exception!(zani._zani_rs, PyFileNotFoundError, pyo3::exceptions::PyFileNotFoundError);

// ==========================================
// THE EXECUTION ENGINE (Releasing the GIL)
// ==========================================

/// Executes the multithreaded Rayon matrix and writes the TSV to disk.
/// This completely bypasses the Python interpreter and drops the GIL.
#[pyclass(name = "Engine", module = "zani._zani_rs")]
pub struct PyEngine {
    pub compression_level: i32,
    pub batch_size: usize,
}

#[pymethods]
impl PyEngine {
    #[new]
    #[pyo3(signature = (compression_level=3, batch_size=10_000))]
    pub fn new(compression_level: i32, batch_size: usize) -> Self {
        Self { compression_level, batch_size }
    }

    pub fn all_vs_all(&self, py: Python, db: &PyDatabase, output_filepath: &str) -> PyResult<()> {
        self.search(py, db, db, output_filepath)
    }

    pub fn search(&self, py: Python, db: &PyDatabase, queries: &PyDatabase, output_filepath: &str) -> PyResult<()> {
        if db.inner.is_empty() || queries.inner.is_empty() {
            return Err(PyValueError::new_err("Cannot run matrix: Database or Queries are empty."));
        }

        let out_path = Path::new(output_filepath);
        let engine = ZaniEngine::new()
            .with_level(self.compression_level)
            .with_batch_size(self.batch_size);

        // RELEASE THE PYTHON GIL!
        let write_result = py.allow_threads(|| {
            let (tx, rx) = std::sync::mpsc::sync_channel(100);
            
            std::thread::scope(|s| {
                s.spawn(|| {
                    engine.query_matrix_batched(&db.inner, &queries.inner, tx);
                });

                io::write_tsv(
                    rx,
                    &db.inner.names,
                    &queries.inner.names,
                    out_path,
                )
            })
        });

        write_result.map_err(|e| {
            PyIOError::new_err(format!("Failed to write TSV output to disk: {}", e))
        })
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
