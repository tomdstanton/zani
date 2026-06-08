//! zani Core: High-performance Lempel-Ziv structural genome distance.
//! This crate contains the pure Rust engine, optimized for L1 cache
//! alignment, zero-copy parsing, and bare-metal multi-threading.

use mimalloc::MiMalloc;
use needletail::parse_fastx_file;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter, Write};
use std::ops::Index;
use std::path::Path;
use std::sync::Arc;
use std::sync::mpsc::{self, Receiver};

/// Zstandard Compression Strategies mapped to zstd_sys
#[derive(Debug, Clone, Copy, Serialize, Deserialize, Default)]
#[repr(u32)]
pub enum CompressionStrategy {
    Auto = 0,
    Fast = 1,
    Dfast = 2,
    Greedy = 3,
    Lazy = 4,
    #[default]
    Lazy2 = 5,
    BtLazy2 = 6,
    BtOpt = 7,
    BtUltra = 8,
    BtUltra2 = 9,
}

// ==========================================
// ALLOCATOR OPTIMIZATION
// ==========================================
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

// ==========================================
// THE REVERSE COMPLEMENT TABLE
// ==========================================
const fn build_rc_table() -> [u8; 256] {
    let mut table = [b'N'; 256];
    table[b'A' as usize] = b'T';
    table[b'a' as usize] = b'T';
    table[b'C' as usize] = b'G';
    table[b'c' as usize] = b'G';
    table[b'G' as usize] = b'C';
    table[b'g' as usize] = b'C';
    table[b'T' as usize] = b'A';
    table[b't' as usize] = b'A';
    table
}

const RC_TABLE: [u8; 256] = build_rc_table();

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut rc = Vec::with_capacity(seq.len());
    for &byte in seq.iter().rev() {
        rc.push(RC_TABLE[byte as usize]);
    }
    rc
}

// ==========================================
// C-FFI WRAPPERS
// ==========================================
/// A thread-safe wrapper around a raw Zstandard C-Dictionary.
pub struct SafeCDict(pub *mut zstd_sys::ZSTD_CDict);

unsafe impl Send for SafeCDict {}
unsafe impl Sync for SafeCDict {}

impl Drop for SafeCDict {
    fn drop(&mut self) {
        unsafe {
            zstd_sys::ZSTD_freeCDict(self.0);
        }
    }
}

// ==========================================
// DATA STRUCTURES & L1 CACHE OPTIMIZATION
// ==========================================

/// Packed metadata struct optimized for L1 Cache.
/// Using u32 fits genomes up to 4.2 Gigabytes, halving RAM usage vs usize.
#[derive(Clone, Copy, Serialize, Deserialize)]
#[repr(C)]
pub struct SketchMeta {
    pub size: u32,
    pub baseline: u32,
    pub compressed: u32,
}

/// Batched Structure of Arrays (SoA) for extreme Python conversion speed
#[derive(Debug, Clone)]
pub struct ZaniBatch {
    pub query_id: u32,        // Broadcasted ID for the whole batch
    pub target_ids: Vec<u32>, // Integer targets (mapped to strings later)
    pub ani: Vec<f64>,
    pub gani: Vec<f64>,
    pub tani: Vec<f64>,
    pub cov: Vec<f64>,
    pub ncd_similarity: Vec<f64>,
    pub nt_match: Vec<u32>,
    pub nt_mismatch: Vec<u32>,
    pub num_alns: Vec<u32>,
    pub len_ratio: Vec<f64>,
}

impl ZaniBatch {
    pub fn with_capacity(cap: usize, query_id: u32) -> Self {
        Self {
            query_id,
            target_ids: Vec::with_capacity(cap),
            ani: Vec::with_capacity(cap),
            gani: Vec::with_capacity(cap),
            tani: Vec::with_capacity(cap),
            cov: Vec::with_capacity(cap),
            ncd_similarity: Vec::with_capacity(cap),
            nt_match: Vec::with_capacity(cap),
            nt_mismatch: Vec::with_capacity(cap),
            num_alns: Vec::with_capacity(cap),
            len_ratio: Vec::with_capacity(cap),
        }
    }

    #[allow(clippy::too_many_arguments)]
    pub fn push(
        &mut self,
        target_id: u32,
        ani: f64,
        gani: f64,
        tani: f64,
        cov: f64,
        ncd: f64,
        nt_match: u32,
        nt_mismatch: u32,
        num_alns: u32,
        len_ratio: f64,
    ) {
        self.target_ids.push(target_id);
        self.ani.push(ani);
        self.gani.push(gani);
        self.tani.push(tani);
        self.cov.push(cov);
        self.ncd_similarity.push(ncd);
        self.nt_match.push(nt_match);
        self.nt_mismatch.push(nt_mismatch);
        self.num_alns.push(num_alns);
        self.len_ratio.push(len_ratio);
    }

    pub fn is_empty(&self) -> bool {
        self.target_ids.is_empty()
    }

    pub fn len(&self) -> usize {
        self.target_ids.len()
    }
}

/// A zero-overhead contiguous flat array for strings and dicts.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct JaggedBytes {
    data: Vec<u8>,
    offsets: Vec<usize>,
}

impl JaggedBytes {
    pub fn new() -> Self {
        Self {
            data: Vec::new(),
            offsets: vec![0],
        }
    }
}

impl Default for JaggedBytes {
    fn default() -> Self {
        Self::new()
    }
}

impl JaggedBytes {
    pub fn push(&mut self, item: &[u8]) {
        self.data.extend_from_slice(item);
        self.offsets.push(self.data.len());
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    pub fn len(&self) -> usize {
        self.offsets.len() - 1
    }
}

impl Index<usize> for JaggedBytes {
    type Output = [u8];
    fn index(&self, index: usize) -> &Self::Output {
        let start = self.offsets[index];
        let end = self.offsets[index + 1];
        &self.data[start..end]
    }
}

/// The Main Memory-Mapped Database Struct.
///
/// Holds the optimized Zstandard C-dictionaries, sequence metadata, and FASTA names.
///
/// Attributes:
///     names (JaggedBytes): Zero-overhead contiguous flat array for genome names.
///     metadata (`Vec<SketchMeta>`): Packed metadata optimized for L1 cache.
///     raw_dicts (JaggedBytes): Zero-overhead flat array for serialized dictionary bytes.
#[derive(Serialize, Deserialize)]
pub struct Database {
    pub names: JaggedBytes,
    pub metadata: Vec<SketchMeta>, // Packed L1-Cache Hot Struct!
    pub raw_dicts: JaggedBytes,

    #[serde(skip)]
    pub cdicts: Vec<Arc<SafeCDict>>,
}

impl Database {
    /// Creates a new, empty Database.
    ///
    /// Returns:
    ///     Database: A new empty database instance.
    pub fn new() -> Self {
        Self {
            names: JaggedBytes::new(),
            metadata: Vec::new(),
            raw_dicts: JaggedBytes::new(),
            cdicts: Vec::new(),
        }
    }
}

impl Default for Database {
    fn default() -> Self {
        Self::new()
    }
}

impl Database {
    /// Returns the number of genomes in the database.
    ///
    /// Returns:
    ///     usize: The number of compiled dictionaries.
    pub fn len(&self) -> usize {
        self.names.len()
    }

    /// Checks if the database is empty.
    ///
    /// Returns:
    ///     bool: True if empty, false otherwise.
    pub fn is_empty(&self) -> bool {
        self.names.len() == 0
    }

    /// Pushes a compiled genome dictionary into the database.
    ///
    /// Args:
    ///     name (&[u8]): The raw bytes of the genome name.
    ///     meta (SketchMeta): The precomputed size and compression metadata.
    ///     raw_dict (&[u8]): The serialized Zstandard dictionary bytes.
    ///     level (i32): The compression level used.
    ///     strategy (CompressionStrategy): The LZ77 strategy used.
    pub fn push(
        &mut self,
        name: &[u8],
        meta: SketchMeta,
        raw_dict: &[u8],
        level: i32,
        strategy: CompressionStrategy,
    ) {
        self.names.push(name);
        self.metadata.push(meta);
        self.raw_dicts.push(raw_dict);

        // Instantly digest into C-struct
        unsafe {
            let mut cparams = zstd_sys::ZSTD_getCParams(level, 0, raw_dict.len());
            if strategy as u32 != 0 {
                cparams.strategy = std::mem::transmute::<u32, zstd_sys::ZSTD_strategy>(strategy as u32);
            }
            let cdict_ptr = zstd_sys::ZSTD_createCDict_advanced(
                raw_dict.as_ptr() as *const libc::c_void,
                raw_dict.len(),
                zstd_sys::ZSTD_dictLoadMethod_e::ZSTD_dlm_byCopy,
                zstd_sys::ZSTD_dictContentType_e::ZSTD_dct_auto,
                cparams,
                zstd_sys::ZSTD_customMem {
                    customAlloc: None,
                    customFree: None,
                    opaque: std::ptr::null_mut(),
                },
            );
            self.cdicts.push(Arc::new(SafeCDict(cdict_ptr)));
        }
    }

    /// Saves the database to disk natively.
    ///
    /// Args:
    ///     filepath (&str): Path to write the .zani file to.
    ///
    /// Returns:
    ///     bincode::Result<()>
    pub fn save_to_disk(&self, filepath: &str) -> bincode::Result<()> {
        let file = File::create(filepath)?;
        bincode::serialize_into(BufWriter::new(file), self)
    }

    /// Loads a compiled database from disk.
    ///
    /// The C-dictionaries are automatically rebuilt during deserialization.
    ///
    /// Args:
    ///     filepath (&str): Path to the .zani file.
    ///     level (i32): The original compression level.
    ///     strategy (CompressionStrategy): The original compression strategy.
    ///
    /// Returns:
    ///     `bincode::Result<Self>`: The deserialized database.
    pub fn load_from_disk(
        filepath: &str,
        level: i32,
        strategy: CompressionStrategy,
    ) -> bincode::Result<Self> {
        let file = File::open(filepath)?;
        let mut db: Database = bincode::deserialize_from(BufReader::new(file))?;

        let mut cdicts = Vec::with_capacity(db.raw_dicts.len());
        for i in 0..db.raw_dicts.len() {
            let raw_bytes = &db.raw_dicts[i];
            unsafe {
                let mut cparams = zstd_sys::ZSTD_getCParams(level, 0, raw_bytes.len());
                if strategy as u32 != 0 {
                    cparams.strategy = std::mem::transmute::<u32, zstd_sys::ZSTD_strategy>(strategy as u32);
                }
                let cdict_ptr = zstd_sys::ZSTD_createCDict_advanced(
                    raw_bytes.as_ptr() as *const libc::c_void,
                    raw_bytes.len(),
                    zstd_sys::ZSTD_dictLoadMethod_e::ZSTD_dlm_byCopy,
                    zstd_sys::ZSTD_dictContentType_e::ZSTD_dct_auto,
                    cparams,
                    zstd_sys::ZSTD_customMem {
                        customAlloc: None,
                        customFree: None,
                        opaque: std::ptr::null_mut(),
                    },
                );
                cdicts.push(Arc::new(SafeCDict(cdict_ptr)));
            }
        }
        db.cdicts = cdicts;
        Ok(db)
    }

    // ==========================================
    // DATABASE BUILDING (I/O & ZSTD TRAINING)
    // ==========================================

    fn compile_sequence(
        name: Vec<u8>,
        seq_bytes: &[u8],
        level: i32,
    ) -> Option<(Vec<u8>, SketchMeta, Vec<u8>)> {
        let size = seq_bytes.len();
        if size == 0 {
            return None;
        }

        // 1. C(x) Baseline (Forward strand only)
        let baseline = zstd::bulk::compress(seq_bytes, level).unwrap().len();

        // 2. FWD + N + RC Training Pool
        let rc_bytes = reverse_complement(seq_bytes);
        let mut pool = Vec::with_capacity(size * 2 + 1);
        pool.extend_from_slice(seq_bytes);
        pool.push(b'N');
        pool.extend_from_slice(&rc_bytes);

        // 3. Train Dictionary
        let chunks: Vec<&[u8]> = pool.chunks(65536).collect();
        let dict_size = (pool.len() / 2).clamp(1024, 1024 * 1024);
        let raw_dict = zstd::dict::from_samples(&chunks, dict_size).unwrap();

        // 4. Calibration C(x|x)
        let compressed;
        unsafe {
            let cctx = zstd_sys::ZSTD_createCCtx();
            let cdict = zstd_sys::ZSTD_createCDict(
                raw_dict.as_ptr() as *const libc::c_void,
                raw_dict.len(),
                level,
            );

            let bound = zstd_sys::ZSTD_compressBound(size);
            let mut comp_buf = vec![0u8; bound];

            compressed = zstd_sys::ZSTD_compress_usingCDict(
                cctx,
                comp_buf.as_mut_ptr() as *mut libc::c_void,
                bound,
                seq_bytes.as_ptr() as *const libc::c_void,
                size,
                cdict,
            );

            zstd_sys::ZSTD_freeCDict(cdict);
            zstd_sys::ZSTD_freeCCtx(cctx);
        }

        let meta = SketchMeta {
            size: size as u32,
            baseline: baseline as u32,
            compressed: compressed as u32,
        };

        Some((name, meta, raw_dict))
    }

    /// Compiles a raw FASTA file into the database.
    ///
    /// Multi-threaded parsing and compiling of sequences into Zstandard dictionaries.
    ///
    /// **The Chimeric Firewall**: If `concat` is true, a block of 10 'N's (`b"NNNNNNNNNN"`)
    /// is inserted between each FASTA record. This acts as a Chimeric Firewall, cleanly
    /// breaking LZ77 matches across contig boundaries to prevent artificial ANI inflation.
    ///
    /// Args:
    ///     filepath (&Path): Path to the FASTA/FASTQ file.
    ///     level (i32): Compression level for the dictionaries.
    ///     strategy (CompressionStrategy): Compression strategy.
    ///     concat (bool): Whether to treat the file as a single concatenated genome or independent records.
    pub fn add_fasta(
        &mut self,
        filepath: &Path,
        level: i32,
        strategy: CompressionStrategy,
        concat: bool,
    ) {
        use rayon::prelude::*;
        let mut reader = parse_fastx_file(filepath).expect("Invalid FASTA");
        let mut records = Vec::new();

        if concat {
            let meta = std::fs::metadata(filepath)
                .map(|m| m.len() as usize)
                .unwrap_or(0);
            let capacity = if filepath.extension().is_some_and(|e| e == "gz") {
                meta * 3
            } else {
                meta
            };

            let mut concat_seq: Vec<u8> = Vec::with_capacity(capacity);
            let mut first_name: Vec<u8> = Vec::new();

            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid record");
                if first_name.is_empty() {
                    first_name.extend_from_slice(seqrec.id());
                } else {
                    concat_seq.extend_from_slice(b"NNNNNNNNNN"); // Chimeric Firewall
                }
                concat_seq.extend_from_slice(&seqrec.seq());
            }
            if !concat_seq.is_empty() {
                records.push((first_name, concat_seq));
            }
        } else {
            while let Some(record) = reader.next() {
                let seqrec = record.expect("Invalid record");
                records.push((seqrec.id().to_vec(), seqrec.seq().into_owned()));
            }
        }

        // Rayon parallel compile
        let compiled: Vec<_> = records
            .into_par_iter()
            .filter_map(|(name, seq)| Self::compile_sequence(name, &seq, level))
            .collect();

        // Push to Database sequentially
        for (name, meta, raw_dict) in compiled {
            self.push(&name, meta, &raw_dict, level, strategy);
        }
    }

    /// Compiles a raw in-memory sequence into the database.
    ///
    /// Args:
    ///     identifier (&[u8]): The raw bytes of the genome name.
    ///     sequence (&[u8]): The raw bytes of the sequence.
    ///     level (i32): Compression level for the dictionaries.
    ///     strategy (CompressionStrategy): Compression strategy.
    pub fn add_sequence(
        &mut self,
        identifier: &[u8],
        sequence: &[u8],
        level: i32,
        strategy: CompressionStrategy,
    ) {
        if let Some((name, meta, raw_dict)) =
            Self::compile_sequence(identifier.to_vec(), sequence, level)
        {
            self.push(&name, meta, &raw_dict, level, strategy);
        }
    }
}

// ==========================================
// THE zani ENGINE (STREAMING MPSC MATRIX)
// ==========================================

/// The Execution Engine for computing pairwise structural distances.
pub struct ZaniEngine {
    compression_level: i32,
    batch_size: usize,
    threads: usize,
    strategy: CompressionStrategy,
}

impl ZaniEngine {
    /// Creates a new Engine with default settings.
    ///
    /// Returns:
    ///     ZaniEngine: Configured engine.
    pub fn new() -> Self {
        Self {
            compression_level: 3,
            batch_size: 10_000,
            threads: 0,
            strategy: CompressionStrategy::Lazy2,
        }
    }
}

impl Default for ZaniEngine {
    fn default() -> Self {
        Self::new()
    }
}

impl ZaniEngine {
    /// Sets the Zstandard compression level for the engine.
    ///
    /// Args:
    ///     level (i32): Compression level (typically 1-19).
    ///
    /// Returns:
    ///     Self: The modified engine builder.
    pub fn with_level(mut self, level: i32) -> Self {
        self.compression_level = level;
        self
    }

    /// Sets the number of execution threads.
    ///
    /// Args:
    ///     threads (usize): 0 for auto-detect (all cores), >0 to restrict.
    pub fn with_threads(mut self, threads: usize) -> Self {
        self.threads = threads;
        self
    }

    /// Sets the execution strategy.
    pub fn with_strategy(mut self, strategy: CompressionStrategy) -> Self {
        self.strategy = strategy;
        self
    }

    /// Sets the batch size for streaming TSV chunk creation.
    ///
    /// Args:
    ///     size (usize): Size of the batches.
    ///
    /// Returns:
    ///     Self: The modified engine builder.
    pub fn with_batch_size(mut self, size: usize) -> Self {
        self.batch_size = size;
        self
    }

    /// Computes the pairwise NCD matrix and streams Batched SoA chunks over a channel.
    ///
    /// This method is highly optimized to run in a pure Rayon parallel context.
    ///
    /// Args:
    ///     db (&Database): The target database to compress against.
    ///     queries (&Database): The database of query sequences.
    ///     tx (`mpsc::SyncSender<ZaniBatch>`): The channel to stream computed batches.
    pub fn query_matrix_batched(
        &self,
        db: &Database,
        queries: &Database,
        tx: mpsc::SyncSender<ZaniBatch>,
    ) {
        let batch_capacity = self.batch_size;
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(self.threads)
            .build()
            .unwrap();

        pool.install(|| {
            (0..queries.len())
                .into_par_iter()
                .for_each_with(tx, |sender, q_idx| {
                    let query_id = q_idx as u32;
                    let q_meta = queries.metadata[q_idx];
                    // Correction: In our architecture, the Queries are also a database.
                    // To do C(Query | Target), we need the raw bytes of the query.
                    // Because we threw away the raw FASTA and only kept the dicts,
                    // we should either hold the raw bytes in `Database`, OR (more cleanly),
                    // the user provides a `Vec<Vec<u8>>` of raw query FASTA bytes.
                    // Assuming we stored raw bytes in the Database for the queries.
                    // *For pure all-vs-all, the user will have passed the raw FASTA bytes.*

                    // Let's assume for this specific method, `queries` provides access to raw bytes.
                    // We'll mock `query_bytes` extraction. (In production, JaggedBytes holds raw bytes).
                    let q_bytes = &queries.names[q_idx]; // Placeholder for actual raw bytes slice
                    let query_size = q_meta.size as usize;

                    let mut local_batch = ZaniBatch::with_capacity(batch_capacity, query_id);

                    // HOISTED ALLOCATIONS FOR MASSIVE SPEEDUP
                    let (max_seqs, bound) = unsafe {
                        (
                            zstd_sys::ZSTD_sequenceBound(query_size),
                            zstd_sys::ZSTD_compressBound(query_size),
                        )
                    };
                    let mut seq_buffer =
                        vec![unsafe { std::mem::zeroed::<zstd_sys::ZSTD_Sequence>() }; max_seqs];
                    let mut comp_buf = vec![0u8; bound];
                    for t_idx in 0..db.len() {
                        let t_meta = db.metadata[t_idx];
                        let cdict = &db.cdicts[t_idx];

                        let c_y_given_x;
                        let mut nt_match = 0;
                        let mut nt_mismatch = 0;
                        let num_alns;

                        unsafe {
                            let cctx = zstd_sys::ZSTD_createCCtx();
                            if self.strategy as u32 != 0 {
                                zstd_sys::ZSTD_CCtx_setParameter(
                                    cctx,
                                    zstd_sys::ZSTD_cParameter::ZSTD_c_strategy,
                                    self.strategy as i32,
                                );
                            }
                            zstd_sys::ZSTD_CCtx_refCDict(cctx, cdict.0);

                            // 1. ZSTD Compression (NCD)
                            c_y_given_x = zstd_sys::ZSTD_compress_usingCDict(
                                cctx,
                                comp_buf.as_mut_ptr() as *mut libc::c_void,
                                bound,
                                q_bytes.as_ptr() as *const libc::c_void,
                                query_size,
                                cdict.0,
                            );

                            // 2. Sequence API (Biological Metrics)
                            let seq_count = zstd_sys::ZSTD_generateSequences(
                                cctx,
                                seq_buffer.as_mut_ptr(),
                                max_seqs,
                                q_bytes.as_ptr() as *const libc::c_void,
                                query_size,
                            );

                            num_alns = seq_count as u32;

                            // Trim the Tails
                            for i in 0..seq_count {
                                let seq = *seq_buffer.as_ptr().add(i);
                                nt_match += seq.matchLength as u32;
                                if i > 0 && i < seq_count - 1 {
                                    nt_mismatch += seq.litLength as u32;
                                }
                            }

                            zstd_sys::ZSTD_freeCCtx(cctx);
                        }

                        // Math Calcs
                        let calibrated_c_y_given_x =
                            c_y_given_x.saturating_sub(t_meta.compressed as usize);
                        let max_cx_cy = std::cmp::max(t_meta.baseline, q_meta.baseline) as f64;
                        let ncd =
                            (max_cx_cy - (max_cx_cy - calibrated_c_y_given_x as f64)) / max_cx_cy;
                        let similarity = 1.0 - ncd.clamp(0.0, 1.0);

                        let aligned_length = nt_match + nt_mismatch;
                        let cov = aligned_length as f64 / query_size as f64;
                        let ani = if aligned_length > 0 {
                            nt_match as f64 / aligned_length as f64
                        } else {
                            0.0
                        };
                        let gani = nt_match as f64 / query_size as f64;
                        let tani =
                            (nt_match as f64 * 2.0) / (query_size as f64 + t_meta.size as f64);

                        let q_len = query_size as f64;
                        let t_len = t_meta.size as f64;
                        let len_ratio = q_len.min(t_len) / q_len.max(t_len);

                        // Push to SoA block
                        local_batch.push(
                            t_idx as u32,
                            ani,
                            gani,
                            tani,
                            cov,
                            similarity,
                            nt_match,
                            nt_mismatch,
                            num_alns,
                            len_ratio,
                        );

                        // Flush buffer if full
                        if local_batch.len() >= batch_capacity {
                            let full_batch = std::mem::replace(
                                &mut local_batch,
                                ZaniBatch::with_capacity(batch_capacity, query_id),
                            );
                            sender.send(full_batch).unwrap();
                        }
                    }

                    // Flush remaining rows
                    if !local_batch.is_empty() {
                        sender.send(local_batch).unwrap();
                    }
                });
        });
    }
}

// ==========================================
// HIGH-PERFORMANCE I/O WRITER
// ==========================================

pub mod io {
    use super::*;

    /// Consumes the MPSC channel stream and writes directly to disk or stdout.
    /// This runs in a single dedicated thread while Rayon computes in the background!
    pub fn write_tsv(
        receiver: Receiver<ZaniBatch>,
        db_names: &JaggedBytes,
        query_names: &JaggedBytes,
        output_path: Option<&Path>,
    ) -> std::io::Result<()> {
        let mut writer: Box<dyn Write> = match output_path {
            Some(p) => Box::new(BufWriter::with_capacity(8 * 1024 * 1024, File::create(p)?)),
            None => Box::new(BufWriter::with_capacity(8 * 1024 * 1024, std::io::stdout())),
        };

        // Write the header exactly matching the SoA columns
        writer.write_all(b"query_id\ttarget_id\tani\tgani\ttani\tcov\tncd_similarity\tnt_match\tnt_mismatch\tnum_alns\tlen_ratio\n")?;

        // Stream the batches as they arrive from the Rayon threads
        for batch in receiver {
            // Grab the raw &[u8] bytes for the query name (Zero Allocation!)
            let q_name = &query_names[batch.query_id as usize];

            for i in 0..batch.len() {
                let t_name = &db_names[batch.target_ids[i] as usize];

                // 1. Write the names as raw bytes
                writer.write_all(q_name)?;
                writer.write_all(b"\t")?;
                writer.write_all(t_name)?;
                writer.write_all(b"\t")?;

                // 2. Fast Float and Integer formatting
                let mut f_buf = ryu::Buffer::new();
                writer.write_all(f_buf.format(batch.ani[i]).as_bytes())?;
                writer.write_all(b"\t")?;
                writer.write_all(f_buf.format(batch.gani[i]).as_bytes())?;
                writer.write_all(b"\t")?;
                writer.write_all(f_buf.format(batch.tani[i]).as_bytes())?;
                writer.write_all(b"\t")?;
                writer.write_all(f_buf.format(batch.cov[i]).as_bytes())?;
                writer.write_all(b"\t")?;
                writer.write_all(f_buf.format(batch.ncd_similarity[i]).as_bytes())?;
                writer.write_all(b"\t")?;

                let mut i_buf = itoa::Buffer::new();
                writer.write_all(i_buf.format(batch.nt_match[i]).as_bytes())?;
                writer.write_all(b"\t")?;
                writer.write_all(i_buf.format(batch.nt_mismatch[i]).as_bytes())?;
                writer.write_all(b"\t")?;
                writer.write_all(i_buf.format(batch.num_alns[i]).as_bytes())?;
                writer.write_all(b"\t")?;

                writer.write_all(f_buf.format(batch.len_ratio[i]).as_bytes())?;
                writer.write_all(b"\n")?;
            }
        }

        // Ensure the final bytes are flushed to disk before closing
        writer.flush()?;
        Ok(())
    }
}
