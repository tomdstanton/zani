use clap::{Parser, Subcommand, ValueEnum};
use log::{info, warn};
use std::path::PathBuf;
use std::time::Instant;

// Import your pure Rust engine from the lib.rs file
use zani::{CompressionStrategy, Database, ZaniEngine, io};

#[derive(ValueEnum, Clone, Debug)]
#[clap(rename_all = "lower")]
enum CliStrategy {
    Auto,
    Fast,
    Dfast,
    Greedy,
    Lazy,
    Lazy2,
    Btlazy2,
    Btopt,
    Btultra,
    Btultra2,
}

impl Into<CompressionStrategy> for CliStrategy {
    fn into(self) -> CompressionStrategy {
        match self {
            CliStrategy::Auto => CompressionStrategy::Auto,
            CliStrategy::Fast => CompressionStrategy::Fast,
            CliStrategy::Dfast => CompressionStrategy::Dfast,
            CliStrategy::Greedy => CompressionStrategy::Greedy,
            CliStrategy::Lazy => CompressionStrategy::Lazy,
            CliStrategy::Lazy2 => CompressionStrategy::Lazy2,
            CliStrategy::Btlazy2 => CompressionStrategy::BtLazy2,
            CliStrategy::Btopt => CompressionStrategy::BtOpt,
            CliStrategy::Btultra => CompressionStrategy::BtUltra,
            CliStrategy::Btultra2 => CompressionStrategy::BtUltra2,
        }
    }
}

/// zani: Zstandard-based Average Nucleotide Identity (ANI).
#[derive(Parser, Debug)]
#[command(author, version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Compute an all-vs-all pairwise distance matrix
    Ava {
        /// Input FASTA/FASTQ files to compare
        #[arg(required = true, num_args = 1..)]
        genomes: Vec<PathBuf>,

        /// Output TSV file path (Defaults to stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Zstandard compression level (1-19). Higher is more accurate but slower.
        #[arg(short, long, default_value_t = 3)]
        level: i32,

        /// Compression match strategy
        #[arg(short, long, value_enum, default_value_t = CliStrategy::Lazy2)]
        strategy: CliStrategy,

        /// Treat multi-record FASTAs as a single concatenated genome
        #[arg(short, long, default_value_t = true)]
        concat: bool,

        /// How many rows to buffer before writing to disk
        #[arg(short, long, default_value_t = 10_000)]
        batch_size: usize,

        /// Number of execution threads (0 = auto)
        #[arg(short, long, default_value_t = 0)]
        threads: usize,

        /// Save the compiled database to disk (.zani file)
        #[arg(short, long)]
        db: Option<PathBuf>,
    },

    /// Compute distances of queries against a compiled target database
    Search {
        /// Input FASTA/FASTQ files to search
        #[arg(required = true, num_args = 1..)]
        genomes: Vec<PathBuf>,

        /// Path to the pre-compiled target database (.zani)
        #[arg(short, long, required = true)]
        db: PathBuf,

        /// Output TSV file path (Defaults to stdout)
        #[arg(short, long)]
        output: Option<PathBuf>,

        /// Zstandard compression level (1-19). Must match database creation!
        #[arg(short, long, default_value_t = 3)]
        level: i32,

        /// Compression match strategy
        #[arg(short, long, value_enum, default_value_t = CliStrategy::Lazy2)]
        strategy: CliStrategy,

        /// Treat multi-record FASTAs as a single concatenated genome
        #[arg(short, long, default_value_t = true)]
        concat: bool,

        /// How many rows to buffer before writing to disk
        #[arg(short, long, default_value_t = 10_000)]
        batch_size: usize,

        /// Number of execution threads (0 = auto)
        #[arg(short, long, default_value_t = 0)]
        threads: usize,
    },

    /// Pre-compile a database of genomes to disk
    Build {
        /// Input FASTA files
        #[arg(required = true, num_args = 1..)]
        genomes: Vec<PathBuf>,

        /// Output database file (.zani)
        #[arg(short, long, required = true)]
        db: PathBuf,

        /// Zstandard compression level (1-19).
        #[arg(short, long, default_value_t = 3)]
        level: i32,

        /// Compression match strategy
        #[arg(short, long, value_enum, default_value_t = CliStrategy::Lazy2)]
        strategy: CliStrategy,

        /// Treat multi-record FASTAs as a single concatenated genome
        #[arg(short, long, default_value_t = true)]
        concat: bool,
    },
}

fn main() -> anyhow::Result<()> {
    // Initialize the logger so we can print beautiful [INFO] tags
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    let cli = Cli::parse();

    match cli.command {
        Commands::Ava {
            genomes,
            output,
            level,
            strategy,
            concat,
            batch_size,
            threads,
            db: db_path,
        } => {
            info!("Starting zani All-vs-All matrix calculation...");
            info!("Compression Level: {}, Strategy: {:?}", level, strategy);
            info!("Input Files: {}", genomes.len());

            let start_time = Instant::now();
            let strat: CompressionStrategy = strategy.into();

            info!("Step 1/3: Compiling genomes into Zstandard dictionaries...");
            let mut db = Database::new();

            for path in &genomes {
                if !path.exists() {
                    warn!("File not found, skipping: {:?}", path);
                    continue;
                }
                db.add_fasta(path, level, strat, concat);
            }

            if db.is_empty() {
                return Err(anyhow::anyhow!("No valid genomes loaded. Exiting."));
            }
            info!("Successfully compiled {} sketches.", db.len());

            if let Some(path) = &db_path {
                info!("Saving compiled database to {:?}...", path);
                if let Some(p_str) = path.to_str() {
                    if let Err(e) = db.save_to_disk(p_str) {
                        warn!("Failed to save database to disk: {}", e);
                    }
                } else {
                    warn!("Invalid path string for database save.");
                }
            }

            info!("Step 2/3: Executing multi-threaded NxN compression matrix...");
            let engine = ZaniEngine::new()
                .with_level(level)
                .with_batch_size(batch_size)
                .with_threads(threads)
                .with_strategy(strat);

            let (tx, rx) = std::sync::mpsc::sync_channel(100);

            std::thread::scope(|s| {
                s.spawn(|| {
                    engine.query_matrix_batched(&db, &db, tx);
                });

                if let Some(p) = &output {
                    info!("Step 3/3: Streaming results to disk at {:?}", p);
                } else {
                    info!("Step 3/3: Streaming results to stdout...");
                }
                io::write_tsv(rx, &db.names, &db.names, output.as_deref()).unwrap();
            });

            let duration = start_time.elapsed();
            info!(
                "Matrix completed successfully in {:.2} seconds!",
                duration.as_secs_f64()
            );
        }

        Commands::Search {
            genomes,
            db: target_db_path,
            output,
            level,
            strategy,
            concat,
            batch_size,
            threads,
        } => {
            info!("Starting zani Search matrix calculation...");

            let start_time = Instant::now();
            let strat: CompressionStrategy = strategy.into();

            info!("Step 1/3: Loading Target Database from disk...");
            let target_db = match target_db_path.to_str() {
                Some(p) => Database::load_from_disk(p, level, strat)?,
                None => return Err(anyhow::anyhow!("Invalid database path.")),
            };

            info!("Step 2/3: Compiling Query genomes into dictionaries...");
            let mut query_db = Database::new();
            for path in &genomes {
                if !path.exists() {
                    warn!("File not found, skipping: {:?}", path);
                    continue;
                }
                query_db.add_fasta(path, level, strat, concat);
            }

            if query_db.is_empty() {
                return Err(anyhow::anyhow!("No valid query genomes loaded. Exiting."));
            }
            info!(
                "Successfully loaded {} targets and compiled {} queries.",
                target_db.len(),
                query_db.len()
            );

            info!("Step 3/3: Executing multi-threaded NxM compression matrix...");
            let engine = ZaniEngine::new()
                .with_level(level)
                .with_batch_size(batch_size)
                .with_threads(threads)
                .with_strategy(strat);

            let (tx, rx) = std::sync::mpsc::sync_channel(100);

            std::thread::scope(|s| {
                s.spawn(|| {
                    engine.query_matrix_batched(&target_db, &query_db, tx);
                });

                if let Some(p) = &output {
                    info!("Streaming results to disk at {:?}", p);
                } else {
                    info!("Streaming results to stdout...");
                }
                io::write_tsv(rx, &target_db.names, &query_db.names, output.as_deref()).unwrap();
            });

            let duration = start_time.elapsed();
            info!(
                "Matrix completed successfully in {:.2} seconds!",
                duration.as_secs_f64()
            );
        }

        Commands::Build {
            genomes,
            db: db_path,
            level,
            strategy,
            concat,
        } => {
            info!("Starting zani database compilation...");
            let strat: CompressionStrategy = strategy.into();

            let start_time = Instant::now();
            let mut db = Database::new();

            for path in &genomes {
                if !path.exists() {
                    warn!("File not found, skipping: {:?}", path);
                    continue;
                }
                db.add_fasta(path, level, strat, concat);
            }

            if db.is_empty() {
                return Err(anyhow::anyhow!("No valid genomes loaded. Exiting."));
            }

            info!("Successfully compiled {} sketches.", db.len());

            if let Some(p_str) = db_path.to_str() {
                db.save_to_disk(p_str)?;
                info!("Database written to: {}", p_str);
            } else {
                return Err(anyhow::anyhow!("Invalid database output path."));
            }

            let duration = start_time.elapsed();
            info!(
                "Build completed successfully in {:.2} seconds!",
                duration.as_secs_f64()
            );
        }
    }

    Ok(())
}
