use clap::{Parser, Subcommand};
use std::path::PathBuf;
use std::time::Instant;
use log::{info, warn};

// Import your pure Rust engine from the lib.rs file
use zani::{ZaniEngine, Database, io};

/// ZANI: High-performance Lempel-Ziv structural genome distance.
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
    AllVsAll {
        /// Input FASTA/FASTQ files to compare
        #[arg(required = true, num_args = 1..)]
        genomes: Vec<PathBuf>,

        /// Output TSV file path
        #[arg(short, long, default_value = "zani_matrix.tsv")]
        output: PathBuf,

        /// Zstandard compression level (1-19). Higher is more accurate but slower.
        #[arg(short, long, default_value_t = 3)]
        level: i32,

        /// Treat multi-record FASTAs as a single concatenated genome
        #[arg(short, long, default_value_t = true)]
        concat: bool,

        /// How many rows to buffer before writing to disk
        #[arg(short, long, default_value_t = 10_000)]
        batch_size: usize,
    },
    
    /// (Future Command Idea) Pre-compile a database of genomes to disk
    Build {
        /// Input FASTA files
        #[arg(required = true)]
        genomes: Vec<PathBuf>,
        
        /// Output database file (.zani)
        #[arg(short, long)]
        output: PathBuf,
    }
}

fn main() -> anyhow::Result<()> {
    // Initialize the logger so we can print beautiful [INFO] tags
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    // Clap parses the command line arguments instantly and handles --help / --version
    let cli = Cli::parse();

    match cli.command {
        Commands::AllVsAll { genomes, output, level, concat, batch_size } => {
            info!("Starting ZANI All-vs-All matrix calculation...");
            info!("Compression Level: {}", level);
            info!("Input Files: {}", genomes.len());

            let start_time = Instant::now();

            // 1. Build the Database (The FWD + N + RC Training Pool)
            info!("Step 1/3: Compiling genomes into Zstandard dictionaries...");
            let mut db = Database::new();
            
            for path in &genomes {
                if !path.exists() {
                    warn!("File not found, skipping: {:?}", path);
                    continue;
                }
                // Call the function from your lib.rs
                db.add_fasta(path, level, concat);
            }

            if db.is_empty() {
                return Err(anyhow::anyhow!("No valid genomes loaded. Exiting."));
            }
            
            info!("Successfully compiled {} sketches.", db.len());
            info!("Step 2/3: Executing multi-threaded NxN compression matrix...");

            // 2. Initialize the Engine
            let engine = ZaniEngine::new()
                .with_level(level)
                .with_batch_size(batch_size);

            // 3. Start the Rayon Matrix computation
            // Note: For All-vs-All, the `queries` is exactly the same as the `db`
            let (tx, rx) = std::sync::mpsc::sync_channel(100);
            
            std::thread::scope(|s| {
                s.spawn(|| {
                    engine.query_matrix_batched(&db, &db, tx);
                });

                // 4. Stream the results directly to disk
                info!("Step 3/3: Streaming results to disk at {:?}", output);
                io::write_tsv(rx, &db.names, &db.names, &output).unwrap();
            });

            let duration = start_time.elapsed();
            info!("Matrix completed successfully in {:.2} seconds!", duration.as_secs_f64());
        }
        
        Commands::Build { genomes: _genomes, output: _output } => {
            // Future implementation: compile DB and call db.save_to_disk(&output)
            info!("Build command not yet implemented.");
        }
    }

    Ok(())
}