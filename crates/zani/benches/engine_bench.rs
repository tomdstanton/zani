use criterion::{criterion_group, criterion_main, Criterion};
use std::path::PathBuf;
use zani::{Database, ZaniEngine};

fn get_test_data_dir() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("../../tests/data");
    path
}

fn bench_database_compilation(c: &mut Criterion) {
    let test_dir = get_test_data_dir();
    let files: Vec<PathBuf> = std::fs::read_dir(test_dir).unwrap()
        .map(|res| res.unwrap().path())
        .filter(|p| p.extension().and_then(|s| s.to_str()) == Some("gz"))
        .collect();

    c.bench_function("db_compilation_5_genomes", |b| {
        b.iter(|| {
            let mut db = Database::new();
            for path in &files {
                db.add_fasta(path, 3, true);
            }
            db
        })
    });
}

fn bench_matrix_execution(c: &mut Criterion) {
    let test_dir = get_test_data_dir();
    let mut db = Database::new();
    
    for entry in std::fs::read_dir(test_dir).unwrap() {
        let entry = entry.unwrap();
        let path = entry.path();
        if path.extension().and_then(|s| s.to_str()) == Some("gz") {
            db.add_fasta(&path, 3, true);
        }
    }

    let engine = ZaniEngine::new().with_level(3).with_batch_size(100);

    c.bench_function("matrix_5x5", |b| {
        b.iter(|| {
            let (tx, rx) = std::sync::mpsc::sync_channel(100);
            std::thread::scope(|s| {
                s.spawn(|| {
                    engine.query_matrix_batched(&db, &db, tx);
                });
                
                // Drain the channel to simulate full matrix traversal
                for _batch in rx { }
            });
        })
    });
}

criterion_group!(benches, bench_database_compilation, bench_matrix_execution);
criterion_main!(benches);
