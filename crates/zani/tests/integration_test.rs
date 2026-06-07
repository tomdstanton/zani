use std::path::PathBuf;
use zani::{Database, ZaniEngine};

fn get_test_data_dir() -> PathBuf {
    let mut path = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    path.push("../../test_data");
    path
}

#[test]
fn test_all_vs_all_matrix() {
    let test_dir = get_test_data_dir();
    let mut db = Database::new();
    
    if test_dir.exists() {
        for entry in std::fs::read_dir(test_dir).unwrap() {
            let entry = entry.unwrap();
            let path = entry.path();
            if path.extension().and_then(|s| s.to_str()) == Some("gz") {
                db.add_fasta(&path, 3, true);
            }
        }
    }
    
    let engine = ZaniEngine::new().with_level(3).with_batch_size(100);
    let (tx, rx) = std::sync::mpsc::sync_channel(100);

    std::thread::scope(|s| {
        s.spawn(|| {
            engine.query_matrix_batched(&db, &db, tx);
        });

        for batch in rx {
            for i in 0..batch.len() {
                if batch.query_id == batch.target_ids[i] {
                    println!("query_id: {}, c_y_given_x: {}, nt_match: {}", batch.query_id, batch.ncd_similarity[i], batch.nt_match[i]);
                }
            }
        }
    });
}
