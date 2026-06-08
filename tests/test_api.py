import glob
import pytest
import zani

def test_database_loading():
    db = zani.Database()
    
    # Load the compressed test genomes
    fastas = glob.glob("tests/data/*.fna.gz")
    assert len(fastas) == 5, "Should have exactly 5 test genomes"
    
    for f in fastas:
        db.add_fasta(f, level=3, concat=True)
        
    assert len(db) == 5, "Database should contain 5 records"

def test_matrix_execution(tmp_path):
    db = zani.Database()
    fastas = glob.glob("tests/data/*.fna.gz")
    
    for f in fastas:
        db.add_fasta(f, level=3, concat=True)
        
    engine = zani.Engine(compression_level=3, batch_size=100)
    
    # Use pytest's tmp_path fixture for safe parallel testing
    output_tsv = str(tmp_path / "test_matrix.tsv")
    
    engine.all_vs_all(db, output_tsv)
    
    # Just verify the output file exists and is populated
    with open(output_tsv, "r") as f:
        lines = f.readlines()
        
    # 5x5 matrix = 25 rows + 1 header
    assert len(lines) == 26, "Matrix TSV should have 26 rows"
