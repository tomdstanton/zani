import glob
import time
import zani

print("1. Initializing Database...")
db = zani.Database()

fastas = glob.glob("../test_genomes/ncbi_dataset/data/*/*.fna")
print(f"Found {len(fastas)} FASTA files.")

start = time.time()
for f in fastas:
    db.add_fasta(f, level=3, concat=True)
print(f"Added genomes in {time.time() - start:.2f}s")
print(f"Database contains {len(db)} records.")

print("\n2. Initializing Engine...")
engine = zani.Engine(compression_level=3, batch_size=10000)

print("\n3. Running All-vs-All Matrix...")
start = time.time()
engine.all_vs_all(db, "py_matrix.tsv")
print(f"Matrix computed and saved in {time.time() - start:.2f}s")
