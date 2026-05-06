from pathlib import Path
from zani import ZaniEngine, Reference, Genome, TsvWriter


def test_engine_with_writer(tmp_path):
    """
    Tests the full ZaniEngine pipeline, writing results with TsvWriter.
    """
    # NOTE: This test uses a hardcoded path. For a real test suite,
    # it's better to use a small, version-controlled dataset.
    genome_dir = Path('/Users/tsta0015/Programming/KpSC_ST23_K1')
    if not genome_dir.exists():
        # Skip the test if the data directory isn't present
        import pytest
        pytest.skip(f"Test data directory not found: {genome_dir}")

    genome_paths = list(genome_dir.glob("*.fasta.gz"))
    output_file = tmp_path / "results.tsv"

    # Use the first genome as the reference
    reference = Reference.from_file(genome_paths[0])

    # Create a generator for the query genomes to save memory
    query_genomes = (Genome.from_file(p) for p in genome_paths[1:])

    with TsvWriter(output_file) as writer, ZaniEngine(reference, max_workers=4) as engine:
        for result in engine(query_genomes):
            writer.write(result)

    # Verify that the output file was created and has content
    assert output_file.exists()
    content = output_file.read_text()
    assert content.startswith("reference\tgenome\tncd")
    assert len(content.strip().split('\n')) > 1  # Check for at least one data row