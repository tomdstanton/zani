"""
Zstandard-based Average Nucleotide Identity (ANI) High-Performance Engine.

This module provides a pure-Rust, GIL-releasing engine for massively parallel
genomic distance calculation using Normalized Compression Distance (NCD).
"""

from typing import Union

class PyFileNotFoundError(Exception):
    """Exception raised when a specified FASTA/FASTQ file is not found."""
    pass

class Database:
    """
    A high-performance Zstandard compression dictionary database.

    The `Database` acts exactly like a native Python object, but holds pure Rust memory.
    It manages the training and storage of Zstandard dictionaries for a collection of genomes.

    Example:
        >>> import zani
        >>> db = zani.Database()
        >>> db.add_fasta("genome.fasta", level=3)
        >>> print(len(db))
        1
    """

    def __init__(self) -> None:
        """
        Initialize an empty Database.
        """
        ...

    def __len__(self) -> int:
        """
        Number of genomes currently loaded in the database.

        Returns:
            int: The number of compiled dictionaries in the database.
        """
        ...

    def add_fasta(self, filepath: str, level: int = 3, strategy: str = "lazy2", concat: bool = True) -> None:
        """
        Reads a FASTA/FASTQ file, compiles the sequences, and adds them to the database.

        Args:
            filepath (str): Path to the `.fna`, `.fasta`, or `.fastq` file.
            level (int, optional): Zstandard compression level. Defaults to 3. (Range: 1-19).
            strategy (str, optional): LZ77 match finding strategy. Defaults to "lazy2".
                Valid options: "auto", "fast", "dfast", "greedy", "lazy", "lazy2",
                "btlazy2", "btopt", "btultra", "btultra2".
            concat (bool, optional): If True, concatenates all sequences in the file
                into a single genome (with a Chimeric Firewall). Defaults to True.

        Raises:
            PyFileNotFoundError: If the specified FASTA file does not exist.
        """
        ...

    def add_sequence(self, identifier: bytes, sequence: bytes, level: int = 3, strategy: str = "lazy2") -> None:
        """
        Compiles a raw in-memory sequence and adds it to the database.

        Args:
            identifier (bytes): The raw bytes of the genome name or identifier.
            sequence (bytes): The raw bytes of the sequence data.
            level (int, optional): Zstandard compression level. Defaults to 3.
            strategy (str, optional): LZ77 match finding strategy. Defaults to "lazy2".
        """
        ...

    def save(self, filepath: str) -> None:
        """
        Saves the compiled database to disk as a binary `.zani` file.

        Args:
            filepath (str): Path to save the `.zani` database.

        Raises:
            IOError: If writing to the disk fails.
        """
        ...

    @staticmethod
    def load(filepath: str, level: int = 3, strategy: str = "lazy2") -> 'Database':
        """
        Load a previously compiled database from disk.

        Args:
            filepath (str): Path to the compiled `.zani` file.
            level (int, optional): The original compression level used during creation. Defaults to 3.
            strategy (str, optional): The original compression strategy. Defaults to "lazy2".

        Returns:
            Database: The loaded database object.

        Raises:
            IOError: If reading from the disk fails.
        """
        ...


class Engine:
    """
    Executes the multithreaded Rayon matrix and writes the TSV to disk.

    This engine completely bypasses the Python interpreter and drops the
    Global Interpreter Lock (GIL) to achieve maximum CPU utilization across all cores.

    Example:
        >>> import zani
        >>> db = zani.Database()
        >>> # ... load genomes ...
        >>> engine = zani.Engine(compression_level=3, batch_size=10_000)
        >>> engine.all_vs_all(db, "results.tsv")
    """

    def __init__(self, compression_level: int = 3, batch_size: int = 10000, threads: int = 0, strategy: str = "lazy2") -> None:
        """
        Initializes the computation Engine.

        Args:
            compression_level (int, optional): Zstandard compression level to use during
                execution. Defaults to 3. (Range: 1-19).
            batch_size (int, optional): Number of rows to process before flushing to disk.
                Defaults to 10000.
            threads (int, optional): Number of worker threads. Defaults to 0 (auto-detect all cores).
            strategy (str, optional): Compression strategy. Defaults to "lazy2".
        """
        ...

    def all_vs_all(self, db: Database, output_filepath: str) -> None:
        """
        Computes the All-vs-All pairwise distance matrix for a database.

        Releases the Python GIL and executes purely in Rust across all CPU cores.

        Args:
            db (Database): The target database to compute against itself.
            output_filepath (str): Path to output the resulting TSV matrix.

        Raises:
            IOError: If writing the TSV fails.
            ValueError: If the database is empty.
        """
        ...

    def search(self, db: Database, queries: Database, output_filepath: str) -> None:
        """
        Computes the structural distances between a query database and a target database.

        Releases the Python GIL and executes purely in Rust across all CPU cores.

        Args:
            db (Database): The reference target database.
            queries (Database): The query database to search against the target.
            output_filepath (str): Path to output the resulting TSV matrix.

        Raises:
            IOError: If writing the TSV fails.
            ValueError: If either database is empty.
        """
        ...
