from io import IOBase
from pathlib import Path
from enum import Enum, auto
from re import compile as re_compile
import threading
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from compression.bz2 import open as bz_open
from compression.gzip import open as gz_open
from compression.lzma import open as lz_open
from compression.zstd import (CompressionParameter, Strategy, ZstdDict, ZstdCompressor, compress, train_dict,
                              open as zst_open)
from typing import IO, Iterable, NamedTuple, Generator
from itertools import islice


class SeqFileType(Enum):
    FASTA = auto()
    FASTQ = auto()
    
    
class SeqFileError(Exception):
    pass


class SeqFile:
    _REGEXES = {
        SeqFileType.FASTA: re_compile(r'\.(?P<ext>f(asta|a|na|fn|as|aa))\.?(?P<comp>(gz|bz2|xz|zst))?$'),
        SeqFileType.FASTQ: re_compile(r'\.(?P<ext>f(astq|q))\.?(?P<comp>(gz|bz2|xz|zst))?$')
    }
    _DELETE_NEWLINES = b'\r\n'
    _OPENERS = {'gz': gz_open, 'bz2': bz_open, 'xz': lz_open, 'zst': zst_open}
    __slots__ = ('_file', '_file_type', '_fh', '_opener', '_is_stream', '_concat', '_base_name')

    def __init__(self, file: str | Path | IO[bytes], file_type: SeqFileType | None = None, concat: bool = True):
        self._file = file
        self._file_type = file_type
        self._fh = None
        self._opener = open
        self._is_stream = isinstance(file, IOBase)
        self._concat = concat
        self._base_name = None

        if self._is_stream:
            if not self._file_type:
                raise SeqFileError('Must provide a SeqFileType when file is a stream')
            self._fh = file
        else:
            self._file = Path(file)
            original_name = self._file.name
            lower_name = original_name.lower()

            # Attempt to automatically determine type and compression
            fasta_match = self._REGEXES[SeqFileType.FASTA].search(lower_name)
            fastq_match = self._REGEXES[SeqFileType.FASTQ].search(lower_name)
            match = fasta_match or fastq_match

            # Extract the pristine base name by slicing off the matched regex
            # e.g., "E_coli_K12.fasta.gz" -> b"E_coli_K12"
            if match:
                self._base_name = original_name[:match.start()].encode('utf-8')
            else:
                self._base_name = self._file.stem.encode('utf-8')

            # Prioritize explicitly passed file_type, fallback to regex detection
            if not self._file_type:
                if fasta_match and fasta_match.group('ext'):
                    self._file_type = SeqFileType.FASTA
                else:
                    raise SeqFileError(f"Could not infer file type from filename: {original_name}")

            # Determine correct file opener based on compression extension
            if match and match.group('comp'):
                self._opener = self._OPENERS[match.group('comp')]

    def __iter__(self) -> Generator[tuple[bytes, bytes], None, None]:
        if self._fh is None:
            raise SeqFileError("File not open. Use the context manager ('with SeqFile(...) as f:').")

        if self._file_type == SeqFileType.FASTA:
            # Pass our pristine base name down to the static parser
            yield from self._parse_fasta(self._fh, self._concat, self._base_name)
        elif self._file_type == SeqFileType.FASTQ:
            yield from self._parse_fastq(self._fh)

    def __enter__(self):
        """Opens the file for reading."""
        if not self._is_stream:
            self._fh = self._opener(self._file, 'rb')
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Closes the file for reading."""
        if not self._is_stream and self._fh is not None:
            self._fh.close()
            self._fh = None

    @staticmethod
    def _parse_fasta(fh: IO[bytes], concat: bool, base_name: bytes | None) -> Generator[
        tuple[bytes, bytes], None, None]:
        content = fh.read()
        if not content:
            return

        records = content.split(b'>')

        # Mode 1: Discrete Pairwise Mode (Yield each record individually)
        if not concat:
            for record in records:
                if not record:
                    continue
                header, _, seq_data = record.partition(b'\n')
                name = header.split(None, 1)[0].rstrip()
                seq = seq_data.translate(None, delete=SeqFile._DELETE_NEWLINES)
                yield name, seq

        # Mode 2: Whole Genome Assembly Mode (Concatenate into a single sequence)
        else:
            first_name = base_name  # Initialize with our stripped filename!
            seq_chunks = []

            for record in records:
                if not record:
                    continue
                header, _, seq_data = record.partition(b'\n')

                # Fallback: If it's a stream, base_name is None, so we grab the first contig header
                if not first_name:
                    first_name = header.split(None, 1)[0].rstrip()

                seq_chunks.append(seq_data.translate(None, delete=SeqFile._DELETE_NEWLINES))

            if first_name and seq_chunks:
                # b"".join is executed in highly optimized C
                yield first_name, b"".join(seq_chunks)

    @staticmethod
    def _parse_fastq(fh: IO[bytes]) -> Generator[tuple[bytes, bytes], None, None]:
        iterator = iter(fh)
        for line in iterator:
            if not line:
                continue
            # 64 corresponds to ASCII '@'
            if line[0] == 64:
                # Split at first whitespace to isolate the ID, stripping off the '@'
                name = line[1:].split(None, 1)[0].rstrip()
                try:
                    seq = next(iterator).rstrip()
                    next(iterator)  # Skip the '+' line
                    next(iterator)  # Skip the quality scores line
                    yield name, seq
                except StopIteration:
                    break  # Handles gracefully truncated FASTQ streams


class Genome:
    """An ultra-lean data container representing the genome sequence of a single sample"""
    __slots__ = ('_name', '_view', '_size')

    def __init__(self, seq: bytes, name: bytes):
        self._name = name
        self._view = memoryview(seq)
        self._size = len(seq)

    @classmethod
    def from_file(cls, file: str | Path, file_type: SeqFileType | None = None, concat: bool = True) -> 'Genome':
        """Creates a Genome instance directly from a sequence file.
        
        Reads the first sequence (or concatenated sequences if concat=True) 
        from the specified file.
        
        Args:
            file (str | Path): The path to the sequence file.
            file_type (SeqFileType | None, optional): The type of sequence file. Defaults to None.
            concat (bool, optional): If True, concatenates all sequences (FASTA only). Defaults to True.
            
        Returns:
            Genome: A new Genome instance.
            
        Raises:
            SeqFileError: If the file is empty or cannot be parsed.
        """
        with SeqFile(file, file_type=file_type, concat=concat) as seq_file:
            for name, seq in seq_file:
                return cls(seq, name)
        raise SeqFileError(f"No sequences found in {file}")

    @property
    def name(self) -> bytes:
        return self._name

    @property
    def view(self) -> memoryview:
        return self._view

    @property
    def size(self) -> int:
        return self._size


class Reference(Genome):
    """
    The ZANI equivalent of a 'sketch'.
    Holds the pre-trained Zstd dictionary and baseline metrics.
    """
    __slots__ = ('_dict', '_baseline_size', '_opts')

    def __init__(self, seq: bytes, name: bytes, opts: dict | None = None):
        super().__init__(seq, name)
        self._opts = opts or self._get_default_dna_options()
        # train_dict expects an iterable of samples, so we wrap our view in a list.
        # 1MB dictionary size is generally optimal for bacterial genomes.
        self._dict = train_dict([self._view], dict_size=1024 * 1024)
        # Calculate the baseline compressed size using the module-level function
        self._baseline_size = len(compress(self._view, options=self._opts))

    @staticmethod
    def _get_default_dna_options() -> dict:
        """Returns default options optimized for a 4-letter alphabet."""
        return {
            CompressionParameter.compression_level: 3,
            CompressionParameter.min_match: 14,  # Ignore random 4-base DNA collisions
            CompressionParameter.hash_log: 15,
            CompressionParameter.strategy: Strategy.lazy2
        }

    @property
    def zstd_dict(self) -> ZstdDict:
        return self._dict

    @property
    def baseline_size(self) -> int:
        return self._baseline_size

    @property
    def opts(self) -> dict:
        return self._opts


class ZaniResult(NamedTuple):
    """Data structure representing the results of a pairwise comparison.

    Attributes:
        reference_name (bytes): The name of the reference genome.
        genome_name (bytes): The name of the queried sample.
        ncd (float): The Normalized Compression Distance ratio.
    """
    reference_name: bytes  # The name of the reference
    genome_name: bytes  # The name of the sample
    ncd: float  # Normalized compression distance


class ZaniEngine:
    """Execution engine for running highly-parallelized ZANI comparisons.

    Utilizes a bounded worker queue and thread-local C-level compression 
    contexts to maximize GIL-released throughput while maintaining a flat memory profile.

    Examples:
        >>> ref = Reference(b"ACGT", b"Reference1")
        >>> engine = ZaniEngine(ref, max_workers=4)
        >>> queries = [Genome(b"ACGA", b"Query1"), Genome(b"TCGT", b"Query2")]
        >>> with engine:
        ...     for result in engine(queries):
        ...         print(result.genome_name, result.ncd)
    """
    __slots__ = ('_ref', '_executor', '_max_workers', '_thread_local')

    def __init__(self, ref: Reference, max_workers: int | None = None):
        """Initializes the ZaniEngine.

        Args:
            ref (Reference): The reference genome sketch to compare against.
            max_workers (int | None, optional): The maximum number of worker threads. 
                If None, defaults to the ThreadPoolExecutor's internal default. Defaults to None.
        """
        self._ref = ref
        self._max_workers = max_workers
        self._executor = None
        # Thread-local storage to cache our C-level ZstdCompressor contexts
        self._thread_local = threading.local()

    def _get_thread_compressor(self) -> ZstdCompressor:
        """
        Retrieves or initializes a reusable C-level compression context
        for the current worker thread.
        """
        if not hasattr(self._thread_local, 'compressor'):
            # Instantiated exactly once per worker thread!
            self._thread_local.compressor = ZstdCompressor(
                options=self._ref.opts,
                zstd_dict=self._ref.zstd_dict
            )
        return self._thread_local.compressor

    @property
    def executor(self) -> ThreadPoolExecutor:
        if self._executor is None:
            self._executor = ThreadPoolExecutor(max_workers=self._max_workers)
        return self._executor

    def __enter__(self):
        executor = self.executor  # Ensure we have an executor
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._executor is not None:
            # Ensures all executing threads complete before the context terminates
            self._executor.shutdown(wait=True, cancel_futures=True)

    def _compress_only(self, query: 'Genome') -> ZaniResult:
        """Optimized for GIL: Math only, GIL is released by zstd."""

        # Grab the thread's persistent C-context
        comp = self._get_thread_compressor()

        # Compress using the persistent workspace
        compressed_bytes = comp.compress(query.view)

        return ZaniResult(
            self._ref.name,
            query.name,
            len(compressed_bytes) / self._ref.baseline_size
        )

    def __call__(self, queries: Iterable[Genome]) -> Generator[ZaniResult, None, None]:
        """Compresses an iterable of Genomes against the reference, yielding results lazily.

        Uses a bounded worker queue to maintain high throughput without exhausting memory.

        Args:
            queries (Iterable[Genome]): An iterable (or generator) yielding Genome objects.

        Yields:
            ZaniResult: The comparison result for each queried Genome.
        """
        queries_iter = iter(queries)

        # Keep workers fed with 2x jobs to avoid thread starvation,
        # whilst preventing memory spikes from loading all Genomes at once.
        max_in_flight = (self._max_workers or 32) * 2

        futures = {self.executor.submit(self._compress_only, query) for query in islice(queries_iter, max_in_flight)}

        while futures:
            done, futures = wait(futures, return_when=FIRST_COMPLETED)

            for future in done:
                yield future.result()

            # Refill the pipeline with the exact number of jobs that just finished
            for query in islice(queries_iter, len(done)):
                futures.add(self.executor.submit(self._compress_only, query))


class TsvWriter:
    """A thread-safe writer for ZaniResult objects to a TSV file.

    This class is designed to be used as a context manager. It handles opening
    and closing the file, and ensures that concurrent writes from multiple
    threads do not corrupt the output file.

    Examples:
        >>> results = [ZaniResult(b"ref", b"g1", 0.5), ZaniResult(b"ref", b"g2", 0.6)]
        >>> with TsvWriter("results.tsv") as writer:
        ...     for result in results:
        ...         writer.write(result)
    """
    __slots__ = ('_file_path', '_fh', '_lock', '_header')

    def __init__(self, file_path: str | Path, header: bool = True):
        """Initializes the TsvWriter.

        Args:
            file_path (str | Path): The path to the output TSV file.
            header (bool, optional): If True, writes a header row. Defaults to True.
        """
        self._file_path = file_path
        self._header = header
        self._fh: IO[str] | None = None
        self._lock = threading.Lock()

    def __enter__(self) -> 'TsvWriter':
        """Opens the file for writing and writes the header."""
        self._fh = open(self._file_path, 'w', encoding='utf-8')
        if self._header:
            self._fh.write("reference\tgenome\tncd\n")
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Closes the file."""
        if self._fh:
            self._fh.close()

    def write(self, result: ZaniResult):
        """Writes a ZaniResult to the file in a thread-safe manner."""
        if self._fh is None:
            raise IOError("TsvWriter is not open. Use a 'with' statement.")

        line = f"{result.reference_name.decode('utf-8')}\t{result.genome_name.decode('utf-8')}\t{result.ncd:.6f}\n"
        with self._lock:
            self._fh.write(line)
