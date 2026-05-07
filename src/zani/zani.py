from io import IOBase
from pathlib import Path
from re import compile as re_compile
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, wait, FIRST_COMPLETED
from compression.bz2 import open as bz_open
from compression.gzip import open as gz_open
from compression.lzma import open as lz_open
from compression.zstd import (CompressionParameter, Strategy, ZstdDict, ZstdCompressor, compress, train_dict,
                              open as zst_open)
from typing import IO, Iterable, NamedTuple, Generator
from itertools import islice


# Classes --------------------------------------------------------------------------------------------------------------
class SeqFileError(Exception):
    pass


class FastaFile:
    _REGEX = re_compile(r'\.(?P<ext>f(asta|a|na|fn|as|aa))\.?(?P<comp>(gz|bz2|xz|zst))?$')
    _DELETE_NEWLINES = b'\r\n'
    _OPENERS = {'gz': gz_open, 'bz2': bz_open, 'xz': lz_open, 'zst': zst_open}
    __slots__ = ('_file', '_fh', '_opener', '_is_stream', '_concat', '_base_name')

    def __init__(self, file: str | Path | IO[bytes], concat: bool = True):
        self._file = file
        self._fh = None
        self._opener = open
        self._is_stream = isinstance(file, IOBase)
        self._concat = concat
        self._base_name = None

        if self._is_stream:
            self._fh = file
        else:
            self._file = Path(file)
            original_name = self._file.name
            lower_name = original_name.lower()

            # Attempt to automatically determine type and compression
            fasta_match = self._REGEX.search(lower_name)
            match = fasta_match
            # Extract the pristine base name by slicing off the matched regex
            # e.g., "E_coli_K12.fasta.gz" -> b"E_coli_K12"
            if match:
                self._base_name = original_name[:match.start()].encode('utf-8')
            else:
                self._base_name = self._file.stem.encode('utf-8')

            # Determine correct file opener based on compression extension
            if match and match.group('comp'):
                self._opener = self._OPENERS[match.group('comp')]

    def __iter__(self) -> Generator[tuple[bytes, bytes], None, None]:
        if self._fh is None:
            raise SeqFileError("File not open. Use the context manager ('with SeqFile(...) as f:').")
        content = self._fh.read()
        if not content:
            return

        # Using an iterator prevents slicing [1:] which would clone the list in RAM
        records = iter(content.split(b'>'))
        next(records, None)  # Swiftly discard the empty chunk before the first '>'

        # Mode 1: Discrete Pairwise Mode (Yield each record individually)
        if not self._concat:
            for record in records:
                if record:
                    header, _, seq_data = record.partition(b'\n')
                    name = header.split(None, 1)[0].rstrip()
                    seq = seq_data.translate(None, delete=FastaFile._DELETE_NEWLINES)
                    yield name, seq

        # Mode 2: Whole Seq Assembly Mode (Concatenate into a single sequence)
        else:
            try:
                first_record = next(records)
            except StopIteration:
                return  # Empty or invalid fasta
                
            header, _, seq_data = first_record.partition(b'\n')
            first_name = self._base_name or header.split(None, 1)[0].rstrip()
            
            # Seed the chunks array with the first sequence
            seq_chunks = [seq_data.translate(None, delete=FastaFile._DELETE_NEWLINES)]
            
            # Use a C-optimized list comprehension to consume the rest!
            # This is measurably faster than a generator expression as it executes entirely in C.
            seq_chunks.extend([
                record.partition(b'\n')[2].translate(None, delete=FastaFile._DELETE_NEWLINES)
                for record in records if record
            ])

            if seq_chunks:
                yield first_name, b"".join(seq_chunks)

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


class SeqRecord:
    """An ultra-lean data container representing the sequence of a single sample"""
    __slots__ = ('_name', '_view', '_size')

    def __init__(self, seq: bytes, name: bytes):
        self._name = name
        self._view = memoryview(seq)
        self._size = len(seq)

    @classmethod
    def from_file(cls, file: str | Path, concat: bool = True) -> 'SeqRecord':
        """Creates a Seq instance directly from a sequence file.
        
        Reads the first sequence (or concatenated sequences if concat=True) 
        from the specified file.
        
        Args:
            file (str | Path): The path to the sequence file.
            concat (bool, optional): If True, concatenates all sequences (FASTA only). Defaults to True.
            
        Returns:
            SeqRecord: A new Seq instance.
            
        Raises:
            SeqFileError: If the file is empty or cannot be parsed.
        """
        with FastaFile(file, concat=concat) as seq_file:
            try:
                name, seq = next(iter(seq_file))
                return cls(seq, name)
            except StopIteration:
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


class Reference(SeqRecord):
    """
    The ZANI equivalent of a 'sketch'.
    Holds the pre-trained Zstd dictionary and baseline metrics.
    """
    __slots__ = ('_dict', '_baseline_size', '_opts')

    def __init__(self, seq: bytes, name: bytes, opts: dict | None = None):
        super().__init__(seq, name)
        self._opts = opts or self._get_default_dna_options()
        
        # Zstandard dictionary training expects multiple independent samples (usually small chunks),
        # not a single massive contiguous block. We can slice the memoryview 
        # without copying bytes to efficiently create these samples!
        chunk_size = 65536  # 64KB chunks
        chunks = [self._view[i:i + chunk_size] for i in range(0, self._size, chunk_size)]
        
        # 1MB dictionary size is optimal for bacterial seq_records. We cap it to half the seq_record 
        # size just in case an unusually small input sequence is provided.
        target_dict_size = min(1024 * 1024, max(self._size // 2, 1024))
        self._dict = train_dict(chunks, dict_size=target_dict_size)
        # Calculate the baseline compressed size using the module-level function
        self._baseline_size = len(compress(self._view, options=self._opts))

    @staticmethod
    def _get_default_dna_options() -> dict:
        """Returns default options optimized for a 4-letter alphabet."""
        return {
            CompressionParameter.compression_level: 3,
            CompressionParameter.min_match: 7,  # Max allowed by zstd; helps ignore short random DNA collisions
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
        reference_name (bytes): The name of the reference seq_record.
        seq_record_name (bytes): The name of the queried sample.
        ncd (float): The Normalized Compression Distance ratio.
    """
    reference_name: bytes  # The name of the reference
    seq_record_name: bytes  # The name of the sample
    ncd: float  # Normalized compression distance

    def write(self, fh: IO[bytes]) -> int:
        return fh.write(b'%b\t%b\t%f\n' % (self.reference_name, self.seq_record_name, self.ncd))


class ZaniEngine:
    """Execution engine for running highly-parallelized ZANI comparisons.

    Utilizes a bounded worker queue and thread-local C-level compression 
    contexts to maximize GIL-released throughput while maintaining a flat memory profile.

    Examples:
        >>> ref = Reference(b"ACGT", b"Reference1")
        >>> engine = ZaniEngine(ref, max_workers=4)
        >>> queries = [SeqRecord(b"ACGA", b"Query1"), SeqRecord(b"TCGT", b"Query2")]
        >>> with engine:
        ...     for result in engine.query(queries):
        ...         print(result.seq_record_name, result.ncd)
    """
    __slots__ = ('_ref', '_executor', '_max_workers', '_thread_local')

    def __init__(self, ref: Reference, max_workers: int | None = None):
        """Initializes the ZaniEngine.

        Args:
            ref (Reference): The reference seq_record sketch to compare against.
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

    def _compress_only(self, query: SeqRecord | FastaFile | str | Path) -> ZaniResult:
        """Executes I/O, parsing, and compression concurrently inside the worker thread."""
        
        # Resolve the input into a Seq object before compressing
        if isinstance(query, FastaFile):
            with query as seq_file:
                try:
                    name, seq = next(iter(seq_file))
                    query = SeqRecord(seq, name)
                except StopIteration:
                    raise SeqFileError(f"No sequences found in {query._file}")

        elif not isinstance(query, SeqRecord):
            query = SeqRecord.from_file(query)

        # Grab the thread's persistent C-context
        comp = self._get_thread_compressor()

        # Compress using the persistent workspace
        compressed_bytes = comp.compress(query.view)

        c_y_given_x = len(compressed_bytes)
        c_x = self._ref.baseline_size
        
        # Estimate the query's baseline compressed size (C(y)) based on uncompressed length ratio.
        # This is incredibly fast and avoids compressing the query twice!
        c_y = c_x * (query.size / self._ref.size)
        
        # Use the standard Normalized Compression Distance (NCD) formula:
        # We approximate C(x,y) as C(x) + C(y|x)
        ncd = (c_x + c_y_given_x - min(c_x, c_y)) / max(c_x, c_y)
        
        # Clamp to [0, 1] to hide tiny compression framing overheads
        ncd = max(0.0, min(1.0, ncd))

        return ZaniResult(self._ref.name, query.name, ncd)

    def query(self, queries: Iterable[SeqRecord | FastaFile | str | Path]) -> Generator[ZaniResult, None, None]:
        """Compresses an iterable of Seqs against the reference, yielding results lazily.

        Uses a bounded worker queue to maintain high throughput without exhausting memory.

        Args:
            queries (Iterable[SeqRecord | FastaFile | str | Path]): An iterable yielding sequences or file paths.

        Yields:
            ZaniResult: The comparison result for each queried Seq.
        """
        queries_iter = iter(queries)

        # Keep workers fed with 2x jobs to avoid thread starvation,
        # whilst preventing memory spikes from loading all Seqs at once.
        max_in_flight = (self._max_workers or 32) * 2

        futures = {self.executor.submit(self._compress_only, query) for query in islice(queries_iter, max_in_flight)}

        while futures:
            done, futures = wait(futures, return_when=FIRST_COMPLETED)

            for future in done:
                yield future.result()

            # Refill the pipeline using set.update and a generator expression (executes loop in C)
            futures.update(
                self.executor.submit(self._compress_only, query) 
                for query in islice(queries_iter, len(done))
            )


class ZaniWriter:
    """A thread-safe writer for ZaniResult objects to a TSV file.

    This class is designed to be used as a context manager. It handles opening
    and closing the file, and ensures that concurrent writes from multiple
    threads do not corrupt the output file.

    Examples:
        >>> results = [ZaniResult(b"ref", b"g1", 0.5), ZaniResult(b"ref", b"g2", 0.6)]
        >>> with ZaniWriter("results.tsv") as writer:
        ...     for result in results:
        ...         writer.write(result)
    """
    _HEADER = b"reference\tseq_record\tncd\n"
    __slots__ = ('_file', '_fh', '_lock', '_header', '_is_stream')

    def __init__(self, file: str | Path | IO[bytes] | None = None, header: bool = True):
        """Initializes the ZaniWriter.

        Args:
            file (str | Path | IO[bytes] | None): The path or stream for output. Defaults to sys.stdout.buffer.
            header (bool, optional): If True, writes a header row. Defaults to True.
        """
        self._file = file if file is not None else sys.stdout.buffer
        self._header = header
        self._is_stream = isinstance(self._file, IOBase)
        self._fh: IO[bytes] | None = self._file if self._is_stream else None
        self._lock = threading.Lock()

    def __enter__(self) -> 'ZaniWriter':
        """Opens the file for writing and writes the header."""
        if not self._is_stream:
            self._fh = open(self._file, mode='wb')
        if self._header and self._fh is not None:
            self._fh.write(self._HEADER)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Closes the file, or flushes if it's a stream."""
        if not self._is_stream and self._fh:
            self._fh.close()
            self._fh = None
        elif self._is_stream and self._fh:
            self._fh.flush()

    def write(self, result: ZaniResult):
        """Writes a ZaniResult to the file in a thread-safe manner."""
        if self._fh is None:
            raise IOError("ZaniWriter is not open. Use a 'with' statement.")
        with self._lock:
            result.write(self._fh)


# Functions ------------------------------------------------------------------------------------------------------------
def run_all_vs_all(seq_records: list[SeqRecord], output_tsv: str | Path, max_workers: int = None):
    """
    Computes a full N x N pairwise distance matrix and streams it to disk.
    Requires ZERO external dependencies.
    """
    with ZaniWriter(output_tsv) as writer:
        for g_ref in seq_records:
            # 1. Initialize the reference (Trains the Zstd dictionary)
            ref = Reference(seq=g_ref.view, name=g_ref.name)

            # 2. Spin up the engine
            with ZaniEngine(ref, max_workers=max_workers) as engine:

                # 3. Stream the queries through the engine directly to disk
                # Because engine.query() is a generator, we never hold more
                # than `max_workers * 2` results in RAM at any given time!
                for result in engine.query(seq_records):
                    writer.write(result)