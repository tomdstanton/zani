from dataclasses import dataclass, field
from io import IOBase
from pathlib import Path
from re import compile as re_compile
from concurrent.futures import ThreadPoolExecutor, as_completed
from compression.bz2 import open as bz_open
from compression.gzip import open as gz_open
from compression.lzma import open as lz_open
from compression.zstd import (CompressionParameter, Strategy, ZstdDict, compress, train_dict, open as zst_open, finalize_dict)
from typing import IO, Iterable, NamedTuple, Generator, Iterator
from itertools import batched
from pickle import dump, load, HIGHEST_PROTOCOL


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
            fasta_match = self._REGEX.search(lower_name)  # Extract the base name by slicing off the matched regex
            match = fasta_match  # e.g., "E_coli_K12.fasta.gz" -> b"E_coli_K12"
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

        if not self._concat:  # Mode 1: Discrete Pairwise Mode (Yield each record individually)
            for record in records:
                if record:
                    header, _, seq_data = record.partition(b'\n')
                    name = header.split(None, 1)[0].rstrip()
                    seq = seq_data.translate(None, delete=FastaFile._DELETE_NEWLINES)
                    yield name, seq

        else:  # Mode 2: Whole Seq Assembly Mode (Concatenate into a single sequence)
            try:
                first_record = next(records)
            except StopIteration:
                return  # Empty or invalid fasta
                
            header, _, seq_data = first_record.partition(b'\n')
            first_name = self._base_name or header.split(None, 1)[0].rstrip()
            
            # Seed the chunks array with the first sequence
            seq_chunks = [seq_data.translate(None, delete=FastaFile._DELETE_NEWLINES)]
            
            # Use a C-optimized list comprehension to consume the rest!
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


@dataclass(slots=True, frozen=True)
class SeqRecord:
    """An ultra-lean data container representing the sequence of a single sample"""
    name: bytes
    view: memoryview
    size: int

    @classmethod
    def from_file(cls, file: str | Path, concat: bool = True) -> 'SeqRecord':
        with FastaFile(file, concat=concat) as seq_file:
            try:
                name, seq = next(iter(seq_file))
                return cls(name, memoryview(seq), len(seq))
            except StopIteration:
                raise SeqFileError(f"No sequences found in {file}")


@dataclass(slots=True, frozen=True)
class Sketch:
    """Ultra-lean pre-computed mathematical state."""
    name: bytes
    size: int
    baseline_size: int
    compressed_size: int
    zstd_dict: ZstdDict


class SketchDatabase:
    """
    A serialized, lock-free collection of pre-compiled ZANI Sketches.
    Safely handles the extraction and re-hydration of C-level ZstdDicts
    for instantaneous database loading.
    """
    __slots__ = ('_sketches',)

    def __init__(self, sketches: Iterable[Sketch]):
        self._sketches = list(sketches)

    def add(self, sketch: Sketch):
        """Appends a pre-compiled Sketch to the database."""
        self._sketches.append(sketch)

    def __iter__(self) -> Iterator[Sketch]:
        return iter(self._sketches)

    def __len__(self) -> int:
        return len(self._sketches)

    def save(self, filepath: Path | str):
        """
        Serializes the database to disk.
        Strips the digested C-structs down to raw bytes for safe pickling.
        """
        # Extract the pure data, utilizing the Python 3.14 .dict_content attribute
        serializable_data = [
            (s.name, s.size, s.baseline_size, s.compressed_size, s.zstd_dict.dict_content)
            for s in self._sketches
        ]

        with open(filepath, 'wb') as f:
            dump(serializable_data, f, protocol=HIGHEST_PROTOCOL)

    @classmethod
    def load(cls, filepath: Path | str) -> 'SketchDatabase':
        """
        Loads a compiled database into RAM instantly.
        """
        with open(filepath, 'rb') as f:
            return cls(
                Sketch(name, size, base, comp, ZstdDict(dict_bytes))
                for name, size, base, comp, dict_bytes in load(f)
            )


class ZaniResult(NamedTuple):
    """Data structure representing the results of a pairwise comparison."""
    reference_name: bytes
    seq_record_name: bytes
    similarity: float

    def write(self, fh: IO[bytes]) -> int:
        return fh.write(b'%b\t%b\t%f\n' % (self.reference_name, self.seq_record_name, self.similarity))


@dataclass(slots=True, frozen=True)
class ZaniOpts:
    """Immutable configuration profile for ZANI compression parameters."""
    level: int
    min_match: int
    hash_log: int
    strategy: int

    @classmethod
    def dna(cls) -> 'ZaniOpts':
        """The standard profile highly optimized for 4-letter bacterial DNA."""
        return cls(
            level=3,
            min_match=7,  # Max allowed by zstd; ignores short random collisions
            hash_log=15,
            strategy=Strategy.lazy2
        )

    @classmethod
    def protein(cls) -> 'ZaniOpts':
        """A hypothetical profile optimized for 20-letter amino acid sequences."""
        return cls(
            level=5,
            min_match=4,
            hash_log=16,
            strategy=Strategy.btlazy2
        )

    @classmethod
    def fast(cls) -> 'ZaniOpts':
        """Optimized for maximum speed on massive databases at the cost of slight accuracy."""
        return cls(
            level=1,
            min_match=7,
            hash_log=14,
            strategy=Strategy.fast
        )

    def as_zstd_dict(self) -> dict:
        """Translates the dataclass into the raw C-enum dictionary required by zstd."""
        return {
            CompressionParameter.compression_level: self.level,
            CompressionParameter.min_match: self.min_match,
            CompressionParameter.hash_log: self.hash_log,
            CompressionParameter.strategy: self.strategy
        }


@dataclass(slots=True, frozen=True)
class ZaniEngine:
    """
    A pure, thread-safe, stateless execution environment.
    """
    opts: ZaniOpts = field(default_factory=ZaniOpts.dna)
    chunk_size: int = 65536

    def sketch(self, record: 'SeqRecord') -> 'Sketch':
        chunks = [record.view[i:i + self.chunk_size] for i in range(0, record.size, self.chunk_size)]
        target_dict_size = min(1024 * 1024, max(record.size // 2, 1024))
        raw_dict = train_dict(chunks, dict_size=target_dict_size)
        final_dict = finalize_dict(
            raw_dict,
            samples=chunks,
            dict_size=target_dict_size,
            level=self.opts.level
        )
        baseline = len(compress(record.view, options=self.opts.as_zstd_dict()))  # Calculate C(x)
        compressed_size = len(compress(  # Calculate C(x|x)
            record.view,
            options=self.opts.as_zstd_dict(),
            zstd_dict=final_dict.as_digested_dict
        ))
        return Sketch(record.name, record.size, baseline, compressed_size, final_dict)

    def compare(self, sketch: 'Sketch', target: 'SeqRecord') -> ZaniResult:
        if sketch.size == 0 or target.size == 0:
            return ZaniResult(sketch.name, target.name, float('nan'))
        compressed_bytes = compress(
            target.view,
            options=self.opts.as_zstd_dict(),
            zstd_dict=sketch.zstd_dict.as_digested_dict
        )
        c_y_given_x = len(compressed_bytes)
        calibrated_c_y_given_x = max(0, c_y_given_x - sketch.compressed_size)  # Zero-point calibration
        c_x = sketch.baseline_size
        c_y = c_x * (target.size / sketch.size)
        ncd = (c_x + calibrated_c_y_given_x - min(c_x, c_y)) / max(c_x, c_y)
        ncd = max(0.0, min(1.0, ncd))
        return ZaniResult(sketch.name, target.name, 1.0 - ncd)


class ZaniPipeline:
    """
    The multithreaded orchestrator for ZANI.
    """
    __slots__ = ('_engine', '_max_workers', '_executor', 'batch_size')

    def __init__(self, engine: ZaniEngine, max_workers: int | None = None, batch_size: int = 1000):
        self._engine = engine
        self._max_workers = max_workers
        self._executor = None
        self.batch_size = batch_size

    def __enter__(self):
        self._executor = ThreadPoolExecutor(max_workers=self._max_workers)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._executor is not None:
            self._executor.shutdown(wait=True)
            self._executor = None

    def __del__(self):
        if self._executor is not None:
            self._executor.shutdown(wait=False)

    @property
    def executor(self) -> ThreadPoolExecutor:
        if self._executor is None:
            self._executor = ThreadPoolExecutor(max_workers=self._max_workers)
        return self._executor

    @staticmethod
    def _all_vs_all_worker(engine: ZaniEngine, ref_record: SeqRecord, all_records: Iterable[SeqRecord]) -> list[ZaniResult]:
        sketch = engine.sketch(ref_record)
        return [engine.compare(sketch, target) for target in all_records]

    def all_vs_all(self, records: Iterable[SeqRecord]) -> Generator[ZaniResult, None, None]:
        records_list = list(records)
        futures = {
            self._executor.submit(self._all_vs_all_worker, self._engine, ref, records_list): ref.name
            for ref in records_list
        }
        for future in as_completed(futures):
            yield from future.result()

    @staticmethod
    def _search_worker(engine: ZaniEngine, sketch: Sketch, targets_chunk: list[SeqRecord]) -> list[ZaniResult]:
        return [engine.compare(sketch, target) for target in targets_chunk]

    def search(self, query: SeqRecord, targets: Iterable[SeqRecord]) -> Generator[ZaniResult, None, None]:
        sketch = self._engine.sketch(query)
        futures = set()
        for targets_chunk in batched(targets, self.batch_size):
            futures.add(self.executor.submit(self._search_worker, self._engine, sketch, list(targets_chunk)))
        for future in as_completed(futures):
            yield from future.result()

    @staticmethod
    def _query_worker(engine: ZaniEngine, database: Iterable[Sketch], queries_chunk: list[SeqRecord]) -> list[ZaniResult]:
        return [engine.compare(sketch, query) for sketch in database for query in queries_chunk]

    def query(self, queries: Iterable[SeqRecord], database: Iterable[Sketch]) -> Generator[ZaniResult, None, None]:
        db_list = list(database)
        futures = {self.executor.submit(self._query_worker, self._engine, db_list, list(queries_chunk))
                   for queries_chunk in batched(queries, self.batch_size)}
        for future in as_completed(futures):
            yield from future.result()

    def sketch(self, records: Iterable[SeqRecord]) -> Generator[Sketch, None, None]:
        futures = [self.executor.submit(self._engine.sketch, record) for record in records]
        for future in as_completed(futures):
            yield future.result()
