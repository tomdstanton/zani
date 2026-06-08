# `zani` 🧬🗜️🤪
*pronounced zany (/ˈzeɪni/)*

[![Release](https://img.shields.io/github/v/release/tomdstanton/zani?style=flat-square)](https://img.shields.io/github/v/release/tomdstanton/zani)
[![License](https://img.shields.io/github/license/tomdstanton/zani?style=flat-square)](https://img.shields.io/github/license/tomdstanton/zani)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19059429.svg?style=flat-square)](https://doi.org/10.5281/zenodo.19059429)
[![PyPI](https://img.shields.io/pypi/v/zani.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/zani)
[![Wheel](https://img.shields.io/pypi/wheel/zani.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/zani/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/zani.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/zani/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/zani.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/zani/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/tomdstanton/zani/)
[![Issues](https://img.shields.io/github/issues/tomdstanton/zani.svg?style=flat-square&maxAge=600)](https://github.com/tomdstanton/zani/issues)

**High-Performance Average Nucleotide Identity (ANI) estimator using Zstandard compression distance.**

## 📖 About
`zani` computes pairwise genomic distances using the Normalized Compression Distance (NCD) metric. 
Originally written in Python, `zani` has been completely rewritten in **Rust** 🦀 for bare-metal multi-threaded performance 🚀. It leverages the blazing-fast **Zstandard (zstd)** compression algorithm to estimate Average Nucleotide Identity (ANI) without the need for expensive sequence alignments or k-mer counting.

### 🧮 The Algorithm
At its core, `zani` treats reference genomes as compression dictionaries. For a given reference genome $x$ and a query genome $y$:

1. **Dictionary Training** 📚: A Zstd dictionary is trained on the reference genome $x$.
2. **Baseline Compression** 📏: We compute $C(x)$, the size of the reference genome compressed with its own dictionary.
3. **Conditional Compression** 🗜️: The query genome $y$ is compressed *using* the dictionary trained on $x$. This yields $C(y|x)$, representing the amount of novel information in $y$ not found in $x$.

### 🔢 The Math
`zani` calculates distance using the standard Normalized Compression Distance (NCD) formula:

$$ NCD(x,y) = \frac{C(x,y) - \min(C(x), C(y))}{\max(C(x), C(y))} $$

To achieve maximum execution speed ⚡, `zani` approximates the joint compression size $C(x,y)$ as:

$$ C(x,y) \approx C(x) + C(y|x) $$

Furthermore, to avoid the performance penalty of compressing the query genome twice to find its baseline $C(y)$, `zani` rapidly estimates $C(y)$ using the ratio of their uncompressed lengths ($|x|$ and $|y|$):

$$ C(y) \approx C(x) \times \frac{|y|}{|x|} $$

### 🧱 The Chimeric Firewall
When concatenating multi-record draft genomes (e.g., thousands of small contigs) into a single sequence, `zani` automatically inserts a "Chimeric Firewall" 🛡️ of exactly 10 'N's (`NNNNNNNNNN`) between each contig boundary. Because Zstandard requires strictly contiguous exact-byte matches, this firewall mathematically guarantees that no LZ77 match can artificially bridge across two independent contigs, perfectly neutralizing artifactual ANI inflation without bloating the NCD baseline scores!

This mathematical approach, combined with zero-copy memoryviews, L1 cache optimization, and thread-safe C-contexts natively in Rust 🦀, allows `zani` to stream thousands of genomes through concurrent worker threads, achieving massive I/O throughput and utilizing 100% of available CPU cores.


## Installation
`zani` is distributed as native pre-compiled binaries for Windows, macOS, and Linux.

```shell
pip install zani
```

## 💻 CLI Usage

`zani` now comes with a compiled Rust CLI for massive throughput:

```shell
❯ zani --help
🧬 zani: Zstandard-based Average Nucleotide Identity (ANI) 🗜️

Usage: zani <COMMAND>

Commands:
  ava     💥 Compute an all-vs-all pairwise distance matrix
  search  🔍 Compute distances of queries against a compiled target database
  build   💾 Pre-compile a database of genomes to disk
  help    Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help
  -V, --version  Print version
```

### 💥 Running an All-vs-All Matrix
```shell
zani ava genome1.fasta genome2.fasta --output zani_matrix.tsv --level 3
```

### 🔍 Searching against a Database
```shell
zani search query1.fasta query2.fasta --db my_database.zani --output results.tsv
```

### 💾 Compiling a Target Database
```shell
zani build genome1.fasta genome2.fasta --db my_database.zani
```


## 🐍 Python API Usage

The Python API has been completely rebuilt as a native PyO3 extension. It acts exactly like native Python objects, but drops the Global Interpreter Lock (GIL) to execute purely in Rust!

> [!NOTE]
> Extensive Python API documentation is available! You can generate the HTML docs locally by running `just doc-api` (outputs to `docs/api`).

```python
import zani

# 1. Initialize an empty native database 🗄️
db = zani.Database()

# 2. Add FASTA files to compile Zstandard dictionaries 🗜️
db.add_fasta("genome_1.fasta", level=3, concat=True)
db.add_fasta("genome_2.fasta", level=3, concat=True)

# 3. Save/Load the compiled database to disk 💾 (Optional)
db.save("my_genomes.zani")

# 4. Initialize the execution Engine 🏎️
engine = zani.Engine(compression_level=3, batch_size=10_000)

# 5. Run the multithreaded matrix! 🚀 (Releases the GIL)
engine.all_vs_all(db, output_filepath="results.tsv")

# Or search queries against a target database! 🔍
engine.search(db, queries=db, output_filepath="results.tsv")
```


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the
[GitHub issue tracker](https://github.com/tomdstanton/zani/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.