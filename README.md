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
[![Downloads](https://img.shields.io/pypi/dm/zani?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/zani)[![License](https://img.shields.io/github/license/tomdstanton/zani)](https://img.shields.io/github/license/tomdstanton/zani)

**Average Nucleotide Identity (ANI) estimator using Zstandard compression distance.**

## About
`zani` computes pairwise genomic distances using the Normalized Compression Distance (NCD) metric.
Inspired by the pioneering work of [LZ-ANI](https://github.com/refresh-bio/LZ-ANI), `zani` leverages the
blazing-fast **Zstandard (zstd)** compression algorithm to estimate Average Nucleotide Identity (ANI) without the need
for expensive sequence alignments or k-mer counting.

### The Algorithm
At its core, `zani` treats reference genomes as compression dictionaries. For a given reference genome $x$ and a query genome $y$:

1. **Dictionary Training**: A Zstd dictionary is trained on the reference genome $x$.
2. **Baseline Compression**: We compute $C(x)$, the size of the reference genome compressed with its own dictionary.
3. **Conditional Compression**: The query genome $y$ is compressed *using* the dictionary trained on $x$. This yields $C(y|x)$, representing the amount of novel information in $y$ not found in $x$.

### The Math
`zani` calculates distance using the standard Normalized Compression Distance (NCD) formula:

$$ NCD(x,y) = \frac{C(x,y) - \min(C(x), C(y))}{\max(C(x), C(y))} $$

To achieve maximum execution speed, `zani` approximates the joint compression size $C(x,y)$ as:

$$ C(x,y) \approx C(x) + C(y|x) $$

Furthermore, to avoid the performance penalty of compressing the query genome twice to find its baseline $C(y)$, `zani` rapidly estimates $C(y)$ using the ratio of their uncompressed lengths ($|x|$ and $|y|$):

$$ C(y) \approx C(x) \times \frac{|y|}{|x|} $$

This mathematical approach, combined with zero-copy memoryviews and thread-local C-contexts, allows `zani` to stream thousands of genomes through concurrent worker threads, achieving massive I/O throughput and utilizing 100% of available CPU cores.


## Installation
`zani` can be installed with `pip`:

```shell
pip install zani
```

## CLI Usage 💻

`zani` has a very basic CLI, use it like so:

```shell
❯ uv run zani -h
usage: zani <genomes ...> [options]

🧬🗜️🤪 Average Nucleotide Identity (ANI) estimator using Zstandard compression distance.

📁:
  Input arguments

  <genomes ...>       Paths to genomes in fasta format; Files may be compressed.
  -a, --allvsall      Run all-vs-all comparison

🛠️:
  Other options

  -t, --max-workers   Maximum number of threads to use for parallelization
  -v, --version       Show version number and exit
  -h, --help          Show this help message and exit
```

## API Usage 💻

```python
from pathlib import Path
from zani import ZaniEngine

genomes = Path('genomes').glob('*.fasta.gz')

with ZaniEngine() as engine:
    for result in engine.query(genomes):
        print(result)
```


## 💭 Feedback

### ⚠️ Issue Tracker

Found a bug ? Have an enhancement request ? Head over to the
[GitHub issue tracker](https://github.com/tomdstanton/zani/issues) if you need to report
or ask something. If you are filing in on a bug, please include as much
information as you can about the issue, and try to recreate the same bug
in a simple, easily reproducible situation.