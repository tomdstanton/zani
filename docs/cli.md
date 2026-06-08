# CLI Reference

```text
Zstandard-based Average Nucleotide Identity (ANI)

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

## `zani ava`

💥 Compute an all-vs-all pairwise distance matrix

```text
💥 Compute an all-vs-all pairwise distance matrix

Usage: zani ava [OPTIONS] <GENOMES>...

Arguments:
  <GENOMES>...  Input FASTA/FASTQ files to compare

Options:
  -o, --output <OUTPUT>          Output TSV file path (Defaults to stdout)
  -l, --level <LEVEL>            Zstandard compression level (1-19). Higher is more accurate but slower [default: 3]
  -s, --strategy <STRATEGY>      Compression match strategy [default: lazy2] [possible values: auto, fast, dfast, greedy, lazy, lazy2, btlazy2, btopt, btultra, btultra2]
  -c, --concat                   Treat multi-record FASTAs as a single concatenated genome
  -b, --batch-size <BATCH_SIZE>  How many rows to buffer before writing to disk [default: 10000]
  -t, --threads <THREADS>        Number of execution threads (0 = auto) [default: 0]
  -d, --db <DB>                  Save the compiled database to disk (.zani file)
  -h, --help                     Print help
  -V, --version                  Print version
```

## `zani search`

🔍 Compute distances of queries against a compiled target database

```text
🔍 Compute distances of queries against a compiled target database

Usage: zani search [OPTIONS] --db <DB> <GENOMES>...

Arguments:
  <GENOMES>...  Input FASTA/FASTQ files to search

Options:
  -d, --db <DB>                  Path to the pre-compiled target database (.zani)
  -o, --output <OUTPUT>          Output TSV file path (Defaults to stdout)
  -l, --level <LEVEL>            Zstandard compression level (1-19). Must match database creation! [default: 3]
  -s, --strategy <STRATEGY>      Compression match strategy [default: lazy2] [possible values: auto, fast, dfast, greedy, lazy, lazy2, btlazy2, btopt, btultra, btultra2]
  -c, --concat                   Treat multi-record FASTAs as a single concatenated genome
  -b, --batch-size <BATCH_SIZE>  How many rows to buffer before writing to disk [default: 10000]
  -t, --threads <THREADS>        Number of execution threads (0 = auto) [default: 0]
  -h, --help                     Print help
  -V, --version                  Print version
```

## `zani build`

💾 Pre-compile a database of genomes to disk

```text
💾 Pre-compile a database of genomes to disk

Usage: zani build [OPTIONS] --db <DB> <GENOMES>...

Arguments:
  <GENOMES>...  Input FASTA files

Options:
  -d, --db <DB>              Output database file (.zani)
  -l, --level <LEVEL>        Zstandard compression level (1-19) [default: 3]
  -s, --strategy <STRATEGY>  Compression match strategy [default: lazy2] [possible values: auto, fast, dfast, greedy, lazy, lazy2, btlazy2, btopt, btultra, btultra2]
  -c, --concat               Treat multi-record FASTAs as a single concatenated genome
  -h, --help                 Print help
  -V, --version              Print version
```

