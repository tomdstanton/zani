import argparse
from zani._version import __version__


# Main CLI entry-point -------------------------------------------------------------------------------------------------
def main():

    # Define args ------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description='🧬🗜️🤪 Average Nucleotide Identity (ANI) estimator using Zstandard compression distance.',
        usage="%(prog)s <genomes ...> [options]", add_help=False, prog=__package__,
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )

    inputs = parser.add_argument_group('📁', 'Input arguments')
    inputs.add_argument('genomes', nargs="+", type=str, metavar='<genomes ...>',
                        help='Paths to genomes in fasta format; Files may be compressed.')

    opts = parser.add_argument_group('🛠️', 'Other options')
    opts.add_argument('-t', '--max-workers', type=int, default=None, metavar='',
                      help='Maximum number of threads to use for parallelization')
    opts.add_argument('-v', '--version', action='version', version=__version__,
                      help='Show version number and exit')
    opts.add_argument('-h', '--help', action='help',
                      help='Show this help message and exit')

    # Parse args -------------------------------------------------------------------------------------------------------
    args = parser.parse_args()

    # Run pipeline -----------------------------------------------------------------------------------------------------
    from zani._engine import ZaniPipeline, ZaniEngine, SeqRecord
    from sys import stdout

    genomes = map(SeqRecord.from_file, args.genomes)
    with ZaniPipeline(ZaniEngine(), args.max_workers) as pipeline:
        for result in pipeline.all_vs_all(genomes):
            result.write(stdout.buffer)
