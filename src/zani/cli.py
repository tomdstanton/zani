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
    inputs.add_argument('-a', '--allvsall', action="store_true",
                        help='Run all-vs-all comparison')

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
    from zani.zani import ZaniEngine
    from sys import stdout

    with ZaniEngine(max_workers=args.max_workers) as engine:
        for result in (engine.all_vs_all if args.allvsall else engine.query)(args.genomes):
            result.write(stdout.buffer)
