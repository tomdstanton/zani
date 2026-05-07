import argparse
from pathlib import Path

from zani._version import __version__


# Main CLI entry-point -------------------------------------------------------------------------------------------------
def main():

    # Define args ------------------------------------------------------------------------------------------------------
    parser = argparse.ArgumentParser(
        description='Graph-aware contextual annotation of targeted genomic features',
        usage="%(prog)s <reference> <genomes ...> [options]", add_help=False, prog=__package__,
        formatter_class=argparse.MetavarTypeHelpFormatter,
    )

    inputs = parser.add_argument_group('📁', 'Input arguments')
    inputs.add_argument('reference', type=Path, metavar='<reference>',
                        help='Path to a reference genome in fasta format; File may be compressed.')
    inputs.add_argument('genomes', type=str, nargs="+", metavar='<genomes ...>',
                        help='Paths to query genomes in fasta format; Files may be compressed.')

    outs = parser.add_argument_group('💾', 'Output arguments')
    outs.add_argument('-o', '--out', type=Path, default=None, metavar='',
                      help='Direct output to file, defaults to stdout')

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
    from zani.zani import ZaniEngine, ZaniWriter, Reference, SeqRecord

    reference = Reference.from_file(args.reference)
    
    # We use map to create a lazy iterator! The main thread acts as a dedicated 
    # I/O & Gzip worker (streaming sequentially off disk), perfectly feeding the 
    # thread pool which acts as dedicated Zstd compression workers.
    queries = map(SeqRecord.from_file, args.genomes)

    with ZaniWriter(args.out) as writer, ZaniEngine(reference, max_workers=args.max_workers) as engine:
        for result in engine.query(queries):
            writer.write(result)
