import argparse
import logging
import os
import shutil
import sys
import warnings
from typing import Optional

from . import __version__
from .logging import logger
from .utils import patch_mp_connection_bpo_17560


def setup_dedup_args(parser: argparse.ArgumentParser, parent: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Helper function to set up a subparser for the `dedup` command.

    Args:
        parser: Argparse parser to add the `dedup` command to
        parent: Argparse parser parent of the newly added subcommand.
            Used to inherit shared commands/flags

    Returns:
        The newly added parser
    """
    parser_dedup = parser.add_parser(
        'dedup',
        description='Deduplicate BAM files',
        help='Deduplicate BAM files',
        parents=[parent],
    )
    parser_dedup._actions[0].help = parser_dedup._actions[0].help.capitalize()

    parser_dedup.add_argument(
        '--input',
        metavar='BAM',
        help='Path to the BAM file',
        type=str,
        required=True
    )
    parser_dedup.add_argument(
        '--output',
        metavar='OUT',
        help='Path to the DEDUP file will be output',
        type=str,
        required=True
    )
    parser_dedup.add_argument(
        '--nasc',
        help=argparse.SUPPRESS,
        action='store_true',
    )

    return parser_dedup


def setup_estimate_args(parser: argparse.ArgumentParser, parent: argparse.ArgumentParser) -> argparse.ArgumentParser:
    """Helper function to set up a subparser for the `estimate` command.

    Args:
        parser: Argparse parser to add the `estimate` command to
        parent: Argparse parser parent of the newly added subcommand.
            Used to inherit shared commands/flags

    Returns:
        The newly added parser
    """
    parser_estimate = parser.add_parser(
        'estimate',
        description='Estimate fraction of labeled RNA',
        help='Estimate fraction of labeled RNA',
        parents=[parent],
    )
    parser_estimate._actions[0].help = parser_estimate._actions[0].help.capitalize()

    parser_estimate.add_argument(
        '--reads',
        help=(
            'Read groups to perform estimation on. '
            'This option can be used multiple times to estimate multiple groups. '
            '(default: all possible reads groups)'
        ),
        action='append',
        choices=['total', 'transcriptome', 'spliced', 'unspliced'],
        default=None
    )

    parser_estimate.add_argument(
        '-o',
        metavar='OUT',
        help='Path to output directory (default: current directory)',
        type=str,
        default='.',
    )
    parser_estimate.add_argument(
        '--gtf',
        metavar='GTF',
        help='Path to GTF file',
        type=str,
        required=True,
    )
    parser_estimate.add_argument(
        '--barcodes',
        metavar='TXT',
        help=(
            'Textfile containing filtered cell barcodes. Only these barcodes will be processed. '
            'This option may be used multiple times when multiple input directories are provided.'
        ),
        action='append',
        default=None,
    )
    parser_estimate.add_argument(
        '--groups',
        metavar='CSV',
        help=(
            'CSV containing cell (barcode) groups, where the first column is the barcode '
            'and the second is the group name the cell belongs to. Estimation will be performed '
            'by aggregating UMIs per group. This option may be used multiple times when multiple '
            'input directories are provided.'
        ),
        action='append',
        default=None,
    )
    parser_estimate.add_argument(
        '--method',
        help=(
            'Correction method to use. '
            'May be `pi_g` to estimate the fraction of labeled RNA for every cell-gene combination, '
            'or `alpha` to use alpha correction as used in the scNT-seq paper. `alpha` is recommended for '
            'UMI-based assays. This option has no effect when used with `--control`. (default: alpha)'
        ),
        choices=['pi_g', 'alpha'],
        default='alpha'
    )
    parser_estimate.add_argument(
        '--ignore-groups-for-est',
        help=(
            'Ignore cell groupings when calculating final estimations for the fraction of labeled RNA. '
            'When `--method pi_g`, groups are ignored when estimating fraction of labeled RNA. '
            'When `--method alpha`, groups are ignored when estimating detection rate. '
            'This option only has an effect when `--groups` is also specified.'
        ),
        action='store_true',
    )
    parser_estimate.add_argument(
        '--genes',
        metavar='TXT',
        help=('Textfile containing list of genes to use. All other genes will be '
              'treated as if they do not exist.'),
        type=str,
    )
    parser_estimate.add_argument(
        '--cell-threshold',
        metavar='COUNT',
        help='A cell must have at least this many reads for correction. (default: 1000)',
        type=int,
        default=1000,
    )
    parser_estimate.add_argument(
        '--cell-gene-threshold',
        metavar='COUNT',
        help=(
            'A cell-gene pair must have at least this many reads for correction. '
            'Only for `--method pi_g`. (default: 16)'
        ),
        type=int,
        default=16
    )
    parser_estimate.add_argument(
        '--gene-names',
        help=('Group counts by gene names instead of gene IDs when generating H5AD file'),
        action='store_true'
    )
    parser_estimate.add_argument(
        '--downsample',
        metavar='NUM',
        help=(
            'Downsample the number of reads (UMIs). If a decimal between 0 and 1 is given, '
            'then the number is interpreted as the proportion of remaining reads. If an integer '
            'is given, the number is interpreted as the absolute number of remaining reads.'
        ),
        type=float,
        default=None
    )
    parser_estimate.add_argument(
        '--downsample-mode',
        metavar='MODE',
        help=(
            'Downsampling mode. Can be one of: `uniform`, `cell`, `group`. If `uniform`, all reads '
            '(UMIs) are downsampled uniformly at random. If `cell`, only cells that have more '
            'reads than the argument to `--downsample` are downsampled to exactly that number. '
            'If `group`, identical to `cell` but per group specified by `--groups`.'
        ),
        type=str,
        choices=['uniform', 'cell', 'group'],
        default='uniform'
    )
    parser_estimate.add_argument(
        '--nasc',
        help=argparse.SUPPRESS,
        action='store_true',
    )
    parser_estimate.add_argument('--seed', help=argparse.SUPPRESS, type=int, default=None)
    parser_estimate.add_argument(
        '--control',
        help=('Indicate this is a control sample, only the background mutation rate '
              'will be estimated.'),
        action='store_true',
    )
    parser_estimate.add_argument(
        '--p-e',
        help='Textfile containing a single number, indicating the estimated background mutation rate',
        type=str,
        default=None,
    )
    parser_estimate.add_argument(
        'count_dirs',
        help=(
            'Path to directory that contains `dynast count` output. When multiple are provided, '
            'the barcodes in each of the count directories are suffixed with `-i` where i is '
            'a 0-indexed integer.'
        ),
        type=str,
        nargs='+',
    )

    return parser_estimate


def parse_dedup(parser: argparse.ArgumentParser, args: argparse.Namespace, temp_dir: Optional[str] = None):
    """Parser for the `dedup` command.

    Args:
        parser: The parser
        args: Command-line arguments dictionary, as parsed by argparse
        temp_dir: Temporary directory
    """
    if not os.path.exists(args.input):
        parser.error(
            f'Input file {args.input} does not exist. '
            'Please check the path and try again.'
        )

    from .dedup import dedup
    dedup(
        args.input,
        args.output,
        nasc=args.nasc,
        n_threads=args.t,
        temp_dir=temp_dir,
    )


def parse_estimate(parser: argparse.ArgumentParser, args: argparse.Namespace, temp_dir: Optional[str] = None):
    """Parser for the `estimate` command.

    Args:
        parser: The parser
        args: Command-line arguments dictionary, as parsed by argparse
        temp_dir: Temporary directory
    """
    # Check control constraints
    if args.control and args.p_e:
        parser.error('`--control` and `--p-e` can not be used together')

    # Check p_e is in correct format (only a single number)
    control_p_e = None
    if args.p_e:
        with open(args.p_e, 'r') as f:
            try:
                control_p_e = float(f.read().strip())
            except ValueError:
                parser.error('`--p-e` must be a textfile containing a single decimal number')

    # group CSV must be proivded when there are multiple input directories
    if len(args.count_dirs) > 1 and not args.groups:
        parser.error(
            '`--group` CSVs must be provided when using multiple input directories. '
            'Otherwise, simply run `dynast estimate` for each directory separately.'
        )

    # If group CSV(s) are provided, the number must match input directories
    if args.groups and len(args.groups) != len(args.count_dirs):
        parser.error('Number of `--group` CSVs must match number of input directories')
    if args.barcodes and len(args.barcodes) != len(args.count_dirs):
        parser.error('Number of `--barcodes` TXTs must match number of input directories')

    # Multiple count dirs can't be used with nasc
    if len(args.count_dirs) > 1 and args.nasc:
        parser.error('`--nasc` does not support multiple input directories')

    # Read genes
    genes = None
    if args.genes:
        if args.by_name:
            logger.warning(
                '`--genes` were provided with `--gene-names`. '
                'Make sure your gene list contains gene names instead of IDs. '
                'IDs should be used for any genes that do not have a name.'
            )
        with open(args.genes, 'r') as f:
            genes = [line.strip() for line in f if not line.isspace()]
        logger.warning(f'Ignoring genes not in the {len(genes)} genes provided by `--genes`')

    barcodes = []
    if args.barcodes:
        for path in args.barcodes:
            with open(path, 'r') as f:
                bcs = set(line.strip() for line in f if not line.isspace())
            barcodes.append(bcs)
        logger.warning(
            f'Ignoring cell barcodes not in the {sum(len(bcs) for bcs in barcodes)} barcodes provided by `--barcodes`'
        )
    else:
        logger.warning('`--barcodes` not provided. All cell barcodes will be processed.')

    # Parse cell groups csv(s)
    groups = []
    if args.groups:
        if args.barcodes:
            logger.warning(
                '`--groups` was provided with `--barcodes`. Only barcodes present in both inputs will be considered.'
            )

        for path in args.groups:
            groups_part = {}
            with open(path, 'r') as f:
                for line in f:
                    if line.isspace():
                        continue
                    barcode, group = line.strip().split(',')

                    if barcode in groups_part:
                        parser.error(f'Found duplicate barcode {barcode} in {path}')

                    groups_part[barcode] = group
            groups.append(groups_part)

    if args.ignore_groups_for_est and not (args.groups):
        parser.error('`--ignore-groups-for-pi` can not be used without `--groups`')

    if not args.reads:
        args.reads = 'complete'

    # Check that --downsample is an integer if --downsample-mode is cell
    if args.downsample and args.downsample_mode in ('cell', 'group'):
        if int(args.downsample) != args.downsample:
            parser.error('`--downsample` must be an integer when using `--downsample-mode cell/group`')
    if args.downsample_mode == 'group' and not args.groups:
        parser.error('`--groups` must be provided when using `--downsample-mode group`')

    from .estimate import estimate
    estimate(
        args.count_dirs,
        args.o,
        args.reads,
        barcodes=args.barcodes,
        groups=groups,
        ignore_groups_for_est=args.ignore_groups_for_est,
        genes=genes,
        downsample=args.downsample,
        downsample_mode=args.downsample_mode,
        cell_threshold=args.cell_threshold,
        cell_gene_threshold=args.cell_gene_threshold,
        control_p_e=control_p_e,
        control=args.control,
        method=args.method,
        n_threads=args.t,
        temp_dir=temp_dir,
        nasc=args.nasc,
        by_name=args.gene_names,
        seed=args.seed,
        genes_filename=args.gtf,
    )


COMMAND_TO_FUNCTION = {
    'dedup': parse_dedup,
    'estimate': parse_estimate,
}


@logger.namespaced('main')
def main():
    parser = argparse.ArgumentParser(
        description=f'{__version__} Complete splicing and labeling quantification from metabolic labeling scRNA-seq'
    )
    parser._actions[0].help = parser._actions[0].help.capitalize()
    subparsers = parser.add_subparsers(
        dest='command',
        metavar='<CMD>',
    )

    # Add common options to this parent parser
    parent = argparse.ArgumentParser(add_help=False)
    parent.add_argument('--tmp', metavar='TMP', help='Override default temporary directory', type=str, default='tmp')
    parent.add_argument('--keep-tmp', help='Do not delete the tmp directory', action='store_true')
    parent.add_argument('--verbose', help='Print debugging information', action='store_true')
    parent.add_argument('-t', metavar='THREADS', help='Number of threads to use (default: 8)', type=int, default=8)

    # Command parsers
    parser_dedup = setup_dedup_args(subparsers, parent)
    parser_estimate = setup_estimate_args(subparsers, parent)
    command_to_parser = {
        'dedup': parser_dedup,
        'estimate': parser_estimate,
    }
    if '--list' in sys.argv:
        print_technologies()

    # Show help when no arguments are given
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    if len(sys.argv) == 2:
        if sys.argv[1] in command_to_parser:
            command_to_parser[sys.argv[1]].print_help(sys.stderr)
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    logger.setLevel(logging.DEBUG if args.verbose else logging.INFO)

    logger.debug('Printing verbose output')
    logger.debug(f'Input args: {args}')
    logger.debug(f'Creating {args.tmp} directory')
    if os.path.exists(args.tmp):
        parser.error(
            f'Temporary directory {args.tmp} already exists. '
            'Is another process running? Please specify a different temporary '
            'directory with the `--tmp` option, or remove the one that already '
            'exists.'
        )
    os.makedirs(args.tmp)
    os.environ['NUMEXPR_MAX_THREADS'] = str(args.t)
    try:
        patch_mp_connection_bpo_17560()  # Monkeypatch for python 3.7
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            COMMAND_TO_FUNCTION[args.command](parser, args, temp_dir=args.tmp)
        logger.info('Done')
    except Exception:
        logger.exception('An exception occurred')
    finally:
        if not args.keep_tmp:
            logger.debug(f'Removing {args.tmp} directory')
            shutil.rmtree(args.tmp, ignore_errors=True)
