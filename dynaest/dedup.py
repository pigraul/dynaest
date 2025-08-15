import multiprocessing
import re
import os
import gc
from functools import partial
from typing import Dict, FrozenSet, List, Optional, Set, Tuple

import ngs_tools as ngs
import numpy as np
import pandas as pd
from typing_extensions import Literal
import pysam

from . import config, utils, constants
from .logging import logger

CONVERSIONS_PARSER = re.compile(
    r'''^
    (?P<read_id>[^,]*),
    (?P<index>[^,]*),
    (?P<contig>[^,]*),
    (?P<genome_i>[^,]*),
    (?P<conversion>[^,]*),
    (?P<quality>[^,]*)\n
    $''', re.VERBOSE
)

ALIGNMENTS_PARSER = re.compile(
    r'''^
    (?P<read_id>[^,]*),
    (?P<index>[^,]*),
    (?P<barcode>[^,]*),
    (?P<umi>[^,]*),
    (?P<GX>[^,]*),
    (?P<A>[^,]*),
    (?P<C>[^,]*),
    (?P<G>[^,]*),
    (?P<T>[^,]*),
    (?P<velocity>[^,]*),
    (?P<transcriptome>[^,]*),
    (?P<score>[^,]*)\n
    $''', re.VERBOSE
)

CONVERSION_IDX = {
    'AC': 0,
    'AG': 1,
    'AT': 2,
    'CA': 3,
    'CG': 4,
    'CT': 5,
    'GA': 6,
    'GC': 7,
    'GT': 8,
    'TA': 9,
    'TC': 10,
    'TG': 11,
}
BASE_IDX = {
    'A': 12,
    'C': 13,
    'G': 14,
    'T': 15,
}
CONVERSION_COMPLEMENT = {
    conversion: ngs.sequence.complement_sequence(conversion, reverse=False)
    for conversion in CONVERSION_IDX.keys()
}
CONVERSION_COLUMNS = sorted(CONVERSION_IDX.keys())
BASE_COLUMNS = sorted(BASE_IDX.keys())
COLUMNS = CONVERSION_COLUMNS + BASE_COLUMNS
CSV_COLUMNS = ['barcode', 'umi', 'GX'] + COLUMNS + ['velocity', 'transcriptome', 'score','st']


def read_alignments(
    alignments_path: str, 
    temp_dir: Optional[str] = None,
    *args, **kwargs
) -> pd.DataFrame:
    """Read counts bam as a pandas dataframe.

    Args:
        counts_path: Path to bam

    Returns:
        Counts dataframe
    """
    dtypes = {
        #'read_id': 'string',
        'barcode': 'category',
        'umi': 'category',
        'GX': 'category',
        'velocity': 'category',
        'transcriptome': bool,
        'score': np.uint16,
        'st': 'category',
        **{column: np.uint8 for column in COLUMNS}
    }

    save = pysam.set_verbosity(0)
    bamfile = pysam.AlignmentFile(alignments_path, "rb", require_index=False)
    pysam.set_verbosity(save)
    batch_size = 1000000  
    first_batch = True
    temp_file = utils.mkstemp(dir=temp_dir)
    temp_records = []

    for idx, read in enumerate(bamfile.fetch(until_eof=True)):
        try:
            snpmatch = re.findall(r'([a-zA-Z]{2})(\d+)', read.get_tag("SC"))
            totmatch = re.findall(r'([a-zA-Z]{1})(\d+)', read.get_tag("TC"))
            if snpmatch and totmatch:
                tmp_sc_dict = {k.upper(): int(v) for k, v in snpmatch}
                sc_dict = {k: v for k, v in tmp_sc_dict.items() if 'N' not in k}
                tc_dict = {k.upper(): int(v) for k, v in totmatch}
                rid_dict = {
                    'barcode': read.get_tag('CB'),
                    'umi': read.get_tag('UB'),
                    'GX': read.get_tag('GX'),
                    'velocity': 'unassigned',
                    'transcriptome': True,
                    'score': read.get_tag('AS'),
                    'st': read.get_tag('ST'),
                }
                merged_dict = {**sc_dict, **tc_dict, **rid_dict}
                temp_records.append(merged_dict)
                if (idx + 1) % batch_size == 0:
                    df = pd.DataFrame(temp_records, columns=CSV_COLUMNS)
                    df.to_csv(temp_file, mode='a', header=first_batch, index=False)
                    first_batch = False
                    temp_records = []
                    gc.collect() 
        except (ValueError, KeyError, TypeError):
            continue

    bamfile.close()
    if temp_records:
        df = pd.DataFrame(temp_records, columns=CSV_COLUMNS)
        df.to_csv(temp_file, mode='a', header=first_batch, index=False)

    return pd.read_csv(temp_file, dtype=dtypes, na_filter=False)


def read_counts(counts_path: str, *args, **kwargs) -> pd.DataFrame:
    """Read counts CSV as a pandas dataframe.

    Any additional arguments and keyword arguments are passed to `pandas.read_csv`.

    Args:
        counts_path: Path to CSV

    Returns:
        Counts dataframe
    """
    dtypes = {
        #'read_id': 'string',
        'barcode': 'category',
        'umi': 'category',
        'GX': 'category',
        'velocity': 'category',
        'transcriptome': bool,
        'score': np.uint16,
        'st': 'category',
        **{column: np.uint8
           for column in COLUMNS}
    }
    return pd.read_csv(counts_path, dtype=dtypes, na_filter=False, *args, **kwargs)


def complement_counts(df_counts: pd.DataFrame) -> pd.DataFrame:
    """Complement the counts in the counts dataframe according to gene strand.

    Args:
        df_counts: Counts dataframe

    Returns:
        counts dataframe with counts complemented for reads mapping to genes on the reverse strand
    """
    # Extract columns that do not need to be complemented
    other_columns = []
    for col in df_counts.columns:
        if col in COLUMNS:
            continue
        other_columns.append(col)
    forward_strand = df_counts['st'] == '+'

    columns = other_columns + COLUMNS
    df_forward = df_counts[forward_strand][columns]
    df_reverse = df_counts[~forward_strand][columns]

    df_reverse.columns = other_columns + CONVERSION_COLUMNS[::-1] + BASE_COLUMNS[::-1]
    df_reverse = df_reverse[columns]

    return pd.concat((df_forward, df_reverse), verify_integrity=True)


def subset_counts(
    df_counts: pd.DataFrame,
    key: Literal['total', 'transcriptome', 'spliced', 'unspliced'],
) -> pd.DataFrame:
    """Subset the given counts DataFrame to only contain reads of the desired key.

    Args:
        df_count: Counts dataframe
        key: Read types to subset

    Returns:s
        Subset dataframe
    """
    if key == 'transcriptome':
        df_counts = df_counts[df_counts['transcriptome']]
    if key in ('spliced', 'unspliced'):
        df_counts = df_counts[df_counts['velocity'] == key]
    return df_counts


def deduplicate_counts(
    df_counts: pd.DataFrame,
    conversions: Optional[FrozenSet[str]] = None,
    use_conversions: bool = True
) -> pd.DataFrame:
    """Deduplicate counts based on barcode, UMI, and gene.

    The order of priority is the following.
    1. If `use_conversions=True`, reads that have at least one such conversion
    2. Reads that align to the transcriptome (exon only)
    3. Reads that have highest alignment score
    4. If `conversions` is provided, reads that have a larger sum of such conversions
       If `conversions` is not provided, reads that have larger sum of all conversions

    Args:
        df_counts: Counts dataframe
        conversions: Conversions to prioritize, defaults to `None`
        use_conversions: Prioritize reads that have conversions first

    Returns:
        Deduplicated counts dataframe
    """
    convs = list(conversions) if conversions is not None else CONVERSION_COLUMNS
    df_counts['conversion_sum'] = df_counts[convs].sum(axis=1)
    # Deduplication priority.
    sort_order = ['transcriptome', 'score', 'conversion_sum']
    to_remove = ['conversion_sum']

    if use_conversions and conversions is not None:
        df_counts['has_conversions'] = df_counts[list(conversions)].sum(axis=1) > 0
        sort_order.insert(0, 'has_conversions')
        to_remove.append('has_conversions')

    # Sort by has desired conversion(s) last, transcriptome last,
    # best alignment last, most conversions last
    df_sorted = df_counts.sort_values(sort_order).drop(columns=to_remove)

    # Restore input dataframe.
    df_counts.drop(columns=to_remove, inplace=True)

    return df_sorted.drop_duplicates(subset=['barcode', 'umi', 'GX'], keep='last').reset_index(drop=True)


def deduplicate_counts_part(
    counter: multiprocessing.Value,
    lock: multiprocessing.Lock,
    split_path: str,
    out_path: str,
    conversions: Optional[FrozenSet[str]],
    use_conversions: bool = True
):
    """Helper function to parallelize :func:`deduplicate_multimappers`.
    """
    deduplicate_counts(
        read_counts(split_path), conversions=conversions, use_conversions=use_conversions
    )[CSV_COLUMNS].to_csv(
        out_path, header=False, index=False
    )
    lock.acquire()
    counter.value += 1
    lock.release()
    return out_path


def split_counts_by_velocity(df_counts: pd.DataFrame) -> Dict[str, pd.DataFrame]:
    """Split the given counts dataframe by the `velocity` column.

    Args:
        df_counts: Counts dataframe

    Returns:
        Dictionary containing `velocity` column values as keys and the subset dataframe as values
    """
    dfs = {}
    for velocity, df_part in df_counts.groupby('velocity', sort=False, observed=True):
        dfs[velocity] = df_part.reset_index(drop=True)
    logger.debug(f'Found the following velocity assignments: {", ".join(dfs.keys())}')
    return dfs


def count_no_conversions(
    alignments_path: str,
    counter: multiprocessing.Value,
    lock: multiprocessing.Lock,
    index: List[Tuple[int, int, int]],
    barcodes: Optional[List[str]] = None,
    temp_dir: Optional[str] = None,
    update_every: int = 10000,
) -> str:
    """Count reads that have no conversion.

    Args:
        alignments_path: Alignments CSV path
        counter: Counter that keeps track of how many reads have been processed
        lock: Semaphore for the `counter` so that multiple processes do not
            modify it at the same time
        index: Index for conversions CSV
        barcodes: List of barcodes to be considered. All barcodes are considered if not provided
        temp_dir: Path to temporary directory
        update_every: Update the counter every this many reads

    Returns:
        Path to temporary counts CSV
    """
    count_path = utils.mkstemp(dir=temp_dir)
    positions = set(tup[2] for tup in index)
    n = 0
    with open(alignments_path, 'r') as f, open(count_path, 'w') as out:
        f.readline()  # header
        while True:
            pos = f.tell()
            line = f.readline()
            if not line:
                break
            n += 1
            if n == update_every:
                lock.acquire()
                counter.value += update_every
                lock.release()
                n = 0
            if pos in positions:
                continue

            groups = ALIGNMENTS_PARSER.match(line).groupdict()
            if barcodes and groups['barcode'] not in barcodes:
                continue
            out.write(
                f'{groups["read_id"]},{groups["barcode"]},{groups["umi"]},{groups["GX"]},'
                f'{",".join(groups.get(key, "0") for key in COLUMNS)},'
                f'{groups["velocity"]},{groups["transcriptome"]},{groups["score"]}\n'
            )
    lock.acquire()
    counter.value += n
    lock.release()

    return count_path


def calculate_mutation_rates(df_counts: pd.DataFrame, rates_path: str, group_by: Optional[List[str]] = None) -> str:
    """Calculate mutation rate for each pair of bases.

    Args:
        df_counts: Counts dataframe, with complemented reverse strand bases
        rates_path: Path to write rates CSV
        group_by: Column(s) to group calculations by, defaults to `None`, which
            combines all rows

    Returns:
        Path to rates CSV
    """
    logger.debug(f'Mutation rates will be grouped by {group_by}')

    if group_by is not None:
        df_sum = df_counts.groupby(group_by, sort=False, observed=True).sum(numeric_only=True).astype(np.uint32)
    else:
        df_sum = pd.DataFrame(df_counts.sum(numeric_only=True).astype(np.uint32)).T

    # Compute mutation rate of each base
    dfs = []
    for base in sorted(BASE_IDX.keys()):
        logger.debug(f'Calculating {base} mutation rate')
        dfs.append(df_sum.filter(regex=f'^{base}[A,C,G,T]$').divide(df_sum[base], axis=0))
    df_rates = pd.concat(dfs, axis=1)
    if group_by is not None:
        df_rates = df_rates.reset_index()
    df_rates.to_csv(rates_path, index=False)
    return rates_path


def count_conversions_part(
    conversions_path: str,
    alignments_path: str,
    counter: multiprocessing.Value,
    lock: multiprocessing.Lock,
    index: List[Tuple[int, int, int]],
    barcodes: Optional[List[str]] = None,
    snps: Optional[Dict[str, Dict[str, Set[int]]]] = None,
    quality: int = 27,
    temp_dir: Optional[str] = None,
    update_every: int = 10000,
) -> str:
    """Count the number of conversions of each read per barcode and gene, along with
    the total nucleotide content of the region each read mapped to, also per barcode
    and gene. This function is used exclusively for multiprocessing.

    Args:
        conversions_path: Path to conversions CSV
        alignments_path: Path to alignments information about reads
        counter: Counter that keeps track of how many reads have been processed
        lock: Semaphore for the `counter` so that multiple processes do not
            modify it at the same time
        index: Index for conversions CSV
        barcodes: List of barcodes to be considered. All barcodes are considered if not provided
        snps: Dictionary of contig as keys and list of genomic positions as
            values that indicate SNP locations
        quality: Only count conversions with PHRED quality greater than this value
        temp_dir: Path to temporary directory, defaults to `None`
        update_every: Update the counter every this many reads

    Returns:
        Path to temporary counts CSV
    """

    def is_snp(gx, conversion, contig, genome_i):
        if not snps:
            return False
        return genome_i in snps.get(conversion, {}).get(contig, set())

    count_path = utils.mkstemp(dir=temp_dir)

    n = 0
    with open(conversions_path, 'r') as f, open(alignments_path, 'r') as f_alignments, open(count_path, 'w') as out:
        for pos, n_lines, pos2 in index:
            f.seek(pos)
            f_alignments.seek(pos2)
            n += 1
            if n == update_every:
                lock.acquire()
                counter.value += update_every
                lock.release()
                n = 0

            alignment = ALIGNMENTS_PARSER.match(f_alignments.readline()).groupdict()
            if barcodes and alignment['barcode'] not in barcodes:
                continue
            counts = [0] * (len(CONVERSION_IDX) + len(BASE_IDX))
            for base, i in BASE_IDX.items():
                counts[i] = alignment[base]
            gx = alignment['GX']
            for _ in range(n_lines):
                groups = CONVERSIONS_PARSER.match(f.readline()).groupdict()
                conversion = groups["conversion"]
                if int(groups['quality']) > quality and not is_snp(gx, conversion, groups['contig'], int(
                        groups['genome_i'])):
                    counts[CONVERSION_IDX[conversion]] += 1

            out.write(
                f'{groups["read_id"]},{alignment["barcode"]},{alignment["umi"]},'
                f'{alignment["GX"]},{",".join(str(c) for c in counts)},'
                f'{alignment["velocity"]},{alignment["transcriptome"]},{alignment["score"]}\n'
            )

        lock.acquire()
        counter.value += n
        lock.release()

    return count_path


def dedup_conversions(
    alignments_path: str,
    counts_dir: str,
    conversions: Optional[FrozenSet[str]] = frozenset({"TC"}),
    dedup_use_conversions: bool = True,
    n_threads: int = 8,
    temp_dir: Optional[str] = None,
    nasc: bool = False,
) -> str:
    """Count the number of conversions of each read per barcode and gene, along with
    the total nucleotide content of the region each read mapped to, also per barcode.
    When a duplicate UMI for a barcode is observed, the read with the greatest
    number of conversions is selected.

    Args:
        alignments_path: Path to alignments information about reads
        counts_dir: Path to write counts CSV
        barcodes: List of barcodes to be considered. All barcodes are considered if not provided
        snps: Dictionary of contig as keys and list of genomic positions as
            values that indicate SNP locations
        conversions: Conversions to prioritize when deduplicating only applicable
            for UMI technologies
        dedup_use_conversions: Prioritize reads that have at least one conversion
            when deduplicating
        quality: Only count conversions with PHRED quality greater than this value
        n_threads: Number of threads
        temp_dir: Path to temporary directory

    Returns:
        Path to counts CSV
    """


    # Filter counts dataframe
    logger.info(f'Loading combined counts from {alignments_path}')
    aln_counts = read_alignments(alignments_path, temp_dir=temp_dir)
    df_counts = complement_counts(aln_counts)
    #umi = all(df_counts['umi'] != 'NA')
    barcode_groupby = df_counts.groupby('barcode', sort=False, observed=True)
    barcode_counts = dict(barcode_groupby.size())
    split_paths = []
    current_split_path = None
    current_split_f = None
    current_split_size = 0
    # Split barcodes into approximately `config.COUNTS_SPLIT_THRESHOLD` bins.
    # Note that a single barcode may have more than this many reads.
    try:
        for barcode, df_counts_barcode in barcode_groupby:
            # Make its own split
            if barcode_counts[barcode] > config.COUNTS_SPLIT_THRESHOLD:
                split_path = utils.mkstemp(dir=temp_dir)
                logger.debug(f'Splitting counts for barcode {barcode} to {split_path}')
                df_counts_barcode.to_csv(split_path, index=False)
                split_paths.append(split_path)
            elif current_split_path is None:
                current_split_path = utils.mkstemp(dir=temp_dir)
                logger.debug(f'Splitting counts for residual barcodes to {current_split_path}')
                current_split_f = open(current_split_path, 'w')
                # Write header
                df_counts_barcode.to_csv(current_split_f, index=False)
                split_paths.append(current_split_path)
            else:
                # Don't write header
                df_counts_barcode.to_csv(current_split_f, index=False, header=False)

                # If we exceeded read threshold, close file & reset.
                current_split_size += df_counts_barcode.shape[0]
                if current_split_size > config.COUNTS_SPLIT_THRESHOLD:
                    current_split_f.close()
                    current_split_path = None
                    current_split_size = 0

    finally:
        if current_split_f is not None:
            current_split_f.close()

    del df_counts

    logger.debug(f'Spawning {n_threads} processes')
    pool, counter, lock = utils.make_pool_with_counter(n_threads)
    paths = [(split_path, utils.mkstemp(dir=temp_dir)) for split_path in split_paths]
    async_result = pool.starmap_async(
        partial(
            deduplicate_counts_part,
            counter,
            lock,
            conversions=conversions,
            use_conversions=dedup_use_conversions,
        ) , paths
    )
    pool.close()

    # Display progres bar
    utils.display_progress_with_counter(counter, len(split_paths), async_result, desc='filtering')
    pool.join()

    # Need to complement counts again (to revert to original)
    counts_path = os.path.join(counts_dir, f'{constants.COUNTS_PREFIX}_{"_".join(list(conversions))}.csv')
    logger.info(f'Counting conversions to {counts_path}')
    with open(counts_path, 'w') as out:
        out.write(f'{",".join(CSV_COLUMNS)}\n')
        for counts_part_path in async_result.get():
            df_part = complement_counts(pd.read_csv(counts_part_path, names=CSV_COLUMNS))
            df_part[CSV_COLUMNS].to_csv(out, header=False, index=False)

    #return counts_path
    # Calculate mutation rates for each group
    if nasc:
        logger.debug(f'Calculate mutation rates')
        df_counts_uncomplemented = read_counts(counts_path)
        df_counts_complemented = complement_counts(df_counts_uncomplemented)
        df_counts = df_counts_uncomplemented if nasc else df_counts_complemented
        rates_all_path = os.path.join(counts_dir, f'{constants.RATES_PREFIX}.csv')
        logger.info(f'Calculating mutation rates for all reads to {rates_all_path}.')
        rates_all_path = calculate_mutation_rates(df_counts, rates_all_path, group_by=['barcode'])



@logger.namespaced('dedup')
def dedup(
    input_path: str,
    output_path: str,
    nasc: bool = False,
    n_threads: int = 8,
    temp_dir: Optional[str] = None
):
    """Main interface for the `dedup` command.

    Args:
        input_path: Path to input BAM/CRAM
        output_path: Path to output BAM
        n_threads: Number of threads, defaults to `8`
        memory: *Suggested* memory to use (this is not guaranteed), in bytes
        temp_dir: Temporary directory
    """
    os.makedirs(output_path, exist_ok=True)
    logger.info(f'Dedup and Calculating mutation rates')
    dedup_conversions(
        input_path, output_path, nasc=nasc, n_threads=n_threads, temp_dir=temp_dir
    )
