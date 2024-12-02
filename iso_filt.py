"""
Filtration of RNA isoforms for selected genes from long-read RNA-seq data.

GitHub: https://github.com/stefmldk/iso-filt
"""

import sys
import os
import gzip
import re
from bisect import bisect

try:
    import pysam
except ModuleNotFoundError as error:
    print('\nIt appears that Pysam is not installed. Please install it before running this script. https://github.com/pysam-developers/pysam\n')
    exit()

__version__ = 1.1
__author__ = "Steffen Møller Bøttger"
__date__ = "12-06-2023"


def get_first_record_gene_annotation_data(gene_names, annotation_file, not_found_exit=True):
    """
    Returns a dictionary with data for the provided gene. If the gene is not found, the script exits unless
    not_found_exit is set to False. This may happen if the provided gene name does not correspond to how it is written
    in the annotation file.

    It returns the first record (first transcript) - not necessarily the one with most exons, but in most cases it will
    span the full gene.

    :param not_found_exit:
    :param annotation_file:
    :param gene_names:
    :return: dict
    """
    gene_data = {}
    for gene_name in gene_names:
        print('Collecting transcript data for gene: ' + gene_name)
        with (gzip.open(annotation_file, 'rt', encoding='utf-8') if annotation_file.endswith('gz') else open(annotation_file, 'r', encoding='utf-8')) as anno_file:
            exons = {}

            search_term = 'gene_id "{}"'.format(gene_name)
            for line in anno_file:
                if search_term in line:
                    line_blocks = line.split('\t')
                    if line_blocks[2] == 'exon':
                        exon_number = re.search(r'\d+', re.search(r'exon_number "\d+"', line).group()).group()
                        if not exon_number in exons:  # This is the control point that ensures the first block only
                            exons[exon_number] = {
                                'start': line_blocks[3],
                                'stop': line_blocks[4],
                            }
                        else:
                            break
                    elif line_blocks[2] == 'gene':
                        gene_data[gene_name] = {
                            'contig': line_blocks[0],
                            'start': line_blocks[3],
                            'stop': line_blocks[4],
                            'strand': line_blocks[6]
                        }
            if not exons:
                if not_found_exit:
                    print('The gene name, {}, was not found in the annotation file. Please make sure it is written correctly'.format(gene_name))
                    exit()
            else:
                gene_data[gene_name]['exons'] = exons
    return gene_data


def get_gene_exons_annotation_data(gene_names, annotation_file, not_found_exit=True):
    """
    Similarly to get_first_record_gene_annotation_data, it returns a dictionary with data for the provided gene. But
    instead of just getting the first record, it gets all exons from all isoforms. The exon numbering is relative to
    genomic exon position (sorted on start position) and may not reflect the numbering inside the isoforms.

    If a gene is not found, the script exits unless not_found_exit is set to False. This may happen if the provided gene
    name does not correspond to how it is written in the annotation file.

    :param gene_names:
    :param annotation_file:
    :param not_found_exit:
    :return: dict
    """
    gene_data = {}
    for gene_name in gene_names:
        print('Collecting exon data for gene: ' + gene_name)
        with (gzip.open(annotation_file, 'rt', encoding='utf-8') if annotation_file.endswith('gz') else open(annotation_file, 'r', encoding='utf-8')) as anno_file:
            exon_set = set()

            search_term = 'gene_id "{}"'.format(gene_name)
            search_term_found = False
            for line in anno_file:
                if search_term in line:
                    search_term_found = True
                    line_blocks = line.split('\t')
                    if line_blocks[2] == 'exon':
                        exon_set.add((int(line_blocks[3]), int(line_blocks[4])))

                    elif line_blocks[2] == 'gene':
                        gene_data[gene_name] = {
                            'contig': line_blocks[0],
                            'start': int(line_blocks[3]),
                            'stop': int(line_blocks[4]),
                            'strand': line_blocks[6]
                        }
                elif search_term_found:  # End of gene annotation block
                    break

            # Sort the exons relative to start coordinate
            sorted_exon_list = sorted(list(exon_set), key=lambda x:x[0])
            if not sorted_exon_list:
                if not_found_exit:
                    print('The gene name, {}, was not found in the annotation file. Please make sure it is written correctly'.format(gene_name))
                    exit()
            else:
                gene_data[gene_name]['exons'] = sorted_exon_list
    return gene_data


def write_bam_file(out_bam_file_name, list_of_reads, template_bam):
    """
    Writes given list of reads to the given out_bam_file_name. The resulting bam file is sorted and indexed.
    :param out_bam_file_name:
    :param list_of_reads:
    :param template_bam:
    :return:
    """

    # Write to temporary, unsorted bam
    template_file = pysam.AlignmentFile(template_bam, "rb")
    bam_file_dir = os.path.dirname(os.path.realpath(template_bam))
    temp_bam_path = bam_file_dir + "/temp_unsorted.bam"

    temp_bam_file = pysam.AlignmentFile(temp_bam_path, "wb", template=template_file)
    for read in list_of_reads:
        temp_bam_file.write(read)
    temp_bam_file.close()

    pysam.sort("-o", out_bam_file_name, temp_bam_path)
    pysam.index(out_bam_file_name)
    os.remove(temp_bam_path)


def get_full_length_rna_reads(bam_file, genes_of_interest, annotation_file):
    """
    Filtering function that filters out reads that are between start and stop exons of provided genes. Reads that extend
    beyond these boundaries by more than 2000 bases in either direction are also removed.
    The function may be useful if for instance standard poly-T based cDNA-synthesis and PCR amplification
    is used as this may lead to a substantial part of the reads being fragmented. Fragmentation is not only due to the
    internal priming but also due to premature termination. Premature termination may happen if the polymerase
    encounters an existing double-strand initiated at an internally primed site or because of RNA degradation.

    :param (str) bam_file: path to the bam file.
    :param (list) genes_of_interest: List of gene names - must match the names as they are written in the annotation file
    :param (bool) gene_dict
    :return: (list) List of reads. They will be a filtered subset of the reads in the input bam file.
    """

    filtered_reads = []

    gene_data = get_first_record_gene_annotation_data(genes_of_interest, annotation_file)

    for gene in gene_data:

        exons = [int(exon) for exon in gene_data[gene]['exons'].keys()]

        forward_first_exon = str(
            max(exons) if gene_data[gene]['strand'] == '-' else min(exons))  # Will be last exon in reverse genes
        forward_last_exon = str(min(exons) if gene_data[gene]['strand'] == '-' else max(exons))

        # Forward refers to the standard p to q direction - i.e. the forward strand
        forward_first_exon_last_base = int(gene_data[gene]['exons'][forward_first_exon][
                                               'stop']) - 2  # Because pysam is 0-based while SAM files are one-based, we subtract (two to be safe)
        forward_last_exon_first_base = int(gene_data[gene]['exons'][forward_last_exon]['start']) + 2  # Same thinking
        before_first_exon_start = int(gene_data[gene]['exons'][forward_first_exon]['start']) - 3000
        before_first_exon_stop = int(gene_data[gene]['exons'][forward_first_exon]['start']) - 2999
        after_last_exon_start = int(gene_data[gene]['exons'][forward_last_exon]['stop']) + 2000
        after_last_exon_stop = int(gene_data[gene]['exons'][forward_last_exon]['stop']) + 2001


        sam_file = pysam.AlignmentFile(bam_file, "rb")
        #
        # print('Creating bam file with reads spanning the full gene')
        before_first_exon_reads = set([read.query_name for read in
                                       sam_file.fetch(contig=gene_data[gene]['contig'], start=before_first_exon_start,
                                                      stop=before_first_exon_stop)])
        after_last_exon_reads = set([read.query_name for read in
                                     sam_file.fetch(contig=gene_data[gene]['contig'], start=after_last_exon_start,
                                                    stop=after_last_exon_stop)])
        start_reads = {read.query_name: read for read in
                       sam_file.fetch(contig=gene_data[gene]['contig'], start=forward_first_exon_last_base - 1,
                                      stop=forward_first_exon_last_base)}  # Query spans only one base - we get all reads containing that base

        filtered_reads_counter = 0
        reads_overlapping_gene = str(len([read for read in sam_file.fetch(contig=gene_data[gene]['contig'], start=forward_first_exon_last_base,
                                   stop=forward_last_exon_first_base)]))
        for read in sam_file.fetch(contig=gene_data[gene]['contig'], start=forward_last_exon_first_base,
                                   stop=forward_last_exon_first_base + 1):  # End reads
            if read.query_name in start_reads and not (
                    read.query_name in before_first_exon_reads or read.query_name in after_last_exon_reads):
                filtered_reads.append(read)
                filtered_reads_counter += 1
        print('Filtered {} full-length-reads from {} reads for the gene: {}'.format(filtered_reads_counter, reads_overlapping_gene, gene))
    return filtered_reads


def get_cigar_iter(cigar_string):
    """
    Given a CIGAR string it returns an iterator with each item being a CIGAR entry, such as 25M, 3I, 6D, etc.
    :param cigar_string:
    :return: Iterator
    """
    start_index = 0
    for i, char in enumerate(cigar_string):
        if not char.isnumeric():
            yield cigar_string[start_index: i + 1]
            start_index = i + 1


def get_read_exons(read):
    """
    Extracts exons (start and stop coordinates) from the CIGAR sequence
    :param read:
    :return: List of (start, stop) tuples for exons - both are inclusive and equals first base and last base.
    """

    # It looks as if minimap only outputs M, D, I, and N (haven't seen any X at mismatches so far - perhaps that takes a stretch of mismatched bases)

    consumes_reference = 'MDX='  # See Sequence Alignment/Map Format Specification. We only update coordinate relative to the reference if entries contain one of these letters
    exons = []
    first_exon_start = read.reference_start + 1  # + 1 due to difference in indexing (samfile pysam)
    exon = [first_exon_start, first_exon_start - 1]  # Start coordinate is after any clippings
    for cigar_entry in get_cigar_iter(read.cigarstring):
        entry_type = cigar_entry[-1]
        if entry_type in consumes_reference:
            exon[1] += int(cigar_entry[:-1])
        elif entry_type == 'N':
            exons.append((exon[0], exon[1]))  # Exon coordinates are both inclusive
            next_exon_start = exon[1] + int(cigar_entry[:-1]) + 1
            exon = [next_exon_start, next_exon_start - 1]
    exons.append((exon[0], exon[1]))
    return exons


def exon_is_contained(query_exon, reference_exon, position_slag=20):
    """
    Determines if a query exon is contained within a reference_exon given a defined slag.
    :param query_exon:
    :param reference_exon:
    :param position_slag:
    :return: Bool
    """

    ref_start = reference_exon[0] - position_slag
    ref_end = reference_exon[1] + position_slag

    return (query_exon[0] >= ref_start) and (query_exon[1] <= ref_end)


def exon_matches_gene_exon(exon, gene_exons, position_slag=20):
    """

    :param exon:
    :param gene_exons:      A list of exon tuples sorted by their start coordinate
    :param position_slag:
    :return: Bool
    """
    first_positions = [x[0] for x in gene_exons]

    """
    If a gene-exon is a "match", the index position will be after the matching gene-exon in cases where the exon's
    start position is >= the start position of the gene-exon. If the exon's start position is smaller, it will be
    before. So, we should test against the gene-exons at index: index_position and index_position - 1 if possible
    (index_position may be larger than the index of the last gene exon).
    
    There is also a theoretical possibility that there are several gene exons with same start coordinate and different
    stop coordinates (the same is possible for the stop coordinate, but we only sort relative to the start coordinate).
    Because of that, we must theoretically test gene_exons at index_position and index_position - i, if the gene-exon at 
    index_position - i has the same start coordinate as the one at index_position. Likewise, if gene_exons at
    index_position + i have the same start coordinate, these must be tested as well. However, we don't do this as we
    assume that there will be other exons without such ambiguities that will be matches making it possible to reach the
    minimal match anyways.
    """
    bisect_position = bisect(first_positions, exon[0])
    index_position = bisect_position if bisect_position < len(first_positions) else len(first_positions) - 1
    before_position = index_position - 1 if index_position > 0 else 0
    return exon_is_contained(exon, gene_exons[index_position], position_slag) or exon_is_contained(exon, gene_exons[before_position], position_slag)


def is_potential_isoform(read, gene_exons, min_match=2, position_slag=20):
    """
    Boolean function that decides if a read is an RNA isoform of a gene. In order to allow fragmented reads (due to
    internal priming) to be identified as isoforms, a number of exons should match - either by start or stop coordinate
    or both. A bit of slag relative start-stop coordinates must be allowed in order to accommodate small alignment
    discrepancies where for instance a base is moved from an acceptor to a donor or vice versa.

    :param read:            A pysam read
    :param gene_exons:      A sorted list of exon tuples for the gene
    :param position_slag:   Allowed mismatch between start/stop coordinates of an exon
    :param min_match:       The number of exon matches to indicate isoform
    :return: Bool           Returns True if matches >= min_match are found, False otherwise
    """

    read_exons = get_read_exons(read)
    matches = 0
    for exon in read_exons:
        if exon_matches_gene_exon(exon, gene_exons, position_slag):
            matches += 1
        if matches == min_match:
            return True
    return False


def get_exon_overlapping_reads(bam_file, genes_of_interest, annotation_file):
    """
    Filters and returns reads that overlap sufficiently with exons in the genes of interest. Sufficiently is defined as at least two exon matches.
    :param (str) bam_file:          Path to bam file
    :param ([]) genes_of_interest:  List of gene names to do filtering for (must match the names in the annotation file)
    :param (str) annotation_file:   Path to annotation file
    :return:
    """

    gene_data = get_gene_exons_annotation_data(genes_of_interest, annotation_file)

    sam_file = pysam.AlignmentFile(bam_file, "rb")

    filtered_reads = []

    for gene in gene_data:
        gene_exons = gene_data[gene]['exons']
        gene_overlapping_reads = sam_file.fetch(contig=gene_data[gene]['contig'], start=int(gene_data[gene]['exons'][0][0]), stop=int(gene_data[gene]['exons'][-1][1]))

        read_count = 0
        filtered_count = 0
        for read in gene_overlapping_reads:
            read_count += 1
            if is_potential_isoform(read, gene_exons):
                filtered_count += 1
                filtered_reads.append(read)

        print( 'Filtered {} reads based on exon match from {} reads for the gene: {}'.format(filtered_count, read_count, gene))

    return filtered_reads


def run_with_snakemake_input(log_file):
    import sys
    orig_stdout = sys.stdout
    orig_stderr = sys.stderr
    log = open(log_file, 'w')
    sys.stdout = log
    sys.stderr = log

    bam_file = snakemake.input['BAM']
    annotation_file = snakemake.input['ANNO']
    genes_of_interest = snakemake.params['GENES']
    full_length_out = snakemake.output['full_length']
    isoform_match_out = snakemake.output['isoform_match']

    full_length_filtered_reads = get_full_length_rna_reads(bam_file, genes_of_interest, annotation_file)
    isoform_filtered_reads = get_exon_overlapping_reads(bam_file, genes_of_interest, annotation_file)

    write_bam_file(full_length_out, full_length_filtered_reads, bam_file)
    write_bam_file(isoform_match_out, isoform_filtered_reads, bam_file)

    sys.stdout = orig_stdout
    sys.stderr = orig_stderr
    log.close()


def main():
    import argparse

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-b', '--bam_file', help='Path to bam file - as relative or absolute path. A corresponding bam index file must be present in the same folder', required=True)
    parser.add_argument('-a', '--annotation_file', help='Path to annotation file - as relative or absolute path', required=True)
    parser.add_argument('-g', '--genes', help='Comma separated (no spaces) names of genes of interest', required=True)
    parser.add_argument('-s', '--sample', help='Sample_name', required=False)

    args = parser.parse_args()
    bam_file = args.bam_file
    annotation_file = args.annotation_file
    genes_of_interest = args.genes.split(',')
    full_length_out = args.sample + '_full_length_filtered_sorted.bam' if args.sample else 'full_length_filtered_sorted.bam'
    isoform_match_out = args.sample + '_isoform_match_filtered_sorted.bam' if args.sample else 'isoform_match_filtered_sorted.bam'


    full_length_filtered_reads = get_full_length_rna_reads(bam_file, genes_of_interest, annotation_file)

    isoform_filtered_reads = get_exon_overlapping_reads(bam_file, genes_of_interest, annotation_file)

    write_bam_file(full_length_out, full_length_filtered_reads, bam_file)
    write_bam_file(isoform_match_out, isoform_filtered_reads, bam_file)


if 'snakemake' in locals() or 'snakemake' in globals():
    log_file = str(snakemake.log[0])
    run_with_snakemake_input(log_file)

else:
    main()
print('RNA-isoform filtering completed.')
