#! /usr/bin/env python
"""
This tool combines results from the dna_pipeline.py and/or rna_pipeline.py
to create an unified table with all the variants filtrated and their epitopes (for each effect).
The table contains useful information for post-analysis.

@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
import statistics
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as np
import math
from scipy import stats
import os
import sys
from collections import defaultdict
from hlapipeline.variants import *

def main(dna_variants,
         dna_names,
         rna_variants,
         rna_names,
         rna_counts,
         cDNA_DICT,
         AA_DICT,
         tumor_coverage,
         tumor_var_depth,
         tumor_var_freq,
         normal_coverage,
         t2n_ratio,
         num_callers,
         num_callers_indel,
         tumor_coverage_rna,
         tumor_var_depth_rna,
         tumor_var_freq_rna,
         num_callers_rna):

    if not dna_variants and not rna_variants:
        sys.stderr.write("Error, no variants given as input (DNA or RNA).\n")
        sys.exit(1)

    # TODO add sanity check for parameters

    AA_seq_dict = dict()
    with open(AA_DICT, "r") as handle:
        for line in handle.readlines():
            tokens = line.split(":")
            AA_seq_dict[tokens[0]] = tokens[1].strip()
    cDNA_seq_dict = dict()
    with open(cDNA_DICT, "r") as handle:
        for line in handle.readlines():
            tokens = line.split(":")
            cDNA_seq_dict[tokens[0]] = tokens[1].strip()

    variant_dict = defaultdict(list)

    if dna_variants and len(dna_variants) > 0 and len(dna_variants) == len(dna_names):
        print('Loading DNA variants..')
        for file, name in zip(dna_variants, dna_names):
            variants = filter_variants_dna(file,
                                           normal_coverage,
                                           tumor_coverage,
                                           tumor_var_depth,
                                           tumor_var_freq,
                                           t2n_ratio,
                                           num_callers,
                                           num_callers_indel,
                                           cDNA_seq_dict,
                                           AA_seq_dict)
            for variant in variants:
                variant_dict[variant.key].append((variant, name))

    if rna_variants and len(rna_variants) > 0 and len(rna_variants) == len(rna_names):
        print('Loading RNA variants..')
        for file, name in zip(rna_variants, rna_names):
            variants = filter_variants_rna(file,
                                           tumor_coverage_rna,
                                           tumor_var_depth_rna,
                                           tumor_var_freq_rna,
                                           num_callers_rna,
                                           cDNA_seq_dict,
                                           AA_seq_dict)
            for variant in variants:
                variant_dict[variant.key].append((variant, name))

    # TODO this could be done more elegantly and efficiently
    counts_dict = defaultdict(lambda: defaultdict(float))
    counts_stats = defaultdict(list)
    if rna_counts and len(rna_counts) > 0 and len(rna_counts) == len(rna_names):
        print('Loading Gene counts..')
        for file, name in zip(rna_counts, rna_names):
            counts_file = open(file)
            counts_file_lines = counts_file.readlines()
            _ = counts_file_lines.pop(0)
            header_counts = counts_file_lines.pop(0).strip().split('\t')
            for line in counts_file_lines:
                columns = line.strip().split('\t')
                gene_id = columns[header_counts.index('Geneid')]
                value = float(columns[-1])
                counts_dict[name][gene_id] = value
            counts_file.close()

        if len(counts_dict) > 0:
            for name, gene_counts in counts_dict.items():
                counts = list(gene_counts.values())
                mean = np.around(statistics.mean(counts), 3)
                percentile = np.around(stats.percentileofscore(counts, 3), 3)
                counts_stats[name] = [mean, percentile]

    print('Creating merged variants..')
    header_final = 'Variant key\tDNA samples (passing)\tNumber of DNA samples (passing)\t' \
                   'DNA samples (failing)\tNumber of DNA samples (failing)\t' \
                   'RNA samples (passing)\tNumber of RNA samples (passing)\t' \
                   'RNA samples (failing)\tRNA of DNA samples (failing)\tEffects\t' \
                   'cDNA change\tAA change\tEpitope creation flags\tWt Epitope\t' \
                   'Mut Epitope\tTranscript\tDNA Callers Sample(Name:NDP;NAD;NVAF;TDP;TAD;TVAF)\t' \
                   'RNA Callers Sample(Name:TDP;TAD;TVAF)\tGeneCount info Sample(gene;exp;mean;percentile)\n'

    final_file = open('overlap_final.txt', 'w')
    final_file.write(header_final)

    final_file_rna = open('overlap_final_rna_unique.txt', 'w')
    final_file_rna.write(header_final)

    final_file_discarded = open('overlap_final_discarded.txt', 'w')
    final_file_discarded.write(header_final)

    final_file_discarded_rna = open('overlap_final_discarded_rna.txt', 'w')
    final_file_discarded_rna.write(header_final)

    for key, value in variant_dict.items():

        # key = variant key
        # value = list of (Variant, sample_name) tuples

        rna_name_pass = [name for variant, name in value if variant.type == 'rna' and variant.status]
        rna_name_fail = [name for variant, name in value if variant.type == 'rna' and not variant.status]
        rna_callers = ';'.join(
            ['{}:({})'.format(name, variant.callers) for variant, name in value if variant.type == 'rna'])
        dna_name_pass = [name for variant, name in value if variant.type == 'dna' and variant.status]
        dna_name_fail = [name for variant, name in value if variant.type == 'dna' and not variant.status]
        dna_callers = ';'.join(
            ['{}:({})'.format(name, variant.callers) for variant, name in value if variant.type == 'dna'])
        num_rna_pass = len(rna_name_pass)
        num_rna_fail = len(rna_name_fail)
        num_dna_pass = len(dna_name_pass)
        num_dna_fail = len(dna_name_fail)

        # All variants share variant key so their epitopes/effects/gene must be the same (we take the first variant)
        epitopes = value[0][0].epitopes
        effects = value[0][0].effects
        gene = value[0][0].geneName

        # Get gene exp. if any (gene should be the same in all the effects)
        gene_locus = []
        if gene is not None:
            for name, gene_counts in counts_dict.items():
                try:
                    gene_count = gene_counts[gene]
                    gene_mean, gene_percentile = counts_stats[name]
                    gene_locus.append('{}:({})'.format(name,
                                                       ';'.join([gene,
                                                                 str(gene_count),
                                                                 str(gene_mean),
                                                                 str(gene_percentile)])))
                except KeyError:
                    gene_locus.append("{}:-".format(name))
        else:
            gene_locus = ["-"]

        for epitope, effect in zip(epitopes, effects):
            to_write = '\t'.join(str(x) for x in [key,
                                                  ';'.join(dna_name_pass), num_dna_pass,
                                                  ';'.join(dna_name_fail), num_dna_fail,
                                                  ';'.join(rna_name_pass), num_rna_pass,
                                                  ';'.join(rna_name_fail), num_rna_fail,
                                                  effect, epitope.dnamut, epitope.aamut, epitope.flags,
                                                  epitope.wtseq, epitope.mutseq, epitope.transcript,
                                                  dna_callers, rna_callers, ';'.join(gene_locus)])
            if num_dna_pass >= 1:
                final_file.write(to_write + '\n')
            elif num_rna_pass >= 1:
                final_file_rna.write(to_write + '\n')
            elif num_dna_fail >= 1:
                final_file_discarded.write(to_write + '\n')
            else:
                final_file_discarded_rna.write(to_write + '\n')

    final_file.close()
    final_file_rna.close()
    final_file_discarded.close()
    final_file_discarded_rna.close()

if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--dna', nargs='+', default=None, required=False,
                        help='List of annotated vcf files with the variants obtained with the DNA pipeline')
    parser.add_argument('--dna-names', nargs='+', default=None, required=False,
                        help='List of names for each DNA sample/file (to include in the report)')
    parser.add_argument('--rna', nargs='+', default=None, required=False,
                        help='List of annotated vcf files with the variants obtained with the RNA pipeline')
    parser.add_argument('--rna-names', nargs='+', default=None, required=False,
                        help='List of names for each RNA sample/file (to include in the report)')
    parser.add_argument('--rna-counts', nargs='+', default=None, required=False,
                        help='List of gene counts files obtained with the RNA pipeline')
    parser.add_argument('--dictAA',
                        help='Path to a dictionary of transcript IDs to peptide sequences', required=True)
    parser.add_argument('--dictcDNA',
                        help='Path to a dictionary of transcript IDs to DNA sequences', required=True)
    parser.add_argument('--filter-dna-tumor-cov', type=int, default=10, required=False, dest='tumor_coverage',
                        help='Filter for DNA variants tumor number of reads (coverage) (DP). Default=10')
    parser.add_argument('--filter-dna-tumor-depth', type=int, default=4, required=False, dest='tumor_var_depth',
                        help='Filter for DNA variants tumor number of allelic reads (AD). Default=4')
    parser.add_argument('--filter-dna-tumor-vaf', type=float, default=7, required=False, dest='tumor_var_freq',
                        help='Filter for DNA variants tumor variant allele frequency (VAF). Default=7')
    parser.add_argument('--filter-dna-normal-cov', type=int, default=10, required=False, dest='normal_coverage',
                        help='Filter for DNA variants normal number of reads (coverage) (DP). Default=10')
    parser.add_argument('--filter-dna-tn-ratio', type=int, default=5, required=False, dest='t2n_ratio',
                        help='Filter for DNA variants tumor-normal VAF ratio. Default=5')
    parser.add_argument('--filter-dna-snv-callers', type=int, default=2, required=False,
                        choices=[1, 2, 3, 4], dest='num_callers',
                        help='Filter for DNA SNPs variants number of callers required. Default=2')
    parser.add_argument('--filter-dna-indel-callers', type=int, default=1, required=False,
                        choices=[1, 2], dest='num_callers_indel',
                        help='Filter for DNA indels variants number of callers required. Default=1')
    parser.add_argument('--filter-rna-tumor-cov', type=int, default=10, required=False,
                        dest='tumor_coverage_rna',
                        help='Filter for RNA variants tumor number of reads (coverage) (DP). Default=10')
    parser.add_argument('--filter-rna-tumor-depth', type=int, default=5, required=False,
                        dest='tumor_var_depth_rna',
                        help='Filter for RNA variants tumor number of allelic reads (AD). Default=5')
    parser.add_argument('--filter-rna-tumor-vaf', type=float, default=3, required=False,
                        dest='tumor_var_freq_rna',
                        help='Filter for RNA variants tumor variant allele frequency (VAF). Default=3')
    parser.add_argument('--filter-rna-callers', type=int, default=1, required=False,
                        choices=[1, 2], dest='num_callers_rna',
                        help='Filter for RNA variants number of callers required. Default=1')

    args = parser.parse_args()
    main(args.dna,
         args.dna_names,
         args.rna,
         args.rna_names,
         args.rna_counts,
         os.path.abspath(args.dictcDNA),
         os.path.abspath(args.dictAA),
         args.tumor_coverage,
         args.tumor_var_depth,
         args.tumor_var_freq,
         args.normal_coverage,
         args.t2n_ratio,
         args.num_callers,
         args.num_callers_indel,
         args.tumor_coverage_rna,
         args.tumor_var_depth_rna,
         args.tumor_var_freq_rna,
         args.num_callers_rna)
