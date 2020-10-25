#! /usr/bin/env python
"""
This tool combines one of several results of the somatic_pipeline.py and/or germline_pipeline.py
to create an unified table with all the variants filtrated and their epitopes.
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

def main(somatic,
         somatic_names,
         germline,
         germline_names,
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
         tumor_coverage_germline,
         tumor_var_depth_germline,
         tumor_var_freq_germline,
         num_callers_germline):

    if not somatic and not germline:
        sys.stderr.write("Error, no variants given as input (somatic or germline).\n")
        sys.exit(1)

    #TODO add sanity check of paramters

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

    if len(somatic) > 0 and len(somatic) == len(somatic_names):
        print('Loading somatic variants..')
        for file, name in zip(somatic, somatic_names):
            variants = filter_variants_somatic(file,
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

    if len(germline) > 0 and len(germline) == len(germline_names):
        print('Loading germline variants..')
        for file, name in zip(germline, germline_names):
            variants = filter_variants_germline(file,
                                                tumor_coverage_germline,
                                                tumor_var_depth_germline,
                                                tumor_var_freq_germline,
                                                num_callers_germline,
                                                cDNA_seq_dict,
                                                AA_seq_dict)
            for variant in variants:
                variant_dict[variant.key].append((variant, name))

    counts_dict = defaultdict(int)
    gene_mean = -1
    gene_percentile = -1
    if os.path.isfile(rna_counts):
        print('Loading Gene counts..')
        counts_file = open(rna_counts)
        counts_file_lines = counts_file.readlines()
        _ = counts_file_lines.pop(0)
        header_counts = counts_file_lines.pop(0).strip().split('\t')
        for line in counts_file_lines:
            columns = line.strip().split('\t')
            gene_id = columns[header_counts.index('Geneid')]
            value = float(columns[-1])
            counts_dict[gene_id] = value
            counts_file.close()
        if len(counts_dict) > 0:
            # Compute mean expression and percentiles
            counts = list(counts_dict.values())
            gene_mean = statistics.mean(counts)
            gene_percentile = np.around(stats.percentileofscore(counts, 3))

    print('Creating merged variants..')
    header_final = 'Variant key\tSomatic samples (passing)\tNumber of Somatic samples (passing)\t' \
                   'Somatic samples (failing)\tNumber of Somatic samples (failing)\t' \
                   'Germline samples (passing)\tNumber of Germline samples (passing)\t' \
                   'Germline samples (failing)\tGermline of Somatic samples (failing)\tEffects\t' \
                   'cDNA change\tAA change\tEpitope creation flags\tWt Epitope\t'\
                   'Mut Epitope\tTranscript\tSomatic Callers(Name:NDP;NAD;NVAF;TDP;TAD;TVAF)\t'\
                   'Germline Callers (Name:TDP;TAD;TVAF)\tGeneCount info (gene;exp;mean;percentile)\n'

    final_file = open('overlap_final.txt', 'w')
    final_file.write(header_final)

    final_file_germline = open('overlap_final_germline_unique.txt', 'w')
    final_file_germline.write(header_final)

    final_file_discarded = open('overlap_final_discarded.txt', 'w')
    final_file_discarded.write(header_final)

    final_file_discarded_germline = open('overlap_final_discarded_germline.txt', 'w')
    final_file_discarded_germline.write(header_final)

    for key, value in variant_dict.items():

        # key = variant key
        # value = list of (Variant, sample_name)

        germline_name_pass = [name for variant, name in value if variant.type == 'germline' and variant.status]
        germline_name_fail = [name for variant, name in value if variant.type == 'germline' and not variant.status]
        germline_callers = ' '.join(['{}-{}'.format(name, variant.callers) for variant, name in value if variant.type == 'germline'])
        somatic_name_pass = [name for variant, name in value if variant.type == 'somatic' and variant.status]
        somatic_name_fail = [name for variant, name in value if variant.type == 'somatic' and not variant.status]
        somatic_callers = ' '.join(['{}-{}'.format(name, variant.callers) for variant, name in value if variant.type == 'somatic'])
        num_germline_pass = len(germline_name_pass)
        num_gemline_fail = len(germline_name_fail)
        num_somatic_pass = len(somatic_name_pass)
        num_somatic_fail = len(somatic_name_fail)

        # All variants shared variant key so their epitopes/effects/gene must be the same (we take the first variant then)
        epitopes, effects, gene = value[0][0].epitopes, value[0][0].effects, value[0][0].ensemblGene

        # Get gene exp. if a ny (gene should be the same in all the effects)
        gene_locus = "-"
        if gene is not None and gene in counts_dict:
            gene_count = counts_dict[gene]
            gene_locus = ';'.join([gene, str(gene_count), str(gene_mean), str(gene_percentile)])

        for epitope, effect in zip(epitopes, effects):
            to_write = '\t'.join(str(x) for x in [key,
                                                  ';'.join(somatic_name_pass), num_somatic_pass,
                                                  ';'.join(somatic_name_fail), num_somatic_fail,
                                                  ';'.join(germline_name_pass), num_germline_pass,
                                                  ';'.join(germline_name_fail), num_gemline_fail,
                                                  effect, epitope.dnamut, epitope.aamut, epitope.flags,
                                                  epitope.wtseq, epitope.mutseq, epitope.transcript,
                                                  somatic_callers, germline_callers, gene_locus])
            if num_somatic_pass >= 1:
                final_file.write(to_write + '\n')
            elif num_germline_pass >= 1:
                final_file_germline.write(to_write + '\n')
            elif num_somatic_fail >= 1:
                final_file_discarded.write(to_write + '\n')
            else:
                final_file_discarded_germline.write(to_write + '\n')

    final_file.close()
    final_file_germline.close()
    final_file_discarded.close()



if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--somatic', nargs='+', default=None, required=False,
                        help='List of annotated vcf files with the variants obtained with the somatic pipeline')
    parser.add_argument('--somatic-names', nargs='+', default=None, required=False,
                        help='List of names for each somatic file (to include in the report)')
    parser.add_argument('--germline', nargs='+', default=None, required=False,
                        help='List of annotated vcf files with the variants obtained with the germline pipeline')
    parser.add_argument('--germline-names', nargs='+', default=None, required=False,
                        help='List of names for each germline file (to include in the report)')
    parser.add_argument('--counts', type=str, default=None, required=False,
                        help='File with RNA gene counts from either the somatic or germline pipeline')
    parser.add_argument('--dictAA',
                        help='Path to a dictionary of transcript IDs to peptide sequences', required=True)
    parser.add_argument('--dictcDNA',
                        help='Path to a dictionary of transcript IDs to DNA sequences', required=True)
    parser.add_argument('--filter-somatic-tumor-cov', type=int, default=10, required=False, dest='tumor_coverage',
                        help='Filter for somatic variants tumor number of reads (coverage) (DP)')
    parser.add_argument('--filter-somatic-tumor-depth', type=int, default=4, required=False, dest='tumor_var_depth',
                        help='Filter for somatic variants tumor number of allelic reads (AD)')
    parser.add_argument('--filter-somatic-tumor-vaf', type=float, default=7, required=False, dest='tumor_var_freq',
                        help='Filter for somatic variants tumor variant allele frecuency (VAF)')
    parser.add_argument('--filter-somatic-normal-cov', type=int, default=10, required=False, dest='normal_coverage',
                        help='Filter for somatic variants normal number of reads (coverage) (DP)')
    parser.add_argument('--filter-somatic-tumor-normal-ratio', type=int, default=5, required=False, dest='t2n_ratio',
                        help='Filter for somatic variants tumor-normal VAF ratio')
    parser.add_argument('--filter-somatic-snv-callers', type=int, default=2, required=False,
                        choices=[1, 2, 3, 4], dest='num_callers',
                        help='Filter for somatic snvs variants number of callers required')
    parser.add_argument('--filter-somatic-indel-callers', type=int, default=1, required=False,
                        choices=[1, 2], dest='num_callers_indel',
                        help='Filter for somatic indels variants number of callers required')
    parser.add_argument('--filter-germline-tumor-cov', type=int, default=5, required=False,
                        dest='tumor_coverage_germline',
                        help='Filter for germline variants tumor number of reads (coverage) (DP)')
    parser.add_argument('--filter-germline-tumor-depth', type=int, default=2, required=False,
                        dest='tumor_var_depth_germline',
                        help='Filter for germline variants tumor number of allelic reads (AD)')
    parser.add_argument('--filter-germline-tumor-vaf', type=float, default=3, required=False,
                        dest='tumor_var_freq_germline',
                        help='Filter for germline variants tumor variant allele frecuency (VAF)')
    parser.add_argument('--filter-germline-callers', type=int, default=1, required=False,
                        choices=[1, 2], dest='num_callers_germline',
                        help='Filter for germline indels variants number of callers required')

    args = parser.parse_args()
    main(args.somatic,
         args.somatic_names,
         args.germline,
         args.germline_names,
         args.counts,
         os.path.abspath(args.dictcDNA),
         os.path.abspath(args.dictAA),
         args.tumor_coverage,
         args.tumor_var_depth,
         args.tumor_var_freq,
         args.normal_coverage,
         args.t2n_ratio,
         args.num_callers,
         args.num_callers_indel,
         args.tumor_coverage_germline,
         args.tumor_var_depth_germline,
         args.tumor_var_freq_germline,
         args.num_callers_germline
         )
