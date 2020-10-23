#! /usr/bin/env python
"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""

#TODO Use a Class to store variants

import statistics
import argparse
import numpy as np
import math
from scipy import stats
import os
import sys
import vcfpy
from collections import defaultdict
import re
from Bio.Seq import translate

def translate_dna(seq):
    return translate(seq, to_stop=True)

def create_epitopes(ref, exonic_func, transcriptID, cDNA_strip, protein_strip, cDNA_seq, AA_seq):
    errors = 'Flags:'
    WT_25mer = '-'
    Mut_25mer = '-'
    position = int(re.findall(r'\d+', protein_strip)[0])
    cDNA_pos = int(re.findall(r'\d+', cDNA_strip)[0])
    protein_seq = AA_seq.get(transcriptID, 'AA_seq not present for this transcript').strip()
    ref_cDNA_seq = cDNA_seq.get(transcriptID, 'cDNA not present for this transcript').strip()
    # Nonsynonymous point mutations
    if exonic_func == 'nonsynonymous_SNV':
        # extract the AA change info
        ref_AA = protein_strip[protein_strip.find('.') + 1]
        var_AA = protein_strip[len(protein_strip) - 1]
        # gather AA seq for transcript
        if protein_seq == 'AA_seq not present for this transcript':
            errors += ' AA_seq not present for this transcript'
        else:
            # Check annotation is correct
            # Create mut and wt epitopes by extracting 12
            # AAs before the mutation, the mutation
            # and then 12 AAs after the mutation
            FASTA_AA = protein_seq[position - 1:position]
            if FASTA_AA == ref_AA:
                if position == 0:
                    errors += ' can not code for this mutated AA_position'
                else:
                    end = position + 12 if position + 12 < len(protein_seq) else None
                    start = position - 13 if position > 13 else 0
                    WT_25mer = protein_seq[start:end]
                    Mut_25mer = protein_seq[start:position - 1] + var_AA + protein_seq[position:end]
                if position == 1:
                    errors += ' mutation occurs in start codon'
            else:
                errors += ' Ref in AA_seq does not match file Ref {} {}'.format(FASTA_AA, ref_AA)
    # frameshift/non-frameshift insertions/deletions/substitutions
    elif 'frameshift' in exonic_func:
        if ref_cDNA_seq == 'cDNA not present for this transcript':
            errors += ' cDNA not present for this transcript'
        else:
            if exonic_func in ['frameshift_deletion', 'nonframeshift_deletion']:
                len_del = len(ref)
                mut_cDNA_left = ref_cDNA_seq[0:cDNA_pos - 1]
                mut_cDNA_right = ref_cDNA_seq[cDNA_pos + len_del - 1:]
                mut_cDNA_seq = mut_cDNA_left + mut_cDNA_right
            elif exonic_func in ['frameshift_insertion', 'nonframeshift_insertion']:
                if 'dup' in cDNA_strip:
                    ins = cDNA_strip[int(cDNA_strip.find('dup')) + 3:]
                elif 'ins' in cDNA_strip:
                    ins = cDNA_strip[int(cDNA_strip.find('ins')) + 3:]
                else:
                    errors += ' could not find mutation in cDNA'
                    ins = ''
                mut_cDNA_left = ref_cDNA_seq[0:cDNA_pos]
                mut_cDNA_right = ref_cDNA_seq[cDNA_pos:]
                mut_cDNA_seq = mut_cDNA_left + ins + mut_cDNA_right
            elif exonic_func in ['frameshift_substitution', 'nonframeshift_substitution']:
                if 'delins' in cDNA_strip:
                    subs = cDNA_strip[int(cDNA_strip.find('delins')) + 6:]
                else:
                    errors += ' could not find mutation in cDNA'
                    subs = ''
                mut_cDNA_left = ref_cDNA_seq[0:cDNA_pos]
                mut_cDNA_right = ref_cDNA_seq[cDNA_pos:]
                mut_cDNA_seq = mut_cDNA_left + subs + mut_cDNA_right
            else:
                mut_cDNA_seq = ref_cDNA_seq
                errors += ' unknown exonic function {}'.format(exonic_func)
            if not protein_strip.startswith('p.'):
                position = 0
            elif protein_strip.startswith('p.X'):
                position = 0
                errors += ' mutation occurs in stop codon'
            # Create mut and wt epitopes by extracting 12
            # AAs before and after the mutation
            if position == 0:
                errors += ' can not code for this mutated AA_position'
            else:
                ref_FASTA = str(translate_dna(ref_cDNA_seq.replace(' ', '')))
                mut_FASTA = str(translate_dna(mut_cDNA_seq.replace(' ', '')))
                end = position + 12 if position + 12 < len(ref_FASTA) else None
                start = position - 13 if position > 13 else 0
                WT_25mer = ref_FASTA[start:end]
                Mut_25mer = mut_FASTA[start:]
            if not ref_cDNA_seq.startswith('ATG'):
                errors += ' no ATG start codon for this transcript cDNA'
            if position == 1:
                errors += ' mutation occurs in start codon'
    elif 'stop' in exonic_func:
        position = ''.join([s for s in protein_strip if s.isdigit()])
        errors += ' stop mutation'
    else:
        errors += ' unknown exonic function {}'.format(exonic_func)
    # NOTE just to make sure that they are not empty
    WT_25mer = "-" if not WT_25mer else WT_25mer
    Mut_25mer = "-" if not Mut_25mer else Mut_25mer
    return(position, errors, WT_25mer, Mut_25mer)

def add_flags(transcript, variant_key, transcript_info, mer_len=25):
    cDNA_flanks = math.floor(mer_len / 2) * 3
    tups = sorted(list(zip(transcript_info[transcript]['position'],
                           transcript_info[transcript]['mutation'],
                           transcript_info[transcript]['aa'],
                           transcript_info[transcript]['status'],
                           transcript_info[transcript]['cDNA'],
                           transcript_info[transcript]['variant_key'])))

    target_pos = [i[0] for i in tups if i[5] == variant_key][0]
    prev_status = None
    prev_position = None
    ss_seen = ''
    flags = ''
    last_variant_to_review = False
    for pos, mut, aa, status, cdna, vkey in tups:
        # Variant in tuple not the one we are looking for
        # does pass filters in at least one sample and we haven't
        # encountered our variant we are reviewing yet
        if vkey != variant_key:
            if not last_variant_to_review:
                status_str = 'Passing' if status else 'Failing'
                if mut == 'stopgain' and ss_seen is None:
                    ss_seen = '{} filters upstream stopgain ({})'.format(status_str, aa)
                elif mut == 'stopgain' and ss_seen is not None:
                    ss_seen += ', {} filters upstream stopgain ({})'.format(status_str, aa)
                if 'frame' in mut and ss_seen is None:
                    ss_seen = '{} filters upstream fs ({})'.format(status_str, aa)
                elif 'frame' in mut and ss_seen is not None:
                    ss_seen += ', {} filters upstream fs ({})'.format(status_str, aa)
                # Create flag if this variant is within distance of note to ours
                if target_pos - pos == 1:
                    flags += ', {} filter dinucleotide change ({})'.format(status_str, cdna)
                elif target_pos - pos <= cDNA_flanks:
                    flags += ', {} filter mutation within {} nt ({})'.format(status_str, cDNA_flanks, aa)
                flags += ', ' + ss_seen
            elif pos > prev_position:
                if pos - prev_position == 1 and status and prev_status:
                    flags += ', Passing filter dinucleotide change ({})'.format(cdna)
                elif pos - prev_position == 1 and not prev_status:
                    flags += ', Failing dinucleotide change ({})'.format(cdna)
                elif pos - prev_position <= cDNA_flanks and prev_status:
                    flags += ', Passing mutation within {} nt ({})'.format(cDNA_flanks, aa)
                elif pos - prev_position <= cDNA_flanks and not prev_status:
                    flags += ', Failing mutation within {} nt ({})'.format(cDNA_flanks, aa)
                flags += ', ' + ss_seen
                # Stop once we get more than 15aa away
                if pos - prev_position > cDNA_flanks:
                    break
        # When we encounter our variant
        # as the 1st variant we want to
        # change the flag that we are reviewing
        # to true so we can wind down
        else:
            last_variant_to_review = True
        prev_position = pos
        prev_status = status
    return flags if flags != '' else '-'

def effects(record, cDNA_seq, AA_seq):
    funcensGene = ''.join(record.INFO['ExonicFunc.ensGene'])
    has_func_ens =  'nonsynonymous' in funcensGene or 'frame' in funcensGene
    funcknownGene = ''.join(record.INFO['ExonicFunc.knownGene'])
    has_func_known = 'nonsynonymous' in funcknownGene or 'frame' in funcknownGene
    funcRefGene = ''.join(record.INFO['ExonicFunc.refGene'])
    has_func_ref = 'nonsynonymous' in funcRefGene or 'frame' in funcRefGene
    already_processed = set()
    records_to_print = []
    dbs_seen_in = defaultdict(list)
    if has_func_ens:
        for mutation in record.INFO['AAChange.ensGene']:
            if len(mutation.split(':')) == 5:
                gene, transcript, exon, mut_dna, mut_aa = mutation.split(':')
                pos, flags, wtmer, mutmer = create_epitopes(record.REF, funcensGene, transcript, mut_dna, mut_aa, cDNA_seq, AA_seq)
                if mutmer not in already_processed and mutmer != "-":
                    already_processed.add(mutmer)
                    records_to_print.append((mut_dna, mut_aa, flags, wtmer, mutmer))
                    dbs_seen_in[mutmer].append('ensGene_{}_{}_{}'.format(gene, transcript, funcensGene))
                elif mutmer != "-":
                    dbs_seen_in[mutmer].append('ensGene_{}_{}_{}'.format(gene, transcript, funcensGene))
    if has_func_known:
        for mutation in record.INFO['AAChange.knownGene']:
            if len(mutation.split(':')) == 5:
                gene, transcript, exon ,mut_dna, mut_aa = mutation.split(':')
                pos, flags, wtmer, mutmer = create_epitopes(record.REF, funcknownGene, transcript, mut_dna, mut_aa, cDNA_seq, AA_seq)
                if mutmer not in already_processed and mutmer != "-":
                    already_processed.add(mutmer)
                    records_to_print.append((mut_dna, mut_aa, flags, wtmer, mutmer))
                    dbs_seen_in[mutmer].append('knownGene_{}_{}_{}'.format(gene, transcript, funcknownGene))
                elif mutmer != "-":
                    dbs_seen_in[mutmer].append('knownGene_{}_{}_{}'.format(gene, transcript, funcknownGene))
    if has_func_ref:
        for mutation in record.INFO['AAChange.refGene']:
            if len(mutation.split(':')) == 5:
                gene, transcript, exon, mut_dna, mut_aa = mutation.split(':')
                pos, flags, wtmer, mutmer = create_epitopes(record.REF, funcRefGene, transcript, mut_dna, mut_aa, cDNA_seq, AA_seq)
                if mutmer not in already_processed and mutmer != "-":
                    already_processed.add(mutmer)
                    records_to_print.append((mut_dna, mut_aa, flags, wtmer, mutmer))
                    dbs_seen_in[mutmer].append('refGene_{}_{}_{}'.format(gene, transcript, funcRefGene))
                elif mutmer != "-":
                    dbs_seen_in[mutmer].append('refGene_{}_{}_{}'.format(gene, transcript, funcRefGene))
    return(records_to_print, dbs_seen_in)

def filter_variants_germline(file, tumor_coverage, tumor_var_depth, tumor_var_freq, num_callers, cDNA_seq, AA_seq):

    variants = list()
    reader = vcfpy.Reader.from_path(file)
    for record in reader:

        funcensGene = ''.join(record.INFO['ExonicFunc.ensGene'])
        has_func_ens = 'nonsynonymous' in funcensGene or 'frame' in funcensGene
        funcknownGene = ''.join(record.INFO['ExonicFunc.knownGene'])
        has_func_known = 'nonsynonymous' in funcknownGene or 'frame' in funcknownGene
        funcRefGene = ''.join(record.INFO['ExonicFunc.refGene'])
        has_func_ref = 'nonsynonymous' in funcRefGene or 'frame' in funcRefGene

        if has_func_ens or has_func_known or has_func_ref:

            called = {x.sample: x.data for x in record.calls if x.called}
            filtered = dict()

            try:
                if 'HaplotypeCaller' in called and 'PASS' in record.FILTER:
                    tumor_DP = int(called['HaplotypeCaller']['DP'])
                    tumor_AD = int(called['HaplotypeCaller']['AD'][1])
                    tumor_VAF = tumor_AD / float(tumor_DP) * 100
                    if tumor_DP >= tumor_coverage and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth:
                        filtered['HaplotypeCaller'] = '{};{};{}'.format(tumor_DP, tumor_AD, tumor_VAF)

                if 'varscan' in called and 'PASS' in record.FILTER:
                    tumor_DP = int(called['varscan']['DP'])
                    tumor_AD = int(called['varscan']['AD'][0])
                    tumor_VAF = float(called['varscan']['FREQ'][0].replace('%', ''))
                    if tumor_DP >= tumor_coverage and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth:
                        filtered['varscan'] = '{};{};{}'.format(tumor_DP, tumor_AD, tumor_VAF)

            except KeyError:
                continue

            is_valid = len(filtered.keys()) >= num_callers
            records_to_print, dbs_seen_in = effects(record, cDNA_seq, AA_seq)
            for r in records_to_print:
                callers = '|'.join(['{}:{}'.format(key, value) for key, value in filtered.items()])
                variant_key = "{}:{} {}>{}".format(record.CHROM, record.POS, record.REF, record.ALT[0].serialize())
                variants.append((variant_key, ';'.join(set(dbs_seen_in[r[-1]])), r, callers, is_valid))

    return variants

def filter_variants_somatic(file, normal_coverage, tumor_coverage, tumor_var_depth,
                            tumor_var_freq, t2n_ratio, num_callers, num_callers_indel, cDNA_seq, AA_seq):

    variants = list()
    reader = vcfpy.Reader.from_path(file)
    for record in reader:

        funcensGene = ''.join(record.INFO['ExonicFunc.ensGene'])
        has_func_ens = 'nonsynonymous' in funcensGene or 'frame' in funcensGene
        funcknownGene = ''.join(record.INFO['ExonicFunc.knownGene'])
        has_func_known = 'nonsynonymous' in funcknownGene or 'frame' in funcknownGene
        funcRefGene = ''.join(record.INFO['ExonicFunc.refGene'])
        has_func_ref = 'nonsynonymous' in funcRefGene or 'frame' in funcRefGene

        if has_func_ens or has_func_known or has_func_ref:

            called = {x.sample: x.data for x in record.calls if x.called}
            filtered = dict()

            try:
                if 'NORMAL.mutect' in called and 'TUMOR.mutect' in called and 'PASS' in record.FILTER:
                    normal_DP = int(called['NORMAL.mutect']['DP'])
                    normal_AD = int(called['NORMAL.mutect']['AD'][1])
                    normal_VAF = float(called['NORMAL.mutect']['AF'][0]) * 100
                    tumor_DP = int(called['TUMOR.mutect']['DP'])
                    tumor_AD = int(called['TUMOR.mutect']['AD'][1])
                    tumor_VAF = float(called['TUMOR.mutect']['AF'][0]) * 100
                    tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                    if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                            and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                            and tumor_normal_ratio >= t2n_ratio:
                        filtered['mutect'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                        normal_AD,
                                                                        normal_VAF,
                                                                        tumor_DP,
                                                                        tumor_AD,
                                                                        tumor_VAF)

                if 'NORMAL.somaticsniper' in called and 'TUMOR.somaticsniper' in called:
                    normal_DP = int(called['NORMAL.somaticsniper']['DP'])
                    normal_AD = sum(called['NORMAL.somaticsniper']['DP4'][2:])
                    normal_VAF = (normal_AD / float(normal_DP)) * 100
                    tumor_DP = int(called['TUMOR.somaticsniper']['DP'])
                    tumor_AD = sum(called['TUMOR.somaticsniper']['DP4'][2:])
                    tumor_VAF = (tumor_AD / float(tumor_DP)) * 100
                    tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                    is_somatic = int(called['TUMOR.somaticsniper']['SS']) == 2
                    if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                            and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                            and tumor_normal_ratio >= t2n_ratio and is_somatic:
                        filtered['somaticsniper'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                               normal_AD,
                                                                               normal_VAF,
                                                                               tumor_DP,
                                                                               tumor_AD,
                                                                               tumor_VAF)

                if ('NORMAL.varscan' in called and 'TUMOR.varscan' in called) \
                        or ('NORMAL.varscan_indel' in called and 'TUMOR.varscan_indel' in called) \
                        and 'PASS' in record.FILTER and 'SOMATIC' in record.INFO:
                    label_index = 'varscan' if 'NORMAL.varscan' in called else 'varscan_indel'
                    normal_DP = int(called['NORMAL.{}'.format(label_index)]['DP'])
                    normal_AD = sum(called['NORMAL.{}'.format(label_index)]['DP4'][2:])
                    normal_VAF = float(called['NORMAL.{}'.format(label_index)]['FREQ'][0].replace('%', ''))
                    tumor_DP = int(called['TUMOR.{}'.format(label_index)]['DP'])
                    tumor_AD = sum(called['TUMOR.{}'.format(label_index)]['DP4'][2:])
                    tumor_VAF = float(called['TUMOR.{}'.format(label_index)]['FREQ'][0].replace('%', ''))
                    tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                    if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                            and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                            and tumor_normal_ratio >= t2n_ratio:
                        filtered[label_index] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                           normal_AD,
                                                                           normal_VAF,
                                                                           tumor_DP,
                                                                           tumor_AD,
                                                                           tumor_VAF)

                if 'NORMAL.strelka' in called and 'TUMOR.strelka' in called and 'PASS' in record.FILTER:
                    ref_index = record.REF + 'U'
                    alt_index = str(record.ALT[0].serialize()) + 'U'
                    # normal_DP = int(called['NORMAL.strelka']['DP'])
                    normal_AD1 = int(called['NORMAL.strelka'][ref_index][0])
                    normal_AD2 = int(called['NORMAL.strelka'][alt_index][0])
                    normal_DP = normal_AD1 + normal_AD2
                    normal_VAF = (normal_AD2 / float(normal_DP)) * 100
                    # tumor_DP = int(called['TUMOR.strelka']['DP'])
                    tumor_AD1 = int(called['TUMOR.strelka'][ref_index][0])
                    tumor_AD2 = int(called['TUMOR.strelka'][alt_index][0])
                    tumor_DP = tumor_AD1 + tumor_AD2
                    tumor_VAF = (tumor_AD2 / float(tumor_DP)) * 100
                    tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                    if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                            and tumor_VAF >= tumor_var_freq and tumor_AD2 >= tumor_var_depth \
                            and tumor_normal_ratio >= t2n_ratio:
                        filtered['strelka'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                         normal_AD,
                                                                         normal_VAF,
                                                                         tumor_DP,
                                                                         tumor_AD,
                                                                         tumor_VAF)

                if 'NORMAL.strelka_indel' in called and 'TUMOR.strelka_indel' in called and 'PASS' in record.FILTER:
                    # normal_DP = int(called['NORMAL.strelka_indel']['DP'])
                    normal_AD1 = int(called['NORMAL.strelka_indel']['TAR'][0])
                    normal_AD2 = int(called['NORMAL.strelka_indel']['TIR'][0])
                    normal_DP = normal_AD1 + normal_AD2
                    normal_VAF = (normal_AD2 / float(normal_DP)) * 100
                    # tumor_DP = int(called['TUMOR.strelka_indel']['DP'])
                    tumor_AD1 = int(called['TUMOR.strelka_indel']['TAR'][0])
                    tumor_AD2 = int(called['TUMOR.strelka_indel']['TIR'][0])
                    tumor_DP = tumor_AD1 + tumor_AD2
                    tumor_VAF = (tumor_AD2 / float(tumor_DP)) * 100
                    tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                    if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                            and tumor_VAF >= tumor_var_freq and tumor_AD2 >= tumor_var_depth \
                            and tumor_normal_ratio >= t2n_ratio:
                        filtered['strelka'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                         normal_AD,
                                                                         normal_VAF,
                                                                         tumor_DP,
                                                                         tumor_AD,
                                                                         tumor_VAF)
            except KeyError:
                continue

            is_valid = len([x for x in filtered.keys() if 'indel' not in x]) >= num_callers \
                       or len([x for x in filtered.keys() if 'indel' in x]) >= num_callers_indel

            records_to_print, dbs_seen_in = effects(record, cDNA_seq, AA_seq)

            for r in records_to_print:
                callers = '|'.join(['{}:{}'.format(key, value) for key, value in filtered.items()])
                variant_key = "{}:{} {}>{}".format(record.CHROM, record.POS, record.REF, record.ALT[0].serialize())
                variants.append((variant_key, ';'.join(set(dbs_seen_in[r[-1]])), r, callers, is_valid))

    return variants

def main(somatic, somatic_names, germline, germline_names, rna_counts, cDNA_DICT, AA_DICT):

    if not somatic and not germline:
        sys.stderr.write("Error, no variants given as input (somatic or germline).\n")
        sys.exit(1)

    normal_coverage = 10
    tumor_coverage = 10
    tumor_var_depth = 4
    tumor_var_freq = 7
    t2n_ratio = 5
    num_callers = 2
    num_callers_indel = 1

    tumor_coverage_germline = 5
    tumor_var_depth_germline = 2
    tumor_var_freq_germline = 3
    num_callers_germline = 1

    variant_dict = defaultdict(lambda: defaultdict(list))

    AA_seq = dict()
    with open(AA_DICT, "r") as handle:
        for line in handle.readlines():
            tokens = line.split(":")
            AA_seq[tokens[0]] = tokens[1].strip()
    cDNA_seq = dict()
    with open(cDNA_DICT, "r") as handle:
        for line in handle.readlines():
            tokens = line.split(":")
            cDNA_seq[tokens[0]] = tokens[1].strip()

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
                                               cDNA_seq,
                                               AA_seq)
            variants.append(name)
            variant_dict[variants[0]]['somatic'].append(variants)
        
<<<<<<< HEAD
    if len(germline) > 0 and len(germline) == len(germline_names):
        print('Loading germline variants..')
        for file, name in zip(germline, germline_names):
            variants = filter_variants_germline(file,
                                                tumor_coverage_germline,
                                                tumor_var_depth_germline,
                                                tumor_var_freq_germline,
                                                num_callers_germline,
                                                cDNA_seq,
                                                AA_seq)
            variants.append(name)
            variant_dict[variants[0]]['germline'].append(variants)

    counts_dict = defaultdict(int)
    gene_mean = -1
    gene_percentile = -1
    if os.path.isfile(rna_counts):
        print('Loading Gene counts..')
        counts_file = open(rna_counts)
=======
    variant_dict = {}

    print('Loading somatic variants..')
    for file in somatic if somatic else []:
        somatic_nonsyn = open(file)
        somatic_nonsyn_lines = somatic_nonsyn.readlines()
        header_somatic = somatic_nonsyn_lines.pop(0).strip().split('\t')
        for line in somatic_nonsyn_lines:
            columns = line.strip().split('\t')
            variant_key = columns[header_somatic.index('VARIANT-KEY')]
            sample = columns[header_somatic.index('SAMPLE_ID')]
            if variant_key not in variant_dict:
                variant_dict[variant_key] = {}
                variant_dict[variant_key]['somatic'] = {}
            if sample not in variant_dict[variant_key]['somatic']:
                variant_dict[variant_key]['somatic'][sample] = {}
            # Compute coverage and pass/fail
            N_cov = int(columns[header_somatic.index('NCOV')])
            T_cov = int(columns[header_somatic.index('TCOV')])
            T_freq = float(columns[header_somatic.index('TVAF')].replace('%', ''))
            N_freq = float(columns[header_somatic.index('NVAF')].replace('%', ''))
            T_reads = int(columns[header_somatic.index('TUMOR_READ2')])
            P_val = columns[header_somatic.index('PVAL')]
            callers = columns[header_somatic.index('CALLERS')]
            ref_gene_mut = columns[header_somatic.index('ExonicFunc.refGene')]
            UCSC_gene_mut = columns[header_somatic.index('ExonicFunc.knownGene')]
            ENS_gene_mut = columns[header_somatic.index('ExonicFunc.ensGene')]
            try:
                no_callers = int(callers.strip().split(':')[0])
            except ValueError:
                no_callers = 0
            cov = '{};{},{},{},{},{},{},{}'.format(sample, T_cov, N_cov, T_freq, N_freq, T_reads, P_val, callers)
            has_frame = 'frame' in ''.join([ref_gene_mut, UCSC_gene_mut, ENS_gene_mut])
            has_cov1 = T_cov >= 10 and N_freq < 1.0 and T_freq >= 7 and T_reads >= 4
            has_cov2 = T_cov >= 10 and N_freq >= 1.0 and T_freq >= 7 and T_reads >= 4 and T_freq / N_freq >= 5
            status = (has_frame and (no_callers >= 1 and (has_cov1 or has_cov2))) \
                     or (no_callers >= 2 and (has_cov1 or has_cov2))
            # Store data, coverage and status
            variant_dict[variant_key]['somatic'][sample]['data'] = columns
            variant_dict[variant_key]['somatic'][sample]['status'] = status
            variant_dict[variant_key]['somatic'][sample]['coverage'] = cov
        somatic_nonsyn.close()
        
    print('Loading germline variants..')
    for file in germline if germline else []:
        germline_nonsyn = open(file)
        germline_nonsyn_lines = germline_nonsyn.readlines()
        header_germline = germline_nonsyn_lines.pop(0).strip().split('\t')
        for line in germline_nonsyn_lines:
            columns = line.strip().split('\t')
            variant_key = columns[header_germline.index('VARIANT-KEY')]
            sample = columns[header_germline.index('SAMPLE_ID')]
            if variant_key not in variant_dict:
                variant_dict[variant_key] = {}
            if 'germline' not in variant_dict[variant_key]:
                variant_dict[variant_key]['germline'] = {}
            if sample not in variant_dict[variant_key]['germline']:
                variant_dict[variant_key]['germline'][sample] = {}
            P_val = columns[header_germline.index('PVAL')]
            callers = columns[header_germline.index('CALLERS')]
            try:
                no_callers = int(callers.strip().split(':')[0])
            except ValueError:
                no_callers = 0
            # Compute coverage and pass/fail
            r1 = int(columns[header_germline.index('TUMOR_READ1')])
            r2 = int(columns[header_germline.index('TUMOR_READ2')])
            rfreq = float(columns[header_germline.index('TVAF')].replace('%', ''))
            rcov = int(columns[header_germline.index('TCOV')])
            cov = '{};{},{},{},{},{},{}'.format(sample, r1, r2, rfreq, rcov, P_val, callers)
            # Storage coverage, data and status
            status = rfreq >= 5 and rcov >= 5 and no_callers >= 2
            variant_dict[variant_key]['germline'][sample]['data'] = columns[0:]
            variant_dict[variant_key]['germline'][sample]['status'] = status
            variant_dict[variant_key]['germline'][sample]['coverage'] = cov
        germline_nonsyn.close()

    print('Loading Epitopes..')
    transcript_dict = {}
    # We assume that the same variant in exactly the same position has the same annotation 
    for file in epitopes:
        epitopes = open(file)
        epitopes_lines = epitopes.readlines()
        header_epitopes = epitopes_lines.pop(0).strip().split('\t')
        for line in epitopes_lines:
            columns = line.strip().split('\t')
            key = columns[header_epitopes.index('VARIANT-KEY')]
            sample = columns[header_epitopes.index('SAMPLE_ID')]
            # Load variant data variant -> mut -> transcript
            if key in variant_dict:
                transcript = columns[header_epitopes.index('Transcript_ID')]
                mut_ep = columns[header_epitopes.index('MUT25MER')]
                function = columns[header_epitopes.index('func_ref_gene')]
                cDNA = columns[header_epitopes.index('NT_CHANGE')]
                AA = columns[header_epitopes.index('AA_CHANGE')]
                cDNAposition = ''.join([s for s in cDNA if s.isdigit()])
                if 'Epitopes' not in variant_dict[key]:
                    variant_dict[key]['Epitopes'] = {}
                if mut_ep not in variant_dict[key]['Epitopes']:
                    variant_dict[key]['Epitopes'][mut_ep] = {}
                variant_dict[key]['Epitopes'][mut_ep][transcript] = columns
                # Load transcript data transcript -> variant
                if transcript not in transcript_dict:
                    transcript_dict[transcript] = {}
                    transcript_dict[transcript]['cDNA'] = []
                    transcript_dict[transcript]['aa'] = []
                    transcript_dict[transcript]['position'] = []
                    transcript_dict[transcript]['mutation'] = []
                    transcript_dict[transcript]['variant_key'] = []
                    transcript_dict[transcript]['status'] = []
                transcript_dict[transcript]['cDNA'].append(cDNA)
                transcript_dict[transcript]['aa'].append(AA)
                transcript_dict[transcript]['position'].append(int(cDNAposition))
                transcript_dict[transcript]['mutation'].append(function)
                transcript_dict[transcript]['variant_key'].append(key)
                status_somatic = False
                if 'somatic' in variant_dict[key] and sample in variant_dict[key]['somatic']:
                    status_somatic = variant_dict[key]['somatic'][sample]['status']
                status_germline = False
                if 'germline' in variant_dict[key] and sample in variant_dict[key]['germline']:
                    status_germline = variant_dict[key]['germline'][sample]['status']
                transcript_dict[transcript]['status'].append(status_somatic or status_germline)
        epitopes.close()

    print('Loading Gene counts..')
    counts_dict = {}
    counts_dict_sample = defaultdict(list)
    for file in rna_counts if rna_counts else []:
        counts_file = open(file)
>>>>>>> 1cb5c5f1bb4f22410b0d6032e8e907ad9fc371f9
        counts_file_lines = counts_file.readlines()
        header_cmd = counts_file_lines.pop(0)
        header_counts = counts_file_lines.pop(0).strip().split('\t')
        for line in counts_file_lines:
            columns = line.strip().split('\t')
            gene_id = columns[header_counts.index('Geneid')]
            value = float(columns[-1])
            counts_dict[gene_id] = value
            counts_file.close()
        if len(counts_dict) > 0:
            # Compute mean expression and percentiles
            gene_mean = statistics.mean(counts_dict.values())
            gene_percentile = np.around(stats.percentileofscore(counts_dict.values(), 3))

    print('Creating final files..')
    header_final = 'Variant key\tSomatic samples (passing)\tNumber of Somatic samples (passing)\t' \
                   'Somatic samples (failing)\tNumber of Somatic samples (failing)\t' \
                   'Germline samples (passing)\tNumber of Germline samples (passing)\t' \
                   'Germline samples (failing)\tGermline of Somatic samples (failing)\tEffects\t' \
                   'cDNA change\tAA change\tError flags\tEpitope creation flags\tWt Epitope\t'\
                   'Mut Epitope\tSomatic Callers(Name:NDP;NAD;NVAF;TDP;TAD;TVAF)\t'\
                   'Germline Callers (Name:TDP;TAD;TVAF)\tGeneCount info (exp, mean, percentile)\n'

    final_file = open('overlap_final.txt', 'w')
    final_file.write(header_final)

    final_file_germline = open('overlap_final_germline_unique.txt', 'w')
    final_file_germline.write(header_final)

    final_file_discarded = open('overlap_final_discarded.txt', 'w')
    final_file_discarded.write(header_final)

    final_file_discarded_germline = open('overlap_final_discarded_germline.txt', 'w')
    final_file_discarded_germline.write(header_final)

    for key, value in variant_dict.items():

        # KEY, EFFECTs, mut_dna, mut_aa, flags, wtmer, mutmer, callers, pass, sample

        germlime_pass = [x[-1] for x in value['germline'] if x[-2]]
        germlime_fail = [x[-1] for x in value['germline'] if not x[-2]]
        num_germline_pass = len(germlime_pass)
        num_gemline_fail = len(germlime_fail)
        has_germline = num_germline_pass > 0 or num_gemline_fail > 0
        germline_callers = value['gemline'][7] if has_germline else '-'

        somatic_pass = [x[-1] for x in value['somatic'] if x[-2]]
        somatic_fail = [x[-1] for x in value['somatic'] if not x[-2]]
        num_somatic_pass = len(somatic_pass)
        num_somatic_fail = len(somatic_fail)
        has_somatic = num_somatic_pass > 0 or num_somatic_pass > 0
        somatic_callers = value['somatic'][7] if has_somatic else '-'

        variant = value['somatic'][0] if num_somatic_pass > 0 else value['germline'][0]

        # TODO make creation_flags work
        creation_flags = "-"

        gene_locus = "-"
        for effect in variant[1]:
            if 'ensGene' in effect:
                gene_count = counts_dict[effect.split("_")[1]]
                gene_locus = ';'.join([gene_count, gene_mean, gene_percentile])
                break

        to_write = '\t'.join(str(x) for x in [key,
                                              ';'.join(somatic_pass), num_somatic_pass,
                                              ';'.join(somatic_fail), num_somatic_fail,
                                              ';'.join(germlime_pass), num_germline_pass,
                                              ';'.join(germlime_fail), num_gemline_fail,
                                              variant[1], variant[2], variant[3], variant[4],
                                              creation_flags, variant[5], variant[6],
                                              somatic_callers, germline_callers, gene_locus])
        if num_somatic_pass >= 1:
            final_file.write(to_write + '\n')
        elif num_germline_pass >= 1:
            final_file_germline.write(to_write + '\n')
        elif has_somatic:
            final_file_discarded.write(to_write + '\n')
        else:
            final_file_discarded_germline.write(to_write + '\n')

    final_file.close()
    final_file_germline.close()
    final_file_discarded.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Script that merges variants to create a final report using the\n'
                                                 'results of the somatic and/or germline pipelines that includes epitopes\n'
                                                 'and other useful information\n'
                                                 'Created by Jose Fernandez <jc.fernandes.navarro@gmail.com>',
                                     prog='merge_results.py',
                                     usage='merge_results.py [options]\n'
                                           '--somatic [somatic variants results files]\n'
                                           '--germine [germine variants results files]\n'
                                           '--counts [rna gene counts results]')
    parser.add_argument('--somatic', nargs='+', default=None, required=False,
                        help='List of annotated vcf files with the variants obtained with the somatic pipeline')
    parser.add_argument('--somatic-names', nargs='+', default=None, required=False,
                        help='List of names for each somatic file (to include in the report)')
    parser.add_argument('--germline', nargs='+', default=None, required=False,
                        help='List of annotated vcf files with the variants obtained with the germline pipeline')
    parser.add_argument('--germline-names', nargs='+', default=None, required=False,
                        help='List of names for each germline file (to include in the report)')
    parser.add_argument('--counts', type=str, default=None, required=False,
                        help='File with gene counts from FeatureCounts')
    parser.add_argument('--dictAA',
                        help='Path to a dictionary of transcript IDs to peptide sequences', required=True)
    parser.add_argument('--dictcDNA',
                        help='Path to a dictionary of transcript IDs to DNA sequences', required=True)

    args = parser.parse_args()
    main(args.somatic, args.somatic_names, args.germline, args.germline_names, args.counts,
         os.path.abspath(args.dictcDNA), os.path.abspath(args.dictAA))
