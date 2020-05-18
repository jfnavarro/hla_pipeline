#! /usr/bin/env python
import statistics
from re import sub
import argparse
import re
import numpy as np
import math
from scipy import stats
from _collections import defaultdict

def add_flags(transcript, variant_key, transcript_info, mer_len=25):
    cDNA_flanks = (math.floor(mer_len / 2)) * 3
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
                elif mut == 'stopgain' and not ss_seen is None:
                    ss_seen += ', {} filters upstream stopgain ({})'.format(status_str, aa)
                if re.search('frame', mut) and ss_seen is None:
                    ss_seen = '{} filters upstream fs ({})'.format(status_str, aa)
                elif re.search('frame', mut) and  not ss_seen is None:
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

def overlap_analysis(exome_variants, exome_epitopes, rna_variants, rna_fpkm):

    variant_dict = {}

    print('Loading exome variants..')
    for file in exome_variants:
        exome_nonsyn = open(file)
        exome_nonsyn_lines = exome_nonsyn.readlines()
        header_exome = exome_nonsyn_lines.pop(0).strip().split('\t')
        for line in exome_nonsyn_lines:
            columns = line.strip().split('\t')
            variant_key = columns[header_exome.index('VARIANT-KEY')]
            sample = columns[header_exome.index('SAMPLE_ID')]
            if variant_key not in variant_dict:
                variant_dict[variant_key] = {}
                variant_dict[variant_key]['Exome'] = {}
            if sample not in variant_dict[variant_key]['Exome']:
                variant_dict[variant_key]['Exome'][sample] = {}
            # Compute coverage and pass/fail
            N_cov = int(columns[header_exome.index('NCOV')])
            T_cov = float(columns[header_exome.index('TCOV')])
            T_freq = float(columns[header_exome.index('TVAF')].replace('%',''))
            N_freq = float(columns[header_exome.index('NVAF')].replace('%',''))
            T_reads = float(columns[header_exome.index('TUMOR_READ2')])
            P_val = columns[header_exome.index('PVAL')]
            callers = columns[header_exome.index('CALLERS')]
            ref_gene_mut = columns[header_exome.index('ExonicFunc.refGene')]
            UCSC_gene_mut = columns[header_exome.index('ExonicFunc.knownGene')]
            ENS_gene_mut = columns[header_exome.index('ExonicFunc.ensGene')]
            try:
                no_callers = int(callers.strip().split(':')[0])
            except ValueError:
                no_callers = 0
            cov = '{};{},{},{},{},{},{},{}'.format(sample, T_cov, N_cov, T_freq, N_freq, T_reads, P_val, callers)
            status = (((re.search(r'frame', ref_gene_mut) or re.search(r'frame', UCSC_gene_mut)
                        or re.search(r'frame', ENS_gene_mut)) and no_callers >= 1) or no_callers >= 2)\
                     and (N_cov > 10 and T_freq >= 7 and T_reads >= 4)
            # Store data, coverage and status
            variant_dict[variant_key]['Exome'][sample]['data'] = columns
            variant_dict[variant_key]['Exome'][sample]['status'] = status
            variant_dict[variant_key]['Exome'][sample]['coverage'] = cov
        exome_nonsyn.close()

    print('Loading epitopes..')
    transcript_dict = {}
    for file in exome_epitopes:
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
                mut_ep = columns[header_epitopes.index('func_ref_gene')]
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
                transcript_dict[transcript]['status'].append(variant_dict[key]['Exome'][sample]['status'])
        epitopes.close()

    print('Loading RNA variants..')
    for file in rna_variants:
        RNA_nonsyn = open(file)
        RNA_nonsyn_lines = RNA_nonsyn.readlines()
        header_rna = RNA_nonsyn_lines.pop(0).strip().split('\t')
        for line in RNA_nonsyn_lines:
            columns = line.strip().split('\t')
            variant_key = columns[header_rna.index('VARIANT-KEY')]
            sample = columns[header_rna.index('SAMPLE_ID')]
            if variant_key not in variant_dict:
                variant_dict[variant_key] = {}
            if 'RNA' not in variant_dict[variant_key]:
                variant_dict[variant_key]['RNA'] = {}
            if sample not in variant_dict[variant_key]['RNA']:
                variant_dict[variant_key]['RNA'][sample] = {}
            # Compute coverage and pass/fail
            r1 = float(columns[header_rna.index('TUMOR_READ1')])
            r2 = float(columns[header_rna.index('TUMOR_READ2')])
            rfreq = float(sub('%', '', columns[header_rna.index('TUMOR_VAR_FREQ')]))
            rcov = r1 + r2
            cov = '{};{},{},{},{}'.format(sample, r1, r2, rfreq, rcov)
            # Storage coverage, data and status
            variant_dict[variant_key]['RNA'][sample]['data'] = columns[0:]
            variant_dict[variant_key]['RNA'][sample]['status'] = rfreq >=5 and rcov >= 5
            variant_dict[variant_key]['RNA'][sample]['coverage'] = cov
        RNA_nonsyn.close()

    print('Loading FPKM data..')
    FPKM_dict = {}
    FPKM_dict_sample = defaultdict(list)
    for file in rna_fpkm:
        FPKM_file = open(file)
        FPKM_file_lines = FPKM_file.readlines()
        header_fpkm = FPKM_file_lines.pop(0).strip().split('\t')
        for line in FPKM_file_lines:
            columns = line.strip().split('\t')
            gene_id = columns[header_fpkm.index('gene_id')]
            locus = columns[header_fpkm.index('locus')]
            fpkm_value = float(columns[header_fpkm.index('FPKM')])
            sample = columns[header_fpkm.index('SAMPLE_ID')]
            #NOTE check for precision errors here
            if fpkm_value != 0:
                if gene_id not in FPKM_dict:
                    FPKM_dict[gene_id] = {}
                if sample not in FPKM_dict[gene_id]:
                    FPKM_dict[gene_id][sample] = {}
                FPKM_dict[gene_id][sample]['locus'] = locus
                FPKM_dict[gene_id][sample]['expression'] = fpkm_value
                FPKM_dict_sample[sample].append(fpkm_value)
        FPKM_file.close()
    # Compute mean expression and percentiles
    FPKM_mean_dict = {}
    FPKM_quan_dict = {}
    #TODO this is slow and ugly, improve!
    for gene, samples in FPKM_dict.items():
        FPKM_mean_dict[gene] = statistics.mean([x['expression'] for x in samples.values()])
        FPKM_quan_dict[gene] = ['{}:{}'.format(sample,
                                               np.around(
                                                   stats.percentileofscore(FPKM_dict_sample[sample], x['expression']), 3))
                                for sample, x in samples.items()]

    print('Creating final files..')
    header_final = 'Variant key\tExomes samples (passing)\tNumber of Exomes samples (passing)\t' \
                   'RNA-seq samples (passing)\tNumber of RNA-seq samples (passing)\t'\
                   'Exomes samples (failing)\tNumber of Exomes samples (failing)\t'\
                   'RNA-seq samples (failing)\tNumber of RNA-seq samples (failing)\t'\
                   'RefGene name\tRefGene mutation\tRefGene AA change\tUCSC gene name\tUCSC mutation\t'\
                   'UCSC AA change\tEnsembl Gene name\tEnsembl mutation\tEnsembl AA change\t'\
                   '1000genome all freq\tdbSNP_ID\tCosmic Info\tGene Name\ttranscript ID\tMutation type\t'\
                   'Exon\tcDNA change\tAA change\tAA position\tEpitope creation flags\tWt Epitope\t'\
                   'Mut Epitope\tExome Coverage info (Sample,Tumor coverage,Normal Coverage,Tumor var freq,'\
                   'Normal var freq,Tumor variant reads,p_value (varscan),callers)\tError flags\t'\
                   'RNA Coverage info (Sample,read1,read2,variant frequency,coverage)\t' \
                   'FPKM info per sample (locus,exp)\tFPKM mean(all samples)\tFPKM percentile (all samples)\n'
    final_file = open('overlap_final.txt', 'w')
    final_file.write(header_final)
    final_file_discarded = open('overlap_final_discarded.txt', 'w')
    final_file_discarded.write(header_final)

    unique_rna = open('overlap_unique_rna.txt', 'w')
    unique_rna.write('Variant key\tPer sample coverage (Sample,read1,read2,variant frequency,coverage)\tGene\t'\
                     'RNA-seq samples (passing)\tNumber of RNA-seq samples (passing)\t'\
                     'RNA-seq samples (failing)\tNumber of RNA-seq samples (failing)\t'\
                     'Mutation type\tFPKM info pers sample (locus,exp)\tFPKM mean(all samples)\tFPKM percentile (all samples)\n')

    for key,value in variant_dict.items():

        rna_cov = ["-;-,-,-,-"]
        rna_samples_pass = []
        rna_samples_fail = []
        if 'RNA' in value:
            rna_cov = '|'.join(x['coverage'] for x in value['RNA'].values())
            rna_samples_pass = [key for key,value in value['RNA'].items() if value['status']]
            rna_samples_fail = [key for key,value in value['RNA'].items() if not value['status']]

        if 'Epitopes' in value and 'Exome' in value:
            # Loop through epitopes for this position and write a line for each individual mut25mer
            for mer in sorted(value['Epitopes'].values()):
                for transcript in sorted(mer.values(), reverse=True):
                    sampleID = transcript[header_epitopes.index('SAMPLE_ID')]
                    exome = value['Exome'][sampleID]
                    exome_samples_pass = [key for key,value in value['Exome'].items() if value['status']]
                    exome_samples_fail = [key for key,value in value['Exome'].items() if not value['status']]
                    exome_cov = '|'.join(x['coverage'] for x in value['Exome'].values())
                    ref_gene_name = exome['data'][header_exome.index('Gene.refGene')]
                    ref_gene_mut = exome['data'][header_exome.index('ExonicFunc.refGene')]
                    ref_gene_change = exome['data'][header_exome.index('AAChange.refGene')]
                    UCSC_gene_name = exome['data'][header_exome.index('Gene.knownGene')]
                    UCSC_gene_mut = exome['data'][header_exome.index('ExonicFunc.knownGene')]
                    UCSC_gene_change = exome['data'][header_exome.index('AAChange.knownGene')]
                    ENS_gene_name = exome['data'][header_exome.index('Gene.ensGene')]
                    ENS_gene_mut = exome['data'][header_exome.index('ExonicFunc.ensGene')]
                    ENS_gene_change = exome['data'][header_exome.index('AAChange.ensGene')]
                    genome_all = exome['data'][header_exome.index('ALL.sites.2015_08')]
                    dbSNP_ID = exome['data'][header_exome.index('snp138NonFlagged')]
                    cosmic = exome['data'][header_exome.index('COSMIC70')]
                    gene_name = transcript[header_epitopes.index('Gene')]
                    transcript_name = transcript[header_epitopes.index('Transcript_ID')]
                    mutation_type = transcript[header_epitopes.index('exonic_func_ref')]
                    exon = transcript[header_epitopes.index('Exon_Numb')]
                    cdna_change = transcript[header_epitopes.index('NT_CHANGE')]
                    aa_change = transcript[header_epitopes.index('AA_CHANGE')]
                    aa_position = transcript[header_epitopes.index('POSITION')]
                    error_flags = transcript[header_epitopes.index('ERRORS')]
                    wt_mer = transcript[header_epitopes.index('WT25MER')]
                    mu_mer = transcript[header_epitopes.index('MUT25MER')]
                    if ENS_gene_name in FPKM_dict:
                        fpkm_info = '|'.join(
                            ['{},{}'.format(x['locus'], x['expression']) for x in FPKM_dict[ENS_gene_name].values()])
                        fpkm_mean = FPKM_mean_dict[ENS_gene_name]
                        percentile = ','.join(FPKM_quan_dict[ENS_gene_name])
                    else:
                        fpkm_info = 'NA'
                        fpkm_mean = 'NA'
                        percentile = 'NA'
                    flags = add_flags(transcript_name, key, transcript_dict, mer_len=25)
                    to_write = '\t'.join(str(x) for x in [key, ','.join(exome_samples_pass), len(exome_samples_pass),
                                                          ','.join(rna_samples_pass), len(rna_samples_pass),
                                                          ','.join(exome_samples_fail), len(exome_samples_fail),
                                                          ','.join(rna_samples_fail), len(rna_samples_fail),
                                                          ref_gene_name, ref_gene_mut, ref_gene_change,
                                                          UCSC_gene_name, UCSC_gene_mut, UCSC_gene_change,
                                                          ENS_gene_name, ENS_gene_mut, ENS_gene_change, genome_all,
                                                          dbSNP_ID, cosmic, gene_name, transcript_name, mutation_type,
                                                          exon, cdna_change, aa_change, aa_position, error_flags, wt_mer,
                                                          mu_mer, exome_cov, flags, rna_cov, fpkm_info, fpkm_mean, percentile])
                    if len(exome_samples_pass) >= 1:
                        final_file.write(to_write + '\n')
                    else:
                        final_file_discarded.write(to_write + '\n')
        elif 'RNA' in value:
            # Assuming the gene name is the same in every sample where this variant was detected
            ENS_gene_name = list(value['RNA'].values())[0]['data'][header_rna.index('Gene')]
            mutation_type = list(value['RNA'].values())[0]['data'][header_rna.index('exonic_func')]
            if ENS_gene_name in FPKM_dict:
                fpkm_info = '|'.join(
                    ['{},{}'.format(x['locus'], x['expression']) for x in FPKM_dict[ENS_gene_name].values()])
                fpkm_mean = FPKM_mean_dict[ENS_gene_name]
                percentile = ','.join(FPKM_quan_dict[ENS_gene_name])
            else:
                fpkm_info = 'NA'
                fpkm_mean = 'NA'
                percentile = 'NA'
            # TODO add more fields to file
            to_write = '\t'.join(str(x) for x in [key, rna_cov, ENS_gene_name,
                                                  ','.join(rna_samples_pass), len(rna_samples_pass),
                                                  ','.join(rna_samples_fail), len(rna_samples_fail),
                                                  mutation_type, fpkm_info, fpkm_mean, percentile])
            unique_rna.write(to_write + '\n')
        else:
            print("Variant {} was only detected in either epitopes or exome, strange!".format(key))

    final_file.close()
    final_file_discarded.close()
    unique_rna.close()

parser = argparse.ArgumentParser(description='Script to aggregate results and create a final report from the exome and rna-seq pipelines results '
                                             '(created by Jose Fernandez <jc.fernandes.navarro@gmail.com>)',
                                 prog='merge_results.py',
                                 usage='merge_results.py [options] '
                                       '--exome [exome results files] '
                                       '--epitote [epitote results files] '
                                       '--rna [rna results files] '
                                       '--fpkm [rna fpkm results]')

parser.add_argument('--exome', nargs='+', default=None, required=True,
                    help='List of files with the results of the exome variants')
parser.add_argument('--epitote', nargs='+', default=None, required=True,
                    help='List of files with the results of the epitotes')
parser.add_argument('--rna', nargs='+', default=None, required=True,
                    help='List of files with the results of the RNA variants')
parser.add_argument('--fpkm', nargs='+', default=None, required=True,
                    help='List of files with the results of the RNA FPKM')

args = parser.parse_args()
overlap_analysis(args.exome, args.epitote, args.rna, args.fpkm)