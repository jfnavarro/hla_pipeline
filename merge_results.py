#! /usr/bin/env python
"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
import statistics
import argparse
import numpy as np
import math
from scipy import stats
from _collections import defaultdict
import sys

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

# TODO the whole approach would be easier by creating a Variant class
# with 2 lists of epitopes (somatic and germline) and a dict of Gene counts
def overlap_analysis(somatic, epitopes, germline, rna_counts):

    if not somatic and not germline:
        sys.stderr.write("Error, no variants given as input (somatic or germline).\n")
        sys.exit(1)
        
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
        counts_file_lines = counts_file.readlines()
        header_counts = counts_file_lines.pop(0).strip().split('\t')
        for line in counts_file_lines:
            columns = line.strip().split('\t')
            gene_id = columns[header_counts.index('Geneid')]
            locus = '{}:{}-{}'.format(columns[header_counts.index('Chr')],
                                      columns[header_counts.index('Start')],
                                      columns[header_counts.index('End')])
            value = float(columns[-1])
            sample = columns[header_counts.index('SAMPLE_ID')]
            # NOTE check for precision errors here
            if value != 0:
                if gene_id not in counts_dict:
                    counts_dict[gene_id] = {}
                if sample not in counts_dict[gene_id]:
                    counts_dict[gene_id][sample] = {}
                counts_dict[gene_id][sample]['locus'] = locus
                counts_dict[gene_id][sample]['expression'] = value
                counts_dict_sample[sample].append(value)
        counts_file.close()
    # Compute mean expression and percentiles
    counts_mean_dict = {}
    counts_quan_dict = {}
    # TODO this is slow and ugly, improve!
    for gene, samples in counts_dict.items():
        counts_mean_dict[gene] = statistics.mean([x['expression'] for x in samples.values()])
        counts_quan_dict[gene] = ['{}:{}'.format(sample,
                                                 np.around(
                                                     stats.percentileofscore(counts_dict_sample[sample], x['expression']), 3))
                                  for sample, x in samples.items()]

    print('Creating final files..')
    header_final = 'Variant key\tSomatic samples (passing)\tNumber of Somatic samples (passing)\t' \
                   'Germline samples (passing)\tNumber of Germline samples (passing)\t'\
                   'Somatic samples (failing)\tNumber of Somatic samples (failing)\t'\
                   'Germline samples (failing)\tNumber of Germline samples (failing)\t'\
                   'RefGene name\tRefGene mutation\tRefGene AA change\tUCSC gene name\tUCSC mutation\t'\
                   'UCSC AA change\tEnsembl Gene name\tEnsembl mutation\tEnsembl AA change\t'\
                   '1000genome all freq\tdbSNP_ID\tCosmic Info\tGene Name\ttranscript ID\tMutation type\t'\
                   'Exon\tcDNA change\tAA change\tAA position\tEpitope creation flags\tWt Epitope\t'\
                   'Mut Epitope\tSomatic Coverage info (Sample,Tumor coverage,Normal Coverage,Tumor var freq,'\
                   'Normal var freq,Tumor variant reads,p_value(varscan),callers)\tError flags\t'\
                   'Germline Coverage info (Sample,read1,read2,variant frequency,coverage,p_value(varscan),callers)\t' \
                   'GeneCounts info per sample (locus,exp)\tGeneCounts mean(all samples)\tGeneCounts percentile (all samples)\n'

    final_file = open('overlap_final.txt', 'w')
    final_file.write(header_final)

    final_file_germline = open('overlap_final_germline_unique.txt', 'w')
    final_file_germline.write(header_final)

    final_file_discarded = open('overlap_final_discarded.txt', 'w')
    final_file_discarded.write(header_final)

    for key, value in variant_dict.items():
        germline_cov = ["-;-,-,-,-,-,-,-"]
        germline_samples_pass = ["-"]
        germline_samples_fail = ["-"]
        has_germline = False
        num_germline_pass = 0
        num_germline_fail = 0
        if 'germline' in value:
            germline_cov = '|'.join(x['coverage'] for x in value['germline'].values())
            germline_samples_pass = [key for key, value in value['germline'].items() if value['status']]
            germline_samples_fail = [key for key, value in value['germline'].items() if not value['status']]
            num_germline_pass = len(germline_samples_pass)
            num_germline_fail = len(germline_samples_fail)
            has_germline = True
            
        somatic_cov = ["-;-,-,-,-,-,-,-"]
        somatic_samples_pass = ["-"]
        somatic_samples_fail = ["-"]
        has_somatic = False
        num_somatic_pass = 0
        num_somatic_fail = 0
        if 'somatic' in value:
            somatic_cov = '|'.join(x['coverage'] for x in value['somatic'].values())
            somatic_samples_pass = [key for key, value in value['somatic'].items() if value['status']]
            somatic_samples_fail = [key for key, value in value['somatic'].items() if not value['status']]
            num_somatic_pass = len(somatic_samples_pass)
            num_somatic_fail = len(somatic_samples_fail)
            has_somatic = True
            
        if not has_somatic and not has_germline:
            print("Variant {} is only present in the epitopes file!".format(key))
            continue
            
        if 'Epitopes' in value:
            # Loop through epitopes for this position and write a line for each individual mut25mer
            for mer in value['Epitopes'].values():
                for transcript in sorted(mer.values(), reverse=True):
                    sampleID = transcript[header_epitopes.index('SAMPLE_ID')]
                    # TODO very ugly hack to distinguish somatic and germline epitopes from the same variant (FIX THIS!!!)
                    if has_somatic and sampleID in value['somatic'] and len(value['somatic'][sampleID]['data']) == 45:
                        data = value['somatic'][sampleID]['data']
                        header = header_somatic
                    elif has_germline and sampleID in value['germline'] and len(value['germline'][sampleID]['data']) == 33:
                        data = value['germline'][sampleID]['data']
                        header = header_germline
                    else:
                        # this should never happen
                        print("Epitope {} was not found in neither somatic or germline variants".format(key))
                        continue
                    # NOTE older versions of the pipeline may produce records with missing values
                    if len(transcript) != len(header_epitopes):
                        print("Epitope {} with incorrect number of columns".format(key))
                        continue
                    # NOTE older versions of the pipeline may produce records with missing values
                    if len(data) != len(header):
                        print("Data in epitope {} with incorrect number of columns".format(key))
                        continue
                    ref_gene_name = data[header.index('Gene.refGene')]
                    ref_gene_mut = data[header.index('ExonicFunc.refGene')]
                    ref_gene_change = data[header.index('AAChange.refGene')]
                    UCSC_gene_name = data[header.index('Gene.knownGene')]
                    UCSC_gene_mut = data[header.index('ExonicFunc.knownGene')]
                    UCSC_gene_change = data[header.index('AAChange.knownGene')]
                    ENS_gene_name = data[header.index('Gene.ensGene')]
                    ENS_gene_mut = data[header.index('ExonicFunc.ensGene')]
                    ENS_gene_change = data[header.index('AAChange.ensGene')]
                    genome_all = data[header.index('ALL.sites.2015_08')]
                    dbSNP_ID = data[header.index('avsnp150')]
                    cosmic = data[header.index('COSMIC70')]
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
                    if ENS_gene_name in counts_dict:
                        counts_info = '|'.join(
                            ['{},{}'.format(x['locus'], x['expression']) for x in counts_dict[ENS_gene_name].values()])
                        counts_mean = counts_mean_dict[ENS_gene_name]
                        percentile = ','.join(counts_quan_dict[ENS_gene_name])
                    else:
                        counts_info = 'NA'
                        counts_mean = 'NA'
                        percentile = 'NA'
                    flags = add_flags(transcript_name, key, transcript_dict, mer_len=25)
                    to_write = '\t'.join(str(x) for x in [key,
                                                          ','.join(somatic_samples_pass), num_somatic_pass,
                                                          ','.join(germline_samples_pass), num_germline_pass,
                                                          ','.join(somatic_samples_fail), num_somatic_fail,
                                                          ','.join(germline_samples_fail), num_germline_fail,
                                                          ref_gene_name, ref_gene_mut, ref_gene_change,
                                                          UCSC_gene_name, UCSC_gene_mut, UCSC_gene_change,
                                                          ENS_gene_name, ENS_gene_mut, ENS_gene_change, genome_all,
                                                          dbSNP_ID, cosmic, gene_name, transcript_name, mutation_type,
                                                          exon, cdna_change, aa_change, aa_position, error_flags, wt_mer,
                                                          mu_mer, somatic_cov, flags, germline_cov, counts_info, counts_mean, percentile])
                    if num_somatic_pass >= 1:
                        final_file.write(to_write + '\n')
                    elif num_germline_pass >= 1:
                        final_file_germline.write(to_write + '\n')
                    else:
                        final_file_discarded.write(to_write + '\n')

    final_file.close()
    final_file_germline.close()
    final_file_discarded.close()

parser = argparse.ArgumentParser(description='Script that merges variants and epitopes to create a final report using the\n'\
                                             'results of the somatic and/or germline pipelines\n'
                                             'Created by Jose Fernandez <jc.fernandes.navarro@gmail.com>',
                                 prog='merge_results.py',
                                 usage='merge_results.py [options]\n'
                                       '--somatic [somatic variants results files]\n'
                                       '--epitope [epitopes results files]\n'
                                       '--germine [germine variants results files]\n'
                                       '--counts [rna gene counts results]')
parser.add_argument('--somatic', nargs='+', default=None, required=False,
                    help='List of files with the variants obtained with the somatic pipeline')
parser.add_argument('--epitope', nargs='+', default=None, required=True,
                    help='List of files with the the epitotes (somatic and/or germline)')
parser.add_argument('--germline', nargs='+', default=None, required=False,
                    help='List of files with the variants obtained with the germline pipeline')
parser.add_argument('--counts', nargs='+', default=None, required=False,
                    help='List of files with the gene counts results of the samples (if any)')

args = parser.parse_args()
overlap_analysis(args.somatic, args.epitope, args.germline, args.counts)
