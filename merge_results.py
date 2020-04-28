#! /usr/bin/env python
import numpy as np
from re import sub
import argparse
import re

def overlap_analysis(exome_variants, exome_epitopes, rna_variants, rna_fpkm):

    print('Loading exome variants..')
    variant_dict = {}
    for file in exome_variants:
        exome_nonsyn = open(file)
        exome_nonsyn_lines = exome_nonsyn.readlines()
        header_exome = exome_nonsyn_lines.pop(0).strip().split('\t')
        for line in exome_nonsyn_lines:
            columns = line.strip().split('\t')
            for n,i in enumerate(columns):
                if i in ['', '-', '.', ' ', 'NULL', 'NA']:
                    columns[n] = '0'
            variant_key = columns[header_exome.index('VARIANT-KEY')]
            sample = columns[header_exome.index('SAMPLE_ID')]
            if variant_key not in variant_dict:
                variant_dict[variant_key] = {}
            if 'Exome' not in variant_dict[variant_key]:
                variant_dict[variant_key]['Exome'] = {}
            if sample not in variant_dict[variant_key]['Exome']:
                variant_dict[variant_key]['Exome'][sample] = []
            variant_dict[variant_key]['Exome'][sample].extend(columns[0:])
        exome_nonsyn.close()

    print('Loading RNA variants..')
    for file in rna_variants:
        RNA_nonsyn = open(file)
        RNA_nonsyn_lines = RNA_nonsyn.readlines()
        header_rna = RNA_nonsyn_lines.pop(0).strip().split('\t')
        for line in RNA_nonsyn_lines:
            columns = line.strip().split('\t')
            for n,i in enumerate(columns):
                if i in ['', '-', '.', ' ', 'NULL', 'NA']:
                    columns[n] = '0'
            variant_key = columns[header_rna.index('VARIANT-KEY')]
            sample = columns[header_rna.index('SAMPLE_ID')]
            if variant_key not in variant_dict:
                variant_dict[variant_key] = {}
            if 'RNA' not in variant_dict[variant_key]:
                variant_dict[variant_key]['RNA'] = {}
            if sample not in variant_dict[variant_key]['RNA']:
                variant_dict[variant_key]['RNA'][sample] = []
            variant_dict[variant_key]['RNA'][sample].extend(columns[0:])
        RNA_nonsyn.close()

    print('Loading FPKM data..')
    FPKM_dict = {}
    for file in rna_fpkm:
        FPKM_file = open(file)
        FPKM_file_lines = FPKM_file.readlines()
        header_fpkm = FPKM_file_lines.pop(0).strip().split('\t')
        for line in FPKM_file_lines:
            columns = line.strip().split('\t')
            gene_id = columns[header_fpkm.index('tracking_id')]
            locus = columns[header_fpkm.index('locus')]
            fpkm_value = columns[header_fpkm.index('FPKM')]
            sample = columns[header_fpkm.index('SAMPLE_ID')]
            if gene_id not in FPKM_dict:
                FPKM_dict[gene_id] = {}
            if sample not in FPKM_dict[gene_id]:
                FPKM_dict[gene_id][sample] = {}
            FPKM_dict[gene_id][sample]['locus'] = locus
            FPKM_dict[gene_id][sample]['expression'] = fpkm_value
        FPKM_file.close()

    print('Loading epitopes..')
    # NOTE these two can be used how many times a given variant occurs in its neighbours (+20 nt or +10 aa)
    aa_info = {}
    nt_info = {}
    for file in exome_epitopes:
        epitopes = open(file)
        epitopes_lines = epitopes.readlines()
        header_epitopes = epitopes_lines.pop(0).strip().split('\t')
        for line in epitopes_lines:
            columns = line.strip().split('\t')
            key = columns[header_epitopes.index('VARIANT-KEY')]
            transcript = columns[header_epitopes.index('Transcript_ID')]
            mut_ep = columns[header_epitopes.index('func_ref_gene')]
            chr = columns[header_epitopes.index('CHR')]
            AA = columns[header_epitopes.index('AA_CHANGE')]
            AAposition = ''.join([s for s in AA if s.isdigit()])
            sample = columns[header_epitopes.index('SAMPLE_ID')]
            start = int(columns[header_epitopes.index('START')])
            end = int(columns[header_epitopes.index('STOP')])

            if variant_key not in variant_dict:
                variant_dict[key] = {}
            if 'Epitopes' not in variant_dict[key]:
                variant_dict[key]['Epitopes'] = {}
            if mut_ep not in variant_dict[key]['Epitopes']:
                variant_dict[key]['Epitopes'][mut_ep] = {}
            if transcript not in variant_dict[key]['Epitopes'][mut_ep]:
                variant_dict[key]['Epitopes'][mut_ep][transcript] = []
            variant_dict[key]['Epitopes'][mut_ep][transcript].extend(columns[0:])

            if sample not in aa_info:
                aa_info[sample] = {}
            if transcript not in aa_info[sample]:
                aa_info[sample][transcript] = []
            aa_info[sample][transcript].append([int(AAposition), variant_key])

            if sample not in nt_info:
                nt_info[sample] = {}
            if chr not in nt_info[sample]:
                nt_info[sample][chr] = []
            nt_info[sample][chr].append([start, end, variant_key])

        epitopes.close()

    print('Creating final files..')
    header_final = 'Variant key\tNumber of Exomes passing filters\tNumber of RNA-seq samples seen in\t'\
                   'Number of Exome samples seen in\tExome samples passing filters\t'\
                   'RNA-seq samples\tExome samples seen in(unfiltered)\t'\
                   'RefGene name\tRefGene mutation\tRefGene AA change\tUCSC gene name\tUCSC mutation\t'\
                   'UCSC AA change\tEnsembl Gene name\tEnsembl mutation\tEnsembl AA change\t'\
                   '1000genome all freq\tdbSNP_ID\tCosmic Info\tGene Name\ttranscript ID\tMutation type\t'\
                   'Exon\tcDNA change\tAA change\tAA position\tEpitope creation flags\tWt Epitope\t'\
                   'Mut Epitope\tExome Coverage info (Sample,Tumor coverage,Normal Coverage,Tumor var freq,'\
                   'Normal var freq,Tumor variant reads,p_value (varscan),callers)\t'\
                   'RNA Coverage info (Sample,read1,read2,variant frequency,coverage)\tFPKM info\tFPKM mean\n'
    final_file = open('overlap_final.txt', 'w')
    final_file.write(header_final)
    final_file_discarded = open('overlap_final_discarded.txt', 'w')
    final_file_discarded.write(header_final)

    unique_rna = open('overlap_unique_rna.txt', 'w')
    unique_rna.write('Variant key\tPer sample coverage\tFPKM info\tFPKM mean\n')

    for key,value in variant_dict.items():
        # Exome
        sample_cov = []
        exome_pass = []
        exome_fail = []
        if 'Exome' in value:
            for sample in sorted(value['Exome'].values()):
                sampleID = sample[header_exome.index('SAMPLE_ID')]
                N_cov = sample[header_exome.index('NCOV')]
                T_cov = sample[header_exome.index('TCOV')]
                T_freq = sample[header_exome.index('TVAF')].replace('%','')
                N_freq = sample[header_exome.index('NVAF')].replace('%','')
                T_reads = sample[header_exome.index('TUMOR_READ2')]
                P_val = sample[header_exome.index('PVAL')]
                callers = sample[header_exome.index('CALLERS')]
                ref_gene_mut = sample[header_exome.index('ExonicFunc.refGene')]
                UCSC_gene_mut = sample[header_exome.index('ExonicFunc.knownGene')]
                ENS_gene_mut = sample[header_exome.index('ExonicFunc.ensGene')]
                try:
                    no_callers = int(callers.strip().split(':')[0])
                except ValueError:
                    no_callers = 0
                cov = '{};{},{},{},{},{},{},{}'.format(sampleID, T_cov, N_cov, T_freq, N_freq, T_reads, P_val, callers)
                sample_cov.append(cov)
                #NOTE coverage filter was already applied when creating the variants
                if ((re.search(r'frame', ref_gene_mut) or re.search(r'frame', UCSC_gene_mut) \
                    or re.search(r'frame', ENS_gene_mut)) and no_callers >= 1) or no_callers >= 2:
                        exome_pass.append(sampleID)
                else:
                    exome_fail.append(sampleID)
        else:
            print("Variant {} is not present in the exome data".format(key))

        # RNA
        sample_covR = []
        if 'RNA' in value:
            for sample in sorted(value['RNA'].values()):
                sampleID = sample[header_rna.index('SAMPLE_ID')]
                r1 = float(sample[header_rna.index('TUMOR_READ1')])
                r2 = float(sample[header_rna.index('TUMOR_READ2')])
                rfreq = sub('%','', sample[header_rna.index('TUMOR_VAR_FREQ')])
                rcov = r1 + r2
                #TODO add a coverage filter
                cov = '{};{},{},{},{}'.format(sampleID, r1, r2, rfreq, rcov)
                sample_covR.append(cov)
        else:
            print("Variant {} is not present in the RNA data".format(key))

        if len(sample_cov) >= 1:
            # Loop through epitopes for this position and write a line for each individual mut25mer
            for mer in sorted(value['Epitopes'].values()):
                # NOTE this assumes that epitoes on the same position in different samples are the same
                for transcript in sorted(mer.values(), reverse=True):
                    sampleID = transcript[header_epitopes.index('SAMPLE_ID')]
                    exome = value['Exome'][sampleID]
                    rna_samples = value['RNA'].keys() if 'RNA' in value else ''
                    num_rna_samples = len(rna_samples)
                    ref_gene_name = exome[header_exome.index('Gene.refGene')]
                    ref_gene_mut = exome[header_exome.index('ExonicFunc.refGene')]
                    ref_gene_change = exome[header_exome.index('AAChange.refGene')]
                    UCSC_gene_name = exome[header_exome.index('Gene.knownGene')]
                    UCSC_gene_mut = exome[header_exome.index('ExonicFunc.knownGene')]
                    UCSC_gene_change = exome[header_exome.index('AAChange.knownGene')]
                    ENS_gene_name = exome[header_exome.index('Gene.ensGene')]
                    ENS_gene_mut = exome[header_exome.index('ExonicFunc.ensGene')]
                    ENS_gene_change = exome[header_exome.index('AAChange.ensGene')]
                    genome_all = exome[header_exome.index('ALL.sites.2015_08')]
                    dbSNP_ID = exome[header_exome.index('snp138NonFlagged')]
                    cosmic = exome[header_exome.index('COSMIC70')]
                    gene_name = transcript[header_epitopes.index('Gene')]
                    transcript_name = transcript[header_epitopes.index('Transcript_ID')]
                    mutation_type = transcript[header_epitopes.index('exonic_func_ref')]
                    exon = transcript[header_epitopes.index('Exon_Numb')]
                    cdna_change = transcript[header_epitopes.index('NT_CHANGE')]
                    aa_change = transcript[header_epitopes.index('AA_CHANGE')]
                    aa_position = transcript[header_epitopes.index('POSITION')]
                    flags = transcript[header_epitopes.index('ERRORS')]
                    wt_mer = transcript[header_epitopes.index('WT25MER')]
                    mu_mer = transcript[header_epitopes.index('MUT25MER')]
                    exome_cov = '|'.join(sample_cov)
                    rna_cov = '|'.join(sample_covR)
                    if ENS_gene_name in FPKM_dict:
                        fpkm_info = '|'.join(
                            ['{},{},{}'.format(s, x['locus'], x['expression']) for s,x in FPKM_dict[ENS_gene_name].items()])
                        fpkm_mean = np.mean(float(x['expression']) for x in FPKM_dict[ENS_gene_name].values())
                    else:
                        fpkm_info = 'NA'
                        fpkm_mean = 'NA'
                    to_write = '\t'.join(str(x) for x in [key, len(exome_pass), num_rna_samples, len(exome_pass) + len(exome_fail),
                                                          ','.join(exome_pass), ','.join(rna_samples), ','.join(exome_fail),
                                                          ref_gene_name, ref_gene_mut, ref_gene_change,
                                                          UCSC_gene_name, UCSC_gene_mut, UCSC_gene_change,
                                                          ENS_gene_name, ENS_gene_mut, ENS_gene_change, genome_all,
                                                          dbSNP_ID, cosmic, gene_name, transcript_name, mutation_type,
                                                          exon, cdna_change, aa_change, aa_position, flags, wt_mer,
                                                          mu_mer, exome_cov, rna_cov, fpkm_info, fpkm_mean])
                    if len(exome_pass) >= 1:
                        final_file.write(to_write + '\n')
                    else:
                        final_file_discarded.write(to_write + '\n')
        elif len(sample_covR) >= 1:
            print('Variant {} was only detected in RNA'.format(key))
            # NOTE assuming same gene for the variant in all samples
            ENS_gene_name = list(value['RNA'].values())[0][header_rna.index('Gene')]
            # TODO add more fields to file
            rna_cov = '|'.join(sample_covR)
            if ENS_gene_name in FPKM_dict:
                fpkm_info = '|'.join(
                    ['{},{},{}'.format(s, x['locus'], x['expression']) for s, x in FPKM_dict[ENS_gene_name].items()])
                fpkm_mean = np.mean(float(x['expression']) for x in FPKM_dict[ENS_gene_name].values())
            else:
                fpkm_info = 'NA'
                fpkm_mean = 'NA'
            to_write = '\t'.join(str(x) for x in [key, rna_cov, fpkm_info, fpkm_mean])
            unique_rna.write(to_write + '\n')


    final_file.close()
    final_file_discarded.close()
    unique_rna.close()

parser = argparse.ArgumentParser(description='Script to aggregate results from Jareds pipeline '
                                             '(adjusted by Jose Fernandez) <jc.fernandes.navarro@gmail.com>',
                                 prog='merge_results.py',
                                 usage='merge_results.py [options] --exome [exome results files] '
                                       '--epitote [epitote results files] --rna [rna results files] --fpkm [rna fpkm results]')

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