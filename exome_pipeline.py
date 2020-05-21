#! /usr/bin/env python
"""
@author: jfnavarro
"""
from re import sub
from hlapipeline.common import *
from hlapipeline.tools import *
from hlapipeline.filters import *
import shutil
import re
import argparse
import multiprocessing
import os
import pandas as pd
from Bio import SeqIO

def HLA_predictionDNA(bamfile, sampleID, outfile, threads):
    OUT_DIR = os.path.abspath(os.path.join('out_hla_', os.path.splitext(outfile)[0]))
    os.makedirs(OUT_DIR, exist_ok=True)
    cmd = '{} --BAM {} --workingDir {} --graph {} --sampleID {}'\
          ' --maxTHREADS {}'.format(HLALA, bamfile, OUT_DIR, 'PRG_MHC_GRCh38_withIMGT', sampleID, threads)
    exec_command(cmd)

    # create a dictionary to store the output for each allele
    hla = pd.read_table(os.path.join(OUT_DIR, sampleID, 'hla', 'R1_bestguess_G.txt'), sep='\t')
    allele_dict = {}
    hla = hla.groupby('Locus')
    for k, g in hla:
        allele = [re.sub('G$', '', a).replace('N', '') for a in g['Allele'].tolist()]
        allele_dict['HLA_' + k] = allele

    # Create formatted output file
    a = open(outfile, 'w')
    for key, value in allele_dict.items():
        a.write('{}\t{}\t{}\n'.format(sampleID, key, '\t'.join(value)))
    a.close()
    
def final_variants(input, output, output_other, vcf_cov_dict, sampleID, tumor_type, header=True):
    nonsyn_snv = open(input)
    nonsyn_snv_lines = nonsyn_snv.readlines()[1:]
    nonsyn_file = open(output, 'w' if header else 'a')
    all_file = open(output_other, 'w' if header else 'a')
    if header:
        header = 'SAMPLE_ID\tCHR\tSTART\tEND\tREF\tALT\tavsnp150\tTUMOR_READ1' \
                 '\tTUMOR_READ2\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tFunc.knownGene\tGene.knownGene' \
                 '\tExonicFunc.knownGene\tAAChange.knownGene\tFunc.ensGene\tGene.ensGene\tExonicFunc.ensGene\tAAChange.ensGene' \
                 '\tALL.sites.2015_08\tEUR.sites.2015_08\tAMR.sites.2015_08\tEAS.sites.2015_08\tAFR.sites.2015_08\tNORMAL_READ1' \
                 '\tNORMAL_READ2\tTRFOR-DP4\tTRREV-DP4\tTVFOR-DP4\tTVREV-DP4\tNRFOR-DP4\tNRREV-DP4\tNVFOR-DP4\tNVAF\tNVREV-DP4' \
                 '\tTVAF\tPVAL\tCALLERS\tTUMOUR\tTCOV\tNCOV\tVARIANT-KEY\tCOSMIC70\n'
        nonsyn_file.write(header)
        all_file.write(header)
    # TODO use header names instead to access the fields
    for line in nonsyn_snv_lines:
        if line.startswith('#'):
            continue
        columns = line.strip().split('\t')
        Chr = columns[0]
        start = columns[1]
        ID = Chr + ':' + start
        end = columns[2]
        ref = columns[3]
        alt = columns[4]
        ref_gene_detail = columns[7]
        known_gene_detail = columns[12]
        ens_gene_detail = columns[17]
        COSMIC = columns[26]
        variant_key = str(Chr) + ':' + str(start) + '-' + str(end) + ' ' + str(ref) + '>' + str(alt)
        if ref_gene_detail != 'NA':
            columns[9] = ref_gene_detail
        if known_gene_detail != 'NA':
            columns[14] = known_gene_detail
        if ens_gene_detail != 'NA':
            columns[19] = ens_gene_detail
        # Can be missing keys if the annotation and the combined variants do not have the same chromosomes
        try:
            p_val = vcf_cov_dict[ID]['pval']
            callers = vcf_cov_dict[ID]['Note']
            trfor = vcf_cov_dict[ID]['trfor']
            trrev = vcf_cov_dict[ID]['trrev']
            tvfor = vcf_cov_dict[ID]['tvfor']
            tvrev = vcf_cov_dict[ID]['tvrev']
            nrfor = vcf_cov_dict[ID]['nrfor']
            nrrev = vcf_cov_dict[ID]['nrrev']
            nvfor = vcf_cov_dict[ID]['nvfor']
            nvrev = vcf_cov_dict[ID]['nvrev']
            tumor_read1 = vcf_cov_dict[ID]['tumor_read1']
            tumor_read2 = vcf_cov_dict[ID]['tumor_read2']
            normal_read1 = vcf_cov_dict[ID]['normal_read1']
            normal_read2 = vcf_cov_dict[ID]['normal_read2']
            tcov = vcf_cov_dict[ID]['tumor_coverage']
            tfreq = vcf_cov_dict[ID]['tumor_freq']
            ncov = vcf_cov_dict[ID]['normal_coverage']
            nfreq = vcf_cov_dict[ID]['normal_freq']
            to_write = '\t'.join([str(x) for x in [sampleID, '\t'.join(columns[0:5]), columns[20], tumor_read1, tumor_read2,
                                                   '\t'.join(columns[5:7]), '\t'.join(columns[8:12]), '\t'.join(columns[13:17]),
                                                   '\t'.join(columns[18:20]), '\t'.join(columns[21:26]), normal_read1, normal_read2,
                                                   trfor, trrev, tvfor, tvrev, nrfor, nrrev, nvfor, nfreq, nvrev, tfreq, p_val,
                                                   callers, tumor_type, tcov, ncov, variant_key, COSMIC]])
            if (re.search(r'nonsynonymous', columns[8]) or re.search(r'frame', columns[8]) or re.search(r'stop', columns[8]) \
                    or re.search(r'nonsynonymous', columns[13]) or re.search(r'frame', columns[13]) or re.search(r'stop', columns[13]) \
                    or re.search(r'nonsynonymous', columns[18]) or re.search(r'frame', columns[18]) or re.search(r'stop', columns[18])):
                nonsyn_file.write(to_write + '\n')
            else:
                all_file.write(to_write + '\n')
        except KeyError:
            print("Missing variant for {}".format(ID))
    nonsyn_file.close()
    all_file.close()
    nonsyn_snv.close()

def exome_pipeline(R1_NORMAL,
                   R2_NORMAL,
                   R1_CANCER,
                   R2_CANCER,
                   IILLUMINA_ADAPTERS,
                   tumor_type,
                   genome,
                   sampleID,
                   THREADS,
                   FASTA_AA_DICT,
                   FASTA_cDNA_DICT,
                   KNOWN_SITE1,
                   KNOWN_SITE2,
                   SNPSITES,
                   GERMLINE,
                   PON,
                   ANNOVAR_DB,
                   ANNOVAR_VERSION,
                   steps):
    print("Exome pipeline")

    # Sample 1 cancer, sample 2 normal
    sample1_ID = sampleID + "_Tumor"
    sample2_ID = sampleID + "_Normal"

    # Create sub-folder to store all results
    os.makedirs('exome', exist_ok=True)
    os.chdir('exome')

    if 'mapping' in steps:
        # TRIMMING
        print('Starting trimming')
        cmd = '{} PE -threads {} -phred33 {} {} R1_normal.fastq.gz R1_normal_unpaired.fastq.gz ' \
              'R2_normal.fastq.gz R2_normal_unpaired.fastq.gz ' \
              'ILLUMINACLIP:{}:2:40:15 HEADCROP:9 CROP:140 SLIDINGWINDOW:4:25 MINLEN:5'.format(TRIPTOMATIC,
                                                                                               THREADS,
                                                                                               R1_NORMAL,
                                                                                               R2_NORMAL,
                                                                                               IILLUMINA_ADAPTERS)
        exec_command(cmd)

        cmd = '{} PE -threads {} -phred33 {} {} R1_cancer.fastq.gz R1_cancer_unpaired.fastq.gz ' \
              'R2_cancer.fastq.gz R2_cancer_unpaired.fastq.gz ' \
              'ILLUMINACLIP:{}:2:40:15 HEADCROP:9 CROP:140 SLIDINGWINDOW:4:25 MINLEN:5'.format(TRIPTOMATIC,
                                                                                               THREADS,
                                                                                               R1_CANCER,
                                                                                               R2_CANCER,
                                                                                               IILLUMINA_ADAPTERS)
        exec_command(cmd)

        # ALIGNMENT
        print('Starting alignment')
        # Normal (paired)
        cmd = '{} -t {} {} R1_normal.fastq.gz R2_normal.fastq.gz | ' \
              '{} sort --threads {} > normal_paired_aligned_sorted.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
        exec_command(cmd)

        # Cancer (paired)
        cmd = '{} -t {} {} R1_cancer.fastq.gz R2_cancer.fastq.gz | ' \
              '{} sort --threads {} > cancer_paired_aligned_sorted.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
        exec_command(cmd)

        # Normal (unpaired R1)
        cmd = '{} -t {} {} R1_normal_unpaired.fastq.gz | ' \
              '{} sort --threads {} > R1_normal_unpaired_aligned_sorted.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
        exec_command(cmd)

        # Cancer (unpaired R1)
        cmd = '{} -t {} {} R1_cancer_unpaired.fastq.gz | ' \
              '{} sort --threads {} > R1_cancer_unpaired_aligned_sorted.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
        exec_command(cmd)

        # Normal (unpaired R2)
        cmd = '{} -t {} {} R2_normal_unpaired.fastq.gz | ' \
              '{} sort --threads {} > R2_normal_unpaired_aligned_sorted.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
        exec_command(cmd)

        # Cancer (unpaired R2)
        cmd = '{} -t {} {} R2_cancer_unpaired.fastq.gz | ' \
              '{} sort --threads {} > R2_cancer_unpaired_aligned_sorted.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
        exec_command(cmd)

        # Merge aligned files
        print('Merging aligned files')
        cmd = '{} merge -f aligned_normal_merged.bam normal_paired_aligned_sorted.bam ' \
              'R1_normal_unpaired_aligned_sorted.bam R2_normal_unpaired_aligned_sorted.bam'.format(SAMTOOLS)
        exec_command(cmd)

        cmd = '{} merge -f aligned_cancer_merged.bam cancer_paired_aligned_sorted.bam ' \
              'R1_cancer_unpaired_aligned_sorted.bam R2_cancer_unpaired_aligned_sorted.bam'.format(SAMTOOLS)
        exec_command(cmd)

    if 'gatk' in steps:
        # Add headers
        print("Adding headers")
        cmd = '{} AddOrReplaceReadGroups I=aligned_cancer_merged.bam O=sample1_header.bam RGID={} RGPL=Illumina RGLB={} RGPU={} RGSM={}' \
              ' RGCN={} RGDS={}'.format(PICARD, sample1_ID, 'DNA', sample1_ID, sample1_ID, 'VHIO', tumor_type)
        exec_command(cmd)
        cmd = '{} AddOrReplaceReadGroups I=aligned_normal_merged.bam O=sample2_header.bam RGID={} RGPL=Illumina RGLB={} RGPU={} RGSM={}' \
              ' RGCN={} RGDS={}'.format(PICARD, sample2_ID, 'DNA', sample2_ID, sample2_ID, 'VHIO', tumor_type)
        exec_command(cmd)

        # Mark duplicates
        print('Marking duplicates')
        # NOTE setting reducers to it works in system that do not allow many files open
        cmd = GATK + ' MarkDuplicatesSpark -I=sample1_header.bam -O=sample1_dedup.bam -M=dedup_sample1.txt'
        exec_command(cmd)
        cmd = GATK + ' MarkDuplicatesSpark -I=sample2_header.bam -O=sample2_dedup.bam -M=dedup_sample2.txt'
        exec_command(cmd)

        # GATK base re-calibration
        print('Starting re-calibration')
        # NOTE BaseRecalibratorSpark needs the system to allow for many open files (ulimit -n)
        cmd = '{} BaseRecalibrator -I sample1_dedup.bam -R {} --known-sites {} --known-sites {}' \
              ' --known-sites {} -O sample1_recal_data.txt'.format(GATK, genome, SNPSITES, KNOWN_SITE1, KNOWN_SITE2)
        exec_command(cmd)
        cmd = '{} BaseRecalibrator -I sample2_dedup.bam -R {} --known-sites {} --known-sites {}' \
              ' --known-sites {} -O sample2_recal_data.txt'.format(GATK, genome, SNPSITES, KNOWN_SITE1, KNOWN_SITE2)
        exec_command(cmd)
        cmd = '{} ApplyBQSR -R {} -I sample1_dedup.bam --bqsr-recal-file sample1_recal_data.txt -O sample1_final.bam'.format(
            GATK, genome)
        exec_command(cmd)
        cmd = '{} ApplyBQSR -R {} -I sample2_dedup.bam --bqsr-recal-file sample2_recal_data.txt -O sample2_final.bam'.format(
            GATK, genome)
        exec_command(cmd)

    if 'hla' in steps:
        # HLA-LA predictions
        print('Performing HLA-LA predictions')
        HLA_predictionDNA('sample1_final.bam', sampleID, 'PRG-HLA-LA_Tumor_output.txt', THREADS)
        HLA_predictionDNA('sample2_final.bam', sampleID, 'PRG-HLA-LA_Normal_output.txt', THREADS)

    if 'variant' in steps:
        # Variant calling (Samtools pile-ups)
        print('Computing pile-ups')
        cmd = '{} mpileup -C50 -B -q 1 -Q 15 -f {} sample1_final.bam > sample1.pileup'.format(SAMTOOLS, genome)
        exec_command(cmd)
        cmd = '{} mpileup -C50 -B -q 1 -Q 15 -f {} sample2_final.bam > sample2.pileup'.format(SAMTOOLS, genome)
        exec_command(cmd)
        print('Pile-ups were computed for tumor and normal samples')

        print('Performing variant calling Mutect2')
        # Variant calling Mutect2
        cmd = '{} Mutect2 -R {} -I sample1_final.bam -I sample2_final.bam -normal {} -O Mutect_unfiltered.vcf' \
              ' --germline-resource {} --panel-of-normals {}'.format(GATK, genome, sample2_ID, GERMLINE, PON)
        exec_command(cmd)
        cmd = '{} FilterMutectCalls -V Mutect_unfiltered.vcf -O Mutect.vcf -R {}'.format(GATK, genome)
        exec_command(cmd)

        # Variant calling Strelka2
        print('Performing variant calling with Strelka2')
        if os.path.isdir('Strelka_output'):
            shutil.rmtree(os.path.abspath('Strelka_output'))
        cmd = '{} --exome --normalBam sample2_final.bam --tumorBam sample1_final.bam --referenceFasta {}' \
              ' --runDir Strelka_output'.format(STRELKA, genome)
        exec_command(cmd)
        cmd = 'Strelka_output/runWorkflow.py -m local -j {}'.format(THREADS)
        exec_command(cmd)

        # Variant calling Somatic Sniper
        print('Performing variant calling with Somatic Sniper')
        cmd = '{} -L -G -F vcf -f {} sample1_final.bam sample2_final.bam SS.vcf'.format(SSNIPER, genome)
        exec_command(cmd)

        # Variant calling VarScan
        print('Performing variant calling with Varscan')
        cmd = VARSCAN + ' somatic sample2.pileup sample1.pileup varscan --tumor-purity .5 --output-vcf 1' \
                        ' --min-coverage 4 --min-var-freq .05 --strand-filter 0'
        exec_command(cmd)

    if 'filter' in steps:
        print('Filtering variants')

        # TODO filters could be done with vcftools and improved
        mutect2_filter('Mutect.vcf', 'mutect_filtered.vcf', sample1_ID, sample2_ID)
        strelka2_filter('Strelka_output/results/variants/somatic.snvs.vcf.gz', 'strelka_filtered.vcf')
        somaticSniper_filter('SS.vcf', 'somaticsniper_filtered.vcf')
        varscan_filter('varscan.snp.vcf', 'varscan_filtered.vcf')
        strelka2_filter_indels('Strelka_output/results/variants/somatic.indels.vcf.gz', 'strelka_indel_filtered.vcf')
        varscan_filter('varscan.indel.vcf', 'varscan_filtered_indel.vcf')

        # Combine with GATK
        print('Combining SNV variants')
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        cmd = '{} -T CombineVariants -R {} -V:varscan varscan_filtered.vcf -V:mutect mutect_filtered.vcf '\
              '-V:strelka strelka_filtered.vcf -V:somaticsniper somaticsniper_filtered.vcf -o combined_calls.vcf '\
              '-genotypeMergeOptions UNIQUIFY'.format(GATK3, genome)
        exec_command(cmd)

        # Annotate with Annovar
        annovardb = '{} -buildver {}'.format(os.path.join(ANNOVAR_PATH, ANNOVAR_DB), ANNOVAR_VERSION)
        print('Running annovar (SNV)')
        cmd = '{} -format vcf4old combined_calls.vcf --withzyg --comment --includeinfo -outfile snp.av'.format(
            os.path.join(ANNOVAR_PATH, 'convert2annovar.pl'))
        exec_command(cmd)
        cmd = '{} snp.av {} -thread {} -out snp.sum -remove -protocol {}'.format(
            os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), annovardb, THREADS,  annovar_anno)
        exec_command(cmd)

        # Combine with GATK
        print('Combining indels variants')
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        cmd = '{} -T CombineVariants -R {} -V:varscan varscan_filtered_indel.vcf ' \
              '-V:strelka strelka_indel_filtered.vcf -o combined_indel_calls.vcf -genotypeMergeOptions UNIQUIFY'.format(
            GATK3, genome)
        exec_command(cmd)

        # Annotate with Annovar
        print('Runnin annovar (indels)')
        cmd = '{} -format vcf4old combined_indel_calls.vcf --withzyg --comment --includeinfo -outfile indel.av'.format(
            os.path.join(ANNOVAR_PATH, 'convert2annovar.pl'))
        exec_command(cmd)
        cmd = '{} indel.av {} -thread {} -out indel.sum -remove -protocol {}'.format(
            os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), annovardb, THREADS, annovar_anno)
        exec_command(cmd)

        # Extract coverage info from vcf file and add to annotation data
        print("Extracting coverage from combined variants (SNV)")
        vcf = open('combined_calls.vcf')
        vcf_cov_dict = {}
        for line in vcf:
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                varscanT = headers.index('TUMOR.varscan')
                varscanN = headers.index('NORMAL.varscan')
                strelkaT = headers.index('TUMOR.strelka')
                strelkaN = headers.index('NORMAL.strelka')
                mutectT = headers.index(sample1_ID + '.mutect')
                mutectN = headers.index(sample2_ID + '.mutect')
                sniperT = headers.index('TUMOR.somaticsniper')
                sniperN = headers.index('NORMAL.somaticsniper')
            if not line.startswith('#'):
                columns = line.split('\t')
                chrm = columns[0]
                pos = columns[1]
                ref = columns[3]
                alt = columns[4]
                info = columns[7]
                DictID = chrm + ':' + pos
                trfor = '.'
                trrev = '.'
                tvfor = '.'
                tvrev = '.'
                nrfor = '.'
                nrrev = '.'
                nvfor = '.'
                nvrev = '.'
                p_val = '.'
                tcov = '.'
                ncov = '.'
                tfreq = '.'
                nfreq = '.'
                tumor_read1 = '.'
                normal_read1 = '.'
                tumor_read2 = '.'
                normal_read2 = '.'
                callers = info.strip().split(';')[-1].replace('set=', '')
                if re.search('Intersection', callers) or re.search('varscan', callers):
                    if callers == 'Intersection':
                        caller_count = 4
                        callers = 'varscan-strelka-mutect-somaticsniper'
                    else:
                        caller_count = callers.count('-') + 1
                    if re.search('SPV=', info):
                        info_split = info.split(';')
                        for x in info_split:
                            if re.search('SPV=', x):
                                p_val = x.replace('SPV=', '')
                    form = columns[8].split(':')
                    DP4 = form.index('DP4')
                    Freq = form.index('FREQ')
                    t_split = columns[varscanT].split(':')
                    n_split = columns[varscanN].split(':')
                    trfor = int(t_split[DP4].split(',')[0])
                    trrev = int(t_split[DP4].split(',')[1])
                    tvfor = int(t_split[DP4].split(',')[2])
                    tvrev = int(t_split[DP4].split(',')[3])
                    tumor_read1 = trfor + trrev
                    tumor_read2 = tvfor + tvrev
                    tcov = trfor + trrev + tvfor + tvrev
                    tfreq = t_split[Freq]
                    nrfor = int(n_split[DP4].split(',')[0])
                    nrrev = int(n_split[DP4].split(',')[1])
                    nvfor = int(n_split[DP4].split(',')[2])
                    nvrev = int(n_split[DP4].split(',')[3])
                    normal_read1 = nrfor + nrrev
                    normal_read2 = nvfor + nvrev
                    ncov = nrfor + nrrev + nvfor + nvrev
                    nfreq = columns[varscanN].split(':')[Freq]
                elif re.search('somaticsniper', callers):
                    caller_count = callers.count('-') + 1
                    form = columns[8].split(':')
                    DP4 = form.index('DP4')
                    t_split = columns[sniperT].split(':')
                    n_split = columns[sniperN].split(':')
                    trfor = int(t_split[DP4].split(',')[0])
                    trrev = int(t_split[DP4].split(',')[1])
                    tvfor = int(t_split[DP4].split(',')[2])
                    tvrev = int(t_split[DP4].split(',')[3])
                    tumor_read1 = trfor + trrev
                    tumor_read2 = tvfor + tvrev
                    tcov = trfor + trrev + tvfor + tvrev
                    nrfor = int(n_split[DP4].split(',')[0])
                    nrrev = int(n_split[DP4].split(',')[1])
                    nvfor = int(n_split[DP4].split(',')[2])
                    nvrev = int(n_split[DP4].split(',')[3])
                    normal_read1 = nrfor + nrrev
                    normal_read2 = nvfor + nvrev
                    ncov = nrfor + nrrev + nvfor + nvrev
                    tfreq = str(round((tumor_read2 / tcov) * 100, 2)) + '%'
                    nfreq = str(round((normal_read2 / ncov) * 100, 2)) + '%'
                elif re.search('strelka', callers):
                    caller_count = callers.count('-') + 1
                    form = columns[8].split(':')
                    AU = form.index('AU')
                    CU = form.index('CU')
                    GU = form.index('GU')
                    TU = form.index('TU')
                    t_split = columns[strelkaT].split(':')
                    n_split = columns[strelkaN].split(':')
                    if ref == 'A':
                        tumor_read1 = int(t_split[AU].split(',')[0])
                        normal_read1 = int(n_split[AU].split(',')[0])
                    elif ref == 'C':
                        tumor_read1 = int(t_split[CU].split(',')[0])
                        normal_read1 = int(n_split[CU].split(',')[0])
                    elif ref == 'G':
                        tumor_read1 = int(t_split[GU].split(',')[0])
                        normal_read1 = int(n_split[GU].split(',')[0])
                    elif ref == 'T':
                        tumor_read1 = int(t_split[TU].split(',')[0])
                        normal_read1 = int(n_split[TU].split(',')[0])
                    if alt == 'A':
                        tumor_read2 = int(t_split[AU].split(',')[0])
                        normal_read2 = int(n_split[AU].split(',')[0])
                    elif alt == 'C':
                        tumor_read2 = int(t_split[CU].split(',')[0])
                        normal_read2 = int(n_split[CU].split(',')[0])
                    elif alt == 'G':
                        tumor_read2 = int(t_split[GU].split(',')[0])
                        normal_read2 = int(n_split[GU].split(',')[0])
                    elif alt == 'T':
                        tumor_read2 = int(t_split[TU].split(',')[0])
                        normal_read2 = int(n_split[TU].split(',')[0])
                    if tumor_read2 != '.':
                        tcov = tumor_read1 + tumor_read2
                        ncov = normal_read1 + normal_read2
                        tfreq = str(round((tumor_read2 / tcov) * 100, 2)) + '%'
                        nfreq = str(round((normal_read2 / ncov) * 100, 2)) + '%'
                elif re.search('mutect', callers):
                    caller_count = callers.count('-') + 1
                    form = columns[8].split(':')
                    if not re.search(',', alt):
                        AD = form.index('AD')
                    t_split = columns[mutectT].split(':')
                    n_split = columns[mutectN].split(':')
                    if not re.search(',', alt):
                        tumor_read1 = int(t_split[AD].split(',')[0])
                        tumor_read2 = int(t_split[AD].split(',')[1])
                        normal_read1 = int(n_split[AD].split(',')[0])
                        normal_read2 = int(n_split[AD].split(',')[1])
                        tcov = tumor_read1 + tumor_read2
                        ncov = normal_read1 + normal_read2
                        tfreq = str(round((tumor_read2 / tcov) * 100, 2)) + '%'
                        nfreq = str(round((normal_read2 / ncov) * 100, 2)) + '%'
                # populate
                vcf_cov_dict[DictID] = {}
                vcf_cov_dict[DictID]['pval'] = p_val
                vcf_cov_dict[DictID]['Note'] = str(caller_count) + ':' + callers
                vcf_cov_dict[DictID]['trfor'] = trfor
                vcf_cov_dict[DictID]['trrev'] = trrev
                vcf_cov_dict[DictID]['tvfor'] = tvfor
                vcf_cov_dict[DictID]['tvrev'] = tvrev
                vcf_cov_dict[DictID]['nrfor'] = nrfor
                vcf_cov_dict[DictID]['nrrev'] = nrrev
                vcf_cov_dict[DictID]['nvfor'] = nvfor
                vcf_cov_dict[DictID]['nvrev'] = nvrev
                vcf_cov_dict[DictID]['tumor_read1'] = tumor_read1
                vcf_cov_dict[DictID]['tumor_read2'] = tumor_read2
                vcf_cov_dict[DictID]['normal_read1'] = normal_read1
                vcf_cov_dict[DictID]['normal_read2'] = normal_read2
                vcf_cov_dict[DictID]['tumor_coverage'] = tcov
                vcf_cov_dict[DictID]['tumor_freq'] = tfreq
                vcf_cov_dict[DictID]['normal_coverage'] = ncov
                vcf_cov_dict[DictID]['normal_freq'] = nfreq
        vcf.close()

        print('Adding annotated SNVs to final report')
        final_variants('snp.sum.{}_multianno.txt'.format(ANNOVAR_VERSION),
                       'nonsyn_SQL_insert.txt', 'all_other_mutations.txt',
                       vcf_cov_dict, sampleID, tumor_type, header=True)

        # Extract coverage info from vcf file
        print("Extracting coverage from combined variants (indels)")
        vcf = open('combined_indel_calls.vcf')
        vcf_cov_dict = {}
        for line in vcf:
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                varscanT = headers.index('TUMOR.varscan')
                varscanN = headers.index('NORMAL.varscan')
                strelkaT = headers.index('TUMOR.strelka')
                strelkaN = headers.index('NORMAL.strelka')
            if not line.startswith('#'):
                columns = line.split('\t')
                chrm = columns[0]
                pos = columns[1]
                info = columns[7]
                DictID = chrm + ':' + pos
                trfor = '.'
                trrev = '.'
                tvfor = '.'
                tvrev = '.'
                nrfor = '.'
                nrrev = '.'
                nvfor = '.'
                nvrev = '.'
                p_val = '.'
                callers = info.strip().split(';')[-1].replace('set=', '')
                if re.search('Intersection', callers) or re.search('varscan', callers):
                    if callers == 'Intersection':
                        caller_count = 2
                        callers = 'varscan-strelka'
                    else:
                        caller_count = callers.count('-') + 1
                    if re.search('SPV=', info):
                        info_split = info.split(';')
                        for x in info_split:
                            if re.search('SPV=', x):
                                p_val = x.replace('SPV=', '')
                    else:
                        p_val = '.'
                    form = columns[8].split(':')
                    DP4 = form.index('DP4')
                    Freq = form.index('FREQ')
                    t_split = columns[varscanT].split(':')
                    n_split = columns[varscanN].split(':')
                    trfor = int(t_split[DP4].split(',')[0])
                    trrev = int(t_split[DP4].split(',')[1])
                    tvfor = int(t_split[DP4].split(',')[2])
                    tvrev = int(t_split[DP4].split(',')[3])
                    tumor_read1 = trfor + trrev
                    tumor_read2 = tvfor + tvrev
                    tcov = trfor + trrev + tvfor + tvrev
                    tfreq = t_split[Freq]
                    nrfor = int(n_split[DP4].split(',')[0])
                    nrrev = int(n_split[DP4].split(',')[1])
                    nvfor = int(n_split[DP4].split(',')[2])
                    nvrev = int(n_split[DP4].split(',')[3])
                    normal_read1 = nrfor + nrrev
                    normal_read2 = nvfor + nvrev
                    ncov = nrfor + nrrev + nvfor + nvrev
                    nfreq = n_split[Freq]
                elif re.search('strelka', callers):
                    caller_count = callers.count('-') + 1
                    form = columns[8].split(':')
                    DP = form.index('DP')
                    TIR = form.index('TIR')
                    t_split = columns[strelkaT].split(':')
                    n_split = columns[strelkaN].split(':')
                    normal_read2 = int(n_split[TIR].split(',')[1])
                    tumor_read2 = int(t_split[TIR].split(',')[1])
                    n_cov = int(n_split[DP])
                    t_cov = int(t_split[DP])
                    tfreq = float((tumor_read2 / t_cov) * 100)
                    if normal_read2 != 0:
                        nfreq = float(normal_read2 / n_cov)
                    else:
                        nfreq = 0
                    tumor_read1 = t_cov - tumor_read2
                    normal_read1 = n_cov - normal_read2
                vcf_cov_dict[DictID] = {}
                vcf_cov_dict[DictID]['pval'] = p_val
                vcf_cov_dict[DictID]['Note'] = str(caller_count) + ':' + callers
                vcf_cov_dict[DictID]['trfor'] = trfor
                vcf_cov_dict[DictID]['trrev'] = trrev
                vcf_cov_dict[DictID]['tvfor'] = tvfor
                vcf_cov_dict[DictID]['tvrev'] = tvrev
                vcf_cov_dict[DictID]['nrfor'] = nrfor
                vcf_cov_dict[DictID]['nrrev'] = nrrev
                vcf_cov_dict[DictID]['nvfor'] = nvfor
                vcf_cov_dict[DictID]['nvrev'] = nvrev
                vcf_cov_dict[DictID]['tumor_read1'] = tumor_read1
                vcf_cov_dict[DictID]['tumor_read2'] = tumor_read2
                vcf_cov_dict[DictID]['normal_read1'] = normal_read1
                vcf_cov_dict[DictID]['normal_read2'] = normal_read2
                vcf_cov_dict[DictID]['tumor_coverage'] = tcov
                vcf_cov_dict[DictID]['tumor_freq'] = tfreq
                vcf_cov_dict[DictID]['normal_coverage'] = ncov
                vcf_cov_dict[DictID]['normal_freq'] = nfreq
        vcf.close()

        print('Adding annotated indels to final report')
        final_variants('indel.sum.{}_multianno.txt'.format(ANNOVAR_VERSION),
                       'nonsyn_SQL_insert.txt', 'all_other_mutations.txt',
                       vcf_cov_dict, sampleID, tumor_type, header=False)

        # Extract peptides
        print("Extracting peptides")
        snv = open('nonsyn_SQL_insert.txt')
        snv_lines = snv.readlines()
        header = snv_lines.pop(0).strip().split('\t')
        epitope_file = open('Formatted_epitope_variant.txt', 'w')
        for line in snv_lines:
            columns = line.strip().split('\t')
            variant_key = columns[header.index('VARIANT-KEY')]
            chrom = columns[header.index('CHR')]
            start = columns[header.index('START')]
            stop = columns[header.index('END')]
            ref = columns[header.index('REF')]
            alt = columns[header.index('ALT')]
            func_ref_gene = columns[header.index('Func.refGene')]
            exonic_func_ref = columns[header.index('ExonicFunc.refGene')]
            AA_change_refGene = columns[header.index('AAChange.refGene')].split(',')
            func_UCSC_gene = columns[header.index('Func.knownGene')]
            exonic_func_UCSC = columns[header.index('ExonicFunc.knownGene')]
            AA_change_UCSCGene = columns[header.index('AAChange.knownGene')].split(',')
            func_ens_gene = columns[header.index('Func.ensGene')]
            exonic_func_ens = columns[header.index('ExonicFunc.ensGene')]
            AA_change_ensGene = columns[header.index('AAChange.ensGene')].split(',')
            if re.search(r'nonsynonymous', exonic_func_ref) or re.search(r'frame', exonic_func_ref):
                for entry in AA_change_refGene:
                    epitope_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sampleID,
                                                                                         variant_key,
                                                                                         chrom,
                                                                                         start,
                                                                                         stop,
                                                                                         ref,
                                                                                         alt,
                                                                                         func_ref_gene,
                                                                                         exonic_func_ref,
                                                                                         sub(':', '\t', entry)))
            if re.search(r'nonsynonymous', exonic_func_UCSC) or re.search(r'frame', exonic_func_UCSC):
                for entry in AA_change_UCSCGene:
                    epitope_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sampleID,
                                                                                         variant_key,
                                                                                         chrom,
                                                                                         start,
                                                                                         stop,
                                                                                         ref,
                                                                                         alt,
                                                                                         func_UCSC_gene,
                                                                                         exonic_func_UCSC,
                                                                                         sub(':', '\t', entry)))
            if re.search(r'nonsynonymous', exonic_func_ens) or re.search(r'frame', exonic_func_ens):
                for entry in AA_change_ensGene:
                    epitope_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sampleID,
                                                                                         variant_key,
                                                                                         chrom,
                                                                                         start,
                                                                                         stop,
                                                                                         ref,
                                                                                         alt,
                                                                                         func_ens_gene,
                                                                                         exonic_func_ens,
                                                                                         sub(':', '\t', entry)))
        epitope_file.close()
        snv.close()

        # Create list of AA and cDNA sequences
        print('Creating epitopes')
        AA_seq = dict()
        with open(FASTA_AA_DICT, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                AA_seq[record.id.split("|")[0]] = record.seq
        cDNA_seq = dict()
        with open(FASTA_cDNA_DICT, "rU") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                cDNA_seq[record.id.split("|")[0]] = record.seq
        epitope_file = open('SQL_Epitopes.txt', 'w')
        input_file = open('Formatted_epitope_variant.txt')
        header = 'SAMPLE_ID\tVARIANT-KEY\tCHR\tSTART\tSTOP\tREF\tALT\tfunc_ref_gene\texonic_func_ref\tGene\t' \
                 'Transcript_ID\tExon_Numb\tNT_CHANGE\tAA_CHANGE\tPOSITION\tERRORS\tWT25MER\tMUT25MER\n'
        epitope_file.write(header)
        for line in input_file:
            columns = line.rstrip('\n').split('\t')
            if len(columns) != 14:
                print("Formatted epitote with wrong number of columns {}".format(','.format(columns)))
                continue
            try:
                ref = columns[5].strip()
                exonic_func = columns[8].strip()
                transcriptID = columns[10].strip()
                cDNA_strip = columns[12].strip()
                protein_strip = columns[13].strip()
                errors = 'Flags:'
                WT_25mer = '-'
                Mut_25mer = '-'
                position = int(re.findall(r'\d+', protein_strip)[0])
                cDNA_pos = int(re.findall(r'\d+', cDNA_strip)[0])
                protein_seq = AA_seq.get(transcriptID, 'AA_seq not present for this transcript').strip()
                ref_cDNA_seq = cDNA_seq.get(transcriptID, 'cDNA not present for this transcript').strip()
                # Nonsynonymous point mutations to 25 mers
                if exonic_func == 'nonsynonymous SNV':
                    # extract the AA change info
                    ref_AA = protein_strip[protein_strip.find('.') + 1]
                    var_AA = protein_strip[len(protein_strip) - 1]
                    # gather AA seq for transcript
                    if protein_seq == 'AA_seq not present for this transcript':
                        errors += 'AA_seq not present for this transcript '
                    else:
                        # check annotation is correct
                        FASTA_AA = protein_seq[position - 1:position]
                        if FASTA_AA == ref_AA:
                            if position >= 13:
                                WT_25mer = protein_seq[position - 13:position + 12]
                                Mut_25mer = protein_seq[position - 13:position - 1] + var_AA + protein_seq[position:position + 12]
                            elif position < 13:
                                WT_25mer = protein_seq[0:position + 12]
                                Mut_25mer = protein_seq[0:position - 1] + var_AA + protein_seq[position:position + 12]
                            if position == 1:
                                errors += ' mutation occurs in start codon'
                        else:
                            errors += ' Ref in AA_seq does not match file Ref'
                # frameshift/non-frameshift insertions/deletions/substitutions
                elif 'frameshift' in exonic_func:
                    if ref_cDNA_seq == 'cDNA not present for this transcript':
                        errors += ' cDNA not present for this transcript'
                    elif exonic_func in ['frameshift deletion', 'nonframeshift deletion']:
                        len_del = len(ref)
                        mut_cDNA_left = ref_cDNA_seq[0:cDNA_pos - 1]
                        mut_cDNA_right = ref_cDNA_seq[cDNA_pos + len_del - 1:]
                        mut_cDNA_seq = mut_cDNA_left + mut_cDNA_right
                    elif exonic_func in ['frameshift insertion', 'nonframeshift insertion']:
                        if re.search(r'dup', cDNA_strip):
                            ins = cDNA_strip[int(cDNA_strip.find('dup')) + 3:]
                        elif re.search(r'ins', cDNA_strip):
                            ins = cDNA_strip[int(cDNA_strip.find('ins')) + 3:]
                        else:
                            errors += ' could not find mutation in cDNA'
                            ins = ''
                        mut_cDNA_left = ref_cDNA_seq[0:cDNA_pos]
                        mut_cDNA_right = ref_cDNA_seq[cDNA_pos:]
                        mut_cDNA_seq = mut_cDNA_left + ins + mut_cDNA_right
                    elif exonic_func in ['frameshift substitution', 'nonframeshift substitution']:
                        if re.search(r'delins', cDNA_strip):
                            subs = cDNA_strip[int(cDNA_strip.find('delins')) + 6:]
                        else:
                            errors += ' could not find mutation in cDNA'
                            subs = ''
                        mut_cDNA_left = ref_cDNA_seq[0:cDNA_pos]
                        mut_cDNA_right = ref_cDNA_seq[cDNA_pos:]
                        mut_cDNA_seq = mut_cDNA_left + subs + mut_cDNA_right
                    else:
                        mut_cDNA_seq = ref_cDNA_seq
                        errors += ' Unknown exonic function'
                    if not protein_strip.startswith('p.'):
                        position = 0
                    elif protein_strip.startswith('p.X'):
                        position = 0
                        errors += ' mutation occurs in stop codon'
                    ref_FASTA = translate_dna(ref_cDNA_seq)
                    mut_FASTA = translate_dna(mut_cDNA_seq)
                    mut_stop = int(mut_FASTA.find('X'))
                    if position >= 13 and mut_stop > 0:
                        WT_25mer = ref_FASTA[position - 13:position + 12].replace('X', '')
                        Mut_25mer = mut_FASTA[position - 13:mut_stop]
                    elif position < 13 and position > 0 and mut_stop > 0:
                        WT_25mer = ref_FASTA[0:position + 12].replace('X', '')
                        Mut_25mer = mut_FASTA[0:mut_stop]
                    elif position >= 13 and mut_stop < 0:
                        WT_25mer = ref_FASTA[position - 13:position + 12].replace('X', '')
                        Mut_25mer = mut_FASTA[position - 13:]
                    elif position < 13 and position > 0 and mut_stop < 0:
                        WT_25mer = ref_FASTA[0:position + 12].replace('X', '')
                        Mut_25mer = mut_FASTA[0:]
                    elif position == 0:
                        errors += ' can not code for this mutated AA_position'
                    if not ref_cDNA_seq.startswith('ATG'):
                        errors += ' No ATG start codon for this transcript cDNA'
                    if position == 1:
                        errors += ' mutation occurs in start codon'
                elif re.search(r'^stop', exonic_func):
                    position = ''.join([s for s in protein_strip if s.isdigit()])
                    errors += ' Stop mutation'
                else:
                    errors += ' Unknown exonic function'
                epitope_file.write('{}\t{}\t{}\t{}\t{}\n'.format('\t'.join(columns[0:]),
                                                                 position,
                                                                 errors,
                                                                 WT_25mer,
                                                                 Mut_25mer))
            except Exception as err:
                print('An unknown error {}\t happened parsing epitope {}'.format(err, ','.join(columns)))
        input_file.close()
        epitope_file.close()

    print("COMPLETED!")

parser = argparse.ArgumentParser(description='Exome variant calling and HLA prediction pipeline '
                                 '(created by Jose Fernandez <jc.fernandes.navarro@gmail.com>)',
                                 prog='exome_pipeline.py',
                                 usage='exome_pipeline.py [options] R1(Normal) R2(Normal) R1(Cancer) R2(Cancer)')
parser.add_argument('R1_NORMAL', help='FASTQ file R1 (Normal)')
parser.add_argument('R2_NORMAL', help='FASTQ file R2 (Normal)')
parser.add_argument('R1_CANCER', help='FASTQ file R1 (Cancer)')
parser.add_argument('R2_CANCER', help='FASTQ file R2 (Cancer)')
parser.add_argument('--adapter',
                    help='Path to the Illumina adapters FASTA file.', required=True)
parser.add_argument('--genome',
                    help='Path to the reference Genome FASTA file (must contain BWA index)', required=True)
parser.add_argument('--sample',
                    help='Name of the sample/experiment. Default is sample', default='sample')
parser.add_argument('--tumor',
                    help='Tumor type. Default is Tumor', default='Tumor')
parser.add_argument('--dir',
                    help='Path to the output folder where output files will be placed', required=True)
parser.add_argument('--known1',
                    help='Path to the file with Mill and 1000G gold standards (GATK bundle)', required=True)
parser.add_argument('--known2',
                    help='Path to the file with 1000G phase indels (GATK bundle)', required=True)
parser.add_argument('--snpsites',
                    help='Path to the file with the SNPs (GATK bundle)', required=True)
parser.add_argument('--germline',
                    help='Path to the file with the germline resources Nomad for Mutect2 (GATK bundle)', required=True)
parser.add_argument('--pon',
                    help='Path to the file with the panel of normals for Mutect2 (GATK bundle)', required=True)
parser.add_argument('--fastaAA',
                    help='Path to the fasta file with the protein sequences (of transcripts)', required=True)
parser.add_argument('--fastacDNA',
                    help='Path to the fasta file with the cDNA sequences (of transcripts)', required=True)
parser.add_argument('--annovar-db',
                    help='String indicated what annovar database to use (default: humandb)',
                    default='humandb', required=False)
parser.add_argument('--annovar-version',
                    help='String indicated what version of the annovar database to use (default: hg19)',
                    default='hg19', required=False)
parser.add_argument('--steps', nargs='+', default=['mapping', 'gatk', 'hla', 'variant', 'filter'],
                    help='Steps to perform in the pipeline', choices=['mapping', 'gatk', 'hla', 'variant', 'filter', "none"])

# Parse arguments
args = parser.parse_args()
DIR = args.dir
R1_NORMAL = os.path.abspath(args.R1_NORMAL)
R2_NORMAL = os.path.abspath(args.R2_NORMAL)
R1_CANCER = os.path.abspath(args.R1_CANCER)
R2_CANCER = os.path.abspath(args.R2_CANCER)
sampleID = args.sample
tumor_type = args.tumor
IILLUMINA_ADAPTERS = os.path.abspath(args.adapter)
GENOME_REF = os.path.abspath(args.genome)
THREADS = multiprocessing.cpu_count() - 1
FASTA_AA_DICT = os.path.abspath(args.fastaAA)
FASTA_cDNA_DICT = os.path.abspath(args.fastacDNA)
KNOWN_SITE1 = os.path.abspath(args.known1)
KNOWN_SITE2 = os.path.abspath(args.known2)
SNPSITES = os.path.abspath(args.snpsites)
GERMLINE = os.path.abspath(args.germline)
PON = os.path.abspath(args.pon)
DNA_STEPS = args.steps
ANNOVAR_DB = args.annovar_db
ANNOVAR_VERSION = args.annovar_version

# Move to output dir
os.makedirs(os.path.abspath(DIR), exist_ok=True)
os.chdir(os.path.abspath(DIR))

exome_pipeline(R1_NORMAL,
               R2_NORMAL,
               R1_CANCER,
               R2_CANCER,
               IILLUMINA_ADAPTERS,
               tumor_type,
               GENOME_REF,
               sampleID,
               THREADS,
               FASTA_AA_DICT,
               FASTA_cDNA_DICT,
               KNOWN_SITE1,
               KNOWN_SITE2,
               SNPSITES,
               GERMLINE,
               PON,
               ANNOVAR_DB,
               ANNOVAR_VERSION,
               DNA_STEPS)
