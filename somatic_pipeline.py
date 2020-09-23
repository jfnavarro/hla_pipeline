#! /usr/bin/env python
"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from hlapipeline.common import *
from hlapipeline.tools import *
from hlapipeline.filters import *
from hlapipeline.epitopes import *
import shutil
import argparse
import os
import sys

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
    # TODO use header names instead of accessing the fields by index
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
        ExonicFunc_refGene = columns[8]
        ExonicFunc_knownGene = columns[13]
        ExonicFunc_ensGene = columns[18]
        variant_key = Chr + ':' + start + '-' + end + ' ' + ref + '>' + alt
        if ref_gene_detail != 'NA':
            columns[9] = ref_gene_detail
        if known_gene_detail != 'NA':
            columns[14] = known_gene_detail
        if ens_gene_detail != 'NA':
            columns[19] = ens_gene_detail
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
            if any(x in ' '.join([ExonicFunc_refGene, ExonicFunc_knownGene, ExonicFunc_ensGene]) for x in ['nonsynonymous', 'frame', 'stop']):
                nonsyn_file.write(to_write + '\n')
            else:
                all_file.write(to_write + '\n')
        except KeyError:
            print("Missing variant for {}".format(ID))
    nonsyn_file.close()
    all_file.close()
    nonsyn_snv.close()

def somatic_pipeline(R1_NORMAL,
                     R2_NORMAL,
                     R1_CANCER,
                     R2_CANCER,
                     tumor_type,
                     genome,
                     genome_star,
                     annotation,
                     sampleID,
                     THREADS,
                     AA_DICT,
                     cDNA_DICT,
                     KNOWN_SITE1,
                     KNOWN_SITE2,
                     SNPSITES,
                     GERMLINE,
                     PON,
                     ANNOVAR_DB,
                     ANNOVAR_VERSION,
                     steps,
                     mode):
    print("Somatic pipeline")

    # Sample 1 cancer, sample 2 normal
    sample1_ID = sampleID + "_Tumor"
    sample2_ID = sampleID + "_Normal"

    # Create sub-folder to store all results
    os.makedirs('somatic', exist_ok=True)
    os.chdir('somatic')

    if 'mapping' in steps:
        # TRIMMING
        print('Starting trimming')
        cmd = '{} --cores {} --paired --basename normal {} {}'.format(TRIMGALORE, THREADS, R1_NORMAL, R2_NORMAL)
        exec_command(cmd)

        cmd = '{} --cores {} --paired --basename cancer {} {}'.format(TRIMGALORE, THREADS, R1_CANCER, R2_CANCER)
        exec_command(cmd)

        # ALIGNMENT
        print('Starting alignment')
        if mode in ['DNA', 'DNA-RNA']:
            # Normal (paired)
            cmd = '{} -t {} {} normal_val_1.fq.gz normal_val_1.fq.gz | ' \
                  '{} sort --threads {} > aligned_normal_merged.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
            exec_command(cmd)
        else:
            # Normal (paired)
            cmd = '{} --genomeDir {} --readFilesIn normal_val_1.fq.gz normal_val_2.fq.gz --outSAMorder Paired' \
                  ' --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {}' \
                  ' --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c' \
                  ' --runThreadN {}'.format(STAR, genome_star, annotation, THREADS)
            exec_command(cmd)
            cmd = 'mv Aligned.sortedByCoord.out.bam aligned_normal_merged.bam'
            exec_command(cmd)

        if mode in ['RNA', 'DNA-RNA']:
            # Cancer (paired)
            cmd = '{} --genomeDir {} --readFilesIn cancer_val_1.fq.gz cancer_val_2.fq.gz --outSAMorder Paired' \
                  ' --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {}' \
                  ' --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c' \
                  ' --runThreadN {}'.format(STAR, genome_star, annotation, THREADS)
            exec_command(cmd)

            cmd = 'mv Aligned.sortedByCoord.out.bam aligned_cancer_merged.bam'
            exec_command(cmd)
        else:
            # Cancer (paired)
            cmd = '{} -t {} {} cancer_val_1.fq.gz cancer_val_2.fq.gz | ' \
                  '{} sort --threads {} > aligned_cancer_merged.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
            exec_command(cmd)

        # Add headers
        print("Adding headers")
        cmd = '{} AddOrReplaceReadGroups I=aligned_cancer_merged.bam O=sample1_header.bam RGID={} RGPL=Illumina ' \
              'RGLB={} RGPU={} RGSM={} RGCN=VHIO RGDS={}'.format(PICARD, sample1_ID, mode, sample1_ID, sample1_ID,
                                                                 tumor_type)
        exec_command(cmd)

        cmd = '{} AddOrReplaceReadGroups I=aligned_normal_merged.bam O=sample2_header.bam RGID={} RGPL=Illumina ' \
              'RGLB={} RGPU={} RGSM={} RGCN=VHIO RGDS={}'.format(PICARD, sample2_ID, mode, sample2_ID, sample2_ID,
                                                                 tumor_type)
        exec_command(cmd)

    if 'gatk' in steps:
        # Mark duplicates
        print('Marking duplicates')
        # NOTE setting reducers to it works in system that do not allow many files open
        cmd = '{} MarkDuplicatesSpark --input sample1_header.bam --output sample1_dedup.bam'.format(GATK)
        exec_command(cmd)
        cmd = '{} MarkDuplicatesSpark --input sample2_header.bam --output sample2_dedup.bam'.format(GATK)
        exec_command(cmd)

        # GATK base re-calibration
        print('Starting re-calibration')
        cmd = '{} BaseRecalibratorSpark --input sample1_dedup.bam --reference {} --known-sites {} --known-sites {}' \
              ' --known-sites {} --output sample1_recal_data.txt'.format(GATK, genome, SNPSITES, KNOWN_SITE1, KNOWN_SITE2)
        exec_command(cmd)

        cmd = '{} BaseRecalibratorSpark --input sample2_dedup.bam --reference {} --known-sites {} --known-sites {}' \
              ' --known-sites {} --output sample2_recal_data.txt'.format(GATK, genome, SNPSITES, KNOWN_SITE1, KNOWN_SITE2)
        exec_command(cmd)

        cmd = '{} ApplyBQSR --reference {} --input sample1_dedup.bam --bqsr-recal-file sample1_recal_data.txt ' \
              '--output sample1_final.bam'.format(GATK, genome)
        exec_command(cmd)

        cmd = '{} ApplyBQSR --reference {} --input sample2_dedup.bam --bqsr-recal-file sample2_recal_data.txt ' \
              '--output sample2_final.bam'.format(GATK, genome)
        exec_command(cmd)

    if 'hla' in steps:
        # HLA-LA predictions
        print('Performing HLA-LA predictions')
        if mode == 'DNA':
            HLA_predictionDNA('sample2_final.bam', sampleID, 'PRG-HLA-LA_Normal_output.txt', THREADS)
            HLA_predictionDNA('sample1_final.bam', sampleID, 'PRG-HLA-LA_Tumor_output.txt', THREADS)
        else:
            HLA_predictionRNA('sample2_final.bam', THREADS)
            HLA_predictionRNA('sample1_final.bam', THREADS)
            
    if 'variant' in steps:
        # Variant calling (Samtools pile-ups)
        print('Computing pile-ups')
        cmd = '{} mpileup -C 50 -B -q 1 -Q 15 -f {} sample1_final.bam > sample1.pileup'.format(SAMTOOLS, genome)
        exec_command(cmd)
        cmd = '{} mpileup -C 50 -B -q 1 -Q 15 -f {} sample2_final.bam > sample2.pileup'.format(SAMTOOLS, genome)
        exec_command(cmd)

        print('Performing variant calling Mutect2')
        # Variant calling Mutect2
        cmd = '{} Mutect2 --reference {} --input sample1_final.bam --input sample2_final.bam --normal-sample {} ' \
              '--output Mutect_unfiltered.vcf --germline-resource {} --panel-of-normals {}'.format(GATK,
                                                                                                   genome,
                                                                                                   sample2_ID,
                                                                                                   GERMLINE,
                                                                                                   PON)
        exec_command(cmd)
        cmd = '{} FilterMutectCalls --variant Mutect_unfiltered.vcf --output Mutect.vcf --reference {}'.format(GATK,
                                                                                                               genome)
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
        print('Performing variant calling with SomaticSniper')
        cmd = '{} -L -G -F vcf -f {} sample1_final.bam sample2_final.bam SS.vcf'.format(SSNIPER, genome)
        exec_command(cmd)

        # Variant calling VarScan
        print('Performing variant calling with VarScan2')
        # TODO improve the filters, enable perhaps the p-value filter
        cmd = '{} somatic sample2.pileup sample1.pileup varscan --tumor-purity 0.5 --output-vcf 1 ' \
              '--min-coverage 4 --min-var-freq 0.05 --strand-filter 0'.format(VARSCAN)
        exec_command(cmd)

        if mode in ['RNA', 'DNA-RNA']:
            # Computing gene counts
            print('Computing gene counts with featureCounts for tumor sample')
            cmd = '{} -T {} --primary --ignoreDup -O -C -t exon ' \
                  '-g gene_name -a {} -o tumor-gene.counts sample1_dedup.bam'.format(FEATURECOUNTS, THREADS, annotation)
            exec_command(cmd)

        if mode in ['RNA', 'RNA-DNA']:
            print('Computing gene counts with featureCounts for normal sample')
            cmd = '{} -T {} --primary --ignoreDup -O -C -t exon ' \
                  '-g gene_name -a {} -o normal-gene.counts sample2_dedup.bam'.format(FEATURECOUNTS, THREADS, annotation)
            exec_command(cmd)

    if 'filter' in steps:
        print('Filtering variants')
        #TODO filters could be performed with VariantFiltration and/or vcftools
        #TODO filters seems arbitrary, they could be improved
        mutect2_filter('Mutect.vcf', 'mutect_filtered.vcf', sample1_ID, sample2_ID)
        strelka2_filter('Strelka_output/results/variants/somatic.snvs.vcf.gz', 'strelka_filtered.vcf')
        somaticSniper_filter('SS.vcf', 'somaticsniper_filtered.vcf')
        varscan_filter('varscan.snp.vcf', 'varscan_filtered.vcf')
        strelka2_filter_indels('Strelka_output/results/variants/somatic.indels.vcf.gz', 'strelka_indel_filtered.vcf')
        varscan_filter('varscan.indel.vcf', 'varscan_filtered_indel.vcf')

        # Combine with GATK
        print('Combining SNP variants')
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        cmd = '{} -T CombineVariants -R {} -V:varscan varscan_filtered.vcf -V:mutect mutect_filtered.vcf '\
              '-V:strelka strelka_filtered.vcf -V:somaticsniper somaticsniper_filtered.vcf -o combined_calls.vcf '\
              '-genotypeMergeOptions UNIQUIFY --num_threads {}'.format(GATK3, genome, THREADS)
        exec_command(cmd)

        # Annotate with Annovar
        annovardb = '{} -buildver {}'.format(os.path.join(ANNOVAR_PATH, ANNOVAR_DB), ANNOVAR_VERSION)
        print('Running annovar (SNP)')
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
                chrm = columns[headers.index('#CHROM')]
                pos = columns[headers.index('POS')]
                ref = columns[headers.index('REF')]
                alt = columns[headers.index('ALT')]
                info = columns[headers.index('INFO')]
                form = columns[headers.index('FORMAT')].split(':')
                DictID = chrm + ':' + pos
                trfor = 0
                trrev = 0
                tvfor = 0
                tvrev = 0
                nrfor = 0
                nrrev = 0
                nvfor = 0
                nvrev = 0
                p_val = -1
                tcov = 0
                ncov = 0
                tfreq = 0
                nfreq = 0
                tumor_read1 = 0
                normal_read1 = 0
                tumor_read2 = 0
                normal_read2 = 0
                callers = info.strip().split(';')[-1].replace('set=', '')
                caller_count = callers.count('-') + 1 if 'Intersection' not in callers else 4
                if 'Intersection' in callers or 'varscan' in callers:
                    if 'SPV=' in info:
                        for x in info.split(';'):
                            if 'SPV=' in x:
                                p_val = x.replace('SPV=', '')
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
                elif 'somaticsniper' in callers:
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
                    if tumor_read2 != 0:
                        tfreq = str(round((tumor_read2 / tcov) * 100, 2)) + '%'
                    if normal_read2 != 0:
                        nfreq = str(round((normal_read2 / ncov) * 100, 2)) + '%'
                elif 'strelka' in callers:
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
                    tcov = tumor_read1 + tumor_read2
                    ncov = normal_read1 + normal_read2
                    if tumor_read2 != 0:
                        tfreq = str(round((tumor_read2 / tcov) * 100, 2)) + '%'
                    if normal_read2 != 0:
                        nfreq = str(round((normal_read2 / ncov) * 100, 2)) + '%'
                elif 'mutect' in callers:
                    t_split = columns[mutectT].split(':')
                    n_split = columns[mutectN].split(':')
                    if ',' not in alt:
                        AD = form.index('AD')
                        tumor_read1 = int(t_split[AD].split(',')[0])
                        tumor_read2 = int(t_split[AD].split(',')[1])
                        normal_read1 = int(n_split[AD].split(',')[0])
                        normal_read2 = int(n_split[AD].split(',')[1])
                        tcov = tumor_read1 + tumor_read2
                        ncov = normal_read1 + normal_read2
                        if tumor_read2 != 0:
                            tfreq = str(round((tumor_read2 / tcov) * 100, 2)) + '%'
                        if normal_read2 != 0:
                            nfreq = str(round((normal_read2 / ncov) * 100, 2)) + '%'
                vcf_cov_dict[DictID] = {}
                vcf_cov_dict[DictID]['pval'] = p_val
                vcf_cov_dict[DictID]['Note'] = str(caller_count) + ":" + callers
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
                chrm = columns[headers.index('#CHROM')]
                pos = columns[headers.index('POS')]
                info = columns[headers.index('INFO')]
                form = columns[headers.index('FORMAT')].split(':')
                DictID = chrm + ':' + pos
                trfor = 0
                trrev = 0
                tvfor = 0
                tvrev = 0
                nrfor = 0
                nrrev = 0
                nvfor = 0
                nvrev = 0
                p_val = -1
                tumor_read1 = 0
                tumor_read2 = 0
                normal_read1 = 0
                normal_read2 = 0
                tcov = 0
                nfreq = 0
                tfreq = 0
                ncov = 0
                callers = info.strip().split(';')[-1].replace('set=', '')
                caller_count = callers.count('-') + 1 if 'Intersection' not in callers else 2
                if 'Intersection' in callers or 'varscan' in callers:
                    if 'SPV=' in info:
                        for x in info.split(';'):
                            if 'SPV=' in x:
                                p_val = x.replace('SPV=', '')
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
                elif 'strelka' in callers:
                    caller_count = callers.count('-') + 1
                    DP = form.index('DP')
                    TIR = form.index('TIR')
                    t_split = columns[strelkaT].split(':')
                    n_split = columns[strelkaN].split(':')
                    normal_read2 = int(n_split[TIR].split(',')[1])
                    tumor_read2 = int(t_split[TIR].split(',')[1])
                    ncov = int(n_split[DP])
                    tcov = int(t_split[DP])
                    if tumor_read2 != 0:
                        tfreq = float((tumor_read2 / tcov) * 100)
                    if normal_read2 != 0:
                        nfreq = float(normal_read2 / ncov)
                    tumor_read1 = tcov - tumor_read2
                    normal_read1 = ncov - normal_read2
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
        extract_peptides('nonsyn_SQL_insert.txt', 'Formatted_epitope_variant.txt', sampleID)

        # Create epitopes
        create_epitopes('Formatted_epitope_variant.txt', 'SQL_Epitopes.txt', AA_DICT, cDNA_DICT)

        if mode in ['RNA', 'DNA-RNA']:
            reformat_gene_counts('tumor-gene.counts', 'Tumor_GeneCounts_SQL_insert.txt', sampleID, tumor_type)

        if mode in ['RNA', 'RNA-DNA']:
            reformat_gene_counts('normal-gene.counts', 'Normal_GeneCounts_SQL_insert.txt', sampleID, tumor_type)

    print("COMPLETED!")

parser = argparse.ArgumentParser(description='DNA/RNA Somatic variant calling and HLA prediction pipeline\n'
                                 'Created by Jose Fernandez <jc.fernandes.navarro@gmail.com>)',
                                 prog='somatic_pipeline.py',
                                 usage='somatic_pipeline.py [options] R1(Normal) R2(Normal) R1(Cancer) R2(Cancer)')
parser.add_argument('R1_NORMAL', help='FASTQ file R1 (Normal)')
parser.add_argument('R2_NORMAL', help='FASTQ file R2 (Normal)')
parser.add_argument('R1_CANCER', help='FASTQ file R1 (Cancer)')
parser.add_argument('R2_CANCER', help='FASTQ file R2 (Cancer)')
parser.add_argument('--genome',
                    help='Path to the reference genome FASTA file (must contain BWA index)', required=True)
parser.add_argument('--genome-star',
                    help='Path to the reference genome STAR index folder (RNA mode)', required=False)
parser.add_argument('--genome-ref',
                    help='Path to the reference genome GTF file (RNA mode)', required=False)
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
parser.add_argument('--dictAA',
                    help='Path to a dictionary of transcript IDs to peptide sequences', required=True)
parser.add_argument('--dictcDNA',
                    help='Path to a dictionary of transcript IDs to DNA sequences', required=True)
parser.add_argument('--annovar-db',
                    help='String indicated which Annovar database to use (default: humandb)',
                    default='humandb', required=False)
parser.add_argument('--annovar-version',
                    help='String indicated which version of the Annovar database to use (default: hg38)',
                    default='hg38', required=False)
parser.add_argument('--threads',
                    help='Number of threads to use in the parallel steps', type=int, default=10, required=False)
parser.add_argument('--steps', nargs='+', default=['mapping', 'gatk', 'hla', 'variant', 'filter'],
                    help='Steps to perform in the pipeline', 
                    choices=['mapping', 'gatk', 'hla', 'variant', 'filter', "none"])
parser.add_argument('--mode', default='DNA',
                    help='Mode to use (DNA (default), RNA, DNA-RNA and RNA-DNA)',
                    choices=['DNA', 'RNA', 'DNA-RNA', 'RNA-DNA'])

# Parse arguments
args = parser.parse_args()
DIR = args.dir
R1_NORMAL = os.path.abspath(args.R1_NORMAL)
R2_NORMAL = os.path.abspath(args.R2_NORMAL)
R1_CANCER = os.path.abspath(args.R1_CANCER)
R2_CANCER = os.path.abspath(args.R2_CANCER)
sampleID = args.sample
tumor_type = args.tumor
GENOME_REF = os.path.abspath(args.genome)
GENOME_REF_STAR = os.path.abspath(args.genome_star) if args.genome_star else None
GENOME_ANNOTATION = os.path.abspath(args.genome_ref) if args.genome_ref else None
THREADS = int(args.threads)
AA_DICT = os.path.abspath(args.dictAA)
cDNA_DICT = os.path.abspath(args.dictcDNA)
KNOWN_SITE1 = os.path.abspath(args.known1)
KNOWN_SITE2 = os.path.abspath(args.known2)
SNPSITES = os.path.abspath(args.snpsites)
GERMLINE = os.path.abspath(args.germline)
PON = os.path.abspath(args.pon)
STEPS = args.steps
ANNOVAR_DB = args.annovar_db
ANNOVAR_VERSION = args.annovar_version
MODE = args.mode
if 'RNA' in MODE and (not GENOME_REF_STAR or not GENOME_ANNOTATION):
    sys.stderr.write("Error, RNA mode but STAR reference or annotation files are missing\n")
    sys.exit(1)
    
# Move to output dir
os.makedirs(os.path.abspath(DIR), exist_ok=True)
os.chdir(os.path.abspath(DIR))

somatic_pipeline(R1_NORMAL,
                 R2_NORMAL,
                 R1_CANCER,
                 R2_CANCER,
                 tumor_type,
                 GENOME_REF,
                 GENOME_REF_STAR,
                 GENOME_ANNOTATION,
                 sampleID,
                 THREADS,
                 AA_DICT,
                 cDNA_DICT,
                 KNOWN_SITE1,
                 KNOWN_SITE2,
                 SNPSITES,
                 GERMLINE,
                 PON,
                 ANNOVAR_DB,
                 ANNOVAR_VERSION,
                 STEPS,
                 MODE)
