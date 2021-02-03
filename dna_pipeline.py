#! /usr/bin/env python
"""
This pipeline computes somatic variants from DNA tumor-normal paired data.

The pipeline trims with trimgalore, aligns with bwa-men,
performs the GATK4 best practices and computes variants with
Mutect2, Strelka2, SomaticSniper and Varscan.
The variants are then combined into one file and annotated with Annovar.

Multiple options are available. To see them type --help

@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from hlapipeline.common import *
from hlapipeline.version import version_number
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os
import sys
import shutil
import glob
import multiprocessing
import logging
import datetime
from hlapipeline.filters import *


def main(R1_NORMAL,
         R2_NORMAL,
         R1_TUMOR,
         R2_TUMOR,
         GENOME,
         SAMPLEID,
         THREADS,
         KNOWN_SITE1,
         KNOWN_SITE2,
         SNPSITES,
         GERMLINE,
         PON,
         INTERVALS,
         ANNOVAR_DB,
         ANNOVAR_VERSION,
         HLA_FASTA,
         KEEP,
         STEPS,
         HLA_NORMAL,
         SPARK):

    # TODO add sanity checks for the parameters

    logging.basicConfig(format='%(asctime)s - %(message)s',
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.DEBUG, filename=SAMPLEID + ".log")
    logger = logging.getLogger(SAMPLEID)

    start_pipeline_time = datetime.datetime.now()

    logger.info('Starting DNA somatic pipeline: {}'.format(start_pipeline_time))
    logger.info('HLA Pipeline version: {}'.format(version_number))
    logger.info('Processing Normal FASTQs {} and {}; and Tumor FASTQs {} and {} with Sample ID {} ' \
                'using genome version {}.'.format(R1_NORMAL, R2_NORMAL, R1_TUMOR, R2_TUMOR, SAMPLEID, ANNOVAR_VERSION))

    # Sample 1 tumor, sample 2 normal
    sample1_ID = SAMPLEID + "_Tumor"
    sample2_ID = SAMPLEID + "_Normal"

    # Create sub-folder to store all results
    os.makedirs('workdir', exist_ok=True)
    os.chdir('workdir')

    if 'mapping' in STEPS:

        start_map_time = datetime.datetime.now()
        logger.info('Starting trimming and mapping step: {}'.format(start_map_time))
        # TRIMMING
        logger.info('Starting trimming')
        cmd = '{} --cores {} --fastqc --paired --basename normal {} {}'.format(TRIMGALORE, THREADS, R1_NORMAL,
                                                                               R2_NORMAL)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} --cores {} --fastqc --paired --basename tumor {} {}'.format(TRIMGALORE, THREADS, R1_TUMOR, R2_TUMOR)
        p2 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()

        # ALIGNMENT
        logger.info('Starting alignment')

        # Normal (paired)
        cmd = '{} -t {} {} -R "@RG\\tID:{}\\tPL:Illumina\\tLB:DNA\\tPU:{}\\tSM:{}\\tCN:{}" normal_val_1.fq.gz normal_val_2.fq.gz | ' \
              '{} sort --threads {} > sample2_header.bam'.format(
            BWA, THREADS, GENOME, sample2_ID, sample2_ID, sample2_ID, sample2_ID, SAMTOOLS, THREADS)
        p1 = exec_command(cmd, detach=True)

        # Tumor (paired)
        cmd = '{} -t {} {} -R "@RG\\tID:{}\\tPL:Illumina\\tLB:DNA\\tPU:{}\\tSM:{}\\tCN:{}" tumor_val_1.fq.gz tumor_val_2.fq.gz | ' \
              '{} sort --threads {} > sample1_header.bam'.format(
            BWA, THREADS, GENOME, sample1_ID, sample1_ID, sample1_ID, sample1_ID, SAMTOOLS, THREADS)
        p2 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()

        if not KEEP:
            if os.path.isfile('normal_val_1.fq.gz'):
                os.remove('normal_val_1.fq.gz')
            if os.path.isfile('normal_val_2.fq.gz'):
                os.remove('normal_val_2.fq.gz')
            if os.path.isfile('tumor_val_1.fq.gz'):
                os.remove('tumor_val_1.fq.gz')
            if os.path.isfile('tumor_val_2.fq.gz'):
                os.remove('tumor_val_2.fq.gz')

        end_map_time = datetime.datetime.now()
        total_map_time = end_map_time - start_map_time
        logger.info('Total trimming and mapping execution time: {}'.format(total_map_time))

    if 'gatk' in STEPS:

        start_gatk_time = datetime.datetime.now()
        logger.info('Starting GATK steps: {}'.format(start_gatk_time))

        # Mark duplicates

        logger.info('Marking Duplicates')

        if SPARK:
            cmd = '{} --java-options "-Xmx32g" MarkDuplicatesSpark -I sample1_header.bam -O sample1_dedup.bam'.format(GATK)
            p1 = exec_command(cmd, detach=True)

            cmd = '{} --java-options "-Xmx32g" MarkDuplicatesSpark -I sample2_header.bam -O sample2_dedup.bam'.format(GATK)
            p2 = exec_command(cmd, detach=True)

            # Wait for the processes to finish in parallel
            p1.wait()
            p2.wait()

        else:
            cmd = '{} --java-options "-Xmx32g" MarkDuplicates -I sample1_header.bam -O sample1_dedup.bam ' \
                      '--CREATE_INDEX true -M sample1_dup_metrics.txt'.format(GATK)
            p1 = exec_command(cmd, detach=True)

            cmd = '{} --java-options "-Xmx32g" MarkDuplicates -I sample2_header.bam -O sample2_dedup.bam ' \
                      '--CREATE_INDEX true -M sample2_dup_metrics.txt'.format(GATK)
            p2 = exec_command(cmd, detach=True)

            # Wait for the processes to finish in parallel
            p1.wait()
            p2.wait()
        
        intervals_cmd = '--intervals {}'.format(INTERVALS) if INTERVALS else ''

        # GATK base re-calibration

        logger.info('Starting re-calibration')

        recal_cmd = 'BaseRecalibratorSpark' if SPARK else 'BaseRecalibrator'

        cmd = '{} --java-options "-Xmx32g" {} --input sample1_dedup.bam --reference {} --known-sites {} --known-sites {}' \
              ' --known-sites {} --output sample1_recal_data.txt {}'.format(GATK, recal_cmd, GENOME, SNPSITES,
                                                                            KNOWN_SITE1, KNOWN_SITE2, intervals_cmd)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} --java-options "-Xmx32g" {} --input sample2_dedup.bam --reference {} --known-sites {} --known-sites {}' \
              ' --known-sites {} --output sample2_recal_data.txt {}'.format(GATK, recal_cmd, GENOME, SNPSITES,
                                                                            KNOWN_SITE1, KNOWN_SITE2, intervals_cmd)
        p2 = exec_command(cmd, detach=True)

        p1.wait()
        cmd = '{} ApplyBQSR --reference {} --input sample1_dedup.bam --bqsr-recal-file sample1_recal_data.txt ' \
              '--output sample1_final.bam'.format(GATK, GENOME)
        exec_command(cmd)

        p2.wait()
        cmd = '{} ApplyBQSR --reference {} --input sample2_dedup.bam --bqsr-recal-file sample2_recal_data.txt ' \
              '--output sample2_final.bam'.format(GATK, GENOME)
        exec_command(cmd)

        # BamQC
        cmd = '{} -bam sample2_final.bam --genome-gc-distr HUMAN -nt {} --java-mem-size=16G ' \
              '-outdir bamQC_Normal -outformat HTML'.format(BAMQC, THREADS)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} -bam sample1_final.bam --genome-gc-distr HUMAN -nt {} --java-mem-size=16G ' \
              '-outdir bamQC_Tumor -outformat HTML'.format(BAMQC, THREADS)
        p2 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()

        if not KEEP:
            if os.path.isfile('sample1_dedup.bam'):
                os.remove('sample1_dedup.bam')
            if os.path.isfile('sample2_dedup.bam'):
                os.remove('sample2_dedup.bam')
            if os.path.isfile('sample1_dedup.bam.bai'):
                os.remove('sample1_dedup.bam.bai')
            if os.path.isfile('sample2_dedup.bam.bai'):
                os.remove('sample2_dedup.bam.bai')
            if os.path.isfile('sample1_recal_data.txt'):
                os.remove('sample1_recal_data.txt')
            if os.path.isfile('sample2_recal_data.txt'):
                os.remove('sample2_recal_data.txt')

        end_gatk_time = datetime.datetime.now()
        total_gatk_time = end_gatk_time - start_gatk_time
        logger.info('Total GATK processing time: {}'.format(total_gatk_time))

    if 'hla' in STEPS:

        start_hla_time = datetime.datetime.now()
        logger.info('Starting HLA prediction: {}'.format(start_hla_time))
        # HLA-LA predictions
        if HLA_NORMAL:
            p1 = multiprocessing.Process(target=HLA_prediction,
                                        args=('sample2_final.bam', THREADS,
                                            'Normal', SAMPLEID, HLA_FASTA, 'dna', KEEP))
            p1.start()

        p2 = multiprocessing.Process(target=HLA_prediction,
                                     args=('sample1_final.bam', THREADS,
                                           'Tumor', SAMPLEID, HLA_FASTA, 'dna', KEEP))
        p2.start()
        p2.join()

        if HLA_NORMAL:
            p1.join()

        end_hla_time = datetime.datetime.now()
        total_hla_time = end_hla_time - start_hla_time
        logger.info('Total HLA prediction time: {}'.format(total_hla_time))

    if 'variant' in STEPS:

        start_variant_time = datetime.datetime.now()
        logger.info('Starting variant calling: {}'.format(start_variant_time))
        logger.info('Performing variant calling Mutect2')

        intervals_cmd = '--intervals {}'.format(INTERVALS) if INTERVALS else ''

        # Variant calling Mutect2
        cmd = '{} Mutect2 --reference {} --input sample1_final.bam --input sample2_final.bam --normal-sample {} ' \
              '--output Mutect_unfiltered.vcf --germline-resource {} --dont-use-soft-clipped-bases ' \
              '--panel-of-normals {} {}'.format(GATK, GENOME, sample2_ID, GERMLINE, PON, intervals_cmd)
        p1 = exec_command(cmd, detach=True)

        # Variant calling Strelka2
        logger.info('Performing variant calling with Strelka2')
        if os.path.isdir('Strelka_output'):
            shutil.rmtree(os.path.abspath('Strelka_output'))
        if INTERVALS:
            if os.path.isfile('intervals.bed.gz'):
                os.remove('intervals.bed.gz')
            if os.path.isfile('intervals.bed.gz.tbi'):
                os.remove('intervals.bed.gz.tbi')
            cmd = 'cp {} intervals.bed && bgzip intervals.bed && tabix intervals.bed.gz'.format(INTERVALS)
            exec_command(cmd)
            intervals_cmd = '--exome --callRegions intervals.bed.gz'
        else:
            intervals_cmd = ''
        cmd = '{} {} --normalBam sample2_final.bam --tumorBam sample1_final.bam --referenceFasta {}' \
              ' --runDir Strelka_output'.format(STRELKA, intervals_cmd, GENOME)
        exec_command(cmd)
        cmd = 'Strelka_output/runWorkflow.py -m local -j {}'.format(THREADS)
        p2 = exec_command(cmd, detach=True)

        # Variant calling Somatic Sniper
        logger.info('Performing variant calling with SomaticSniper')
        cmd = '{} -Q 15 -L -G -F vcf -f {} sample1_final.bam sample2_final.bam SS.vcf'.format(SSNIPER, GENOME)
        p3 = exec_command(cmd, detach=True)

        # Variant calling (Samtools pile-ups)
        logger.info('Computing pile-ups')
        intervals_cmd = '--positions {}'.format(INTERVALS) if INTERVALS else ''
        cmd = '{} mpileup -C 50 -B -q 1 -Q 15 {} -f {} sample1_final.bam > sample1.pileup'.format(SAMTOOLS,
                                                                                                  intervals_cmd,
                                                                                                  GENOME)
        p4 = exec_command(cmd, detach=True)

        cmd = '{} mpileup -C 50 -B -q 1 -Q 15 {} -f {} sample2_final.bam > sample2.pileup'.format(SAMTOOLS,
                                                                                                  intervals_cmd,
                                                                                                  GENOME)
        p5 = exec_command(cmd, detach=True)

        # Variant calling VarScan
        p4.wait()
        p5.wait()
        logger.info('Performing variant calling with VarScan2')
        cmd = '{} somatic sample2.pileup sample1.pileup varscan --tumor-purity .5 --output-vcf 1 ' \
              '--min-coverage 4 --min-var-freq .05 --min-reads 2 --strand-filter 1'.format(VARSCAN)
        p6 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()
        p3.wait()
        p6.wait()

        if not KEEP:
            if os.path.isfile('sample1.pileup'):
                os.remove('sample1.pileup')
            if os.path.isfile('sample2.pileup'):
                os.remove('sample2.pileup')
        
        end_variant_time = datetime.datetime.now()
        total_variant_time = end_variant_time - start_variant_time
        logger.info('Total variant calling processing time: {}'.format(total_variant_time))

    if 'filter' in STEPS:
        start_filter_time = datetime.datetime.now()
        logger.info('Starting variant filtering and annotation: {}'.format(start_filter_time))
        logger.info('Filtering variants')
        cmd = '{} FilterMutectCalls --variant Mutect_unfiltered.vcf --stats Mutect_unfiltered.vcf.stats ' \
              '--output Mutect.vcf --reference {}'.format(GATK, GENOME)
        exec_command(cmd)

        mutect2_filter('Mutect.vcf', 'mutect_filtered.vcf', sample1_ID, sample2_ID)
        strelka2_filter('Strelka_output/results/variants/somatic.snvs.vcf.gz', 'strelka_filtered.vcf')
        somaticSniper_filter('SS.vcf', 'somaticsniper_filtered.vcf')
        varscan_filter('varscan.snp.vcf', 'varscan_filtered.vcf')
        strelka2_filter_indels('Strelka_output/results/variants/somatic.indels.vcf.gz', 'strelka_indel_filtered.vcf')
        varscan_filter('varscan.indel.vcf', 'varscan_filtered_indel.vcf')

        # Combine with GATK
        logger.info('Combining variants')
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        # TODO replace this with jacquard merge
        cmd = '{} -T CombineVariants -R {} -V:varscan_indel varscan_filtered_indel.vcf -V:varscan varscan_filtered.vcf ' \
              '-V:mutect mutect_filtered.vcf -V:strelka_indel strelka_indel_filtered.vcf -V:strelka strelka_filtered.vcf ' \
              '-V:somaticsniper somaticsniper_filtered.vcf -o combined_calls.vcf ' \
              '-genotypeMergeOptions UNIQUIFY --num_threads {}'.format(GATK3, GENOME, THREADS)
        exec_command(cmd)

        # Annotate with Annovar
        logger.info('Annotating variants')
        annotate_variants('combined_calls.vcf', 'annotated', ANNOVAR_DB, ANNOVAR_VERSION, THREADS)
        # Replace UTF-8 code to equivalent characters
        cmd = "sed -i -e 's/{}{}/-/g' -e 's/{}{}/:/g' annotated.{}_multianno.vcf".format("\\", "\\x3b", "\\", "\\x3d",
                                                                                         ANNOVAR_VERSION)
        exec_command(cmd)

        # Summary of basic statistic of annotated VCF file
        annotated_vcf = "annotated.{}_multianno.vcf".format(ANNOVAR_VERSION)
        vcf_stats(annotated_vcf, SAMPLEID)

        end_filer_time = datetime.datetime.now()
        total_filter_time = end_filer_time - start_filter_time
        logger.info('Total filtering and annotation time: {}'.format(total_filter_time))

        # Moving result files to output
        if os.path.isfile('combined_calls.vcf'):
            shutil.move('combined_calls.vcf', '../combined_calls.vcf')
        if os.path.isfile('{}.relatedness2'.format(SAMPLEID)):
            shutil.move('{}.relatedness2'.format(SAMPLEID), '../{}.relatedness2'.format(SAMPLEID))
        if os.path.isfile('{}.TsTv.summary'.format(SAMPLEID)):
            shutil.move('{}.TsTv.summary'.format(SAMPLEID), '../{}.TsTv.summary'.format(SAMPLEID))
        if os.path.isfile('{}.vchk'.format(SAMPLEID)):
            shutil.move('{}.vchk'.format(SAMPLEID), '../{}.vchk'.format(SAMPLEID))
        if os.path.isfile('annotated.{}_multianno.vcf'.format(ANNOVAR_VERSION)):
            shutil.move('annotated.{}_multianno.vcf'.format(ANNOVAR_VERSION),
                        '../annotated.{}_multianno.vcf'.format(ANNOVAR_VERSION))
        if os.path.isfile('Tumor_{}_hla_genotype_result.tsv'.format(SAMPLEID)):
            shutil.move('Tumor_{}_hla_genotype_result.tsv'.format(SAMPLEID),
                        '../Tumor_hla_genotype.tsv')
        if os.path.isfile('Normal_{}_hla_genotype_result.tsv'.format(SAMPLEID)):
            shutil.move('Normal_{}_hla_genotype_result.tsv'.format(SAMPLEID),
                        '../Normal_hla_genotype.tsv')
        if os.path.isfile('sample1_final.bam'):
            shutil.move('sample1_final.bam', '../tumor_final.bam')
        if os.path.isfile('sample2_final.bam'):
            shutil.move('sample2_final.bam', '../normal_final.bam')
        if os.path.isdir('../{}_bamQCNormal'.format(SAMPLEID)):
            shutil.rmtree(os.path.abspath('../{}_bamQCNormal'.format(SAMPLEID)))
        if os.path.isdir('bamQC_Normal'):
            shutil.move('bamQC_Normal', '../{}_bamQCNormal'.format(SAMPLEID))
        if os.path.isdir('../{}_bamQCTumor'.format(SAMPLEID)):
            shutil.rmtree(os.path.abspath('../{}_bamQCTumor'.format(SAMPLEID)))
        if os.path.isdir('bamQC_Tumor'):
            shutil.move('bamQC_Tumor', '../{}_bamQCTumor'.format(SAMPLEID))
        for file in glob.glob('*_fastqc*'):
            shutil.move(file, '../{}_{}'.format(SAMPLEID, file))
        for file in glob.glob('*_trimming_report*'):
            shutil.move(file, '../{}_{}'.format(SAMPLEID, file))

    end_pipeline_time = datetime.datetime.now()
    total_pipeline_time = end_pipeline_time - start_pipeline_time
    logger.info('Total pipeline execution time: {}'.format(total_pipeline_time))
    
    logger.info('COMPLETED!')


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('R1_NORMAL', help='FASTQ file R1 (Normal)')
    parser.add_argument('R2_NORMAL', help='FASTQ file R2 (Normal)')
    parser.add_argument('R1_TUMOR', help='FASTQ file R1 (Tumor)')
    parser.add_argument('R2_TUMOR', help='FASTQ file R2 (Tumor)')
    parser.add_argument('--genome', type=str, required=True,
                        help='Path to the reference genome FASTA file (must contain BWA index)')
    parser.add_argument('--sample', type=str,
                        help='Name of the sample/experiment. Default is sample', default='sample')
    parser.add_argument('--outdir', type=str, required=True,
                        help='Path to the output folder where output files will be placed')
    parser.add_argument('--known1', type=str, required=True,
                        help='Path to the file with Mill and 1000G gold standards (GATK bundle)')
    parser.add_argument('--known2', type=str, required=True,
                        help='Path to the file with 1000G phase indels (GATK bundle)')
    parser.add_argument('--snpsites', type=str, required=True,
                        help='Path to the file with the SNPs (GATK bundle dbSNP)')
    parser.add_argument('--germline', type=str, required=True,
                        help='Path to the file with the germline resources Nomad for Mutect2 (GATK bundle)')
    parser.add_argument('--pon', type=str, required=True,
                        help='Path to the file with the panel of normals for Mutect2 (GATK bundle)')
    parser.add_argument('--intervals', type=str, default=None, required=False,
                        help='Path to the file with the intervals to operate in BaseRecalibrator and Mutect2 (BED)')
    parser.add_argument('--annovar-db', type=str, default='humandb', required=False,
                        help='String indicated which Annovar database to use (default: humandb)')
    parser.add_argument('--annovar-version', type=str, default='hg38', required=False,
                        help='String indicated which version of the Annovar database to use (default: hg38)')
    parser.add_argument("--hla-fasta", type=str, default=None, required=True,
                        help="Path to the HLA reference fasta file located for Optype.")
    parser.add_argument('--threads',
                        help='Number of threads to use in the parallel steps', type=int, default=10, required=False)
    parser.add_argument('--steps', nargs='+', default=['mapping', 'gatk', 'hla', 'variant', 'filter'],
                        help='Steps to perform in the pipeline',
                        choices=['mapping', 'gatk', 'hla', 'variant', 'filter'])
    parser.add_argument('--keep-intermediate', default=False, action='store_true', required=False,
                        help='Avoid intermediate files from being removed.')
    parser.add_argument('--normal-hla', default=False, action='store_true', required=False,
                        help='Perform HLA typing also in normal sample.')
    parser.add_argument('--use-gatk-spark', default=False, action='store_true', required=False,
                        help='Enable the use of MarkDuplicatesSpark and BaseRecalibratorSpark.')

    # Parse arguments
    args = parser.parse_args()
    DIR = args.outdir
    R1_NORMAL = os.path.abspath(args.R1_NORMAL)
    R2_NORMAL = os.path.abspath(args.R2_NORMAL)
    R1_TUMOR = os.path.abspath(args.R1_TUMOR)
    R2_TUMOR = os.path.abspath(args.R2_TUMOR)
    SAMPLEID = args.sample
    GENOME_REF = os.path.abspath(args.genome)
    THREADS = int(args.threads)
    KNOWN_SITE1 = os.path.abspath(args.known1)
    KNOWN_SITE2 = os.path.abspath(args.known2)
    SNPSITES = os.path.abspath(args.snpsites)
    GERMLINE = os.path.abspath(args.germline)
    PON = os.path.abspath(args.pon)
    INTERVALS = os.path.abspath(args.intervals) if args.intervals else None
    STEPS = args.steps
    ANNOVAR_DB = args.annovar_db
    ANNOVAR_VERSION = args.annovar_version
    HLA_FASTA = os.path.abspath(args.hla_fasta)
    KEEP = args.keep_intermediate
    HLA_NORMAL = args.normal_hla
    SPARK =  args.use_gatk_spark

    # Move to output dir
    os.makedirs(os.path.abspath(DIR), exist_ok=True)
    os.chdir(os.path.abspath(DIR))

    main(R1_NORMAL,
         R2_NORMAL,
         R1_TUMOR,
         R2_TUMOR,
         GENOME_REF,
         SAMPLEID,
         THREADS,
         KNOWN_SITE1,
         KNOWN_SITE2,
         SNPSITES,
         GERMLINE,
         PON,
         INTERVALS,
         ANNOVAR_DB,
         ANNOVAR_VERSION,
         HLA_FASTA,
         KEEP,
         STEPS,
         HLA_NORMAL,
         SPARK)
