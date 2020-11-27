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
import shutil
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import os
import sys
import shutil
import glob
import multiprocessing
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
         STEPS):

    # TODO add sanity checks for the parameters
    # TODO better log info

    print("DNA somatic pipeline")

    # Sample 1 tumor, sample 2 normal
    sample1_ID = SAMPLEID + "_Tumor"
    sample2_ID = SAMPLEID + "_Normal"

    # Create sub-folder to store all results
    os.makedirs('workdir', exist_ok=True)
    os.chdir('workdir')

    if 'mapping' in STEPS:
        # TRIMMING
        print('Starting trimming')
        cmd = '{} --cores {} --fastqc --paired --basename normal {} {}'.format(TRIMGALORE, THREADS, R1_NORMAL, R2_NORMAL)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} --cores {} --fastqc --paired --basename tumor {} {}'.format(TRIMGALORE, THREADS, R1_TUMOR, R2_TUMOR)
        p2 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()

        # ALIGNMENT
        print('Starting alignment')

        # Normal (paired)
        cmd = '{} -t {} {} -R @RG\\tID:{}\\tPL:Illumina\\tLB:DNA\\tPU:{}\\tSM:{}\\tCN:{} normal_val_1.fq.gz normal_val_2.fq.gz | ' \
              '{} sort --threads {} > sample2_header.bam'.format(
                    BWA, THREADS, GENOME, sample1_ID, sample1_ID, sample1_ID, sample1_ID, SAMTOOLS, THREADS)
        p1 = exec_command(cmd, detach=True)

        # Tumor (paired)
        cmd = '{} -t {} {} -R "@RG\\tID:{}\\tPL:Illumina\\tLB:DNA\\tPU:{}\\tSM:{}\\tCN:{}" tumor_val_1.fq.gz tumor_val_2.fq.gz | ' \
              '{} sort --threads {} > sample1_header.bam'.format(
                    BWA, THREADS, GENOME, sample1_ID, sample1_ID, sample1_ID, sample1_ID, SAMTOOLS, THREADS)
        p2 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()

    if 'gatk' in STEPS:
        # Mark duplicates
        print('Marking duplicates')
        cmd = '{} MarkDuplicatesSpark --input sample1_header.bam --output sample1_dedup.bam'.format(GATK)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} MarkDuplicatesSpark --input sample2_header.bam --output sample2_dedup.bam'.format(GATK)
        p2 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()

        intervals_cmd = '--intervals {}'.format(INTERVALS) if INTERVALS else ''

        # GATK base re-calibration
        print('Starting re-calibration')
        cmd = '{} BaseRecalibratorSpark --input sample1_dedup.bam --reference {} --known-sites {} --known-sites {}' \
              ' --known-sites {} --output sample1_recal_data.txt {}'.format(GATK, GENOME, SNPSITES,
                                                                            KNOWN_SITE1, KNOWN_SITE2, intervals_cmd)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} BaseRecalibratorSpark --input sample2_dedup.bam --reference {} --known-sites {} --known-sites {}' \
              ' --known-sites {} --output sample2_recal_data.txt {}'.format(GATK, GENOME, SNPSITES,
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
        cmd = '{} -bam sample2_final.bam --genome-gc-distr HUMAN -nt {} ' \
              '-outdir bamQC_Normal -outformat HTML'.format(BAMQC, THREADS)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} -bam sample1_final.bam --genome-gc-distr HUMAN -nt {} ' \
              '-outdir bamQC_Tumor -outformat HTML'.format(BAMQC, THREADS)
        p2 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()

    if 'hla' in STEPS:
        # HLA-LA predictions
        print('Performing HLA-LA predictions')
        p1 = multiprocessing.Process(target=HLA_prediction,
                                    args=("sample2_final.bam", THREADS,
                                    "Normal", SAMPLEID, HLA_FASTA))
        p1.start()
        p2 = multiprocessing.Process(target=HLA_prediction,
                                    args=("sample1_final.bam", THREADS,
                                    "Tumor", SAMPLEID, HLA_FASTA))
        p2.start()

        # Wait for the processes to finish in parallel
        p1.join()
        p2.join()

    if 'variant' in STEPS:
        print('Performing variant calling Mutect2')
        intervals_cmd = '--intervals {}'.format(INTERVALS) if INTERVALS else ''
        # Variant calling Mutect2
        cmd = '{} Mutect2 --reference {} --input sample1_final.bam --input sample2_final.bam --normal-sample {} ' \
              '--output Mutect_unfiltered.vcf --germline-resource {} --dont-use-soft-clipped-bases ' \
              '--panel-of-normals {} {}'.format(GATK, GENOME, sample2_ID, GERMLINE, PON, intervals_cmd)
        p1 = exec_command(cmd, detach=True)

        # Variant calling Strelka2
        print('Performing variant calling with Strelka2')
        if os.path.isdir('Strelka_output'):
            shutil.rmtree(os.path.abspath('Strelka_output'))
        if INTERVALS:
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
        print('Performing variant calling with SomaticSniper')
        cmd = '{} -Q 15 -L -G -F vcf -f {} sample1_final.bam sample2_final.bam SS.vcf'.format(SSNIPER, GENOME)
        p3 = exec_command(cmd, detach=True)

        # Variant calling (Samtools pile-ups)
        print('Computing pile-ups')
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
        print('Performing variant calling with VarScan2')
        cmd = '{} somatic sample2.pileup sample1.pileup varscan --tumor-purity .5 --output-vcf 1 ' \
              '--min-coverage 4 --min-var-freq .05 --min-reads 2 --strand-filter 1'.format(VARSCAN)
        p6 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()
        p3.wait()
        p6.wait()

    if 'filter' in STEPS:
        print('Filtering variants')
        cmd = '{} FilterMutectCalls --variant Mutect_unfiltered.vcf --output Mutect.vcf --reference {}'.format(GATK,
                                                                                                               GENOME)
        exec_command(cmd)
        mutect2_filter('Mutect.vcf', 'mutect_filtered.vcf', sample1_ID, sample2_ID)
        strelka2_filter('Strelka_output/results/variants/somatic.snvs.vcf.gz', 'strelka_filtered.vcf')
        somaticSniper_filter('SS.vcf', 'somaticsniper_filtered.vcf')
        varscan_filter('varscan.snp.vcf', 'varscan_filtered.vcf')
        strelka2_filter_indels('Strelka_output/results/variants/somatic.indels.vcf.gz', 'strelka_indel_filtered.vcf')
        varscan_filter('varscan.indel.vcf', 'varscan_filtered_indel.vcf')

        # Combine with GATK
        print('Combining variants')
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        # TODO replace this with jacquard merge
        cmd = '{} -T CombineVariants -R {} -V:varscan_indel varscan_filtered_indel.vcf -V:varscan varscan_filtered.vcf ' \
              '-V:mutect mutect_filtered.vcf -V:strelka_indel strelka_indel_filtered.vcf -V:strelka strelka_filtered.vcf ' \
              '-V:somaticsniper somaticsniper_filtered.vcf -o combined_calls.vcf '\
              '-genotypeMergeOptions UNIQUIFY --num_threads {}'.format(GATK3, GENOME, THREADS)
        exec_command(cmd)

        # Annotate with Annovar
        print('Annotating variants')
        annotate_variants('combined_calls.vcf', 'annotated', ANNOVAR_DB, ANNOVAR_VERSION, THREADS)

        # Moving result files to output
        shutil.move('combined_calls.vcf', '../combined_calls.vcf')
        shutil.move('annotated.{}_multianno.vcf'.format(ANNOVAR_VERSION),
                    '../annotated.{}_multianno.vcf'.format(ANNOVAR_VERSION))
        shutil.move('Tumor_{}_hla_genotype_results.tsv'.format(SAMPLEID),
                    '../Tumor_hla_genotype.tsv')
        shutil.move('Normal_{}_hla_genotype_results.tsv'.format(SAMPLEID),
                    '../Normal_hla_genotype.tsv')
        shutil.move('sample1_final.bam', '../tumor_final.bam')
        shutil.move('sample2_final.bam', '../normal_final.bam')
        if os.path.isdir('../{}_bamQCNormal'.format(SAMPLEID)):
            shutil.rmtree(os.path.abspath('../{}_bamQCNormal'.format(SAMPLEID)))
        shutil.move('bamQC_Normal', '../{}_bamQCNormal'.format(SAMPLEID))
        if os.path.isdir('../{}_bamQCTumor'.format(SAMPLEID)):
            shutil.rmtree(os.path.abspath('../{}_bamQCTumor'.format(SAMPLEID)))
        shutil.move('bamQC_Tumor', '../{}_bamQCTumor'.format(SAMPLEID))
        for file in glob.glob('*_fastqc*'):
            shutil.move(file, '../{}_{}'.format(SAMPLEID, file))
        for file in glob.glob('*_trimming_report*'):
            shutil.move(file, '../{}_{}'.format(SAMPLEID, file))

    print('COMPLETED!')

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
    parser.add_argument("--hla-fasta",type=str,default="hla_reference_rna.fasta",required=False,
                        help="Path to the dna hla reference fasta file located in shared. (default: hla_reference_dna.fasta)")
    parser.add_argument('--threads',
                        help='Number of threads to use in the parallel steps', type=int, default=10, required=False)
    parser.add_argument('--steps', nargs='+', default=['mapping', 'gatk', 'hla', 'variant', 'filter'],
                        help='Steps to perform in the pipeline',
                        choices=['mapping', 'gatk', 'hla', 'variant', 'filter'])

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
    HLA_FASTA = args.hla_fasta

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
         STEPS)