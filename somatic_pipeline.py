#! /usr/bin/env python
"""
This pipeline computes somatic variants from DNA or RNA tumor-normal paired data.

The pipline trims with trimgalore, aligns with STAR or bwa-men,
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
from hlapipeline.filters import *

def main(R1_NORMAL,
         R2_NORMAL,
         R1_CANCER,
         R2_CANCER,
         genome,
         genome_star,
         annotation,
         sampleID,
         THREADS,
         KNOWN_SITE1,
         KNOWN_SITE2,
         SNPSITES,
         GERMLINE,
         PON,
         ANNOVAR_DB,
         ANNOVAR_VERSION,
         steps,
         mode):

    # TODO add sanity checks for the parameters
    # TODO better log info
    # TODO remove temp files
    # TODO put output files somewhere else

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

        # Normal (paired)
        cmd = '{} -t {} {} normal_val_1.fq.gz normal_val_2.fq.gz | ' \
              '{} sort --threads {} > aligned_normal_merged.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
        exec_command(cmd)

        if mode == 'RNA':
            # Cancer (paired)
            cmd = '{} --genomeDir {} --readFilesIn cancer_val_1.fq.gz cancer_val_2.fq.gz --outSAMorder Paired' \
                  ' --twopassMode None --outSAMunmapped None --outSAMtype BAM SortedByCoordinate' \
                  ' --readFilesCommand gunzip -c --runThreadN {}'.format(STAR, genome_star, THREADS)
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
        cmd = '{} AddOrReplaceReadGroups -I aligned_cancer_merged.bam -O sample1_header.bam -RGID {} -RGPL Illumina ' \
              '-RGLB {} -RGPU {} -RGSM {} -RGCN VHIO'.format(PICARD, sample1_ID, mode, sample1_ID, sample1_ID)
        exec_command(cmd)

        cmd = '{} AddOrReplaceReadGroups -I aligned_normal_merged.bam -O sample2_header.bam -RGID {} -RGPL Illumina ' \
              '-RGLB {} -RGPU {} -RGSM {} -RGCN VHIO'.format(PICARD, sample2_ID, mode, sample2_ID, sample2_ID)
        exec_command(cmd)

    if 'gatk' in steps:
        # Mark duplicates
        print('Marking duplicates')
        # NOTE setting reducers to it works in system that do not allow many files open
        cmd = '{} MarkDuplicatesSpark --input sample1_header.bam --output sample1_dedup.bam'.format(GATK)
        exec_command(cmd)
        cmd = '{} MarkDuplicatesSpark --input sample2_header.bam --output sample2_dedup.bam'.format(GATK)
        exec_command(cmd)

        # Split N and cigards for RNA data
        if mode == 'RNA':
            cmd = '{} SplitNCigarReads --reference {} --input sample1_dedup.bam --output sample1_split.bam'.format(GATK,
                                                                                                                   genome)
            exec_command(cmd)

            cmd = 'rm -rf sample1_dedup* && mv sample1_split.bam sample1_dedup.bam && mv sample1_split.bai sample1_dedup.bai'
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
        HLA_predictionDNA('sample2_final.bam', sampleID, 'PRG-HLA-LA_Normal_output.txt', THREADS)
        if mode == 'DNA':
            HLA_predictionDNA('sample1_final.bam', sampleID, 'PRG-HLA-LA_Tumor_output.txt', THREADS)
        else:
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
        cmd = '{} somatic sample2.pileup sample1.pileup varscan --tumor-purity .5 --output-vcf 1 ' \
              '--min-coverage 4 --min-var-freq .05 --min-reads 2 --strand-filter 1'.format(VARSCAN)
        exec_command(cmd)

        if mode in ['RNA']:
            # Computing gene counts
            print('Computing gene counts with featureCounts for tumor sample')
            cmd = '{} -T {} --primary --ignoreDup -O -C -t exon ' \
                  '-g gene_name -a {} -o tumor-gene.counts sample1_dedup.bam'.format(FEATURECOUNTS, THREADS, annotation)
            exec_command(cmd)

    if 'filter' in steps:
        print('Filtering variants')
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
              '-genotypeMergeOptions UNIQUIFY --num_threads {}'.format(GATK3, genome, THREADS)
        exec_command(cmd)

        # Annotate with Annovar
        annovardb = '{} -buildver {}'.format(os.path.join(ANNOVAR_PATH, ANNOVAR_DB), ANNOVAR_VERSION)
        print('Running Annovar')
        cmd = '{} combined_calls.vcf {} -thread {} -out annotated -vcfinput -remove -protocol {}'.format(
            os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), annovardb, THREADS,  annovar_anno)
        exec_command(cmd)

    print("COMPLETED!")

if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
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
                        help='DNA if tumor sample is from DNA or else RNA if tumor sample is from RNA [DNA].',
                        choices=['DNA', 'RNA'])

    # Parse arguments
    args = parser.parse_args()
    DIR = args.dir
    R1_NORMAL = os.path.abspath(args.R1_NORMAL)
    R2_NORMAL = os.path.abspath(args.R2_NORMAL)
    R1_CANCER = os.path.abspath(args.R1_CANCER)
    R2_CANCER = os.path.abspath(args.R2_CANCER)
    sampleID = args.sample
    GENOME_REF = os.path.abspath(args.genome)
    GENOME_REF_STAR = os.path.abspath(args.genome_star) if args.genome_star else None
    GENOME_ANNOTATION = os.path.abspath(args.genome_ref) if args.genome_ref else None
    THREADS = int(args.threads)
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

    main(R1_NORMAL,
         R2_NORMAL,
         R1_CANCER,
         R2_CANCER,
         GENOME_REF,
         GENOME_REF_STAR,
         GENOME_ANNOTATION,
         sampleID,
         THREADS,
         KNOWN_SITE1,
         KNOWN_SITE2,
         SNPSITES,
         GERMLINE,
         PON,
         ANNOVAR_DB,
         ANNOVAR_VERSION,
         STEPS,
         MODE)
