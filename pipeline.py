#! /usr/bin/env python

import argparse
from common import *
from Full_exome_pipeline import *
from HLA_two_sample import *
import multiprocessing

parser=argparse.ArgumentParser(description='Jared pipeline (adjusted by Jose Fernandez <jc.fernandes.navarro@gmail.com>',
                               prog='pipeline.py',
                               usage='pipeline.py [options] R1(Normal) R2(Normal) R1(Cancer) R2(Cancer)')
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
                    help='Tumor type. Default is NA', default='NA')
parser.add_argument('--dir',
                    help='Path to the output file', required=True)
parser.add_argument('--known1',
                    help='Path to the file with Mill and 1000G gold standards (GATK bundle)', required=True)
parser.add_argument('--known2',
                    help='Path to the file with 1000G phase indels (GATK bundle)', required=True)
parser.add_argument('--snpsites',
                    help='Path to the file with the SNPs (GATK buldle)', required=True)
parser.add_argument('--germline',
                    help='Path to the file with the germline resources Nomad (GATK buldle)', required=True)
parser.add_argument('--fastaAA',
                    help='Path to the file with the dictionary of FASTA to AA', required=True)
parser.add_argument('--fastacDNA',
                    help='Path to the file with the dictionary of FASTA to cDNA', required=True)

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
FASTA_cDNA = os.path.abspath(args.fastacDNA)
KNOWN_SITE1 = os.path.abspath(args.known1)
KNOWN_SITE2 = os.path.abspath(args.known2)
SNPSITES = os.path.abspath(args.snpsites)
GERMLINE = os.path.abspath(args.germline)

# Recommend to install with Anaconda
TRIPTOMATIC = 'trimmomatic'
BWA = 'bwa mem'

# Move to output dir
os.makedirs(os.path.abspath(DIR), exist_ok=True)
os.chdir(os.path.abspath(DIR))

# TRIMMING
print('Starting trimming')

cmd = '{} PE -threads {} -phred33 {} {} R1_normal.fastq.gz R1_normal_unpaired.fastq.gz '\
      'R2_normal.fastq.gz R2_normal_unpaired.fastq.gz '\
      'ILLUMINACLIP:{}:2:40:15 HEADCROP:9 CROP:140 SLIDINGWINDOW:4:25 MINLEN:5'.format(TRIPTOMATIC,
                                                                                       THREADS,
                                                                                       R1_NORMAL,
                                                                                       R2_NORMAL,
                                                                                       IILLUMINA_ADAPTERS)
exec_command(cmd)

cmd = '{} PE -threads {} -phred33 {} {} R1_cancer.fastq.gz R1_cancer_unpaired.fastq.gz '\
      'R2_cancer.fastq.gz R2_cancer_unpaired.fastq.gz '\
      'ILLUMINACLIP:{}:2:40:15 HEADCROP:9 CROP:140 SLIDINGWINDOW:4:25 MINLEN:5'.format(TRIPTOMATIC,
                                                                                        THREADS,
                                                                                        R1_CANCER,
                                                                                        R2_CANCER,
                                                                                        IILLUMINA_ADAPTERS)
exec_command(cmd)

print('Trimming of tumor and normal samples completed.')

# ALIGNMENT
print('Starting alignment')

# Normal (paired)
cmd = '{} -t {} {} R1_normal.fastq.gz R2_normal.fastq.gz | {} view -bS > normal_paired_aligned.bam'.format(BWA,
                                                                                                           THREADS,
                                                                                                           GENOME_REF,
                                                                                                           SAMTOOLS)
exec_command(cmd)

cmd = '{} sort --threads {} normal_paired_aligned.bam > normal_paired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                 THREADS)
exec_command(cmd)

# Cancer (paired)
cmd = '{} -t {} {} R1_cancer.fastq.gz R2_cancer.fastq.gz | {} view -bS > cancer_paired_aligned.bam'.format(BWA,
                                                                                                           THREADS,
                                                                                                           GENOME_REF,
                                                                                                           SAMTOOLS)
exec_command(cmd)

cmd = '{} sort --threads {} cancer_paired_aligned.bam > cancer_paired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                 THREADS)
exec_command(cmd)

# Normal (unpaired R1)
cmd = '{} -t {} {} R1_normal_unpaired.fastq.gz | {} view -bS > R1_normal_unpaired_aligned.bam'.format(BWA,
                                                                                                      THREADS,
                                                                                                      GENOME_REF,
                                                                                                      SAMTOOLS)
exec_command(cmd)

cmd = '{} sort --threads {} R1_normal_unpaired_aligned.bam > R1_normal_unpaired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                           THREADS)
exec_command(cmd)

# Cancer (unpaired R1)
cmd = '{} -t {} {} R1_cancer_unpaired.fastq.gz | {} view -bS > R1_cancer_unpaired_aligned.bam'.format(BWA,
                                                                                                      THREADS,
                                                                                                      GENOME_REF,
                                                                                                      SAMTOOLS)
exec_command(cmd)

cmd = '{} sort --threads {} R1_cancer_unpaired_aligned.bam > R1_cancer_unpaired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                           THREADS)
exec_command(cmd)

# Normal (unpaired R2)
cmd = '{} -t {} {} R2_normal_unpaired.fastq.gz | {} view -bS > R2_normal_unpaired_aligned.bam'.format(BWA,
                                                                                                      THREADS,
                                                                                                      GENOME_REF,
                                                                                                      SAMTOOLS)
exec_command(cmd)

cmd = '{} sort --threads {} R2_normal_unpaired_aligned.bam > R2_normal_unpaired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                           THREADS)
exec_command(cmd)

# Cancer (unpaired R2)
cmd = '{} -t {} {} R2_cancer_unpaired.fastq.gz | {} view -bS > R2_cancer_unpaired_aligned.bam'.format(BWA,
                                                                                                      THREADS,
                                                                                                      GENOME_REF,
                                                                                                      SAMTOOLS)
exec_command(cmd)

cmd = '{} sort --threads {} R2_cancer_unpaired_aligned.bam > R2_cancer_unpaired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                           THREADS)
exec_command(cmd)

print('Aligment of tumor and normal samples completed.')

# Merge aligned files
print('Merging aligned files')

cmd = '{} merge -f aligned_normal_merged.bam normal_paired_aligned_sorted.bam '\
       'R1_normal_unpaired_aligned_sorted.bam R2_normal_unpaired_aligned_sorted.bam'.format(SAMTOOLS)
exec_command(cmd)

cmd = '{} merge -f aligned_cancer_merged.bam cancer_paired_aligned_sorted.bam '\
       'R1_cancer_unpaired_aligned_sorted.bam R2_cancer_unpaired_aligned_sorted.bam'.format(SAMTOOLS)
exec_command(cmd)

print('Merging of tumor and normal aligned samples completed.')

# Final pìpeline
Full_exome_pipeline('aligned_cancer_merged.bam',
                    'aligned_normal_merged.bam',
                    tumor_type,
                    GENOME_REF,
                    sampleID,
                    THREADS,
                    FASTA_AA_DICT,
                    FASTA_cDNA,
                    KNOWN_SITE1,
                    KNOWN_SITE2,
                    SNPSITES,
                    GERMLINE)
#HLA_pipeline(loc, sample1, sample2, THREADS)


