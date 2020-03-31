#! /usr/bin/env python

import argparse
import os
import multiprocessing
from Full_exome_pipeline import *
from HLA_two_sample import *

parser=argparse.ArgumentParser(description='Jared pipeline (adjusted by Jose Fernandez <jc.fernandes.navarro@gmail.com>',
                               prog='pipeline.py',
                               usage='pipeline.py [options] R1(Normal) R2(Normal) R1(Cancer) R2(Cancer)')
parser.add_argument('R1_NORMAL', help='FASTQ file R1 (Normal)')
parser.add_argument('R2_NORMAL', help='FASTQ file R2 (Normal)')
parser.add_argument('R1_CANCER', help='FASTQ file R1 (Cancer)')
parser.add_argument('R2_CANCER', help='FASTQ file R2 (Cancer)')
parser.add_argument('-a', '--adapter',
                    help='Path to the Illumina adapters FASTA file. Default TruSeq2-PE.fa', default='TruSeq2-PE.fa')
parser.add_argument('-s', '--sample',
                    help='Name of the sample/experiment. Default is sample', default='sample')
parser.add_argument('-t', '--tumor',
                    help='Tumor type. Default is NA', default='NA')
parser.add_argument('-D', '--dir',
                    help='Path to a location to store all results. Default is ./data', default='./data')

# Parse arguments
args = parser.parse_args()
DIR = args.dir
R1_NORMAL = args.R1_NORMAL
R2_NORMAL = args.R2_NORMAL
R1_CANCER = args.R1_CANCER
R2_CANCER = args.R2_CANCER
sampleID = args.sample
tumor_type = args.tumor
IILLUMINA_ADAPTERS = os.path.abspath(args.adapter)
THREADS = multiprocessing.cpu_count() - 1

# Assumed to be in ~/shared/ for convenience
# BWA index must be present here
GENOME_REF = "~/shared/hg19.fa"

# These must be installed in the system or in PATH
TRIPTOMATIC = 'java -jar ~/shared/trimmomatic.jar'
BWA = '~/shared/bwa mem'
SAMTOOLS = 'samtools'

# Move to output dir
os.makedirs(os.path.abspath(DIR), exist_ok=True)
os.chdir(os.path.abspath(DIR))

# TRIMMING
print('Starting trimming')
cmd1 = '{} PE -threads {} -phred33 {} {} R1_normal.fastq.gz R1_normal_unpaired.fastq.gz '\
       'R2_normal.fastq.gz R2_normal_unpaired.fastq.gz '\
       'ILLUMINACLIP:{}:2:40:15 HEADCROP:9 CROP:140 SLIDINGWINDOW:4:25 MINLEN:5'.format(TRIPTOMATIC,
                                                                                        THREADS,
                                                                                        R1_NORMAL,
                                                                                        R2_NORMAL,
                                                                                        IILLUMINA_ADAPTERS)
cmd2 = '{} PE -threads {} -phred33 {} {} R1_cancer.fastq.gz R1_cancer_unpaired.fastq.gz '\
       'R2_cancer.fastq.gz R2_cancer_unpaired.fastq.gz '\
       'ILLUMINACLIP:{}:2:40:15 HEADCROP:9 CROP:140 SLIDINGWINDOW:4:25 MINLEN:5'.format(TRIPTOMATIC,
                                                                                        THREADS,
                                                                                        R1_CANCER,
                                                                                        R2_CANCER,
                                                                                        IILLUMINA_ADAPTERS)
p1 = subprocess.Popen(cmd1, shell=True)
p2 = subprocess.Popen(cmd2, shell=True)
p1.wait()
p2.wait()
print('Trimming of tumor and normal samples completed.')

# ALIGNMENT

# Normal (paired)
cmd1 = '{} -H -t {} {} R1_normal.fastq.gz R2_normal.fastq.gz | {} view -bS > normal_paired_aligned.bam'.format(BWA,
                                                                                                               THREADS,
                                                                                                               GENOME_REF,
                                                                                                               SAMTOOLS)
cm2 = '{} sort -@ {} normal_paired_aligned.bam > normal_paired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                          THREADS)
p1 = subprocess.Popen(cmd1 + '; wait ; ' + cmd2, shell=True)
p1.wait()
# Cancer (paired)
cmd1 = '{} -H -t {} {} R1_cancer.fastq.gz R2_cancer.fastq.gz | {} view -bS > cancer_paired_aligned.bam'.format(BWA,
                                                                                                               THREADS,
                                                                                                               GENOME_REF,
                                                                                                               SAMTOOLS)
cm2 = '{} sort -@ {} cancer_paired_aligned.bam > cancer_paired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                          THREADS)
p1 = subprocess.Popen(cmd1 + '; wait ; ' + cmd2, shell=True)
p1.wait()

# Normal (unpaired R1)
cmd1 = '{} -H -t {} {} R1_normal_unpaired.fastq.gz | {} view -bS > R1_normal_unpaired_aligned.bam'.format(BWA,
                                                                                                          THREADS,
                                                                                                          GENOME_REF,
                                                                                                          SAMTOOLS)
cm2 = '{} sort -@ {} R1_normal_unpaired_aligned.bam > R1_normal_unpaired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                    THREADS)
p1 = subprocess.Popen(cmd1 + '; wait ; ' + cmd2, shell=True)
p1.wait()
# Cancer (unpaired R1)
cmd1 = '{} -H -t {} {} R1_cancer_unpaired.fastq.gz | {} view -bS > R1_cancer_unpaired_aligned.bam'.format(BWA,
                                                                                                          THREADS,
                                                                                                          GENOME_REF,
                                                                                                          SAMTOOLS)
cm2 = '{} sort -@ {} R1_cancer_unpaired_aligned.bam > R1_cancer_unpaired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                    THREADS)
p1 = subprocess.Popen(cmd1 + '; wait ; ' + cmd2, shell=True)
p1.wait()

# Normal (unpaired R2)
cmd1 = '{} -H -t {} {} R2_normal_unpaired.fastq.gz | {} view -bS > R2_normal_unpaired_aligned.bam'.format(BWA,
                                                                                                          THREADS,
                                                                                                          GENOME_REF,
                                                                                                          SAMTOOLS)
cm2 = '{} sort -@ {} R2_normal_unpaired_aligned.bam > R2_normal_unpaired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                    THREADS)
p1 = subprocess.Popen(cmd1 + '; wait ; ' + cmd2, shell=True)
p1.wait()
# Cancer (unpaired R2)
cmd1 = '{} -H -t {} {} R2_cancer_unpaired.fastq.gz | {} view -bS > R2_cancer_unpaired_aligned.bam'.format(BWA,
                                                                                                          THREADS,
                                                                                                          GENOME_REF,
                                                                                                          SAMTOOLS)
cm2 = '{} sort -@ {} R2_cancer_unpaired_aligned.bam > R2_cancer_unpaired_aligned_sorted.bam'.format(SAMTOOLS,
                                                                                                    THREADS)
p1 = subprocess.Popen(cmd1 + '; wait ; ' + cmd2, shell=True)
p1.wait()

print('Aligment of tumor and normal samples completed.')

# Merge aligned files
print('Merging aligned files')
cmd1 = '{} merge aligned_normal_merged.bam normal_paired_aligned_sorted.bam '\
       'R1_normal_unpaired_aligned_sorted.bam R2_normal_unpaired_aligned_sorted.bam'

cmd2 = '{} merge aligned_cancer_merged.bam cancer_paired_aligned_sorted.bam '\
       'R1_cancer_unpaired_aligned_sorted.bam R2_cancer_unpaired_aligned_sorted.bam'
p1 = subprocess.Popen(cmd1, shell=True)
p2 = subprocess.Popen(cmd2, shell=True)
p1.wait()
p2.wait()
print('Merging of tumor and normal aligned samples completed.')

# Final pìpeline
#Full_exome_pipeline("aligned_cancer_merged.bam", "aligned_normal_merged.bam ", tumor_type, GENOME_REF, sampleID)
#HLA_pipeline(loc, sample1, sample2, THREADS)


