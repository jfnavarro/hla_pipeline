#! /usr/bin/env python

import argparse
from Full_exome_pipeline import *
from RNA_seq_pipeline import *
import multiprocessing

parser = argparse.ArgumentParser(description='Jareds pipeline (adjusted by Jose Fernandez <jc.fernandes.navarro@gmail.com>',
                                 prog='pipeline.py',
                                 usage='pipeline.py [options] R1(Normal) R2(Normal) R1(Cancer) R2(Cancer) R1(RNA) R2(RNA)')
parser.add_argument('R1_NORMAL', help='FASTQ file R1 (Normal)')
parser.add_argument('R2_NORMAL', help='FASTQ file R2 (Normal)')
parser.add_argument('R1_CANCER', help='FASTQ file R1 (Cancer)')
parser.add_argument('R2_CANCER', help='FASTQ file R2 (Cancer)')
parser.add_argument('R1_RNA', help='FASTQ file R1 (RNA)')
parser.add_argument('R2_RNA', help='FASTQ file R2 (RNA)')
parser.add_argument('--adapter',
                    help='Path to the Illumina adapters FASTA file.', required=True)
parser.add_argument('--genome',
                    help='Path to the reference Genome FASTA file (must contain BWA index)', required=True)
parser.add_argument('--genome-star',
                    help='Path to the reference Genome STAR index folder', required=True)
parser.add_argument('--genome-ref',
                    help='Path to the reference Genome GTF file', required=True)
parser.add_argument('--sample',
                    help='Name of the sample/experiment. Default is sample', default='sample')
parser.add_argument('--tumor',
                    help='Tumor type. Default is Tumor', default='Tumor')
parser.add_argument('--dir',
                    help='Path to the output files', required=True)
parser.add_argument('--known1',
                    help='Path to the file with Mill and 1000G gold standards (GATK bundle)', required=True)
parser.add_argument('--known2',
                    help='Path to the file with 1000G phase indels (GATK bundle)', required=True)
parser.add_argument('--snpsites',
                    help='Path to the file with the SNPs (GATK buldle)', required=True)
parser.add_argument('--germline',
                    help='Path to the file with the germline resources Nomad (GATK buldle)', required=True)
parser.add_argument('--pon',
                    help='Path to the file with the panel of normals (GATK buldle)', required=True)
parser.add_argument('--fastaAA',
                    help='Path to the file with the dictionary of FASTA to AA', required=True)
parser.add_argument('--fastacDNA',
                    help='Path to the file with the dictionary of FASTA to cDNA', required=True)
parser.add_argument('--dna-steps', nargs='+', default=['mapping', 'gatk', 'hla', 'variant', 'filter'],
                    help='Steps to apply in the DNA pipeline', choices=['mapping', 'gatk', 'hla', 'variant', 'filter', "none"])
parser.add_argument('--rna-steps', nargs='+', default=['mapping', 'gatk', 'hla', 'variant', 'filter'],
                    help='Steps to apply in the RNA pipeline', choices=['mapping', 'gatk', 'hla', 'variant', 'filter', "none"])

# Parse arguments
args = parser.parse_args()
DIR = args.dir
R1_NORMAL = os.path.abspath(args.R1_NORMAL)
R2_NORMAL = os.path.abspath(args.R2_NORMAL)
R1_CANCER = os.path.abspath(args.R1_CANCER)
R2_CANCER = os.path.abspath(args.R2_CANCER)
R1_RNA = os.path.abspath(args.R1_RNA)
R2_RNA = os.path.abspath(args.R2_RNA)
sampleID = args.sample
tumor_type = args.tumor
IILLUMINA_ADAPTERS = os.path.abspath(args.adapter)
GENOME_REF = os.path.abspath(args.genome)
GENOME_REF_STAR = os.path.abspath(args.genome_star)
GENOME_ANNOTATION = os.path.abspath(args.genome_ref)
THREADS = multiprocessing.cpu_count() - 1
FASTA_AA_DICT = os.path.abspath(args.fastaAA)
FASTA_cDNA_DICT = os.path.abspath(args.fastacDNA)
KNOWN_SITE1 = os.path.abspath(args.known1)
KNOWN_SITE2 = os.path.abspath(args.known2)
SNPSITES = os.path.abspath(args.snpsites)
GERMLINE = os.path.abspath(args.germline)
PON = os.path.abspath(args.pon)
DNA_STEPS = args.dna_steps
RNA_STEPS = args.rna_steps

# Move to output dir
os.makedirs(os.path.abspath(DIR), exist_ok=True)
os.chdir(os.path.abspath(DIR))

# Exome pìpeline
Full_exome_pipeline(R1_NORMAL,
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
                    DNA_STEPS)
os.chdir('..')

# RNA pìpeline
RNA_seq_pipeline(R1_RNA,
                 R2_RNA,
                 sampleID,
                 GENOME_REF,
                 GENOME_REF_STAR,
                 GENOME_ANNOTATION,
                 tumor_type,
                 SNPSITES,
                 KNOWN_SITE1,
                 KNOWN_SITE2,
                 THREADS,
                 RNA_STEPS)
os.chdir('..')


