#! /usr/bin/env python
"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from hlapipeline.common import *
from hlapipeline.tools import *
from hlapipeline.epitopes import *
import os
import argparse

def RNAseq_pipeline(sample1,
                    sample2,
                    sampleID,
                    genome,
                    genome_star,
                    annotation,
                    tumor_type,
                    SNPSITES,
                    KNOWN_SITE1,
                    KNOWN_SITE2,
                    FASTA_AA_DICT,
                    FASTA_cDNA_DICT,
                    THREADS,
                    ANNOVAR_DB,
                    ANNOVAR_VERSION,
                    steps):
    print("RNA-seq pipeline")

    # Create sub-folder to store all results
    os.makedirs('rna', exist_ok=True)
    os.chdir('rna')

    if 'mapping' in steps:
        print('Trimming reads')
        cmd = '{} --paired --basename sample {} {}'.format(TRIMGALORE, sample1, sample2)
        exec_command(cmd)

        print('Aligning with STAR')
        #cmd = '{} --genomeDir {} --readFilesIn sample_val_1.fq.gz sample_val_2.fq.gz --outSAMmultNmax 1 --outSAMorder Paired' \
        #      ' --outSAMprimaryFlag OneBestScore --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {} --outFilterIntronMotifs' \
        #      ' RemoveNoncanonical --outFilterType Normal --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c' \
        #      ' --runThreadN {} --outFilterMultimapNmax 20'.format(STAR, genome_star, annotation, THREADS)

        cmd = '{} --genomeDir {} --readFilesIn sample_val_1.fq.gz sample_val_2.fq.gz --outSAMorder Paired' \
              ' --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {}' \
              ' --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c' \
              ' --runThreadN {}'.format(STAR, genome_star, annotation, THREADS)

        exec_command(cmd)

        # Add headers
        print("Adding headers")
        cmd = '{} AddOrReplaceReadGroups I=Aligned.sortedByCoord.out.bam O=sample_header.bam SO=coordinate RGID=RNA RGLB={}' \
              ' RGPL=Illumina RGPU=PU{} RGSM=rna-seq_{} Create_Index=true Validation_Stringency=SILENT'.format(PICARD,
                                                                                                               'VHIO',
                                                                                                               'VHIO',
                                                                                                               sampleID)
        exec_command(cmd)

    if 'hla' in steps:
        print('Predicting HLAs')
        HLA_predictionRNA('sample_header.bam', THREADS)

    if 'gatk' in steps:
        # Mark duplicates
        print('Marking duplicates')
        cmd = GATK + ' MarkDuplicatesSpark --VALIDATION_STRINGENCY SILENT -I=sample_header.bam ' \
                     '-O=sample_dedup.bam -M=dedup_sample.txt'
        exec_command(cmd)

        # Split N and cigars
        print('Splitting NCigar Reads')
        cmd = '{} SplitNCigarReads -R {} -I sample_dedup.bam -O sample_split.bam'.format(GATK, genome)
        exec_command(cmd)

        # GATK base re-calibration
        print('Starting re-calibration')
        # NOTE BaseRecalibratorSpark needs the system to allow for many open files (ulimit -n)
        cmd = '{} BaseRecalibrator --use-original-qualities -I sample_dedup.bam -R {} --known-sites {} --known-sites {}' \
              ' --known-sites {} -O sample_recal_data.txt'.format(GATK, genome, SNPSITES, KNOWN_SITE1, KNOWN_SITE2)
        exec_command(cmd)
        cmd = '{} ApplyBQSR --use-original-qualities --add-output-sam-program-record -R {} -I sample_dedup.bam ' \
              '--bqsr-recal-file sample_recal_data.txt -O sample_final.bam'.format(GATK, genome)
        exec_command(cmd)

    if 'variant' in steps:
        # Variant calling (Samtools pile-ups)
        print('Computing pile-ups')
        cmd = '{} mpileup -C50 -B -q 1 -Q 15 -f {} sample_final.bam > sample.pileup'.format(SAMTOOLS, genome)
        exec_command(cmd)

        # Variant calling VarScan
        print('Variant calling with varscan')
        cmd = '{} mpileup2cns sample.pileup varscan --variants 0 --min-coverage 2 --min-reads2 1 --output-vcf 1 ' \
              '--min-var-freq .01 --p-value 0.99 > varscan.vcf'.format(VARSCAN)
        exec_command(cmd)

        # Variant calling (HaplotypeCaller)
        print('Variant calling with HaplotypeCaller')
        cmd = '{} HaplotypeCaller -R {} -I sample_final.bam -O haplotype_caller.vcf ' \
              '-dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20 ' \
              '--dbsnp {}'.format(GATK, genome, SNPSITES)
        exec_command(cmd)

        # Computing gene counts
        print('Computing gene counts with featureCounts')
        cmd = '{} -T {} --primary --ignoreDup -O -C -t exon ' \
              '-g gene_name -a {} -o gene.counts sample_dedup.bam'.format(FEATURECOUNTS, THREADS, annotation)
        exec_command(cmd)

    if 'filter' in steps:
        # Filtering variants (HaplotypeCaller)
        print("Filtering HaplotypeCaller variants")
        cmd = '{} VariantFiltration --R {} --V haplotype_caller.vcf --window 35 --cluster 3 --filter-name "FS" ' \
              '--filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" -O haplotype_caller_filtered.vcf'.format(GATK, genome)
        exec_command(cmd)

        # Filtering variants (VarScan)
        print("Filtering VarScan variants")
        cmd = '{} VariantFiltration --R {} --V varscan.vcf --window 35 --cluster 3 --filter-name "FS" ' \
              '--filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" -O varscan_filtered.vcf'.format(GATK, genome)
        exec_command(cmd)

        # Combine with GATK
        print('Combining variants')
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        cmd = '{} -T CombineVariants -R {} -V:varscan varscan_filtered.vcf ' \
              '-V:HaplotypeCaller haplotype_caller_filtered.vcf -o combined_calls.vcf '\
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

        # TODO extract peptides and epitopes

        # Extract peptides
        #extract_peptides('nonsyn_SQL_insert.txt', 'Formatted_epitope_variant.txt', sampleID)

        # Create epitopes
        #create_epitopes('Formatted_epitope_variant.txt', 'SQL_Epitopes.txt', FASTA_AA_DICT, FASTA_cDNA_DICT)

        # Reformat Gene counts file
        print('Creating Gene counts info file')
        counts_file = open('gene.counts')
        lines = counts_file.readlines()
        if lines[0].startswith("#"):
            lines.pop(0)
        secondline = lines.pop(0)
        counts_out = open('GeneCounts_SQL_insert.txt', 'w')
        header = 'SAMPLE_ID\tTUMOUR\t' + secondline
        counts_out.write(header)
        for line in lines:
            counts_out.write('{}\t{}\t{}'.format(sampleID, tumor_type, line))
        counts_out.close()
        counts_file.close()

    print("COMPLETED!")


parser = argparse.ArgumentParser(description='RNA-seq variant calling and HLA prediction pipeline\n'
                                             'Created by Jose Fernandez <jc.fernandes.navarro@gmail.com>)',
                                 prog='rnaseq_pipeline.py',
                                 usage='rnaseq_pipeline.py [options] R1(RNA) R2(RNA)')
parser.add_argument('R1_RNA', help='FASTQ file R1 (RNA)')
parser.add_argument('R2_RNA', help='FASTQ file R2 (RNA)')
parser.add_argument('--genome',
                    help='Path to the reference Genome FASTA file', required=True)
parser.add_argument('--genome-star',
                    help='Path to the reference Genome STAR index folder', required=True)
parser.add_argument('--genome-ref',
                    help='Path to the reference Genome GTF file', required=True)
parser.add_argument('--sample',
                    help='Name of the sample/experiment. Default is sample', default='sample')
parser.add_argument('--tumor',
                    help='Tumor type. Default is Tumor', default='Tumor')
parser.add_argument('--dir',
                    help='Path to the folder where to put output files', required=True)
parser.add_argument('--known1',
                    help='Path to the file with Mill and 1000G gold standards (GATK bundle)', required=True)
parser.add_argument('--known2',
                    help='Path to the file with 1000G phase indels (GATK bundle)', required=True)
parser.add_argument('--snpsites',
                    help='Path to the file with the SNPs (GATK buldle)', required=True)
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
parser.add_argument('--threads',
                    help='Number of threads to use in the parallel steps', type=int, default=10, required=False)
parser.add_argument('--steps', nargs='+', default=['mapping', 'gatk', 'hla', 'variant', 'filter'],
                    help='Steps to perform in the pipeline',
                    choices=['mapping', 'gatk', 'hla', 'variant', 'filter', "none"])

# Parse arguments
args = parser.parse_args()
DIR = args.dir
R1_RNA = os.path.abspath(args.R1_RNA)
R2_RNA = os.path.abspath(args.R2_RNA)
sampleID = args.sample
tumor_type = args.tumor
GENOME_REF = os.path.abspath(args.genome)
GENOME_REF_STAR = os.path.abspath(args.genome_star)
GENOME_ANNOTATION = os.path.abspath(args.genome_ref)
THREADS = int(args.threads)
KNOWN_SITE1 = os.path.abspath(args.known1)
KNOWN_SITE2 = os.path.abspath(args.known2)
SNPSITES = os.path.abspath(args.snpsites)
FASTA_AA_DICT = os.path.abspath(args.fastaAA)
FASTA_cDNA_DICT = os.path.abspath(args.fastacDNA)
RNA_STEPS = args.steps
ANNOVAR_DB = args.annovar_db
ANNOVAR_VERSION = args.annovar_version

# Move to output dir
os.makedirs(os.path.abspath(DIR), exist_ok=True)
os.chdir(os.path.abspath(DIR))

RNAseq_pipeline(R1_RNA,
                R2_RNA,
                sampleID,
                GENOME_REF,
                GENOME_REF_STAR,
                GENOME_ANNOTATION,
                tumor_type,
                SNPSITES,
                KNOWN_SITE1,
                KNOWN_SITE2,
                FASTA_AA_DICT,
                FASTA_cDNA_DICT,
                THREADS,
                ANNOVAR_DB,
                ANNOVAR_VERSION,
                RNA_STEPS)
