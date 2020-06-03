#! /usr/bin/env python
"""
@author: jfnavarro
"""
from hlapipeline.common import exec_command
from hlapipeline.tools import *
from hlapipeline.epitopes import *
import re
import os
import argparse
import multiprocessing

def HLA_predictionRNA(sample, threads):
    cmd = '{} extract --threads {} --paired {}'.format(ARCASHLA, threads, sample)
    exec_command(cmd)

    clean_name = os.path.splitext(sample)[0]

    cmd = '{} genotype --threads {} {}.extracted.1.fq.gz {}.extracted.2.fq.gz'.format(ARCASHLA,
                                                                                      threads,
                                                                                      clean_name,
                                                                                      clean_name)
    exec_command(cmd)

    cmd = '{} partial --threads {} -G {}.genotype.json {}.extracted.1.fq.gz {}.extracted.2.fq.gz'.format(ARCASHLA,
                                                                                                         threads,
                                                                                                         clean_name,
                                                                                                         clean_name,
                                                                                                         clean_name)
    exec_command(cmd)

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
        print('Trimmimg reads')
        cmd = '{} --paired --basename sample {} {}'.format(TRIMGALORE, sample1, sample2)
        exec_command(cmd)

        print('Aligining with STAR')
        cmd = '{} --genomeDir {} --readFilesIn sample_val_1.fq.gz sample_val_2.fq.gz --outSAMmultNmax 1 --outSAMorder Paired'\
              ' --outSAMprimaryFlag OneBestScore --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {} --outFilterIntronMotifs'\
              ' RemoveNoncanonical --outFilterType Normal --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c'\
              ' --runThreadN {} --outFilterMultimapNmax 20'.format(STAR, genome_star, annotation, THREADS)
        exec_command(cmd)

        # Add headers
        print("Adding headers")
        cmd = '{} AddOrReplaceReadGroups I=Aligned.sortedByCoord.out.bam O=sample_header.bam SO=coordinate RGID=SB_RNA-seq RGLB={}'\
              ' RGPL=Illumina RGPU=PU{} RGSM=rna-seq_{} Create_Index=true Validation_Stringency=SILENT'.format(PICARD,
                                                                                                               'VHIO',
                                                                                                               'VHIO',
                                                                                                               sampleID)
        exec_command(cmd)

    if 'hla' in steps:
        print('Predicting HLA')
        HLA_predictionRNA('sample_header.bam', THREADS)

    if 'gatk' in steps:
        # Mark duplicates
        print('Marking duplicates')
        cmd = GATK + ' MarkDuplicatesSpark -I=sample_header.bam -O=sample_dedup.bam -M=dedup_sample.txt'
        exec_command(cmd)

        # Split N and cigars
        print('Splitting NCigar Reads')
        cmd = '{} SplitNCigarReads -R {} -I sample_dedup.bam -O sample_split.bam'.format(GATK, genome)
        exec_command(cmd)

        # GATK base re-calibration
        print('Starting re-calibration')
        # NOTE BaseRecalibratorSpark needs the system to allow for many open files (ulimit -n)
        cmd = '{} BaseRecalibrator -I sample_dedup.bam -R {} --known-sites {} --known-sites {}'\
              ' --known-sites {} -O sample_recal_data.txt'.format(GATK, genome, SNPSITES, KNOWN_SITE1, KNOWN_SITE2)
        exec_command(cmd)
        cmd = '{} ApplyBQSR -R {} -I sample_dedup.bam --bqsr-recal-file sample_recal_data.txt -O sample_final.bam'.format(GATK, genome)
        exec_command(cmd)

    if 'variant' in steps:
        # Variant calling (Samtools pile-ups)
        print('Computing pile-ups')
        cmd = '{} mpileup -C50 -B -q 1 -Q 15 -f {} sample_final.bam > sample.pileup'.format(SAMTOOLS, genome)
        exec_command(cmd)

        # Variant calling VarScan
        print('Variant calling with varscan')
        cmd = VARSCAN + ' mpileup2cns sample.pileup varscan --variants 0 --min-coverage 2 --min-reads2 1 --output-vcf 1'\
              + ' --min-var-freq .01 --p-value 0.99 > varscan.vcf'
        exec_command(cmd)
        cmd = VARSCAN + ' mpileup2cns sample.pileup varscan --variants 0 --min-coverage 2 --min-reads2 1'\
              + ' --min-var-freq .01 --p-value 0.99 > varscan.pileup'
        exec_command(cmd)

        # Running Cufflinks (output is genes.fpkm_tracking)
        print('Running Cufflinks')
        cmd = '{} -p {} -G {} --library-type fr-firststrand sample_dedup.bam'.format(CUFFLINKS, THREADS, annotation)
        exec_command(cmd)

    if 'filter' in steps:
        # TODO apply a filter here with vcftools on varscan.vcf

        # Run annovar to annotate variants
        print('Running annovar')
        annovardb = '{} -buildver {}'.format(os.path.join(ANNOVAR_PATH, ANNOVAR_DB), ANNOVAR_VERSION)
        cmd = '{} -format vcf4 varscan.vcf --comment --includeinfo -outfile snp.av'.format(
            os.path.join(ANNOVAR_PATH, 'convert2annovar.pl'))
        exec_command(cmd)
        cmd = '{} snp.av {} -thread {} -out snp.sum -remove -protocol {}'.format(
            os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), annovardb, THREADS, annovar_anno)
        exec_command(cmd)

        print('Formatting Varscan variants')
        snv = open('snp.sum.{}_multianno.txt'.format(ANNOVAR_VERSION))
        snv_lines = snv.readlines()
        header = snv_lines.pop(0).strip().split('\t')
        insert_file = open('SQL_variant_input.txt', 'w')
        for line in snv_lines:
            if line.startswith('#'):
                continue
            columns = line.rstrip('\n').split('\t')
            Chr = columns[header.index('Chr')]
            start = columns[header.index('Start')]
            end = columns[header.index('End')]
            ref = columns[header.index('Ref')]
            alt = columns[header.index('Alt')]
            func_ref_gene = columns[header.index('Func.refGene')]
            gene_ref_gene = columns[header.index('Gene.refGene')]
            ref_gene_detail = columns[header.index('GeneDetail.refGene')]
            exonic_func_ref = columns[header.index('ExonicFunc.refGene')]
            AA_change_refGene = columns[header.index('AAChange.refGene')]
            func_known_gene = columns[header.index('Func.knownGene')]
            gene_known_gene = columns[header.index('Gene.knownGene')]
            known_gene_detail = columns[header.index('GeneDetail.knownGene')]
            exonic_known_ref = columns[header.index('ExonicFunc.knownGene')]
            AA_change_knownGene = columns[header.index('AAChange.knownGene')]
            func_ens_gene = columns[header.index('Func.ensGene')]
            gene_ens_gene = columns[header.index('Gene.ensGene')]
            ens_gene_detail = columns[header.index('GeneDetail.ensGene')]
            exonic_ens_ref = columns[header.index('ExonicFunc.ensGene')]
            AA_change_ensGene = columns[header.index('AAChange.ensGene')]
            avsnp150 = columns[header.index('avsnp150')]
            apr_all = columns[header.index('ALL.sites.2015_08')]
            apr_eur = columns[header.index('EUR.sites.2015_08')]
            apr_amr = columns[header.index('AMR.sites.2015_08')]
            apr_asn = columns[header.index('EAS.sites.2015_08')]
            apr_afr = columns[header.index('AFR.sites.2015_08')]
            gDNA = 'chr' + Chr + ':' + start
            variant_key = Chr + ':' + start + '-' + end + ' ' + ref + '>' + alt
            if ref_gene_detail != 'NA':
                AA_change_refGene = ref_gene_detail
            if known_gene_detail != 'NA':
                AA_change_knownGene = known_gene_detail
            if ens_gene_detail != 'NA':
                AA_change_ensGene = ens_gene_detail
            to_write = '\t'.join([str(x) for x in [gDNA, sampleID, Chr, start, end, ref, alt, avsnp150,
                                                   func_ref_gene, gene_ref_gene, exonic_func_ref, AA_change_refGene,
                                                   func_known_gene, gene_known_gene, exonic_known_ref, AA_change_knownGene,
                                                   func_ens_gene, gene_ens_gene, exonic_ens_ref, AA_change_ensGene,
                                                   apr_all, apr_eur, apr_amr, apr_asn, apr_afr, variant_key]])
            insert_file.write(to_write + "\n")
        insert_file.close()
        snv.close()

        # This will format the  coverage information into a format that will be joined with the variant information
        print('Formatting varscan pile-ups')
        pileup = open('varscan.pileup')
        pileup_lines = pileup.readlines()
        header = pileup_lines.pop(0).strip().split('\t')
        insert_file = open('SQL_coverage_input.txt', 'w')
        for line in pileup_lines:
            columns = line.rstrip('\n').split('\t')
            Chr = columns[header.index('Chrom')]
            start = columns[header.index('Position')]
            ref = columns[header.index('Ref')]
            alt = columns[header.index('Var')]
            columns2 = columns[header.index('Cons:Cov:Reads1:Reads2:Freq:P-value')].split(':')
            cons = columns2[0]
            cov = columns2[1]
            read1 = columns2[2]
            read2 = columns2[3]
            freq = columns2[4]
            p_val = columns2[5]
            columns3 = columns[header.index('StrandFilter:R1+:R1-:R2+:R2-:pval')].split(':')
            r1_plus = columns3[1]
            r1_minus = columns3[2]
            r2_plus = columns3[3]
            r2_minus = columns3[4]
            p_val2 = columns3[5]
            gDNA = 'chr' + Chr + ':' + start
            join_key = gDNA
            if re.search(r'-', alt):
                join_key = 'chr' + Chr + ':' + str(int(start) + 1)
            to_write = '\t'.join([str(x) for x in [join_key, sampleID, Chr, start, ref, alt, cons, cov, read1, read2,
                                                   freq, p_val, r1_plus, r1_minus, r2_plus, r2_minus, p_val2]])
            insert_file.write(to_write + "\n")
        pileup.close()
        insert_file.close()

        print('Sorting and joining files')
        cmd = 'sort -k1b,1 SQL_coverage_input.txt > join_coverage_sort.txt'
        exec_command(cmd)
        cmd = 'sort -k1b,1 SQL_variant_input.txt > join_variants_sort.txt'
        exec_command(cmd)
        cmd = 'join -t$\'\\t\' -a1 -e"-" join_variants_sort.txt join_coverage_sort.txt > joined_coverage_variants.txt'
        exec_command(cmd)

        # This will extract just the nonsynonymous mutations:
        print('Creating final variants')
        joined_variants = open('joined_coverage_variants.txt')
        nonsyn_file = open('nonsyn_SQL_insert.txt', 'w')
        all_file = open('all_other_mutations.txt', 'w')
        header = 'SAMPLE_ID\tCHR\tSTART\tEND\tREF\tALT\tsnp138JJG\tTUMOR_READ1\tTUMOR_READ2' \
                 '\tTUMOR_VAR_FREQ\tAPR_ALL\tFunc_refGene\tGene_refGene\tExonicFunc_refGene\tAAChange_refGene\tFunc_knownGene' \
                 '\tGene_knownGene\tExonicFunc_knownGene\tAAChange_knownGene\tFunc_ensGene\tGene\texonic_func\tNT-AA_CHANGE' \
                 '\tAPR_EUR\tAPR_AMR2\tAPR_ASN3\tAPR_AFR4\tREAD1_PLUS\tREAD1_MINUS\tREAD2_PLUS\tREAD2_MINUS' \
                 '\ttumour_type\tsource_of_RNA_used_for_sequencing\tVARIANT-KEY\n'
        nonsyn_file.write(header)
        all_file.write(header)
        for line in joined_variants:
            columns = line.rstrip('\n').split('\t')
            # Very ugly way to check that the variant was called by the two methods (TODO improve)
            if len(columns) < 41:
                continue
            # TODO use header names instead
            Chr = columns[2]
            Start = columns[3]
            End = columns[4]
            Ref = columns[5]
            Alt = columns[6]
            snp138JJG = columns[7]
            Func_refGene = columns[8]
            Gene_refGene = columns[9]
            ExonicFunc_refGene = columns[10]
            AAChange_refGene = columns[11]
            Func_knownGene = columns[12]
            Gene_knownGene = columns[13]
            ExonicFunc_knownGene = columns[14]
            AAChange_knownGene = columns[15]
            Func_ensGene = columns[16]
            Gene_ensGene = columns[17]
            ExonicFunc_ensGene = columns[18]
            AAChange_ensGene = columns[19]
            apr_all = columns[20]
            apr_eur = columns[21]
            apr_amr = columns[22]
            apr_asn = columns[23]
            apr_afr = columns[24]
            variant_key = columns[25]
            source_of_RNA_used_for_sequencing = columns[32]
            tumor_reads1 = columns[33]
            tumor_reads2 = columns[34]
            tumor_var_freq = columns[35].replace('%', '')
            read1_plus = columns[37]
            read1_minus = columns[38]
            read2_plus = columns[39]
            read2_minus = columns[40]
            to_write = '\t'.join([str(x) for x in [sampleID, Chr, Start, End, Ref, Alt, snp138JJG, tumor_reads1,
                                                   tumor_reads2, tumor_var_freq, apr_all, Func_refGene, Gene_refGene,
                                                   ExonicFunc_refGene, AAChange_refGene, Func_knownGene, Gene_knownGene,
                                                   ExonicFunc_knownGene, AAChange_knownGene, Func_ensGene, Gene_ensGene,
                                                   ExonicFunc_ensGene, AAChange_ensGene, apr_eur, apr_amr, apr_asn,
                                                   apr_afr, read1_plus, read1_minus, read2_plus, read2_minus, tumor_type,
                                                   source_of_RNA_used_for_sequencing, variant_key]])
            if (re.search(r'nonsynonymous', ExonicFunc_refGene) or re.search(r'frame', ExonicFunc_refGene) or re.search(r'stop', ExonicFunc_refGene)\
                    or re.search(r'nonsynonymous', ExonicFunc_knownGene) or re.search(r'frame', ExonicFunc_knownGene) or re.search(r'stop', ExonicFunc_knownGene)\
                    or re.search(r'nonsynonymous', ExonicFunc_ensGene) or re.search(r'frame', ExonicFunc_ensGene) or re.search(r'stop', ExonicFunc_ensGene)):
                nonsyn_file.write(to_write + '\n')
            else:
                all_file.write(to_write + '\n')
        joined_variants.close()
        nonsyn_file.close()
        all_file.close()

        # Extract peptides
        extract_peptides('nonsyn_SQL_insert.txt', 'Formatted_epitope_variant.txt', sampleID)

        # Create epitopes
        create_epitopes('Formatted_epitope_variant.txt', 'SQL_Epitopes.txt', FASTA_AA_DICT, FASTA_cDNA_DICT)
                        
        # Reformat FPKM file
        print('Creating FPKM info file')
        fpkm = open('genes.fpkm_tracking')
        fpkm_lines = fpkm.readlines()
        firstline = fpkm_lines.pop(0)
        FPKM_ins = open('FPKM_SQL_insert.txt', 'w')
        header = 'SAMPLE_ID\tTUMOUR\t' + firstline
        FPKM_ins.write(header)
        for line in fpkm_lines:
            FPKM_ins.write('{}\t{}\t{}'.format(sampleID, tumor_type, line))
        FPKM_ins.close()
        fpkm.close()

    print("COMPLETED!")

parser = argparse.ArgumentParser(description='RNA-seq variant calling and HLA prediction pipeline'
                                 '(created by Jose Fernandez <jc.fernandes.navarro@gmail.com>)',
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
                    help='Steps to perform in the pipeline', choices=['mapping', 'gatk', 'hla', 'variant', 'filter', "none"])

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
