#! /usr/bin/env python
"""
This pipeline computes germline variants from DNA or RNA data.

The pipline trims with trimgalore, aligns with STAR or bwa-men,
performs the GATK4 best practices and computes variants with
HaplotypeCaller and Varscan. The variants are then combined into
one file and annotated with Annovar.

Multiple options are available. To see them type --help

@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from hlapipeline.common import *
from hlapipeline.tools import *
import os
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def main(sample1,
         sample2,
         sampleID,
         genome,
         genome_star,
         annotation,
         SNPSITES,
         KNOWN_SITE1,
         KNOWN_SITE2,
         THREADS,
         ANNOVAR_DB,
         ANNOVAR_VERSION,
         steps,
         mode):

    # TODO add sanity checks for the parameters
    # TODO better log info
    # TODO remove temp files
    # TODO put output files somewhere else

    print("Germline pipeline")

    # Create sub-folder to store all results
    os.makedirs('germline', exist_ok=True)
    os.chdir('germline')

    if 'mapping' in steps:
        print('Trimming reads')
        cmd = '{} --cores {} --paired --basename sample {} {}'.format(TRIMGALORE, THREADS, sample1, sample2)
        exec_command(cmd)

        # ALIGNMENT
        print('Starting alignment')
        if mode == 'DNA':
            cmd = '{} -t {} {} sample_val_1.fq.gz sample_val_2.fq.gz | ' \
                  '{} sort --threads {} > trimmed_paired_aligned_sorted.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
            exec_command(cmd)
        else:
            cmd = '{} --genomeDir {} --readFilesIn sample_val_1.fq.gz sample_val_2.fq.gz --outSAMorder Paired' \
                  ' --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {}' \
                  ' --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c' \
                  ' --runThreadN {}'.format(STAR, genome_star, annotation, THREADS)
            exec_command(cmd)

            cmd = 'mv Aligned.sortedByCoord.out.bam trimmed_paired_aligned_sorted.bam'
            exec_command(cmd)

        # Add headers
        print("Adding headers")
        cmd = '{} AddOrReplaceReadGroups I=trimmed_paired_aligned_sorted.bam O=sample_header.bam SO=coordinate RGID={} RGLB=VHIO' \
              ' RGPL=Illumina RGPU=VHIO RGSM={} Create_Index=true Validation_Stringency=SILENT'.format(PICARD, mode, sampleID)
        exec_command(cmd)

    if 'gatk' in steps:
        # Mark duplicates
        print('Marking duplicates')
        cmd = '{} MarkDuplicatesSpark --input sample_header.bam --output sample_dedup.bam'.format(GATK)
        exec_command(cmd)

        if mode == 'RNA':
            # Split N and cigars
            print('Splitting NCigar Reads')
            cmd = '{} SplitNCigarReads --reference {} --input sample_dedup.bam --output sample_split.bam'.format(GATK, genome)
            exec_command(cmd)

            cmd = 'rm -rf sample_dedup* && mv sample_split.bam sample_dedup.bam && mv sample_split.bai sample_dedup.bai'
            exec_command(cmd)

        # GATK base re-calibration
        print('Starting re-calibration')
        cmd = '{} BaseRecalibratorSpark --use-original-qualities --input sample_dedup.bam --reference {} --known-sites {} ' \
              '--known-sites {} --known-sites {} --output sample_recal_data.txt'.format(GATK,
                                                                                        genome,
                                                                                        SNPSITES,
                                                                                        KNOWN_SITE1,
                                                                                        KNOWN_SITE2)
        exec_command(cmd)
        cmd = '{} ApplyBQSR --use-original-qualities --add-output-sam-program-record --reference {} --input sample_dedup.bam ' \
              '--bqsr-recal-file sample_recal_data.txt --output sample_final.bam'.format(GATK, genome)
        exec_command(cmd)

    if 'hla' in steps:
        print('Predicting HLAs')
        if mode == 'DNA':
            HLA_predictionDNA('sample_final.bam', sampleID, 'PRG-HLA-LA_output.txt', THREADS)
        else:
            HLA_predictionRNA('sample_final.bam', THREADS)

    if 'variant' in steps:
        # Variant calling (Samtools pile-ups)
        print('Computing pile-ups')
        cmd = '{} mpileup -C 50 -B -q 1 -Q 15 -f {} sample_final.bam > sample.pileup'.format(SAMTOOLS, genome)
        exec_command(cmd)

        # Variant calling VarScan
        print('Variant calling with VarScan2')
        cmd = '{} mpileup2cns sample.pileup --variants 0 --min-coverage 2 --min-reads2 1 --output-vcf 1 ' \
              '--min-var-freq 0.01 --min-avg-qual 15 --p-value 0.99 --strand-filter 1 > varscan.vcf'.format(VARSCAN)
        exec_command(cmd)

        # Variant calling (HaplotypeCaller)
        print('Variant calling with HaplotypeCaller')
        cmd = '{} HaplotypeCaller --reference {} --input sample_final.bam --output haplotypecaller.vcf ' \
              '--dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20 ' \
              '--dbsnp {}'.format(GATK, genome, SNPSITES)
        exec_command(cmd)

        if mode == 'RNA':
            # Computing gene counts
            print('Computing gene counts with featureCounts')
            cmd = '{} -T {} --primary --ignoreDup -O -C -t exon ' \
                  '-g gene_name -a {} -o gene.counts sample_dedup.bam'.format(FEATURECOUNTS, THREADS, annotation)
            exec_command(cmd)

    if 'filter' in steps:
        # Filtering variants (HaplotypeCaller)
        print("Filtering HaplotypeCaller variants")
        cmd = '{} VariantFiltration --reference {} --variant haplotype_caller.vcf --window 35 --cluster 3 --filter-name "FS" ' \
              '--filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" --output haplotype_caller_filtered.vcf'.format(GATK, genome)
        exec_command(cmd)

        cmd = 'sed -i \'s/{}/HaplotypeCaller/g\' haplotype_caller_filtered.vcf'.format(sampleID)
        exec_command(cmd)

        # NOTE replacing IUPAC codes from varscan output (TEMP HACK till the bug is fixed in VarScan)
        # NOTE this will also skip variants whose REF and ALT fields are identical (another bug in VarScan)
        cmd = 'awk \'{if ($1 ~ /#/) {print} else if ($4 != $5) {gsub(/W|K|B|Y|D|H|V|R|S|M/,"N",$4); OFS="\t"; print}}\' ' \
              'varscan.vcf > varscan_filtered.vcf'
        exec_command(cmd)

        cmd = 'sed -i \'s/Sample1.varscan/varscan/g\' varscan_filtered.vcf'.format(sampleID)
        exec_command(cmd)

        # Combine with GATK
        print('Combining variants')
        # TODO replace this with jacquard merge
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        cmd = '{} -T CombineVariants -R {} -V:varscan varscan_filtered.vcf ' \
              '-V:HaplotypeCaller haplotype_caller_filtered.vcf -o combined_calls.vcf '\
              '-genotypeMergeOptions UNIQUIFY --num_threads {}'.format(GATK3, genome, THREADS)
        exec_command(cmd)

        # Annotate with Annovar
        annovardb = '{} -buildver {}'.format(os.path.join(ANNOVAR_PATH, ANNOVAR_DB), ANNOVAR_VERSION)
        print('Running annovar (SNV)')
        cmd = '{} combined_calls.vcf {} -thread {} -out annotated -vcfinput -remove -protocol {}'.format(
            os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), annovardb, THREADS,  annovar_anno)
        exec_command(cmd)

    print("COMPLETED!")

if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('R1', help='FASTQ file R1 (DNA or RNA)')
    parser.add_argument('R2', help='FASTQ file R2 (DNA or RNA)')
    parser.add_argument('--genome',
                        help='Path to the reference genome FASTA file', required=True)
    parser.add_argument('--genome-star',
                        help='Path to the reference genome STAR index folder (when in RNA mode)', required=False)
    parser.add_argument('--genome-ref',
                        help='Path to the reference genome GTF file (when in RNA mode)', required=False)
    parser.add_argument('--sample',
                        help='Name of the sample/experiment. Default is sample', default='sample')
    parser.add_argument('--dir',
                        help='Path to the folder where to put output files', required=True)
    parser.add_argument('--known1',
                        help='Path to the file with Mill and 1000G gold standards (GATK bundle)', required=True)
    parser.add_argument('--known2',
                        help='Path to the file with 1000G phase indels (GATK bundle)', required=True)
    parser.add_argument('--snpsites',
                        help='Path to the file with the SNP sites (GATK bundle)', required=True)
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
    parser.add_argument('--mode', default='RNA',
                        help='Mode to use (RNA (default) or DNA)',
                        choices=['DNA', 'RNA'])

    # Parse arguments
    args = parser.parse_args()
    DIR = args.dir
    R1_RNA = os.path.abspath(args.R1)
    R2_RNA = os.path.abspath(args.R2)
    sampleID = args.sample
    GENOME_REF = os.path.abspath(args.genome)
    GENOME_REF_STAR = os.path.abspath(args.genome_star)
    GENOME_ANNOTATION = os.path.abspath(args.genome_ref)
    THREADS = int(args.threads)
    KNOWN_SITE1 = os.path.abspath(args.known1)
    KNOWN_SITE2 = os.path.abspath(args.known2)
    SNPSITES = os.path.abspath(args.snpsites)
    STEPS = args.steps
    ANNOVAR_DB = args.annovar_db
    ANNOVAR_VERSION = args.annovar_version
    MODE = args.mode
    if MODE == "RNA" and (not GENOME_REF_STAR or not GENOME_ANNOTATION):
        sys.stderr.write("Error, RNA mode but STAR reference or annotation files are missing\n")
        sys.exit(1)

    # Move to output dir
    os.makedirs(os.path.abspath(DIR), exist_ok=True)
    os.chdir(os.path.abspath(DIR))

    main(R1_RNA,
         R2_RNA,
         sampleID,
         GENOME_REF,
         GENOME_REF_STAR,
         GENOME_ANNOTATION,
         SNPSITES,
         KNOWN_SITE1,
         KNOWN_SITE2,
         THREADS,
         ANNOVAR_DB,
         ANNOVAR_VERSION,
         STEPS,
         MODE)
