#! /usr/bin/env python
"""
This pipeline computes somatic variants from RNA data.
The pipeline trims with trimgalore, aligns with STAR,
performs the GATK4 best practices and computes variants with
HaplotypeCaller and Varscan. The variants are then combined into
one file and annotated with Annovar. Gene counts are also
computed with featureCounts.
Multiple options are available. To see them type --help
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from hlapipeline.common import *
from hlapipeline.tools import *
import os
import shutil
import glob
from argparse import ArgumentParser, RawDescriptionHelpFormatter

def main(R1,
         R2,
         SAMPLEID,
         GENOME,
         GENOME_STAR,
         ANNOTATION,
         SNPSITES,
         KNOWN_SITE1,
         KNOWN_SITE2,
         THREADS,
         ANNOVAR_DB,
         ANNOVAR_VERSION,
         STEPS,
         HLA_FASTA):

    # TODO add sanity checks for the parameters
    # TODO better log info

    print("RNA somatic pipeline")

    # Create sub-folder to store all results
    os.makedirs('workdir', exist_ok=True)
    os.chdir('workdir')

    if 'mapping' in STEPS:
        print('Trimming reads')
        cmd = '{} --cores {} --fastqc --paired --basename sample {} {}'.format(TRIMGALORE, THREADS, R1, R2)
        exec_command(cmd)

        # ALIGNMENT
        print('Starting alignment')
        cmd = '{} --genomeDir {} --readFilesIn sample_val_1.fq.gz sample_val_2.fq.gz --outSAMorder Paired' \
              ' --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {}' \
              ' --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c' \
              ' --runThreadN {}'.format(STAR, GENOME_STAR, ANNOTATION, THREADS)
        exec_command(cmd)

        # Add headers
        print("Adding headers")
        cmd = '{} AddOrReplaceReadGroups --INPUT Aligned.sortedByCoord.out.bam --OUTPUT sample_header.bam ' \
              '--SORT_ORDER coordinate --RGID {} --RGPL Illumina --RGLB DNA --RGPU {} --RGSM {} --RGCN {} ' \
              '--CREATE_INDEX true --VALIDATION_STRINGENCY SILENT'.format(PICARD, SAMPLEID, SAMPLEID, SAMPLEID, SAMPLEID)
        exec_command(cmd)

    if 'gatk' in STEPS:
        # Mark duplicates
        print('Marking duplicates')
        cmd = '{} MarkDuplicatesSpark --input sample_header.bam --output sample_dedup.bam'.format(GATK)
        exec_command(cmd)

        # Split N and cigars
        print('Splitting NCigar Reads')
        cmd = '{} SplitNCigarReads --reference {} --input sample_dedup.bam --output sample_split.bam'.format(GATK, GENOME)
        exec_command(cmd)

        # GATK base re-calibration
        print('Starting re-calibration')
        cmd = '{} BaseRecalibratorSpark --use-original-qualities --input sample_split.bam --reference {} --known-sites {} ' \
              '--known-sites {} --known-sites {} --output sample_recal_data.txt'.format(GATK,
                                                                                        GENOME,
                                                                                        SNPSITES,
                                                                                        KNOWN_SITE1,
                                                                                        KNOWN_SITE2)
        exec_command(cmd)
        cmd = '{} ApplyBQSR --use-original-qualities --add-output-sam-program-record --reference {} --input sample_split.bam ' \
              '--bqsr-recal-file sample_recal_data.txt --output sample_final.bam'.format(GATK, GENOME)
        exec_command(cmd)

        # BamQC
        cmd = '{} -bam sample_final.bam -gtf {} --paired -outdir bamQCRNA ' \
              '--java-mem-size=16000M -outformat HTML'.format(BAMQCRNA, ANNOTATION)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} -bam sample_final.bam --genome-gc-distr HUMAN -nt {} ' \
              '-outdir bamQC -outformat HTML'.format(BAMQC, THREADS)
        p2 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()

    if 'hla' in STEPS:
        print('Predicting HLAs')
        HLA_prediction('sample_final.bam', THREADS, "rna", SAMPLEID, HLA_FASTA)

    if 'variant' in STEPS:
        # Variant calling (Samtools pile-ups)
        print('Computing pile-ups')
        cmd = '{} mpileup -C 50 -B -q 1 -Q 15 -f {} sample_final.bam > sample.pileup'.format(SAMTOOLS, GENOME)
        p1 = exec_command(cmd, detach=True)

        # Variant calling (HaplotypeCaller)
        print('Variant calling with HaplotypeCaller')
        cmd = '{} HaplotypeCaller --reference {} --input sample_final.bam --output haplotypecaller.vcf ' \
              '--dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20 ' \
              '--dbsnp {}'.format(GATK, GENOME, SNPSITES)
        p2 = exec_command(cmd, detach=True)

        # Computing gene counts
        print('Computing gene counts with featureCounts')
        cmd = '{} -T {} --primary --ignoreDup -O -C -t exon ' \
              '-g gene_name -a {} -o gene.counts sample_dedup.bam'.format(FEATURECOUNTS, THREADS, ANNOTATION)
        p3 = exec_command(cmd, detach=True)

        # Variant calling VarScan
        p1.wait()
        print('Variant calling with VarScan2')
        cmd = '{} mpileup2cns sample.pileup --variants 0 --min-coverage 2 --min-reads2 1 --output-vcf 1 ' \
              '--min-var-freq 0.01 --min-avg-qual 15 --p-value 0.99 --strand-filter 1 > varscan.vcf'.format(VARSCAN)
        p4 = exec_command(cmd, detach=True)

        # Wait for processes to finish
        p2.wait()
        p3.wait()
        p4.wait()

    if 'filter' in STEPS:
        # Filtering variants (HaplotypeCaller)
        print("Filtering HaplotypeCaller variants")
        cmd = '{} VariantFiltration --reference {} --variant haplotypecaller.vcf --window 35 --cluster 3 --filter-name "FS" ' \
              '--filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" --output haplotype_caller_filtered.vcf'.format(GATK, GENOME)
        exec_command(cmd)

        # NOTE replacing IUPAC codes from VCF
        # NOTE this will also skip variants whose REF and ALT fields are identical
        cmd = 'awk \'{if ($1 ~ /#/) {print} else if ($4 != $5) {gsub(/W|K|B|Y|D|H|V|R|S|M/,"N",$4); OFS="\t"; print}}\' ' \
              'varscan.vcf > varscan_filtered.vcf'
        exec_command(cmd)

        # Combine with GATK
        print('Combining variants')
        # TODO replace this with jacquard merge
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        cmd = '{} -T CombineVariants -R {} -V:varscan varscan_filtered.vcf ' \
              '-V:HaplotypeCaller haplotype_caller_filtered.vcf -o combined_calls.vcf '\
              '-genotypeMergeOptions UNIQUIFY --num_threads {}'.format(GATK3, GENOME, THREADS)
        exec_command(cmd)

        # Replace name of the caller in the VCF file
        cmd = 'sed -i \'s/{}.HaplotypeCaller/HaplotypeCaller/g\' combined_calls.vcf'.format(SAMPLEID)
        exec_command(cmd)

        # Replace name of the caller in the VCF file
        cmd = 'sed -i \'s/Sample1.varscan/varscan/g\' combined_calls.vcf'
        exec_command(cmd)

        # Annotate with Annovar
        print('Annotating variants')
        annotate_variants('combined_calls.vcf', 'annotated', ANNOVAR_DB, ANNOVAR_VERSION, THREADS)

        # Moving result files to output
        shutil.move('combined_calls.vcf', '../combined_calls.vcf')
        shutil.move('annotated.{}_multianno.vcf'.format(ANNOVAR_VERSION),
                    '../annotated.{}_multianno.vcf'.format(ANNOVAR_VERSION))
        shutil.move('rna_{}_hla_genotype_results.tsv'.format(SAMPLEID),
                    '../hla_genotype.tsv')
        shutil.move('gene.counts', '../gene.counts')
        shutil.move('gene.counts.summary', '../{}_gene.counts.summary'.format(SAMPLEID))
        shutil.move('Log.final.out', '../{}_Log.final.out'.format(SAMPLEID))
        shutil.move('sample_final.bam', '../sample_final.bam')
        if os.path.isdir('../{}_bamQC'.format(SAMPLEID)):
            shutil.rmtree(os.path.abspath('../{}_bamQC'.format(SAMPLEID)))
        shutil.move('bamQC', '../{}_bamQC'.format(SAMPLEID))
        if os.path.isdir('../{}_bamQCRNA'.format(SAMPLEID)):
            shutil.rmtree(os.path.abspath('../{}_bamQCRNA'.format(SAMPLEID)))
        shutil.move('bamQCRNA', '../{}_bamQCRNA'.format(SAMPLEID))
        for file in glob.glob('*_fastqc*'):
            shutil.move(file, '../{}_{}'.format(SAMPLEID, file))
        for file in glob.glob('*_trimming_report*'):
            shutil.move(file, '../{}_{}'.format(SAMPLEID, file))

    print("COMPLETED!")

if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('R1', help='FASTQ file R1 (RNA)')
    parser.add_argument('R2', help='FASTQ file R2 (RNA)')
    parser.add_argument('--genome',
                        help='Path to the reference genome FASTA file', required=True)
    parser.add_argument('--genome-star',
                        help='Path to the reference genome STAR index folder', required=False)
    parser.add_argument('--genome-ref',
                        help='Path to the reference genome GTF file', required=False)
    parser.add_argument('--sample',
                        help='Name of the sample/experiment. Default is sample', default='sample')
    parser.add_argument('--outdir',
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
                        choices=['mapping', 'gatk', 'hla', 'variant', 'filter'])
    parser.add_argument("--hla-fasta", type=str, default="hla_reference_rna.fasta", required=False, 
                        help="Path to the rna hla reference fasta file located in shared. (default: hla_reference_rna.fasta)")

    # Parse arguments
    args = parser.parse_args()
    DIR = args.outdir
    R1 = os.path.abspath(args.R1)
    R2 = os.path.abspath(args.R2)
    SAMPLEID = args.sample
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
    HLA_FASTA = args.hla_fasta

    # Move to output dir
    os.makedirs(os.path.abspath(DIR), exist_ok=True)
    os.chdir(os.path.abspath(DIR))

    main(R1,
         R2,
         SAMPLEID,
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
         HLA_FASTA)