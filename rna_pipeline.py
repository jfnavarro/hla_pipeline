#! /usr/bin/env python
"""
This pipeline computes somatic variants from RNA data.

The pipeline trims with trimgalore, aligns with STAR,
performs the GATK4 best practices and computes variants with
HaplotypeCaller and Varscan. The variants are then combined into
one file and annotated with VEP. Gene counts are also
computed with featureCounts. HLA typing is performed with Optitype.

Multiple options are available. To see them type --help

@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from hlapipeline.common import *
from hlapipeline.tools import *
from hlapipeline.version import version_number
from pathlib import Path
import os
import shutil
import glob
import logging
import datetime
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
         ASSEMBLY,
         VERSION,
         CACHEDIR,
         STEPS,
         HLA_FASTA,
         KEEP,
         SPARK):
    # TODO add sanity checks for the parameters
    if 'filter' in STEPS:
        if not CACHEDIR:
            if not Path.home().joinpath('.vep').exists():
                raise Exception('Cache directory doesn\'t exist at default location, please provide a valid path.')
        elif not Path(CACHEDIR).exists():
            raise Exception('The cache directory provided doesn\'t exist. Please provide a different one.')

    logging.basicConfig(format='%(asctime)s - %(message)s', 
                        datefmt='%d-%b-%y %H:%M:%S',
                        level=logging.DEBUG, filename=SAMPLEID + ".log")
    logger = logging.getLogger(SAMPLEID)

    start_pipeline_time = datetime.datetime.now()

    logger.info('Starting RNA somatic pipeline: {}'.format(start_pipeline_time))
    logger.info('HLA Pipeline version: {}'.format(version_number))
    logger.info('Processing FASTQs {} and {} with Sample ID {} using ' \
                'genome version {}.'.format(R1, R2, SAMPLEID, ASSEMBLY))

    # Create sub-folder to store all results
    os.makedirs('workdir', exist_ok=True)
    os.chdir('workdir')

    if 'mapping' in STEPS:
        start_map_time = datetime.datetime.now()
        logger.info('Starting trimming and mapping step: {}'.format(start_map_time))
        logger.info('Trimming reads')
        cmd = '{} --cores {} --fastqc --paired --basename sample {} {}'.format(TRIMGALORE, THREADS, R1, R2)
        exec_command(cmd)

        # ALIGNMENT
        logger.info('Starting alignment')
        cmd = '{} --genomeDir {} --readFilesIn sample_val_1.fq.gz sample_val_2.fq.gz --outSAMorder Paired' \
              ' --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {}' \
              ' --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c' \
              ' --runThreadN {}'.format(STAR, GENOME_STAR, ANNOTATION, THREADS)
        exec_command(cmd)

        # Add headers
        logger.info("Adding headers")
        cmd = '{} AddOrReplaceReadGroups --INPUT Aligned.sortedByCoord.out.bam --OUTPUT sample_header.bam ' \
              '--SORT_ORDER coordinate --RGID {} --RGPL Illumina --RGLB DNA --RGPU {} --RGSM {} --RGCN {} ' \
              '--CREATE_INDEX true --VALIDATION_STRINGENCY SILENT'.format(PICARD, SAMPLEID, SAMPLEID, SAMPLEID,
                                                                          SAMPLEID)
        exec_command(cmd)

        if not KEEP:
            if os.path.isfile('Aligned.sortedByCoord.out.bam'):
                os.remove('Aligned.sortedByCoord.out.bam')
            if os.path.isfile('sample_val_1.fq.gz'):
                os.remove('sample_val_1.fq.gz')
            if os.path.isfile('sample_val_2.fq.gz'):
                os.remove('sample_val_2.fq.gz')
        
        end_map_time = datetime.datetime.now()
        total_map_time = end_map_time - start_map_time
        logger.info('Total trimming and mapping execution time: {}'.format(total_map_time))

    if 'gatk' in STEPS:

        start_gatk_time = datetime.datetime.now()
        logger.info('Starting GATK steps: {}'.format(start_gatk_time))

        # Mark duplicates
        logger.info('Marking duplicates')

        if SPARK:
            cmd = '{} --java-options "-Xmx32g" MarkDuplicatesSpark -I sample_header.bam -O sample_dedup.bam'.format(GATK)

        else:
            cmd = '{} --java-options "-Xmx32g" MarkDuplicates -I sample_header.bam -O sample_dedup.bam ' \
                    '--CREATE_INDEX true -M sample_dup_metrics.txt'.format(GATK)

        exec_command(cmd)

        # Split N and cigars
        logger.info('Splitting NCigar Reads')
        cmd = '{} SplitNCigarReads --reference {} --input sample_dedup.bam --output sample_split.bam'.format(GATK,
                                                                                                             GENOME)
        exec_command(cmd)

        # GATK base re-calibration
        logger.info('Starting re-calibration')

        recal_cmd = 'BaseRecalibratorSpark' if SPARK else 'BaseRecalibrator'

        cmd = '{} --java-options "-Xmx32g" {} --use-original-qualities --input sample_split.bam --reference {} --known-sites {} ' \
              '--known-sites {} --known-sites {} --output sample_recal_data.txt'.format(GATK,
                                                                                        recal_cmd,
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
              '--java-mem-size=16G -outformat HTML'.format(BAMQCRNA, ANNOTATION)
        p1 = exec_command(cmd, detach=True)

        cmd = '{} -bam sample_final.bam --genome-gc-distr HUMAN -nt {} ' \
              '--java-mem-size=16G -outdir bamQC -outformat HTML'.format(BAMQC, THREADS)
        p2 = exec_command(cmd, detach=True)

        # Wait for the processes to finish in parallel
        p1.wait()
        p2.wait()

        if not KEEP:
            if os.path.isfile('sample_recal_data.txt'):
                os.remove('sample_recal_data.txt')
            if os.path.isfile('sample_split.bam'):
                os.remove('sample_split.bam')
            if os.path.isfile('sample_split.bai'):
                os.remove('sample_split.bai')

        end_gatk_time = datetime.datetime.now()
        total_gatk_time = end_gatk_time - start_gatk_time
        logger.info('Total GATK processing time: {}'.format(total_gatk_time))

    if 'hla' in STEPS:

        start_hla_time = datetime.datetime.now()
        logger.info('Starting HLA prediction: {}'.format(start_hla_time))

        HLA_prediction('sample_final.bam', THREADS, 'rna', SAMPLEID, HLA_FASTA, 'rna', KEEP)

        end_hla_time = datetime.datetime.now()
        total_hla_time = end_hla_time - start_hla_time
        logger.info('Total HLA prediction time: {}'.format(total_hla_time))

    if 'variant' in STEPS:

        start_variant_time = datetime.datetime.now()
        logger.info('Starting variant calling: {}'.format(start_variant_time))

        # Variant calling (Samtools pile-ups)
        logger.info('Computing pile-ups')
        cmd = '{} mpileup -C 50 -B -q 1 -Q 15 -f {} sample_final.bam > sample.pileup'.format(SAMTOOLS, GENOME)
        p1 = exec_command(cmd, detach=True)

        # Variant calling (HaplotypeCaller)
        logger.info('Variant calling with HaplotypeCaller')
        cmd = '{} HaplotypeCaller --reference {} --input sample_final.bam --output haplotypecaller.vcf ' \
              '--dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20 ' \
              '--dbsnp {}'.format(GATK, GENOME, SNPSITES)
        p2 = exec_command(cmd, detach=True)

        # Computing gene counts
        logger.info('Computing gene counts with featureCounts')
        cmd = '{} -T {} --primary --ignoreDup -O -C -t exon ' \
              '-g gene_name -a {} -o gene.counts sample_dedup.bam'.format(FEATURECOUNTS, THREADS, ANNOTATION)
        p3 = exec_command(cmd, detach=True)

        # Variant calling VarScan
        p1.wait()
        logger.info('Variant calling with VarScan2')
        cmd = '{} mpileup2cns sample.pileup --variants 0 --min-coverage 2 --min-reads2 1 --output-vcf 1 ' \
              '--min-var-freq 0.01 --min-avg-qual 15 --p-value 0.99 --strand-filter 1 > varscan.vcf'.format(VARSCAN)
        p4 = exec_command(cmd, detach=True)

        # Wait for processes to finish
        p2.wait()
        p3.wait()
        p4.wait()

        if not KEEP:
            if os.path.isfile('sample.pileup'):
                os.remove('sample.pileup')
            if os.path.isfile('sample_dedup.bam'):
                os.remove('sample_dedup.bam')
            if os.path.isfile('sample_dedup.bam.bai'):
                os.remove('sample_dedup.bam.bai')

        end_variant_time = datetime.datetime.now()
        total_variant_time = end_variant_time - start_variant_time
        logger.info('Total variant calling processing time: {}'.format(total_variant_time))

    if 'filter' in STEPS:

        start_filter_time = datetime.datetime.now()
        logger.info('Starting variant filtering and annotation: {}'.format(start_filter_time))
        
        # Filtering variants (HaplotypeCaller)
        logger.info("Filtering HaplotypeCaller variants")
        cmd = '{} VariantFiltration --reference {} --variant haplotypecaller.vcf --window 35 --cluster 3 --filter-name "FS" ' \
              '--filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" --output haplotype_caller_filtered.vcf'.format(
            GATK, GENOME)
        exec_command(cmd)

        # NOTE replacing IUPAC codes from VCF
        # NOTE this will also skip variants whose REF and ALT fields are identical
        cmd = 'awk \'{if ($1 ~ /#/) {print} else if ($4 != $5) {gsub(/W|K|B|Y|D|H|V|R|S|M/,"N",$4); OFS="\t"; print}}\' ' \
              'varscan.vcf > varscan_filtered.vcf'
        exec_command(cmd)

        # Combine with GATK
        logger.info('Combining variants')
        # TODO replace this with jacquard merge
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        cmd = '{} -T CombineVariants -R {} -V:varscan varscan_filtered.vcf ' \
              '-V:HaplotypeCaller haplotype_caller_filtered.vcf -o combined_calls.vcf ' \
              '-genotypeMergeOptions UNIQUIFY --num_threads {}'.format(GATK3, GENOME, THREADS)
        exec_command(cmd)

        # Replace name of the caller in the VCF file
        cmd = 'sed -i \'s/{}.HaplotypeCaller/HaplotypeCaller/g\' combined_calls.vcf'.format(SAMPLEID)
        exec_command(cmd)

        # Replace name of the caller in the VCF file
        cmd = 'sed -i \'s/Sample1.varscan/varscan/g\' combined_calls.vcf'
        exec_command(cmd)

        # Annotate with VEP

        logger.info('Annotating variants')
        annotate_variants('combined_calls.vcf', ASSEMBLY, VERSION, THREADS, GENOME_REF, CACHEDIR)

        # Summary of basic statistic of the annotated VCF file
        annotated_vcf = "annotated.{}_multianno.vcf".format(ASSEMBLY)
        vcf_stats(annotated_vcf, SAMPLEID)

        end_filer_time = datetime.datetime.now()
        total_filter_time = end_filer_time - start_filter_time
        logger.info('Total filtering and annotation time: {}'.format(total_filter_time))

        # Moving result files to output
        if os.path.isfile('{}.relatedness2'.format(SAMPLEID)):
            shutil.move('{}.relatedness2'.format(SAMPLEID), '../{}.relatedness2'.format(SAMPLEID))
        if os.path.isfile('{}.TsTv.summary'.format(SAMPLEID)):
            shutil.move('{}.TsTv.summary'.format(SAMPLEID), '../{}.TsTv.summary'.format(SAMPLEID))
        if os.path.isfile('{}.vchk'.format(SAMPLEID)):
            shutil.move('{}.vchk'.format(SAMPLEID), '../{}.vchk'.format(SAMPLEID))
        if os.path.isfile('combined_calls.vcf'):
            shutil.move('combined_calls.vcf', '../combined_calls.vcf')
        if os.path.isfile('annotated.{}_multianno.vcf'.format(ASSEMBLY)):
            shutil.move('annotated.{}_multianno.vcf'.format(ASSEMBLY),
                        '../annotated.{}_multianno.vcf'.format(ASSEMBLY))
        if os.path.isfile('rna_{}_hla_genotype_result.tsv'.format(SAMPLEID)):
            shutil.move('rna_{}_hla_genotype_result.tsv'.format(SAMPLEID), '../hla_genotype.tsv')
        if os.path.isfile('gene.counts'):
            shutil.move('gene.counts', '../gene.counts')
        if os.path.isfile('gene.counts.summary'):
            shutil.move('gene.counts.summary', '../{}_gene.counts.summary'.format(SAMPLEID))
        if os.path.isfile('Log.final.out'):
            shutil.move('Log.final.out', '../{}_Log.final.out'.format(SAMPLEID))
        if os.path.isfile('sample_final.bam'):
            shutil.move('sample_final.bam', '../sample_final.bam')
        if os.path.isfile('sample_final.bai'):
            shutil.move('sample_final.bai', '../sample_final.bai')
        if os.path.isdir('../{}_bamQC'.format(SAMPLEID)):
            shutil.rmtree(os.path.abspath('../{}_bamQC'.format(SAMPLEID)))
        if os.path.isdir('bamQC'):
            shutil.move('bamQC', '../{}_bamQC'.format(SAMPLEID))
        if os.path.isdir('../{}_bamQCRNA'.format(SAMPLEID)):
            shutil.rmtree(os.path.abspath('../{}_bamQCRNA'.format(SAMPLEID)))
        if os.path.isdir('bamQCRNA'):
            shutil.move('bamQCRNA', '../{}_bamQCRNA'.format(SAMPLEID))
        for f in glob.glob('*_fastqc*'):
            shutil.move(f, '../{}_{}'.format(SAMPLEID, f))
        for f in glob.glob('*_trimming_report*'):
            shutil.move(f, '../{}_{}'.format(SAMPLEID, f))

        end_pipeline_time = datetime.datetime.now()
        total_pipeline_time = end_pipeline_time - start_pipeline_time
        logger.info('Total pipeline execution time: {}'.format(total_pipeline_time))

        logger.info("COMPLETED!")


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('R1', help='FASTQ file R1 (RNA)')
    parser.add_argument('R2', help='FASTQ file R2 (RNA)')
    parser.add_argument('--genome',
                        help='Path to the reference genome FASTA file', required=True)
    parser.add_argument('--genome-star',
                        help='Path to the reference genome STAR index folder', required=False, default='')
    parser.add_argument('--genome-ref',
                        help='Path to the reference genome GTF file', required=False, default='')
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
    parser.add_argument('--vep-db', type=str, default='GRCh38', required=False,
                        help='Genome assembly version to be used in VEP (default: GRCh38)')
    parser.add_argument('--vep-version', type=str, default='102', required=False,
                        help='Ensembl version number to be used in VEP (default: 102)')
    parser.add_argument('--vep-dir', type=str, default=None, required=False,
                        help='Path to the VEP cache directory (default: $HOME/.vep)')
    parser.add_argument("--hla-fasta", type=str, default=None, required=True,
                        help="Path to the HLA reference FASTA file to be used in OptiType (HLA)")
    parser.add_argument('--threads',
                        help='Number of threads to use in the parallel steps', type=int, default=10, required=False)
    parser.add_argument('--steps', nargs='+', default=['mapping', 'gatk', 'hla', 'variant', 'filter'],
                        help='Steps to perform in the pipeline',
                        choices=['mapping', 'gatk', 'hla', 'variant', 'filter'])
    parser.add_argument('--keep-intermediate', default=False, action='store_true', required=False,
                        help='Do not remove temporary files')
    parser.add_argument('--use-gatk-spark', default=False, action='store_true', required=False,
                        help='Enable the use of Spark in MarkDuplicates and BaseRecalibrator (GATK)')

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
    ASSEMBLY = args.vep_db
    VERSION = args.vep_version
    CACHEDIR = os.path.abspath(args.vep_dir) if args.vep_dir else None
    HLA_FASTA = os.path.abspath(args.hla_fasta)
    KEEP = args.keep_intermediate
    SPARK = args.use_gatk_spark

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
         ASSEMBLY,
         VERSION,
         CACHEDIR,
         STEPS,
         HLA_FASTA,
         KEEP,
         SPARK)
