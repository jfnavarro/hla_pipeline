"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com
"""
import subprocess
import sys
from hlapipeline.tools import *
import os
import logging


def exec_command(cmd, detach=False):
    logger = logging.getLogger()
    logger.info(cmd)
    if detach:
        return subprocess.Popen(cmd, shell=True)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = p.communicate()
    if p.returncode != 0:
        for line in output.decode("utf-8").split("\n") if output else "":
            logger.error(line.rstrip())
        for line in error.decode("utf-8").split("\n") if error else "":
            logger.error(line.rstrip())
        sys.exit(-1)


def HLA_prediction(inputbam, threads, origin, sample, fasta, nacid, KEEP):
    """
    Performs HLA typing with OptiType for either RNA or DNA data.
    :param inputbam: BAM file with aligned reads
    :param threads: the number of threads
    :param origin: prefix for the output hla genotype file
    :param sample: prefix for the sample name
    :param fasta: HLA reference fasta file
    :param KEEP: do not discard temp files when True
    """

    # We do not need many threads for Samtools view
    SAMTOOLS_THREADS = max(int(threads / 2), 1)

    # TODO use os.makedirs instead
    cmd = 'mkdir -p {}/index'.format(os.getcwd())
    exec_command(cmd)

    cmd = '{} {} -o {}/index/hla_reference'.format(YARAI, fasta, os.getcwd())
    exec_command(cmd)

    cmd = '{} view -@ {} -h -f 0x40 {} > {}_output_1.bam'.format(SAMTOOLS, SAMTOOLS_THREADS, inputbam, origin)
    p1 = exec_command(cmd, detach=True)

    cmd = '{} view -@ {} -h -f 0x80 {} > {}_output_2.bam'.format(SAMTOOLS, SAMTOOLS_THREADS, inputbam, origin)
    p2 = exec_command(cmd, detach=True)

    p1.wait()
    p2.wait()

    cmd = '{} bam2fq {}_output_1.bam > {}_output_1.fastq'.format(SAMTOOLS, origin, origin)
    p1 = exec_command(cmd, detach=True)

    cmd = '{} bam2fq {}_output_2.bam > {}_output_2.fastq'.format(SAMTOOLS, origin, origin)
    p2 = exec_command(cmd, detach=True)

    p1.wait()
    p2.wait()

    if not KEEP:
        if os.path.isfile('{}_output_1.bam'.format(origin)):
            os.remove('{}_output_1.bam'.format(origin))
        if os.path.isfile('{}_output_2.bam'.format(origin)):
            os.remove('{}_output_2.bam'.format(origin))

    cmd = '{} -e 3 -t {} -f bam {}/index/hla_reference {}_output_1.fastq {}_output_2.fastq > {}_output.bam'.format(
        YARAM, threads, os.getcwd(), origin, origin, origin)
    exec_command(cmd)

    if not KEEP:
        if os.path.isfile('{}_output_1.fastq'.format(origin)):
            os.remove('{}_output_1.fastq'.format(origin))
        if os.path.isfile('{}_output_2.fastq'.format(origin)):
            os.remove('{}_output_2.fastq'.format(origin))

    cmd = '{} view -@ {} -h -F 4 -f 0x40 -b1 {}_output.bam > {}_mapped_1.bam'.format(SAMTOOLS,
                                                                                     SAMTOOLS_THREADS,
                                                                                     origin,
                                                                                     origin)
    p1 = exec_command(cmd, detach=True)

    cmd = '{} view -@ {} -h -F 4 -f 0x80 -b1 {}_output.bam > {}_mapped_2.bam'.format(SAMTOOLS,
                                                                                     SAMTOOLS_THREADS,
                                                                                     origin,
                                                                                     origin)
    p2 = exec_command(cmd, detach=True)

    p1.wait()
    p2.wait()

    cmd = '{} --input {}_mapped_1.bam {}_mapped_2.bam --{} --prefix {}_{}_hla_genotype --outdir {}'.format(
        OPTITYPE, origin, origin, nacid, origin, sample, os.getcwd())
    exec_command(cmd)

    if not KEEP:
        if os.path.isfile('{}_output.bam'.format(origin)):
            os.remove('{}_output.bam'.format(origin))
        if os.path.isfile('{}_mapped_1.bam'.format(origin)):
            os.remove('{}_mapped_1.bam'.format(origin))
        if os.path.isfile('{}_mapped_2.bam'.format(origin)):
            os.remove('{}_mapped_2.bam'.format(origin))


def annotate_variants(input, db, version, threads, fasta, cache):
    """
    Annotate a VCF using VEP
    :param input: the VCF file
    :param output: the annotated VCF file
    :param db: the genome assembly (GRCh37, GRCh38)
    :param version: the ensembl version (75, 102, etc..)
    :param threads: the number of threads to use
    :param cache: the location of the VEP cache, can be None (default location)
    """
    cache_cmd = '--dir_cache {}'.format(cache) if cache is not None else ''
    cmd = '{} -i {} --fork {} -o annotated.{}_multianno.vcf --fasta {} --format vcf --vcf --assembly {} '\
        '--cache_version {} --species homo_sapiens {} {}'.format(VEP, input, threads, db, fasta, db, version, VEP_OPTIONS, cache_cmd)
    exec_command(cmd)


def vcf_stats(annotated_VCF, sampleID):
    """
    Performs summary of basic statistics of annotated VCF file using vcftools and bcftools
    :param annotated_VCF: annotated VCF file
    :param sampleID: the ID to give to the sample
    """

    # VCFtools: pairwise individual relatedness using relatedness2 method
    cmd = '{} --vcf {} --relatedness2 --out {}'.format(VCFTOOLS, annotated_VCF, sampleID)
    exec_command(cmd)

    # VCFtools: summary of all Transitions and Transversions
    cmd = '{} --vcf {} --TsTv-summary --out {}'.format(VCFTOOLS, annotated_VCF, sampleID)
    exec_command(cmd)

    # bcftools multiple stats
    cmd = '{} -c {} > {}.gz'.format(BGZIP, annotated_VCF, annotated_VCF)
    exec_command(cmd)

    cmd = '{} -p vcf {}.gz'.format(TABIX, annotated_VCF)
    exec_command(cmd)

    cmd = '{} stats {}.gz > {}.vchk'.format(BCFTOOLS, annotated_VCF, sampleID)
    exec_command(cmd)
