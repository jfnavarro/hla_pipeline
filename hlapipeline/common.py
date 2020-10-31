"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com
"""
import subprocess
import sys
from hlapipeline.tools import *
import pandas as pd
import re
import os


def exec_command(cmd):
    print(cmd)
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    output, error = p.communicate()
    if p.returncode != 0:
        for line in output.decode("utf-8").split("\n") if output else "":
            print(line.rstrip())
        for line in error.decode("utf-8").split("\n") if error else "":
            print(line.rstrip())
        sys.exit(-1)


def HLA_predictionDNA(inputbam, sampleID, outfile, threads):
    """
    Performs HLA typing with HLA-LA for aligned DNA data.
    The output of HLA is reformatted to make it easier to use.
    :param inputbam: aligned DNA data
    :param sampleID: the ID to give to the sample
    :param outfile: the name of the reformatted output file
    :param threads: the number of threads to use
    """
    # HLA-LA only supports alphanumeric
    sampleID_clean = ''.join(ch for ch in sampleID if ch.isalnum())
    OUT_DIR = os.path.abspath(os.path.join('out_hla_', os.path.splitext(outfile)[0]))
    os.makedirs(OUT_DIR, exist_ok=True)
    cmd = '{} --BAM {} --workingDir {} --graph {} --sampleID {}' \
          ' --maxTHREADS {}'.format(HLALA, inputbam, OUT_DIR, 'PRG_MHC_GRCh38_withIMGT', sampleID_clean, threads)
    exec_command(cmd)

    # create a dictionary to store the output for each allele
    hla = pd.read_table(os.path.join(OUT_DIR, sampleID_clean, 'hla', 'R1_bestguess_G.txt'), sep='\t')
    allele_dict = {}
    hla = hla.groupby('Locus')
    for k, g in hla:
        allele = [re.sub('G$', '', a).replace('N', '') for a in g['Allele'].tolist()]
        allele_dict['HLA_' + k] = allele

    # Create formatted output file
    a = open(outfile, 'w')
    for key, value in allele_dict.items():
        a.write('{}\t{}\t{}\n'.format(sampleID, key, '\t'.join(value)))
    a.close()


def HLA_predictionRNA(inputbam, threads):
    """
    Performs HLA typing with arcasHLA for RNA data.
    :param inputbam: BAM file with aligned reads
    :param threads: the number of threads
    """
    cmd = '{} extract --threads {} --paired {}'.format(ARCASHLA, threads, inputbam)
    exec_command(cmd)

    clean_name = os.path.splitext(inputbam)[0]

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


def annotate_variants(input, output, db, version, threads):
    """
    Annotate a VCF using Annovar
    :param input: the VCF file
    :param output: the annotated VCF file
    :param db: the species (humandb or mousedb)
    :param version: the version (hg19 or hg38)
    :param threads: the number of threads to use
    """
    annovardb = '{} -buildver {}'.format(os.path.join(ANNOVAR_PATH, db), version)
    cmd = '{} {} {} -thread {} -out output -vcfinput -remove -protocol {}'.format(
        os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), annovardb, input, threads, output, ANNOVAR_ANNO)
    exec_command(cmd)
