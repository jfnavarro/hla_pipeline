"""
@author: jfnavarro
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

def HLA_predictionDNA(bamfile, sampleID, outfile, threads):
    OUT_DIR = os.path.abspath(os.path.join('out_hla_', os.path.splitext(outfile)[0]))
    os.makedirs(OUT_DIR, exist_ok=True)
    cmd = '{} --BAM {} --workingDir {} --graph {} --sampleID {}'\
          ' --maxTHREADS {}'.format(HLALA, bamfile, OUT_DIR, 'PRG_MHC_GRCh38_withIMGT', sampleID, threads)
    exec_command(cmd)

    # create a dictionary to store the output for each allele
    hla = pd.read_table(os.path.join(OUT_DIR, sampleID, 'hla', 'R1_bestguess_G.txt'), sep='\t')
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

def annotate_varinats(input, output, ANNOVAR_DB, ANNOVAR_VERSION, annovar_anno, THREADS):
    # Annotate with Annovar
    annovardb = '{} -buildver {}'.format(os.path.join(ANNOVAR_PATH, ANNOVAR_DB), ANNOVAR_VERSION)
    cmd = '{} {} {} -thread {} -out output -vcfinput -remove -protocol {}'.format(
        os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), annovardb, input, THREADS, output, annovar_anno)
    exec_command(cmd)