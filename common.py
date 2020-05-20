import subprocess
import os
import sys
import re
import pandas as pd
from Bio.Seq import translate

PICARD = 'picard'
GATK = 'gatk'
GATK3 = 'java -jar ' + os.path.join(os.environ['GATK3_PATH'])
VARSCAN = 'varscan'
STRELKA = os.path.join(os.environ['STRELKA_PATH'], 'bin', 'configureStrelkaSomaticWorkflow.py')
SAMTOOLS = 'samtools'
SSNIPER = 'bam-somaticsniper'
HLALA = 'HLA-LA.pl'
SAMTOOLS = 'samtools'
CUFFLINKS = 'cufflinks'
STAR = 'STAR'
TRIMGALORE = 'trim_galore'
TRIPTOMATIC = 'trimmomatic'
BWA = 'bwa mem'
ARCASHLA = os.path.join(os.environ['ARCASHLA_PATH'], 'arcasHLA')
VCFTOOLS = 'vcftools'

# ANNOVAR location must be in $ANNOVAR
ANNOVAR_PATH = os.environ['ANNOVAR_PATH']
annovar_anno = 'refGene,knownGene,ensGene,avsnp150,ALL.sites.2015_08,EUR.sites.2015_08,AMR.sites.2015_08,EAS.sites.2015_08,AFR.sites.2015_08,cosmic70 -operation g,g,g,f,f,f,f,f,f,f -nastring NA'

def index_column_substring(your_list, substring):
    for i, s in enumerate(your_list):
        if substring in s:
            return i
    return -1

def translate_dna(seq):
    return translate(seq)

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

def HLA_LA(bamfile, sampleID, outfile, threads):
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
    for key,value in allele_dict.items():
        a.write('{}\t{}\t{}\n'.format(sampleID, key, '\t'.join(value)))
    a.close()
