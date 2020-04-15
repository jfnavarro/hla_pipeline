import subprocess
import os
import sys
import re
import datetime

PICARD = 'picard'
GATK = 'gatk'
GATK3 = 'java -jar ' + os.path.join(os.environ['GATK3_PATH'])
VARSCAN = 'varscan'
STRELKA = os.path.join(os.environ['STRELKA_PATH'], 'bin', 'configureStrelkaSomaticWorkflow.py')
SAMTOOLS = 'samtools'
SSNIPER = 'bam-somaticsniper'
HLA = "HLA-LA.pl"
SAMTOOLS = 'samtools'
CUFFLINKS = 'cufflinks'
HLA = 'python PHLAT.py'
HLA_INDEX = '~/shared/index4phlat'
HLA_PATH = "~/phlat"
BOWTIE2 = 'bowtie2'
STAR = 'STAR'
TRIMGALORE = 'trim_galore'

# ANNOVAR location must be in $ANNOVAR
ANNOVAR_PATH = os.environ['ANNOVAR_PATH']
annovar_db = os.path.join(ANNOVAR_PATH, 'humandb') + ' -buildver hg19'
annovar_anno = 'refGene,knownGene,ensGene,snp138NonFlagged,ALL.sites.2015_08,EUR.sites.2015_08,AMR.sites.2015_08,EAS.sites.2015_08,AFR.sites.2015_08,cosmic70 -operation g,g,g,f,f,f,f,f,f,f -nastring NA'

MRN = "MRN"
SEQ_CENTER = "VHIO"
LIBRARY = "Library"
SOURCE = "Source"
SAMPLE_NOTE = "Sample"
RESECTION_DATE = "Na"
RUN_DATE = "Na"
SEQUENCER = "Na"
KIT = "Na"
NOTE = "Na"
INDEX = "Na"

FASTA_AA = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def index_column_substring(your_list, substring):
    for i, s in enumerate(your_list):
        if substring in s:
            return i
    return -1

def translate_dna(seq):
    rna = ""
    # Generate the RNA string
    for i in seq:
        # Replace all occurrences of T with U
        if i == "T":
            rna += "U"
        else:
            rna += i
    return rna

def translate_dna_to_protein(seq):
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i:i + 3]
            protein += FASTA_AA[codon]
    return protein

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

# This function will predict HLAs using PHLAT one sample
def HLA_prediction(sample1, sample2, threads, outfile):
    cmd = "{} -1 {} -2 {} -index {} -b2url {} -tag normal -e {} -o . -p {}".format(HLA,
                                                                                   sample1,
                                                                                   sample2,
                                                                                   BOWTIE2,
                                                                                   HLA_PATH,
                                                                                   threads)
    exec_command(cmd)

def HLA_PRG(bamfile, sampleID, outfile, threads):
    OUT_DIR = "out_hla"
    cmd = HLA + ' --BAM {} --workingDir {} --graph {} --sampleID {}'\
          + ' --maxTHREADS {}'.format(bamfile, OUT_DIR, 'PRG_MHC_GRCh38_withIMGT', sampleID, threads)
    exec_command(cmd)

    # create a dictionary to store the output for each allele
    hla = pd.read_table(os.path.join(OUT_DIR, sampleID, 'hla', 'R1_bestguess_G.txt'), sep='\t')
    allele_dict = {}
    hla = hla.groupby('Locus')
    for k, g in hla:
        allele = [re.sub('G$', '', a).replace('N', '') for a in g['Allele'].tolist()]
        allele_dict['HLA_' + k] = allele

    # Create formatted output file
    today = datetime.now().strftime('%B_%d_%Y')
    a = open(outfile, 'w')
    for x in sorted(allele_dict):
        a.write('{}\t{}\tExome {}\t{}\t{}\t{}\tPRG-HLA-LA\t-\t-\n'.format(MRN,
                                                                          sampleID,
                                                                          today,
                                                                          x,
                                                                          allele_dict[x][0],
                                                                          allele_dict[x][1]))
    a.close()