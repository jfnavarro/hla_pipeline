"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com
"""
import os

# TODO Some tools are expected to be in an Anaconda environment, make it more generic

PICARD = 'picard'
GATK = 'gatk'
GATK3 = 'java -jar ' + os.path.join(os.environ['GATK3_PATH'])
VARSCAN = 'varscan'
STRELKA = os.path.join(os.environ['STRELKA_PATH'], 'bin', 'configureStrelkaSomaticWorkflow.py')
SAMTOOLS = 'samtools'
SSNIPER = 'bam-somaticsniper'
HLALA = 'HLA-LA.pl'
SAMTOOLS = 'samtools'
STAR = 'STAR'
TRIMGALORE = 'trim_galore'
BWA = 'bwa mem'
ARCASHLA = os.path.join(os.environ['ARCASHLA_PATH'], 'arcasHLA')
FEATURECOUNTS = 'featureCounts'

# ANNOVAR location must be in $ANNOVAR
ANNOVAR_PATH = os.environ['ANNOVAR_PATH']
ANNOVAR_ANNO = 'refGene,knownGene,ensGene,avsnp150,ALL.sites.2015_08,EUR.sites.2015_08,'\
    'AMR.sites.2015_08,EAS.sites.2015_08,AFR.sites.2015_08,cosmic70 -operation g,g,g,f,f,f,f,f,f,f -nastring .'