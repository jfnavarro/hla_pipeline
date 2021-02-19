"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com
"""
import os

# TODO Some tools are expected to be in an Anaconda environment, make it more generic

PICARD = 'picard'
GATK = 'gatk'
GATK3 = 'java -jar ' + os.path.join(os.environ['GATK3_PATH'])
VARSCAN = 'varscan'
STRELKA = os.path.join(os.environ['STRELKA_PATH'], 'bin', 'configureStrelkaSomaticWorkflow.py')
SAMTOOLS = 'samtools'
SSNIPER = 'bam-somaticsniper'
SAMTOOLS = 'samtools'
STAR = 'STAR'
TRIMGALORE = 'trim_galore'
BWA = 'bwa mem'
FEATURECOUNTS = 'featureCounts'
BAMQC = 'qualimap bamqc'
BAMQCRNA = 'qualimap rnaseq'
OPTITYPE = 'OptiTypePipeline.py'
YARAI = 'yara_indexer'
YARAM = 'yara_mapper'
# ANNOVAR location must be in $ANNOVAR_PATH
ANNOVAR_PATH = os.environ['ANNOVAR_PATH']
ANNOVAR_ANNO = 'ensGene,avsnp150,gnomad211_exome,cosmic70 -operation g,f,f,f -nastring .'
VCFTOOLS = 'vcftools'
BGZIP = 'bgzip'
TABIX = 'tabix'
BCFTOOLS = 'bcftools'
