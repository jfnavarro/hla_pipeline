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
VEP = 'vep'
VEP_OPTIONS = '--no_escape --offline --no_stats --af_gnomad --hgvs --cache --fields "Allele,Consequence,SYMBOL,Gene,Feature_type,Feature,BIOTYPE,EXON,INTRON,HGVSc,HGVSp,' \
    'cDNA_position,CDS_position,Protein_position,Existing_variation,FLAGS,gnomAD_AF"' # Could add --pick
VCFTOOLS = 'vcftools'
BGZIP = 'bgzip'
TABIX = 'tabix'
BCFTOOLS = 'bcftools'
