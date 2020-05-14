# hla_pipeline
This a private repository that contains a simplified and more generic version 
of Jared's pipeline to process WES (tumor and norma) and RNA-seq data (normal).

There are 2 pipelines and 2 tools:

**exome_pipeline.py** processes the WES data and generates a list of unified
filtered and annotated variants and their epitotes. The variant callers
are Mutect2, Strelka, Varscan and SomaticSniper and both indels and snvs are
kept. Annotation is done using annovar. 
The pipeline uses bwa to align, trimmomatic to quality trim and GATK4
best practices. The filters are applied based on allele frecuency and reads
frecuency. The pipeline will also perform HLA predictions with HLA-LA.

**rnaseq_pipeline.py** processes the RNA-seq data and generates a list of unified
annotated variants (unfiltered) and also a list of FPKM values. The variant
callers used are Varscan and mpileup. Annmotation is done with annovar.
The pipeline uses trim-galore to trim,
STAR to align and GATK4 best practices. The FPKM values are computed with
cufflinks. The pipeline will also perform HLA predictions with arcasHLA.

**merge_resuls.py** combines results from 1 or several runs of the exome and rna-seq
pipelines in order to generate an unified Excel sheet with useful information where
variants are filtered by occurrences (2 or more callers). 

**mhc_predict.py** can take the file generated with merge_results.py and the HLA files
generated in the exome and rna-seq pipeline and output a list of predicted MHCs
using the most common alleles for each HLA class 1. 

Each tool/pipeline uses a command line interface with parameters which
can be shown and described with --help.

See INSTALL.txt for installation instructions. 

See RUN.txt for running instructions and explanation on the different
files generated. 

Contact: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>


