# DNA (WES and WGS) and RNA-Seq variant calling pipeline with HLA and MHC predictions
This a private repository that contains a simplified, more generic, up-to-date
and optimized version of a pipeline provided by Jared Gartner <jared.gartner@nhi.gov>
The pipeline can process WES/WGS (tumor and norma) and RNA-seq data (tumor)
and generate a list of variants, HLAs and MHCs. 

The DNA and RNA pipelines use the GATK4 best-practices.
All the tools and pipelines are fully parametrised and optimized for speed. 
The pipeline generates several intermediate files as well as files containing
the discarded variants.

There are 2 pipelines and 2 tools:

**exome_pipeline.py** processes the DNA data and generates a list of unified
filtered and annotated variants and their epitotes. The variant callers
are Mutect2, Strelka, Varscan and SomaticSniper and both indels and snvs are
kept. Annotation is done using Annovar. 
The pipeline uses bwa to align, trimmomatic to quality trim and GATK4
best practices. The filters are applied based on allele frecuency and reads
frecuency. The pipeline will also perform HLA predictions with HLA-LA.

**rnaseq_pipeline.py** processes the RNA-seq data and generates a list of unified
annotated variants (unfiltered) and also a list of FPKM values. The variant
callers used are Varscan and mpileup. Annmotation is done with Annovar.
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

See RUN.txt for running instructions.

Main output files:

**exome_pipeline.py** 
- nonsyn_SQL_insert.txt (nonsynonymous variants filtered and annotated)
- all_other_mutations.txt (synonymous variants filtered and annotated)
- SQL_Epitopes.txt (all the epitotes for the differents exons in each variant including the reference and mutated sequences)

**rnaseq_pipeline.py** 
- nonsyn_SQL_insert.txt (nonsynonymous variants filtered and annotated)
- all_other_mutations.txt (synonymous variants filtered and annotated)
- FPKM_SQL_insert.txt (all the FPKM values of the genes found)

**merge_resuls.py** 
- overlap_final.txt (all the AND and RNA variants/epitotes collapsed and filtered with useful information)
- overlap_final_discarded.txt (all the discarded AND and RNA variants/epitotes collapsed with useful information)

**mhc_predict.py** 
- predictions.csv (all the prediction MHCs)

Contact: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>


