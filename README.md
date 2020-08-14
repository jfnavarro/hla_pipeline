# DNA (WES) and RNA-Seq variant calling pipelines with HLA and MHC predictions
This a private repository that contains a simplified, more generic, up-to-date
and optimized version of a pipeline provided by Jared Gartner <jared.gartner@nhi.gov>
The pipeline can process WES (tumor and normal) and RNA-seq data (tumor)
and generate a list of variants, HLAs and neo-antigens with affinity scores. 

The DNA and RNA pipelines use the GATK4 best-practices.
All the tools and pipelines are fully parametrised and optimized for speed. 
The pipeline generates several intermediate files as well as files containing
the discarded variants.

There are 2 pipelines and 2 tools:

**exome_pipeline.py** processes the WES data and generates a list of unified
filtered and annotated somatic variants and their epitopes. The variant callers
are Mutect2, Strelka, Varscan and SomaticSniper and both indels and SNPs are
reported. Annotation is performed using Annovar. 
The pipeline uses bwa-men to align, trimmomatic to quality trim and GATK4
best practices. The filters are applied based on allele frequency and reads
frequency. The pipeline will also perform HLA predictions with HLA-LA.
This pipeline can be used with RNA-seq data as tumor input using the --mode
flag. If this is the case the STAR index and annotation reference (GTF) must be provided
as arguments.

**rnaseq_pipeline.py** processes the RNA-seq data and generates a list of unified
annotated germline variants (weak filtered), epitopes and also a list of gene counts values. 
The variant callers used are Varscan and HaplotypeCaller. Annmotation is performed with Annovar.
The pipeline uses trim-galore to trim, STAR to align and GATK4 best practices. 
The gene counts values are computed with featureCounts.
The pipeline will also perform HLA predictions with arcasHLA.

**merge_resuls.py** combines results from 1 or several runs of the exome and/or rna-seq
pipelines in order to generate an unified file with useful information where
variants are filtered by certain criterias. 

**mhc_predict.py** can take the file generated with merge_results.py and the HLA files
generated in the exome and/or rna-seq pipelines and then output a list of predicted neo-antigens.
Variants are filtered by certain criterias and only the most common alleles for each HLA class 1
are used. 

Each tool/pipeline uses a command line interface with parameters which
can be shown and described with --help.

See INSTALL.txt for installation instructions. 

See RUN.txt for running instructions.

Main output files:

**exome_pipeline.py** 
- nonsyn_SQL_insert.txt (nonsynonymous variants filtered and annotated)
- all_other_mutations.txt (synonymous variants filtered and annotated)
- SQL_Epitopes.txt (all the epitopes for the different isoforms in each variant including the reference and mutated sequences)

**rnaseq_pipeline.py** 
- nonsyn_SQL_insert.txt (nonsynonymous variants filtered and annotated)
- all_other_mutations.txt (synonymous variants filtered and annotated)
- GeneCounts_SQL_insert.txt (all the counts of the genes found)
- SQL_Epitopes.txt (all the epitopes for the different isoforms in each variant including the reference and mutated sequences)

**merge_resuls.py** 
- overlap_final.txt (all the AND and RNA variants/epitopes collapsed and filtered with useful information)
- overlap_final_discarded.txt (all the discarded AND and RNA variants/epitopes collapsed with useful information)

**mhc_predict.py** 
- predictions_mut.csv (all the mutated peptides predictions)
- predictions_wt.csv (all the WT peptides predictions)

Contact: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>


