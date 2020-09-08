# DNA and RNA-Seq variant calling pipelines with HLA and MHC predictions
This a private repository that contains a simplified, more generic, up-to-date
and optimized version of a pipeline provided by Jared Gartner <jared.gartner@nhi.gov>
The pipeline can process DNA/RNA in somatic (tumor and normal) and germline (tumor) modes
and generate a list of variants, HLAs and neo-antigens with affinity scores. 

The pipelines use the GATK4 best-practices.
All the tools and pipelines are fully parametrised and optimized for speed. 
The pipeline generates several intermediate files as well as files containing
the discarded variants.

There are 2 pipelines and 2 tools:

**somatic_pipeline.py** processes DNA/RNA data and generates a list of unified
filtered and annotated somatic variants, their epitopes and also a list of gene counts values (when in RNA mode). 
The variant callers are Mutect2, Strelka, Varscan and SomaticSniper and both indels and SNPs are
reported. Annotation is performed using Annovar. 
The pipeline uses trim-galore to trim, STAR (RNA) or bwa-men (DNA) to align and GATK4 best practices. 
The pipeline also performs HLA predictions with HLA-LA (DNA) and arcasHLA (RNA).
The gene counts values are computed with featureCounts (when in RNA mode).
When in RNA mode the STAR index and annotation reference (GTF) must be provided as arguments.

**germline_pipeline.py** processes DNA/RNA data and generates a list of unified
annotated germline variants (weak filtered), their epitopes and also a list of gene counts values (when in RNA mode). 
The variant callers used are Varscan and HaplotypeCaller. Annotation is performed with Annovar.
The pipeline uses trim-galore to trim, STAR (RNA) or bwa-men (DNA) to align and GATK4 best practices. 
The pipeline also performs HLA predictions with HLA-LA (DNA) and arcasHLA (RNA).
The gene counts values are computed with featureCounts (when in RNA mode).
When in RNA mode the STAR index and annotation reference (GTF) must be provided as arguments.

**merge_resuls.py** combines results from 1 or several runs of the somatic and germline
pipelines in order to generate an unified file with useful information where
variants are filtered by certain criteria. 

**mhc_predict.py** can take the file generated with merge_results.py and the HLA files
generated in the germline or somatic pipelines and then output a list of predicted neo-antigens.
Variants are filtered by certain criterias and only the most common alleles for each HLA class 1
are used. 

Each tool/pipeline uses a command line interface with parameters which
can be shown and described with --help.

See INSTALL.txt for installation instructions. 

See RUN.txt for running instructions.

Main output files:

**somatic_pipeline.py** 
- nonsyn_SQL_insert.txt (nonsynonymous variants filtered and annotated)
- all_other_mutations.txt (synonymous variants filtered and annotated)
- Tumor_GeneCounts_SQL_insert.txt (all the counts of the genes found in the tumor sample if RNA mode is used)
- Normal_GeneCounts_SQL_insert.txt (all the counts of the genes found in the normal sample if RNA mode is used)
- SQL_Epitopes.txt (all the epitopes for the different isoforms in each variant including the reference and mutated sequences)
- HLA predictions (Tumor and normal)

**germline_pipeline.py** 
- nonsyn_SQL_insert.txt (nonsynonymous variants filtered and annotated)
- all_other_mutations.txt (synonymous variants filtered and annotated)
- GeneCounts_SQL_insert.txt (all the counts of the genes found if RNA mode is used)
- SQL_Epitopes.txt (all the epitopes for the different isoforms in each variant including the reference and mutated sequences)
- HLA predictions

**merge_resuls.py** 
- overlap_final.txt (all the AND and RNA variants/epitopes collapsed and filtered with useful information)
- overlap_final_discarded.txt (all the discarded AND and RNA variants/epitopes collapsed with useful information)

**mhc_predict.py** 
- predictions_mut.csv (all the mutated peptides predictions)
- predictions_wt.csv (all the WT peptides predictions)

Contact: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>


