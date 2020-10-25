# DNA and RNA NGS variant calling pipelines with HLA and MHC predictions
The pipeline can process DNA/RNA data in somatic (tumor and normal) and germline (tumor) modes
and generate a list of variants, HLAs and neo-antigens with affinity scores. 

The somatic and germline pipelines make use of the latest GATK4 best-practices.
All the tools and pipelines are fully parametrised and optimized for speed. 

There are 2 pipelines and 2 tools:

**somatic_pipeline.py** processes DNA/RNA data and generates a list of unified
filtered and annotated somatic variants and also a list of gene counts values (when in RNA mode). 
The variant callers are Mutect2, Strelka, Varscan and SomaticSniper and both indels and SNPs are
reported. Annotation is performed using Annovar. 
The pipeline uses trim-galore to trim, STAR (RNA) or bwa-men (DNA) to align and follows GATK4 best practices. 
The pipeline also performs HLA predictions with HLA-LA (DNA) and arcasHLA (when in RNA).
The gene counts values are computed with featureCounts (when in RNA mode).
When in RNA mode the STAR index and annotation reference (GTF) must be provided as arguments.

**germline_pipeline.py** processes DNA/RNA data and generates a list of unified
annotated germline variants (weak filtered) and also a list of gene counts values (when in RNA mode). 
The variant callers used are Varscan and HaplotypeCaller. Annotation is performed with Annovar.
The pipeline uses trim-galore to trim, STAR (RNA) or bwa-men (DNA) to align and follows GATK4 best practices. 
The pipeline also performs HLA predictions with HLA-LA (DNA) and arcasHLA (RNA).
The gene counts values are computed with featureCounts (when in RNA mode).
When in RNA mode the STAR index and annotation reference (GTF) must be provided as arguments.

**merge_resuls.py** combines results from 1 or several runs of the somatic and germline
pipelines in order to generate an unified table with useful information where
variants are filtered by certain criteria and epitopes are created for each of the variants
effects. 

**mhc_predict.py** can take the file generated with merge_results.py and the HLA files
generated in the germline and/or somatic pipelines and then output a list of predicted neo-antigens.
Variants are filtered by certain criterias and only the most common alleles for each HLA class 1
are used. 

Each tool/pipeline uses a command line interface with parameters which
can be shown and described with --help.

See INSTALL.txt for installation instructions. 

See RUN.txt for running instructions.

## cDNA and Peptides dictionaries
merge_results.py requires two dictionaries, one mapping transcript ids to DNA sequences and another
one mapping transcript ids to peptide sequences. The format is the following for both files:

TRANSCRIPT_ID:SEQUENCE 

To build these dictionaries you can use as reference the Jupyter Notebooks located in dictionaries

## Output (important files)

**somatic_pipeline.py** 
- annotated.hg38_multianno.vcf (annotated and combined somatic variants)
- gene.counts (all the counts of the genes found in the tumor sample if RNA mode is used)
- HLA predictions (Tumor and normal)

**germline_pipeline.py** 
- annotated.hg38_multianno.vcf (annotated and combined germline variants)
- gene.counts (all the counts of the genes found in the tumor sample if RNA mode is used)
- HLA predictions

**merge_resuls.py** 
- overlap_final.txt (all the somatic and germline variants collapsed and filtered with useful information and epitopes)
- overlap_final_unique_germline.txt (all the unique germline variants collapsed and filtered with useful information and epitopes)
- overlap_final_discarded.txt (all the discarded somatic and germline variants collapsed with useful information and epitopes)

**mhc_predict.py** 
- predictions_mut.csv (all the mutated peptides predictions)
- predictions_wt.csv (all the WT peptides predictions)

## Contact
Contact: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>


