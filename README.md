# DNA and RNA variant calling pipelines with HLA typing and Neoantigen predictions
This toolset can process DNA (tumor and normal) and RNA (tumor) sequencing data 
and generate a list of somatic variants, HLAs and neoantigens with affinity scores. 

<p align="center">
<img src="diagram.png" height="800">
</p>


The DNA and RNA pipelines make use of the latest GATK4 best-practices.
All the tools and pipelines are fully parametrised and optimized for speed. 

There are 2 pipelines and 2 tools:

**dna_pipeline.py** processes DNA data and generates a list of unified
filtered and annotated somatic variants. 
The variant callers are Mutect2, Strelka2, Varscan and SomaticSniper and both indels and SNPs are
reported. Annotation is performed using Annovar. 
The pipeline uses trim-galore to trim, bwa-men to align and follows GATK4 best practices. 
The pipeline also performs HLA predictions with OptiType (tumor and normal).
QC is performed with FastQC and BamQC.

**rna_pipeline.py** processes RNA data and generates a list of unified
annotated somatic variants (weak filtered) and also a list of gene counts values. 
The variant callers used are Varscan and HaplotypeCaller. Annotation is performed with Annovar.
The pipeline uses trim-galore to trim, STAR to align and follows GATK4 best practices. 
The pipeline also performs HLA predictions with OptiType.
The gene counts values are computed with featureCounts.
QC is performed with FastQC and BamQC.

**merge_resuls.py** combines results from 1 or several runs of the DNA and RNA
pipelines in order to generate an unified table with useful information where
variants are filtered by certain criteria (defined by the user) and epitopes 
are created for each of the variants somatic effects. The user can define
the values of the filters for both dna and rna variants. 

**mhc_predict.py** can take the file generated with merge_results.py and the HLA files
generated in the DNA and/or RNA pipelines and then generate a list of predicted neo-antigens.
Variants are filtered by certain criteria and only the most common alleles for each HLA class 1
are used. 

Each tool/pipeline uses a command line interface with parameters which
can be shown and described with --help.


## Requirements
We strongly recommend to use Anaconda or Miniconda, otherwise you may need to create aliases
for some tools as expected in the file hlapipeline/tools.py. 

See environment.yml for a list of required packages. 

## Install
See INSTALL.txt for installation instructions. 

See REFERENCES.txt for instructions to download the references needed
to run the pipelines. 

## How to run
See RUN.txt for a running example.

It is recommended to use a Linux machine with at least 40 cores, 64GB of RAM
and 500GB of disk space. 

## Output (important files)

**dna_pipeline.py** 
- annotated.hgXX_multianno.vcf (annotated and combined somatic variants)
- HLA predictions DNA (Tumor_hla_genotype.tsv and Normal_hla_genotype.tsv)

Other files:
  - combined_calls.vcf
  - tumor_final.bam
  - normal_final.bam
  - fastqc files
  - cutadapt stats
  - vcf stats
  - bamQC_Normal folder
  - bamQC_Tumor folder
  
**rna_pipeline.py** 
- annotated.hgXX_multianno.vcf (annotated and combined germline variants)
- gene.counts (gene counts from featureCounts)
- HLA predictions (hla_genotype.tsv)

Other files:
  - combined_calls.vcf
  - sample_final.bam
  - fastqc files
  - cutadapt stats
  - STAR log
  - vcf stats
  - bamQC folder
  - bamQCRNA folder
  
**merge_results.py** 
- overlap_final.txt (all the DNA and RNA variants collapsed and filtered with useful information and epitopes)
- overlap_final_unique_rna.txt (all the RNA variants collapsed and filtered with useful information and epitopes)
- overlap_final_discarded.txt (all the discarded DNA variants collapsed with useful information and epitopes)
- overlap_final_discarded_rna.txt (all the discarded RNA variants collapsed with useful information and epitopes)
- gene.counts.final (same file as gene.counts with three additional columns for TPM, RPKM and Percentile)

**mhc_predict.py** 
- predictions_mut.csv (all the mutated peptides predictions)
- predictions_wt.csv (all the WT peptides predictions)

Other files:
  - protein_sequences_mu.fasta
  - protein_sequences_wt.fasta
  
## Contact
Contact: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>


