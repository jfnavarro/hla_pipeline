#! /usr/bin/env python
"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from hlapipeline.common import *
from hlapipeline.tools import *
from hlapipeline.epitopes import *
import os
import argparse

def final_variants(input, output, output_other, vcf_cov_dict, sampleID, tumor_type):
    snv = open(input)
    snv_lines = snv.readlines()
    header = snv_lines.pop(0).strip().split('\t')
    nonsyn_file = open(output, 'w')
    all_file = open(output_other, 'w')
    header_str = 'SAMPLE_ID\tCHR\tSTART\tEND\tREF\tALT\tavsnp150\tTUMOR_READ1\tTUMOR_READ2' \
                 '\tFunc.refGene\tGene.refGene\tExonicFunc.refGene\tAAChange.refGene\tFunc.knownGene' \
                 '\tGene.knownGene\tExonicFunc.knownGene\tAAChange.knownGene\tFunc.ensGene\tGene.ensGene' \
                 '\tExonicFunc.ensGene\tAAChange.ensGene\tALL.sites.2015_08\tEUR.sites.2015_08\tAMR.sites.2015_08' \
                 '\tEAS.sites.2015_08\tAFR.sites.2015_08\tTVAF\tTCOV\tPVAL\tCALLERS\tTUMOUR\tVARIANT-KEY\tCOSMIC70\n'
    nonsyn_file.write(header_str)
    all_file.write(header_str)
    for line in snv_lines:
        if line.startswith('#'):
            continue
        columns = line.rstrip('\n').split('\t')
        Chr = columns[header.index('Chr')]
        start = columns[header.index('Start')]
        end = columns[header.index('End')]
        ref = columns[header.index('Ref')]
        alt = columns[header.index('Alt')]
        func_ref_gene = columns[header.index('Func.refGene')]
        gene_ref_gene = columns[header.index('Gene.refGene')]
        ref_gene_detail = columns[header.index('GeneDetail.refGene')]
        exonic_func_ref = columns[header.index('ExonicFunc.refGene')]
        AA_change_refGene = columns[header.index('AAChange.refGene')]
        func_known_gene = columns[header.index('Func.knownGene')]
        gene_known_gene = columns[header.index('Gene.knownGene')]
        known_gene_detail = columns[header.index('GeneDetail.knownGene')]
        exonic_known_ref = columns[header.index('ExonicFunc.knownGene')]
        AA_change_knownGene = columns[header.index('AAChange.knownGene')]
        func_ens_gene = columns[header.index('Func.ensGene')]
        gene_ens_gene = columns[header.index('Gene.ensGene')]
        ens_gene_detail = columns[header.index('GeneDetail.ensGene')]
        exonic_ens_ref = columns[header.index('ExonicFunc.ensGene')]
        AA_change_ensGene = columns[header.index('AAChange.ensGene')]
        avsnp150 = columns[header.index('avsnp150')]
        apr_all = columns[header.index('ALL.sites.2015_08')]
        apr_eur = columns[header.index('EUR.sites.2015_08')]
        apr_amr = columns[header.index('AMR.sites.2015_08')]
        apr_asn = columns[header.index('EAS.sites.2015_08')]
        apr_afr = columns[header.index('AFR.sites.2015_08')]
        cosmic = columns[header.index('cosmic70')]
        ID = Chr + ':' + start
        variant_key = Chr + ':' + start + '-' + end + ' ' + ref + '>' + alt
        if ref_gene_detail != '.':
            AA_change_refGene = ref_gene_detail
        if known_gene_detail != '.':
            AA_change_knownGene = known_gene_detail
        if ens_gene_detail != '.':
            AA_change_ensGene = ens_gene_detail
        try:
            tumor_reads1 = vcf_cov_dict[ID]['read1']
            tumor_reads2 = vcf_cov_dict[ID]['read2']
            tumor_var_freq = vcf_cov_dict[ID]['freq']
            cov = vcf_cov_dict[ID]['cov']
            pval = vcf_cov_dict[ID]['pval']
            callers = vcf_cov_dict[ID]['Note']
            to_write = '\t'.join([str(x) for x in [sampleID, Chr, start, end, ref, alt, avsnp150, tumor_reads1,
                                                   tumor_reads2, func_ref_gene, gene_ref_gene, exonic_func_ref,
                                                   AA_change_refGene, func_known_gene, gene_known_gene, exonic_known_ref,
                                                   AA_change_knownGene, func_ens_gene, gene_ens_gene, exonic_ens_ref,
                                                   AA_change_ensGene, apr_all, apr_eur, apr_amr, apr_asn, apr_afr,
                                                   tumor_var_freq, cov, pval, callers, tumor_type, variant_key, cosmic]])
            if any(x in ' '.join([exonic_func_ref, exonic_known_ref, exonic_ens_ref]) for x in
                   ['nonsynonymous', 'frame', 'stop']):
                nonsyn_file.write(to_write + '\n')
            else:
                all_file.write(to_write + '\n')
        except KeyError:
            print("Missing variant for {}".format(ID))
    nonsyn_file.close()
    all_file.close()
    snv.close()

def germline_pipeline(sample1,
                      sample2,
                      sampleID,
                      genome,
                      genome_star,
                      annotation,
                      tumor_type,
                      SNPSITES,
                      KNOWN_SITE1,
                      KNOWN_SITE2,
                      AA_DICT,
                      cDNA_DICT,
                      THREADS,
                      ANNOVAR_DB,
                      ANNOVAR_VERSION,
                      steps,
                      mode):
    print("Germline pipeline")

    # Create sub-folder to store all results
    os.makedirs('germline', exist_ok=True)
    os.chdir('germline')

    if 'mapping' in steps:
        print('Trimming reads')
        cmd = '{} --cores {} --paired --basename sample {} {}'.format(TRIMGALORE, THREADS, sample1, sample2)
        exec_command(cmd)

        # ALIGNMENT
        print('Starting alignment')
        if mode == 'DNA':
            cmd = '{} -t {} {} sample_val_1.fq.gz sample_val_2.fq.gz | ' \
                  '{} sort --threads {} > trimmed_paired_aligned_sorted.bam'.format(BWA, THREADS, genome, SAMTOOLS, THREADS)
            exec_command(cmd)
        else:
            cmd = '{} --genomeDir {} --readFilesIn sample_val_1.fq.gz sample_val_2.fq.gz --outSAMorder Paired' \
                  ' --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {}' \
                  ' --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c' \
                  ' --runThreadN {}'.format(STAR, genome_star, annotation, THREADS)
            exec_command(cmd)

            cmd = 'mv Aligned.sortedByCoord.out.bam trimmed_paired_aligned_sorted.bam'
            exec_command(cmd)

        # Add headers
        print("Adding headers")
        cmd = '{} AddOrReplaceReadGroups I=trimmed_paired_aligned_sorted.bam O=sample_header.bam SO=coordinate RGID={} RGLB=VHIO' \
              ' RGPL=Illumina RGPU=VHIO RGSM={} Create_Index=true Validation_Stringency=SILENT'.format(PICARD, mode, sampleID)
        exec_command(cmd)

    if 'gatk' in steps:
        # Mark duplicates
        print('Marking duplicates')
        cmd = '{} MarkDuplicatesSpark --input sample_header.bam --output sample_dedup.bam'.format(GATK)
        exec_command(cmd)

        if mode == 'RNA':
            # Split N and cigars
            print('Splitting NCigar Reads')
            cmd = '{} SplitNCigarReads --reference {} --input sample_dedup.bam --output sample_split.bam'.format(GATK, genome)
            exec_command(cmd)

            cmd = 'mv sample_split.bam sample_dedup.bam'
            exec_command(cmd)

        # GATK base re-calibration
        print('Starting re-calibration')
        cmd = '{} BaseRecalibratorSpark --use-original-qualities --input sample_dedup.bam --reference {} --known-sites {} ' \
              '--known-sites {} --known-sites {} --output sample_recal_data.txt'.format(GATK,
                                                                                        genome,
                                                                                        SNPSITES,
                                                                                        KNOWN_SITE1,
                                                                                        KNOWN_SITE2)
        exec_command(cmd)
        cmd = '{} ApplyBQSR --use-original-qualities --add-output-sam-program-record --reference {} --input sample_dedup.bam ' \
              '--bqsr-recal-file sample_recal_data.txt --output sample_final.bam'.format(GATK, genome)
        exec_command(cmd)

    if 'hla' in steps:
        print('Predicting HLAs')
        if mode == 'DNA':
            HLA_predictionDNA('sample_final.bam', sampleID, 'PRG-HLA-LA_output.txt', THREADS)
        else:
            HLA_predictionRNA('sample_final.bam', THREADS)

    if 'variant' in steps:
        # Variant calling (Samtools pile-ups)
        print('Computing pile-ups')
        cmd = '{} mpileup -C 50 -B -q 1 -Q 15 -f {} sample_final.bam > sample.pileup'.format(SAMTOOLS, genome)
        exec_command(cmd)

        # Variant calling VarScan
        print('Variant calling with VarScan2')
        cmd = '{} mpileup2cns sample.pileup --variants 0 --min-coverage 2 --min-reads2 1 --output-vcf 1 ' \
              '--min-var-freq 0.01 --min-avg-qual 15 --p-value 0.99 --strand-filter 1 > varscan.vcf'.format(VARSCAN)
        exec_command(cmd)

        # Variant calling (HaplotypeCaller)
        print('Variant calling with HaplotypeCaller')
        cmd = '{} HaplotypeCaller --reference {} --input sample_final.bam --output haplotype_caller.vcf ' \
              '--dont-use-soft-clipped-bases --standard-min-confidence-threshold-for-calling 20 ' \
              '--dbsnp {}'.format(GATK, genome, SNPSITES)
        exec_command(cmd)

        if mode == 'RNA':
            # Computing gene counts
            print('Computing gene counts with featureCounts')
            cmd = '{} -T {} --primary --ignoreDup -O -C -t exon ' \
                  '-g gene_name -a {} -o gene.counts sample_dedup.bam'.format(FEATURECOUNTS, THREADS, annotation)
            exec_command(cmd)

    if 'filter' in steps:
        # Filtering variants (HaplotypeCaller)
        print("Filtering HaplotypeCaller variants")
        cmd = '{} VariantFiltration --reference {} --variant haplotype_caller.vcf --window 35 --cluster 3 --filter-name "FS" ' \
              '--filter "FS > 30.0" --filter-name "QD" --filter "QD < 2.0" --output haplotype_caller_filtered.vcf'.format(GATK,
                                                                                                                          genome)
        exec_command(cmd)

        # TODO add a filter to VarScan2 variants

        # NOTE replacing IUPAC codes from varscan output (TEMP HACK till the bug is fixed in VarScan)
        # NOTE this will also skip variants whose REF and ALT fields are identical (another bug in VarScan)
        cmd = 'awk \'{if ($1 ~ /#/) {print} else if ($4 != $5) {gsub(/W|K|B|Y|D|H|V|R|S|M/,"N",$4); OFS="\t"; print}}\' ' \
              'varscan.vcf > varscan_filtered.vcf'
        exec_command(cmd)

        # Combine with GATK
        print('Combining variants')
        # CombineVariants is not available in GATK 4 so we need to use the 3.8 version
        cmd = '{} -T CombineVariants -R {} -V:varscan varscan_filtered.vcf ' \
              '-V:HaplotypeCaller haplotype_caller_filtered.vcf -o combined_calls.vcf '\
              '-genotypeMergeOptions UNIQUIFY --num_threads {}'.format(GATK3, genome, THREADS)
        exec_command(cmd)

        # Annotate with Annovar
        annovardb = '{} -buildver {}'.format(os.path.join(ANNOVAR_PATH, ANNOVAR_DB), ANNOVAR_VERSION)
        print('Running annovar (SNV)')
        cmd = '{} combined_calls.vcf {} -thread {} -out snp.sum -vcfinput -remove -protocol {}'.format(
            os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), annovardb, THREADS,  annovar_anno)
        exec_command(cmd)

        # Extract coverage info from vcf file
        print("Extracting coverage from combined variants")
        vcf = open('combined_calls.vcf')
        vcf_cov_dict = {}
        for line in vcf:
            if line.startswith('#CHROM'):
                headers = line.strip().split('\t')
                varscanT = headers.index('Sample1.varscan')
                HaplotypeCallerT = headers.index('{}.HaplotypeCaller'.format(sampleID))
            if not line.startswith('#'):
                columns = line.split('\t')
                chrm = columns[headers.index('#CHROM')]
                pos = columns[headers.index('POS')]
                ref = columns[headers.index('REF')]
                alt = columns[headers.index('ALT')]
                info = columns[headers.index('INFO')]
                form = columns[headers.index('FORMAT')].split(':')
                if len(ref) != len(alt):
                    pos = str(int(pos) + 1)
                ID = chrm + ':' + pos
                pval = -1
                freq = 0
                read1 = 0
                read2 = 0
                cov = 0
                callers = info.strip().split(';')[-1].replace('set=', '')
                caller_count = callers.count('-') + 1 if 'Intersection' not in callers else 2
                vcf_cov_dict[ID] = {}
                vcf_cov_dict[ID]['Note'] = str(caller_count) + ':' + callers
                if 'Intersection' in callers or 'varscan' in callers:
                    t_split = columns[varscanT].split(':')
                    pval = t_split[form.index('PVAL')]
                    freq = t_split[form.index('FREQ')]
                    read1 = t_split[form.index('RDF')]
                    read2 = t_split[form.index('RDR')]
                    cov = t_split[form.index('RD')]
                elif 'HaplotypeCaller' in callers:
                    t_split = columns[HaplotypeCallerT].split(':')
                    if ',' not in alt:
                        AD = form.index('AD')
                        tumor_read1 = int(t_split[AD].split(',')[0])
                        tumor_read2 = int(t_split[AD].split(',')[1])
                        cov = tumor_read1 + tumor_read2
                        if tumor_read2 != 0:
                            freq = str(round((tumor_read2 / cov) * 100, 2)) + '%'
                vcf_cov_dict[ID]['pval'] = pval
                vcf_cov_dict[ID]['freq'] = freq
                vcf_cov_dict[ID]['read1'] = read1
                vcf_cov_dict[ID]['read2'] = read2
                vcf_cov_dict[ID]['cov'] = cov
        vcf.close()

        print('Adding annotated variants to final report')
        final_variants('snp.sum.{}_multianno.txt'.format(ANNOVAR_VERSION),
                       'nonsyn_SQL_insert.txt', 'all_other_mutations.txt',
                       vcf_cov_dict, sampleID, tumor_type)

        # Extract peptides
        extract_peptides('nonsyn_SQL_insert.txt', 'Formatted_epitope_variant.txt', sampleID)

        # Create epitopes
        create_epitopes('Formatted_epitope_variant.txt', 'SQL_Epitopes.txt', AA_DICT, cDNA_DICT)

        if mode == 'RNA':
            reformat_gene_counts('gene.counts', 'GeneCounts_SQL_insert.txt', sampleID, tumor_type)

    print("COMPLETED!")


parser = argparse.ArgumentParser(description='RNA-seq variant calling and HLA prediction pipeline\n'
                                             'Created by Jose Fernandez <jc.fernandes.navarro@gmail.com>)',
                                 prog='germline_pipeline.py',
                                 usage='germline_pipeline.py [options] R1(RNA) R2(RNA)')
parser.add_argument('R1_RNA', help='FASTQ file R1 (RNA)')
parser.add_argument('R2_RNA', help='FASTQ file R2 (RNA)')
parser.add_argument('--genome',
                    help='Path to the reference genome FASTA file', required=True)
parser.add_argument('--genome-star',
                    help='Path to the reference genome STAR index folder (when in RNA mode)', required=False)
parser.add_argument('--genome-ref',
                    help='Path to the reference genome GTF file (when in RNA mode)', required=False)
parser.add_argument('--sample',
                    help='Name of the sample/experiment. Default is sample', default='sample')
parser.add_argument('--tumor',
                    help='Tumor type. Default is Tumor', default='Tumor')
parser.add_argument('--dir',
                    help='Path to the folder where to put output files', required=True)
parser.add_argument('--known1',
                    help='Path to the file with Mill and 1000G gold standards (GATK bundle)', required=True)
parser.add_argument('--known2',
                    help='Path to the file with 1000G phase indels (GATK bundle)', required=True)
parser.add_argument('--snpsites',
                    help='Path to the file with the SNP sites (GATK bundle)', required=True)
parser.add_argument('--dictAA',
                    help='Path to a dictionary of transcript IDs to peptide sequences', required=True)
parser.add_argument('--dictcDNA',
                    help='Path to a dictionary of transcript IDs to DNA sequences', required=True)
parser.add_argument('--annovar-db',
                    help='String indicated which Annovar database to use (default: humandb)',
                    default='humandb', required=False)
parser.add_argument('--annovar-version',
                    help='String indicated which version of the Annovar database to use (default: hg38)',
                    default='hg38', required=False)
parser.add_argument('--threads',
                    help='Number of threads to use in the parallel steps', type=int, default=10, required=False)
parser.add_argument('--steps', nargs='+', default=['mapping', 'gatk', 'hla', 'variant', 'filter'],
                    help='Steps to perform in the pipeline',
                    choices=['mapping', 'gatk', 'hla', 'variant', 'filter', "none"])
parser.add_argument('--mode', default='RNA',
                    help='Mode to use (RNA (default) or DNA)',
                    choices=['DNA', 'RNA'])

# Parse arguments
args = parser.parse_args()
DIR = args.dir
R1_RNA = os.path.abspath(args.R1_RNA)
R2_RNA = os.path.abspath(args.R2_RNA)
sampleID = args.sample
tumor_type = args.tumor
GENOME_REF = os.path.abspath(args.genome)
GENOME_REF_STAR = os.path.abspath(args.genome_star)
GENOME_ANNOTATION = os.path.abspath(args.genome_ref)
THREADS = int(args.threads)
KNOWN_SITE1 = os.path.abspath(args.known1)
KNOWN_SITE2 = os.path.abspath(args.known2)
SNPSITES = os.path.abspath(args.snpsites)
AA_DICT = os.path.abspath(args.dictAA)
cDNA_DICT = os.path.abspath(args.dictcDNA)
STEPS = args.steps
ANNOVAR_DB = args.annovar_db
ANNOVAR_VERSION = args.annovar_version
MODE = args.mode
if MODE == "RNA" and (not GENOME_REF_STAR or not GENOME_ANNOTATION):
    sys.stderr.write("Error, RNA mode but STAR reference or annotation files are missing\n")
    sys.exit(1)

# Move to output dir
os.makedirs(os.path.abspath(DIR), exist_ok=True)
os.chdir(os.path.abspath(DIR))

germline_pipeline(R1_RNA,
                  R2_RNA,
                  sampleID,
                  GENOME_REF,
                  GENOME_REF_STAR,
                  GENOME_ANNOTATION,
                  tumor_type,
                  SNPSITES,
                  KNOWN_SITE1,
                  KNOWN_SITE2,
                  AA_DICT,
                  cDNA_DICT,
                  THREADS,
                  ANNOVAR_DB,
                  ANNOVAR_VERSION,
                  STEPS,
                  MODE)
