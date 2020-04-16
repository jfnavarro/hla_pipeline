from common import *
import re
import datetime
import os

def RNA_seq_pipeline(sample1, sample2, sampleID, genome, genome_star, annotation, tumor_type, SNPSITES, KNOWN_SITE1, KNOWN_SITE2, THREADS):
    print("RNA-seq pipeline")

    # Create sub-folder to store all results
    os.makedirs('rna', exist_ok=True)
    os.chdir('rna')

    print('Trimmimg reads')
    cmd = '{} --paired --basename sample {} {}'.format(TRIMGALORE, sample1, sample2)
    exec_command(cmd)
    print('Trimming completed.')

    print('Aligining with STAR')
    cmd = '{} --genomeDir {} --readFilesIn sample_val_1.fq.gz sample_val_2.fq.gz --outSAMmultNmax 1 --outSAMorder Paired'\
          ' --outSAMprimaryFlag OneBestScore --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {} --outFilterIntronMotifs'\
          ' RemoveNoncanonical --outFilterType Normal --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c'\
          ' --runThreadN {} --outFilterMultimapNmax 20'.format(STAR, genome_star, annotation, max(int(THREADS / 2), 1))
    exec_command(cmd)
    print('Aligment completed.')

    # Add headers
    print("Adding headers")
    cmd = '{} AddOrReplaceReadGroups I=Aligned.sortedByCoord.out.bam O=sample_header.bam SO=coordinate RGID=SB_RNA-seq RGLB={}'\
          ' RGPL=Illumina RGPU=PU{} RGSM=rna-seq_{} Create_Index=true Validation_Stringency=SILENT'.format(PICARD,
                                                                                                           SEQ_CENTER,
                                                                                                           SEQ_CENTER,
                                                                                                           sampleID)
    exec_command(cmd)
    print('Headers added.')

    # Mark duplicates
    print('Marking duplicates')
    cmd = GATK + ' MarkDuplicatesSpark -I=sample_header.bam -O=sample_dedup.bam -M=dedup_sample.txt'
    exec_command(cmd)
    print('Duplicates marked')

    # Split N and cigars
    cmd = '{} SplitNCigarReads -R {} -I sample_dedup.bam -O sample_split.bam'.format(GATK, genome)
    exec_command(cmd)

    # GATK base re-calibration
    print('Starting re-calibration')
    cmd = '{} BaseRecalibratorSpark -I sample_dedup.bam -R {} --known-sites {} --known-sites {}'\
          ' --known-sites {} -O sample_recal_data.txt'.format(GATK, genome, SNPSITES, KNOWN_SITE1, KNOWN_SITE2)
    exec_command(cmd)
    cmd = '{} ApplyBQSR -R {} -I sample_dedup.bam --bqsr-recal-file sample_recal_data.txt -O sample_final.bam'.format(GATK, genome)
    exec_command(cmd)
    print('Re-calibration was performed on the tumor and normal samples.')

    # Variant calling (Samtools pile-ups)
    print('Computing pile-ups')
    cmd = '{} mpileup -C50 -B -q 1 -Q 15 -f {} sample_final.bam > sample.pileup'.format(SAMTOOLS, genome)
    exec_command(cmd)
    print('Pile-ups computed.')

    # Variant calling VarScan
    print('Variant calling with varscan')
    cmd = VARSCAN + ' mpileup2cns sample.pileup varscan --variants 0 --min-coverage 2 --min-reads2 1 --output-vcf 1'\
                    + ' --min-var-freq .01 --p-value 0.99 > varscan.vcf'
    exec_command(cmd)
    print('Variant calling completed.')

    # Run annovar to annotate variants
    print('Running annovar')
    # TODO ensure GHRC37 works with annovar (hg19)
    cmd = '{} -format vcf4 varscan.vcf --comment --includeinfo -outfile snp.av'.format(
        os.path.join(ANNOVAR_PATH, 'convert2annovar.pl'))
    exec_command(cmd)
    cmd = '{} snp.av {} -thread {} -out snp.sum -remove -protocol {}'.format(
        os.path.join(ANNOVAR_PATH, 'table_annovar.pl'), annovar_db, THREADS, annovar_anno)
    exec_command(cmd)

    # Running Cufflinks (output is genes.fpkm_tracking)
    print('Running Cufflinks')
    cmd = '{} -p {} -G {} --library-type fr-firststrand sample_dedup.bam'.format(CUFFLINKS, THREADS, annotation)
    exec_command(cmd)
    print('Cufflinks completed')

    print('Filtering Varscan')
    snv = open('snp.sum.hg19_multianno.txt')
    insert_file = open('SQL_variant_input.txt', 'w')
    nonsyn_file = open('nonsyn_SQL_insert.txt', 'w')
    date = datetime.datetime.now().replace(microsecond=0)
    for line in snv:
        if line.startswith('#') or line.startswith("Chr"):
            continue
        columns = line.rstrip('\n').split('\t')
        Chr = columns[0]
        start = columns[1]
        end = columns[2]
        ref = columns[3]
        alt = columns[4]
        func_ref_gene = columns[5]
        gene_ref_gene = columns[6]
        ref_gene_detail = columns[7]
        exonic_func_ref = columns[8]
        AA_change_refGene = columns[9]
        func_known_gene = columns[10]
        gene_known_gene = columns[11]
        known_gene_detail = columns[12]
        exonic_known_ref = columns[13]
        AA_change_knownGene = columns[14]
        func_ens_gene = columns[15]
        gene_ens_gene = columns[16]
        ens_gene_detail = columns[17]
        exonic_ens_ref = columns[18]
        AA_change_ensGene = columns[19]
        snp138NonFlagged = columns[20]
        apr_all = columns[21]
        apr_eur = columns[22]
        apr_amr = columns[23]
        apr_asn = columns[24]
        apr_afr = columns[25]
        mrn = MRN
        seq_center = SEQ_CENTER
        sampleID = sampleID
        sample_gDNA = sampleID + ' chr' + Chr + ':' + start
        gDNA = 'chr' + Chr + ':' + start
        tumor_type = tumor_type
        source_of_DNA = SOURCE
        sample_note = SAMPLE_NOTE
        sample_center = SEQ_CENTER
        variant_key = Chr + ':' + start + '-' + end + ' ' + ref + '>' + alt
        if ref_gene_detail != 'NA':
            AA_change_refGene = ref_gene_detail
        if known_gene_detail != 'NA':
            AA_change_knownGene = known_gene_detail
        if ens_gene_detail != 'NA':
            AA_change_ensGgene = ens_gene_detail
        txt = str(gDNA) + "\t" + str(mrn) + "\t" + str(seq_center) + "\t" + str(sampleID) + "\t" + str(Chr) + "\t" + str(start) + "\t"\
                          + str(end) + "\t" + str(ref) + "\t" + str(alt) + "\t" + str(snp138NonFlagged) + "\t" + str(func_ref_gene)\
                          + "\t" + str(gene_ref_gene) + "\t" + str(exonic_func_ref) + "\t" + str(AA_change_refGene) + "\t" + str(func_known_gene)\
                          + "\t" + str(gene_known_gene) + "\t" + str(exonic_known_ref) + "\t" + str(AA_change_knownGene) + "\t"\
                          + str(func_ens_gene) + "\t" + str(gene_ens_gene) + "\t" + str(exonic_ens_ref) + "\t" + str(AA_change_ensGene)\
                          + "\t" + str(apr_all) + "\t" + str(apr_eur) + "\t" + str(apr_amr) + "\t" + str(apr_asn) + "\t" + str(apr_afr)\
                          + "\t" + str(sample_gDNA) + "\t" + str(gDNA) + "\t" + str(sample_center) + "\t" + str(tumor_type)\
                          + "\t" + str(sample_note) + '_' + str(date) + "\t" + str(source_of_DNA) + "\t" + str(variant_key) + "\n"
        insert_file.write(txt)
        if (re.search(r'nonsynonymous', ExonicFunc_refGene) or re.search(r'frame', ExonicFunc_refGene)or re.search(r'stop', ExonicFunc_refGene)\
                or re.search(r'nonsynonymous', ExonicFunc_knownGene)or re.search(r'frame', ExonicFunc_knownGene) or re.search(r'stop', ExonicFunc_knownGene)\
                or re.search(r'nonsynonymous', ExonicFunc_ensGene) or re.search(r'frame', ExonicFunc_ensGene)or re.search(r'stop', ExonicFunc_ensGene)):
            nonsyn_file.write(txt)
    insert_file.close()
    snv.close()
    nonsyn_file.close()

    # Reformat FPKM file
    print('Creating FPKM insert file.')
    fpkm = open('genes.fpkm_tracking')
    FPKM_ins = open('FPKM_SQL_insert.txt', 'w')
    first_line = True
    for line in fpkm:
        if first_line:
            first_line = False
            continue
        FPKM_ins.write(str(MRN) + '\t' + str(SEQ_CENTER) + '\t' + str(sampleID) + '\t' + str(SOURCE) + '\t' + str(tumor_type)\
                       + '\t' + str(SAMPLE_NOTE) + '_' + date + '\t' + str(sampleID) + ' ' + str(SEQ_CENTER) + '\t' + line)
    FPKM_ins.close()
    fpkm.close()