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
    print('Re-calibration was performed.')

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
    cmd = VARSCAN + ' mpileup2cns sample.pileup varscan --variants 0 --min-coverage 2 --min-reads2 1'\
                    + ' --min-var-freq .01 --p-value 0.99 > varscan.pileup'
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
        insert_file.write(str(gDNA) + "\t" + str(mrn) + "\t" + str(seq_center) + "\t" + str(sampleID) + "\t" + str(Chr) + "\t" + str(start) + "\t"\
                          + str(end) + "\t" + str(ref) + "\t" + str(alt) + "\t" + str(snp138NonFlagged) + "\t" + str(func_ref_gene)\
                          + "\t" + str(gene_ref_gene) + "\t" + str(exonic_func_ref) + "\t" + str(AA_change_refGene) + "\t" + str(func_known_gene)\
                          + "\t" + str(gene_known_gene) + "\t" + str(exonic_known_ref) + "\t" + str(AA_change_knownGene) + "\t"\
                          + str(func_ens_gene) + "\t" + str(gene_ens_gene) + "\t" + str(exonic_ens_ref) + "\t" + str(AA_change_ensGene)\
                          + "\t" + str(apr_all) + "\t" + str(apr_eur) + "\t" + str(apr_amr) + "\t" + str(apr_asn) + "\t" + str(apr_afr)\
                          + "\t" + str(sample_gDNA) + "\t" + str(gDNA) + "\t" + str(sample_center) + "\t" + str(tumor_type)\
                          + "\t" + str(sample_note) + '_' + str(date) + "\t" + str(source_of_DNA) + "\t" + str(variant_key) + "\n")
    insert_file.close()
    snv.close()

    # This will format the  coverage information into a format that will be joined with the variant information
    print('Filtering Pile-ups')
    pileup = open('varscan.pileup')
    insert_file = open('SQL_coverage_input.txt', 'w')
    date = datetime.datetime.now().replace(microsecond=0)
    first_line = True
    for line in pileup:
        if first_line:
            first_line = False
            continue
        columns = line.rstrip('\n').split('\t')
        Chr = columns[0]
        start = columns[1]
        ref = columns[2]
        alt = columns[3]
        columns2 = columns[4].rstrip('\n').split(':')
        cons = columns2[0]
        cov = columns2[1]
        read1 = columns2[2]
        read2 = columns2[3]
        freq = columns2[4]
        p_val = columns2[5]
        columns3 = columns[5].rstrip('\n').split(':')
        length = len(columns3)
        r1_plus = columns3[1]
        r1_minus = columns3[2]
        r2_plus = columns3[3]
        r2_minus = columns3[4]
        p_val2 = columns3[5]
        mrn = MRN
        seq_center = SEQ_CENTER
        sampleID = sampleID
        sample_gDNA = sampleID + ' chr' + Chr + ':' + start
        gDNA = 'chr' + Chr + ':' + start
        join_key = gDNA
        if re.search(r'-', alt):
            join_key = 'chr' + Chr + ':'+  str(int(start) + 1)
        tumor_type = tumor_type
        source_of_DNA = SOURCE
        sample_note = SAMPLE_NOTE
        sample_center = SEQ_CENTER
        insert_file.write(str(join_key) + "\t" + str(mrn) + "\t" + str(sampleID) + "\t" + str(seq_center) + "\t" + str(Chr)\
                          + "\t" + str(start) + "\t" + str(ref) + "\t" + str(alt) + "\t" + str(cons) + "\t" + str(cov)\
                          + "\t" + str(read1) + "\t" + str(read2) + "\t" + str(freq) + "\t" + str(p_val) + "\t" + str(r1_plus)\
                          + "\t" + str(r1_minus) + "\t" + str(r2_plus) + "\t" + str(r2_minus) + "\t" + str(p_val2)\
                          + "\t" + str(sample_center) + "\t" + str(sample_gDNA) + "\t" + str(gDNA) + "\t" + str(tumor_type)\
                          + "\t" +str(sample_note) + '_' + str(date) + "\t" + str(source_of_DNA) + "\n")
    pileup.close()
    insert_file.close()

    print('Sorting and joining files')
    cmd = 'sort -k1b,1 SQL_coverage_input.txt > join_coverage_sort.txt'
    exec_command(cmd)
    cmd = 'sort -k1b,1 SQL_variant_input.txt > join_variants_sort.txt'
    exec_command(cmd)
    cmd = 'join -t$\'\\t\' -a1 -e"-" join_variants_sort.txt join_coverage_sort.txt > joined_coverage_variants.txt'
    exec_command(cmd)

    # This will extract just the nonsynonymous mutations:
    print('Creating nonsynonymous file for insert into database')
    joined_variants = open('joined_coverage_variants.txt')
    nonsyn_file = open('nonsyn_SQL_insert.txt', 'w')
    for line in joined_variants:
        columns = line.rstrip('\n').split('\t')
        if (len(columns) < 51):
            continue
        mrn = columns[1]
        seq_center = columns[2]
        sampleID = columns[3]
        Chr = columns[4]
        Start = columns[5]
        End = columns[6]
        Ref = columns[7]
        Alt = columns[8]
        snp138JJG = columns[9]
        Func_refGene = columns[10]
        Gene_refGene = columns[11]
        ExonicFunc_refGene = columns[12]
        AAChange_refGene = columns[13]
        Func_knownGene = columns[14]
        Gene_knownGene = columns[15]
        ExonicFunc_knownGene = columns[16]
        AAChange_knownGene = columns[17]
        Func_ensGene = columns[18]
        Gene_ensGene = columns[19]
        ExonicFunc_ensGene = columns[20]
        AAChange_ensGene = columns[21]
        apr_all = columns[22]
        apr_eur = columns[23]
        apr_amr = columns[24]
        apr_asn = columns[25]
        apr_afr = columns[26]
        sample_gDNA = columns[27]
        gDNA = columns[28]
        sample_center = columns[29]
        tumor_type = columns[30]
        Note = columns[31]
        source_of_RNA_used_for_sequencing = columns[32]
        tumor_reads1 = columns[43]
        tumor_reads2 = columns[44]
        tumor_var_freq = columns[45].replace('%','')
        read1_plus = columns[47]
        read1_minus = columns[48]
        read2_plus = columns[49]
        read2_minus = columns[50]
        variant_key = columns[33]
        if (re.search(r'nonsynonymous', ExonicFunc_refGene)or re.search(r'frame', ExonicFunc_refGene)or re.search(r'stop', ExonicFunc_refGene)\
                or re.search(r'nonsynonymous', ExonicFunc_knownGene)or re.search(r'frame', ExonicFunc_knownGene) or re.search(r'stop', ExonicFunc_knownGene)\
                or re.search(r'nonsynonymous', ExonicFunc_ensGene) or re.search(r'frame', ExonicFunc_ensGene)or re.search(r'stop', ExonicFunc_ensGene)):
            nonsyn_file.write(str(mrn) + "\t" + str(seq_center) + "\t" + str(sampleID) + "\t" + str(Chr) + "\t" + str(Start) + "\t" + str(End)\
                              + "\t" + str(Ref) + "\t" + str(Alt) + "\t" + "\t" + "\t" + str(snp138JJG) + "\t" + str(tumor_reads1)\
                              + "\t" + str(tumor_reads2) + "\t" + str(tumor_var_freq) + "\t" + str(apr_all) + "\t" + str(Func_refGene)\
                              + "\t" + str(Gene_refGene) + "\t" + str(ExonicFunc_refGene) + "\t" + str(AAChange_refGene) + "\t" + str(Func_knownGene)\
                              + "\t" + str(Gene_knownGene) + "\t" + str(ExonicFunc_knownGene) + "\t" + str(AAChange_knownGene) + "\t" + str(Func_ensGene)\
                              + "\t" + str(Gene_ensGene) + "\t" + str(ExonicFunc_ensGene) + "\t" + str(AAChange_ensGene) + "\t" + str(apr_eur)\
                              + "\t" + str(apr_amr) + "\t" + str(apr_asn) + "\t" + str(apr_afr) + "\t" + "\t" + "\t" + "\t" + str(read1_plus)\
                              + "\t" + str(read1_minus) + "\t" + str(read2_plus) + "\t" + str(read2_minus) + "\t" + "\t" + "\t" + str(sample_gDNA)\
                              + "\t" +str(gDNA) + "\t" + str(sample_center) + "\t" + str(tumor_type) + "\t" + str(sample_note) + '_' + str(date)\
                              + "\t" + str(source_of_RNA_used_for_sequencing) + "\t" + str(variant_key) + "\n")
    joined_variants.close()
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