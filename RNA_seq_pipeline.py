from common import *
import re
import datetime
import os

def RNA_seq_pipeline(sample1,
                     sample2,
                     sampleID,
                     genome,
                     genome_star,
                     annotation,
                     tumor_type,
                     SNPSITES,
                     KNOWN_SITE1,
                     KNOWN_SITE2,
                     THREADS,
                     steps):
    print("RNA-seq pipeline")

    # Create sub-folder to store all results
    os.makedirs('rna', exist_ok=True)
    os.chdir('rna')

    if 'mapping' in steps:
        print('Trimmimg reads')
        cmd = '{} --paired --basename sample {} {}'.format(TRIMGALORE, sample1, sample2)
        exec_command(cmd)

        print('Aligining with STAR')
        cmd = '{} --genomeDir {} --readFilesIn sample_val_1.fq.gz sample_val_2.fq.gz --outSAMmultNmax 1 --outSAMorder Paired'\
              ' --outSAMprimaryFlag OneBestScore --twopassMode Basic --outSAMunmapped None --sjdbGTFfile {} --outFilterIntronMotifs'\
              ' RemoveNoncanonical --outFilterType Normal --outSAMtype BAM SortedByCoordinate --readFilesCommand gunzip -c'\
              ' --runThreadN {} --outFilterMultimapNmax 20'.format(STAR, genome_star, annotation, THREADS)
        exec_command(cmd)

        # Add headers
        print("Adding headers")
        cmd = '{} AddOrReplaceReadGroups I=Aligned.sortedByCoord.out.bam O=sample_header.bam SO=coordinate RGID=SB_RNA-seq RGLB={}'\
              ' RGPL=Illumina RGPU=PU{} RGSM=rna-seq_{} Create_Index=true Validation_Stringency=SILENT'.format(PICARD,
                                                                                                               SEQ_CENTER,
                                                                                                               SEQ_CENTER,
                                                                                                               sampleID)
        exec_command(cmd)

    if 'hla' in steps:
        print('Predicting HLA')
        HLA_predictionRNA('sample_header.bam', THREADS)

    if 'gatk' in steps:
        # Mark duplicates
        print('Marking duplicates')
        cmd = GATK + ' MarkDuplicatesSpark -I=sample_header.bam -O=sample_dedup.bam -M=dedup_sample.txt'
        exec_command(cmd)

        # Split N and cigars
        print('Splitting NCirgar Reads')
        cmd = '{} SplitNCigarReads -R {} -I sample_dedup.bam -O sample_split.bam'.format(GATK, genome)
        exec_command(cmd)

        # GATK base re-calibration
        print('Starting re-calibration')
        # NOTE BaseRecalibratorSpark needs the system to allow for many open files (ulimit -n)
        cmd = '{} BaseRecalibrator -I sample_dedup.bam -R {} --known-sites {} --known-sites {}'\
              ' --known-sites {} -O sample_recal_data.txt'.format(GATK, genome, SNPSITES, KNOWN_SITE1, KNOWN_SITE2)
        exec_command(cmd)
        cmd = '{} ApplyBQSR -R {} -I sample_dedup.bam --bqsr-recal-file sample_recal_data.txt -O sample_final.bam'.format(GATK, genome)
        exec_command(cmd)

    if 'variant' in steps:
        # Variant calling (Samtools pile-ups)
        print('Computing pile-ups')
        cmd = '{} mpileup -C50 -B -q 1 -Q 15 -f {} sample_final.bam > sample.pileup'.format(SAMTOOLS, genome)
        exec_command(cmd)

        # Variant calling VarScan
        print('Variant calling with varscan')
        cmd = VARSCAN + ' mpileup2cns sample.pileup varscan --variants 0 --min-coverage 2 --min-reads2 1 --output-vcf 1'\
                        + ' --min-var-freq .01 --p-value 0.99 > varscan.vcf'
        exec_command(cmd)
        cmd = VARSCAN + ' mpileup2cns sample.pileup varscan --variants 0 --min-coverage 2 --min-reads2 1'\
                        + ' --min-var-freq .01 --p-value 0.99 > varscan.pileup'
        exec_command(cmd)

        # TODO apply a filter to the variants

        # Run annovar to annotate variants
        print('Running annovar')
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

    if 'filter' in steps:
        print('Formatting Varscan variants')
        snv = open('snp.sum.hg19_multianno.txt')
        snv_lines = snv.readlines()
        header = snv_lines.pop(0).strip().split('\t')
        insert_file = open('SQL_variant_input.txt', 'w')
        date = datetime.datetime.now().replace(microsecond=0)
        # TODO remove unnecessary fields
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
            gene_known_gene = columns[header.index('Gene.known')]
            known_gene_detail = columns[header.index('GeneDetail.knownGene')]
            exonic_known_ref = columns[header.index('ExonicFunc.knownGene')]
            AA_change_knownGene = columns[header.index('AAChange.knownGene')]
            func_ens_gene = columns[header.index('Func.ensGene')]
            gene_ens_gene = columns[header.index('Gene.ensGene')]
            ens_gene_detail = columns[header.index('GeneDetail.ensGene')]
            exonic_ens_ref = columns[header.index('ExonicFunc.ensGene')]
            AA_change_ensGene = columns[header.index('AAChange.ensGene')]
            snp138NonFlagged = columns[header.index('snp138NonFlagged')]
            apr_all = columns[header.index('ALL.sites.2015_08')]
            apr_eur = columns[header.index('EUR.sites.2015_08')]
            apr_amr = columns[header.index('AMR.sites.2015_08')]
            apr_asn = columns[header.index('EAS.sites.2015_08')]
            apr_afr = columns[header.index('AFR.sites.2015_08')]
            sample_gDNA = sampleID + ' chr' + Chr + ':' + start
            gDNA = 'chr' + Chr + ':' + start
            variant_key = Chr + ':' + start + '-' + end + ' ' + ref + '>' + alt
            if ref_gene_detail != 'NA':
                AA_change_refGene = ref_gene_detail
            if known_gene_detail != 'NA':
                AA_change_knownGene = known_gene_detail
            if ens_gene_detail != 'NA':
                AA_change_ensGgene = ens_gene_detail
            insert_file.write(str(gDNA) + "\t" + str(MRN) + "\t" + str(SEQ_CENTER) + "\t" + str(sampleID) + "\t" + str(Chr) + "\t" + str(start)\
                              + "\t" + str(end) + "\t" + str(ref) + "\t" + str(alt) + "\t" + str(snp138NonFlagged) + "\t" + str(func_ref_gene)\
                              + "\t" + str(gene_ref_gene) + "\t" + str(exonic_func_ref) + "\t" + str(AA_change_refGene) + "\t" + str(func_known_gene)\
                              + "\t" + str(gene_known_gene) + "\t" + str(exonic_known_ref) + "\t" + str(AA_change_knownGene)\
                              + "\t" + str(func_ens_gene) + "\t" + str(gene_ens_gene) + "\t" + str(exonic_ens_ref) + "\t" + str(AA_change_ensGene)\
                              + "\t" + str(apr_all) + "\t" + str(apr_eur) + "\t" + str(apr_amr) + "\t" + str(apr_asn) + "\t" + str(apr_afr)\
                              + "\t" + str(sample_gDNA) + "\t" + str(gDNA) + "\t" + str(SEQ_CENTER) + "\t" + str(tumor_type)\
                              + "\t" + str(SAMPLE_NOTE) + '_' + str(date) + "\t" + str(SOURCE) + "\t" + str(variant_key) + "\n")
        insert_file.close()
        snv.close()

        # This will format the  coverage information into a format that will be joined with the variant information
        print('Formatting varscan pile-ups')
        pileup = open('varscan.pileup')
        pileup_lines = pileup.readlines()
        header = pileup_lines.pop(0).strip().split('\t')
        insert_file = open('SQL_coverage_input.txt', 'w')
        date = datetime.datetime.now().replace(microsecond=0)
        # TODO remove unnecessary fields
        for line in pileup_lines:
            columns = line.rstrip('\n').split('\t')
            Chr = columns[header.index('Chrom')]
            start = columns[header.index('Position')]
            ref = columns[header.index('Ref')]
            alt = columns[header.index('Var')]
            columns2 = columns[header.index('Cons:Cov:Reads1:Reads2:Freq:P-value')].split(':')
            cons = columns2[0]
            cov = columns2[1]
            read1 = columns2[2]
            read2 = columns2[3]
            freq = columns2[4]
            p_val = columns2[5]
            columns3 = columns[header.index('StrandFilter:R1+:R1-:R2+:R2-:pval')].split(':')
            r1_plus = columns3[1]
            r1_minus = columns3[2]
            r2_plus = columns3[3]
            r2_minus = columns3[4]
            p_val2 = columns3[5]
            sample_gDNA = sampleID + ' chr' + Chr + ':' + start
            gDNA = 'chr' + Chr + ':' + start
            join_key = gDNA
            if re.search(r'-', alt):
                join_key = 'chr' + Chr + ':'+  str(int(start) + 1)
            insert_file.write(str(join_key) + "\t" + str(MRN) + "\t" + str(sampleID) + "\t" + str(SEQ_CENTER) + "\t" + str(Chr)\
                              + "\t" + str(start) + "\t" + str(ref) + "\t" + str(alt) + "\t" + str(cons) + "\t" + str(cov)\
                              + "\t" + str(read1) + "\t" + str(read2) + "\t" + str(freq) + "\t" + str(p_val) + "\t" + str(r1_plus)\
                              + "\t" + str(r1_minus) + "\t" + str(r2_plus) + "\t" + str(r2_minus) + "\t" + str(p_val2)\
                              + "\t" + str(SEQ_CENTER) + "\t" + str(sample_gDNA) + "\t" + str(gDNA) + "\t" + str(tumor_type)\
                              + "\t" + str(SAMPLE_NOTE) + '_' + str(date) + "\t" + str(SOURCE) + "\n")
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
        print('Creating final variants')
        joined_variants = open('joined_coverage_variants.txt')
        nonsyn_file = open('nonsyn_SQL_insert.txt', 'w')
        all_file = open('all_other_mutations.txt', 'w')
        header = 'NAME\tSEQ_CENTER\tSAMPLE_ID\tCHR\tSTART\tEND\tREF\tALT\tColumna1\tColumna2\tsnp138JJG\tTUMOR_READ1\tTUMOR_READ2' \
                 '\tTUMOR_VAR_FREQ\tAPR_ALL\tFunc_refGene\tGene_refGene\tExonicFunc_refGene\tAAChange_refGene\tFunc_knownGene' \
                 '\tGene_knownGene\tExonicFunc_knownGene\tAAChange_knownGene\tFunc_ensGene\tGene\texonic_func\tNT-AA_CHANGE' \
                 '\tAPR_AMR2\tAPR_ASN3\tAPR_AFR4\tColumna3\tColumna4\tColumna5\tREAD1_PLUS\tREAD1_MINUS\tREAD2_PLUS\tREAD2_MINUS' \
                 '\tColumna6\tColumna7\tSAMPLE_ID_CHR:START\tCHR:START\tSAMPLE_CENTER\tTUMOR\tSAMPLE_NOTE_DATE' \
                 '\tsource_of_RNA_used_for_sequencing\tVARIANT-KEY\n'
        nonsyn_file.write(header)
        all_file.write(header)
        # TODO remove unnecessary fields
        for line in joined_variants:
            columns = line.rstrip('\n').split('\t')
            mrn = columns[1]
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
            to_write = str(mrn) + "\t" + str(SEQ_CENTER) + "\t" + str(sampleID) + "\t" + str(Chr) + "\t" + str(Start) + "\t" + str(End)\
                       + "\t" + str(Ref) + "\t" + str(Alt) + "\t" + "\t" + "\t" + str(snp138JJG) + "\t" + str(tumor_reads1)\
                       + "\t" + str(tumor_reads2) + "\t" + str(tumor_var_freq) + "\t" + str(apr_all) + "\t" + str(Func_refGene)\
                       + "\t" + str(Gene_refGene) + "\t" + str(ExonicFunc_refGene) + "\t" + str(AAChange_refGene) + "\t" + str(Func_knownGene)\
                       + "\t" + str(Gene_knownGene) + "\t" + str(ExonicFunc_knownGene) + "\t" + str(AAChange_knownGene) + "\t" + str(Func_ensGene)\
                       + "\t" + str(Gene_ensGene) + "\t" + str(ExonicFunc_ensGene) + "\t" + str(AAChange_ensGene) + "\t" + str(apr_eur)\
                       + "\t" + str(apr_amr) + "\t" + str(apr_asn) + "\t" + str(apr_afr) + "\t" + "\t" + "\t" + "\t" + str(read1_plus)\
                       + "\t" + str(read1_minus) + "\t" + str(read2_plus) + "\t" + str(read2_minus) + "\t" + "\t" + "\t" + str(sample_gDNA)\
                       + "\t" +str(gDNA) + "\t" + str(SEQ_CENTER) + "\t" + str(tumor_type) + "\t" + str(SAMPLE_NOTE) + '_' + str(date)\
                       + "\t" + str(source_of_RNA_used_for_sequencing) + "\t" + str(variant_key) + "\n"
            if (re.search(r'nonsynonymous', ExonicFunc_refGene)or re.search(r'frame', ExonicFunc_refGene)or re.search(r'stop', ExonicFunc_refGene)\
                    or re.search(r'nonsynonymous', ExonicFunc_knownGene)or re.search(r'frame', ExonicFunc_knownGene) or re.search(r'stop', ExonicFunc_knownGene)\
                    or re.search(r'nonsynonymous', ExonicFunc_ensGene) or re.search(r'frame', ExonicFunc_ensGene)or re.search(r'stop', ExonicFunc_ensGene)):
                nonsyn_file.write(to_write)
            else:
                all_file.write(to_write)
        joined_variants.close()
        nonsyn_file.close()
        all_file.close()

        # Reformat FPKM file
        print('Creating FPKM info file')
        fpkm = open('genes.fpkm_tracking')
        fpkm_lines = fpkm.readlines()
        FPKM_ins = open('FPKM_SQL_insert.txt', 'w')
        firstline = fpkm.pop(0)
        date = datetime.datetime.now().replace(microsecond=0)
        header = 'NAME\tSEQ_CENTER\tSAMPLEID\tSOURCE\tTUMOUR\tSAMPLE_NOTE_DATE\tSAMPLEID\tSEQ_CENTER\t' + firstline
        FPKM_ins.write(header)
        # TODO remove unnecessary fields
        for line in fpkm_lines:
            FPKM_ins.write(str(MRN) + '\t' + str(SEQ_CENTER) + '\t' + str(sampleID) + '\t' + str(SOURCE) + '\t' + str(tumor_type)\
                           + '\t' + str(SAMPLE_NOTE) + '_' + str(date) + '\t' + str(sampleID) + ' ' + str(SEQ_CENTER) + '\t' + line)
        FPKM_ins.close()
        fpkm.close()

    print("COMPLETED!")