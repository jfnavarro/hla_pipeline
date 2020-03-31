import os
import subprocess
import pandas as pd
import re
import datetime
import multiprocessing
import time
import translate_dna

# All these programs must be in system or in PATH
PICARD = 'java -Xmx32g picard.jar'
SAMTOOLS  = 'samtools'
GATK = 'java -Xmx64g GenomeAnalysisTK.jar'
MUTECT = 'muTect'
STRELKA = 'configureStrelkaWorkflow.pl'
SSNIPER = 'bam-somaticsniper'
VARSCAN = 'java -Xmx32g -jar VarScan.jar'
HLA_PRG = "HLA-PRG-LA.pl"
#ANOVAR location must be in PATH
annovar_db = 'humandb -buildver hg19'
annovar_anno = 'refGene,knownGene,ensGene,snp138NonFlagged,1000g2012apr_all,1000g2012apr_eur,1000g2012apr_amr,1000g2012apr_asn,1000g2012apr_afr,cosmic70 -operation g,g,g,f,f,f,f,f,f,f -nastring NA'
# All these files will be assumed to be in HOME/shared for convenience
SNP_ANNO = "~/shared/snp.sum.hg19_multianno.txt"
INDEL_ANNO = "~/shared/indel.sum.hg19_multianno.txt'"
FASTA_AA_DICT = "~/shared/FASTA_AA.dict"
FASTA_cDNA = "~/shared/FASTA_cDNA.dict"
known_site1 = '~/shared/Mills_and_1000G_gold_standard.indels.b37.vcf'
known_site2 = '~/shared/1000G_phase1.indels.b37.vcf'
snpsites = '~/shared/dbsnp_139.b37.vcf.gz'

def index_column_substring(your_list, substring):
	for i, s in enumerate(your_list):
		if substring in s:
			return i
	return -1

def HLA_PRG(bamfile, sampleID):
	cpu_count = multiprocessing.cpu_count() - 1
	dirID = "out_HLA"
	cmd1 = HLA_PRG + ' --BAM {} --graph PRG_MHC_GRCh38_withIMGT --sampleID {} --maxThreads {}'.format(bamfile, dirID, cpu_count)
	p = subprocess.Popen(cmd1, shell = True)
	p.wait()

    # create a dictionary to store the output for each allele
	hla = pd.read_table(dirID + '/hla/R1_bestguess_G.txt', sep = '\t')
	allele_dict = {}
	hla = hla.groupby('Locus')
	for k,g in hla:
		allele = [re.sub('G$','', a).replace('N','') for a in g['Allele'].tolist()]
		allele_dict['HLA_' + k] = allele

    # Create formatted output file
	today = datetime.now().strftime('%B_%d_%Y')
	a = open('PRG-HLA-LA_output.txt', 'w')
	for x in sorted(allele_dict):
		a = open('PRG-HLA-LA_output.txt','a')
		a.write('{}\t{}\tExome {}\t{}\t{}\t{}\tPRG-HLA-LA\t-\t-\n'.format("mrn",sampleID,today,x,allele_dict[x][0],allele_dict[x][1]))
		a.close()

# Sample 2 cancer, sample 1 normal
# Genome must be hg19
def Full_exome_pipeline(sample1, sample2, tumor_type, genome, sampleID):

	WORKING_DIR = os.path.abspath(os.getcwd())
	os.chdir(WORKING_DIR)

	# Mark duplicates
	PICARD_TMP = os.path.join(WORKING_DIR, 'picard_tmp')
	os.makedirs(PICARD_TMP, exist_ok=True)
	cmdT_mark = PICARD + ' MarkDuplicates I=' + sample1 + 'O=sample1_dedup.bam M=dedup_sample1.txt TMP_DIR=' + PICARD_TMP
	cmdN_mark = PICARD + ' MarkDuplicates I=' + sample2 + 'O=sample2_dedup.bam M=dedup_sample2.txt TMP_DIR=' + PICARD_TMP
	print('Marking duplicates')
	p1 = subprocess.Popen(cmdT_mark, shell=True)
	p2 = subprocess.Popen(cmdN_mark, shell=True)
	p1.wait()
	p2.wait()
	print('Tumor and normal bam files had their optical and PCR duplicates marked.')

	# Add headers
	cmdT_RG = PICARD + ' AddOrReplaceReadGroups I=sample1_dedup.bam O=sample1_final.bam RGID=' + sample1 +\
			  ' RGPL=Illumina RGLB=library RGPU=' + sample1 + ' RGSM=' + sample1 + ' RGCN=VHIO RGDS=' + tumor_type + ' TMP_DIR=' + PICARD_TMP
	cmdN_RG = PICARD + ' AddOrReplaceReadGroups I=sample2_dedup.bam O=sample2_final.bam RGID=' + sample2 +\
			  ' RGPL=Illumina RGLB=library RGPU=' + sample2 + ' RGSM=' + sample2 + ' RGCN=VHIO RGDS=' + tumor_type + ' TMP_DIR=' + PICARD_TMP
	print('Adding headers')
	p1 = subprocess.Popen(cmdT_RG, shell=True)
	p2 = subprocess.Popen(cmdN_RG, shell=True)
	p1.wait()
	p2.wait()
	print('Tumor and normal bam files had read group information added.')

	# Create index
	cmdT_index = SAMTOOLS + ' index sample1_final.bam sample1_final.bai'
	cmdN_index = SAMTOOLS + ' index sample2_final.bam sample2_final.bai'
	print('Creating indexes')
	p1 = subprocess.Popen(cmdT_index, shell=True)
	p2 = subprocess.Popen(cmdN_index, shell=True)
	p1.wait()
	p2.wait()
	print('Tumor and normal final bam files were indexed.')

	# GATK create intervals
	cmdT_int = GATK + ' -T RealignerTargetCreator -I sample1_final.bam -R ' + genome +\
			   ' -known ' + known_site1 + ' -known ' + known_site2 + ' -o sample1.intervals &> sample1_target.log'
	cmdN_int = GATK + ' -T RealignerTargetCreator -I sample2_final.bam -R ' + genome +\
			   ' -known ' + known_site1 + ' -known ' + known_site2 + ' -o sample2.intervals &> sample2_target.log'
	print('Creating realignment target intervals')
	p1 = subprocess.Popen(cmdT_int, shell=True)
	p2 = subprocess.Popen(cmdN_int, shell=True)
	p1.wait()
	p2.wait()
	print('Indel target intervals were created for the tumor and normal samples.')

	# GATK re-align
	cmdT_realign = GATK + ' -T IndelRealigner -R ' + genome + ' -I sample1_final.bam -targetIntervals sample1.intervals ' +\
				   '-known ' + known_site1 + ' -known ' + known_site2 + ' -o sample1_realign.bam &> sample1_realign.log'
	cmdN_realign = GATK + ' -T IndelRealigner -R ' + genome + ' -I sample2_final.bam -targetIntervals sample2.intervals ' +\
				   '-known ' +known_site1 + ' -known ' + known_site2 + ' -o sample2_realign.bam &> sample2_realign.log'
	print('Starting realignments')
	p1 = subprocess.Popen(cmdT_realign, shell=True)
	p2 = subprocess.Popen(cmdN_realign, shell=True)
	p1.wait()
	p2.wait()
	print('Realignment around indels was performed on the tumor and normal samples.')

	# GATK base re-calibration
	print('Starting re-calibration')
	cmd1 = GATK + ' -T BaseRecalibrator -I sample1_realign.bam' + ' -R ' + genome + ' -knownSites ' + snpsites +\
		   ' -knownSites ' + known_site1 + ' -knownSites ' + known_site2 + ' -o sample1_recal_data.txt'
	cmd2 = GATK + ' -T BaseRecalibrator -I sample2_realign.bam' + ' -R ' + genome + ' -knownSites ' + snpsites +\
		   ' -knownSites ' + known_site1 + ' -knownSites ' + known_site2 + ' -o sample2_recal_data.txt '
	cmd3 = GATK + ' -T PrintReads -R ' + genome + ' -I sample1_realign.bam ' + ' -BQSR sample1_recal_data.txt -o sample1_recal.bam'
	cmd4 = GATK + ' -T PrintReads -R ' + genome + ' -I sample2_realign.bam ' + ' -BQSR sample2_recal_data.txt -o sample2_recal.bam'
	print('Starting recalibration...')
	p1 = subprocess.Popen(cmd1 + '; wait ; ' + cmd3, shell=True)
	p2 = subprocess.Popen(cmd2 + '; wait ; ' + cmd4, shell=True)
	p1.wait()
	p2.wait()
	print('Recalibration was performed on the tumor and normal samples.')

	# HLA predictions
	print('Performing HLA predictions')
	HLA_PRG('sample1_recal.bam', sampleID)
	HLA_PRG('sample2_recal.bam', sampleID)
	print('HLA predictions completed for tumor and normal samples')

	# Variant calling (SAMTOOLS PILEUP)
	print('Computing pile-ups')
	PILEUP_DIR = os.path.join(WORKING_DIR, "pileups")
	os.makedirs(PILEUP_DIR, exist_ok=True)
	cmd1 = SAMTOOLS + ' mpileup -C50 -B -q 1 -Q 15 -f ' + genome + ' sample1_recal.bam' + ' > ' + os.path.join(PILEUP_DIR, 'sample1.pileup')
	cmd2 = SAMTOOLS + ' mpileup -C50 -B -q 1 -Q 15 -f ' + genome + ' sample2_recal.bam' + ' > ' + os.path.join(PILEUP_DIR, 'sample2.pileup')
	p1 = subprocess.Popen(cmd1, shell=True)
	p2 = subprocess.Popen(cmd2, shell=True)
	print('Pile-ups were computed for tumor and normal samples')

	# Variant calling Mutect
	cmd_mutect = MUTECT + ' --memory 16g -T MuTect --reference_sequence ' + genome + ' --dbsnp ' + snpsites +\
				 ' --input_file:normal sample2_recal.bam --input_file:tumor sample1_recal.bam --vcf Mutect1.vcf'
	p3 = subprocess.Popen(cmd_mutect, shell=True)

	# Variant calling Strelka
	cmd_Strelka = STRELKA + ' --normal=sample2_recal.bam --tumor=sample1_recal.bam --ref=' +\
				  genome + ' --output-dir=Strelka_output ; cd Strelka_output ; make -j 4 ; cd..'
	p4 = subprocess.Popen(cmd_Strelka, shell=True)

	# Variant calling Somatic Sniper
	cmd_SomaticSniper = SSNIPER + ' -L -G -F vcf -f '+ genome + ' sample1_recal.bam sample2_recal.bam SS.vcf'
	p6 = subprocess.Popen(cmd_SomaticSniper, shell=True)

	# Wait till pileups are completed
	while p1.poll() is None or p2.poll() is None:
		time.sleep(300)

	# Variant calling VarScan
	cmd_varscan = VARSCAN + ' somatic ' + os.path.join(PILEUP_DIR, 'sample2.pileup') +\
				  ' ' + os.path.join(PILEUP_DIR, 'sample1.pileup') + ' . ' +\
				  '--tumor-purity .5 --output-vcf 1 --min-reads2 2 --min-coverage 4 --min-var-freq .05 --strand-filter 0'
	p7 = subprocess.Popen(cmd_varscan, shell=True)

	print('Performing variant calling')
	p3.wait()
	p4.wait()
	p6.wait()
	p7.wait()
	print('Done calling with Varscan, Mutect, SomaticSniper & Strelka.')

	##############################
	##Filtering Mutect snv calls##
	##############################
	print("Filtering Mutect SNV")
	filtered_vcf = open('mutect_filtered.vcf','w')
	filtered_vcf.close()
	vcf = open('Mutect1.vcf')
	for line in vcf:
		if line.startswith('#') and not line.startswith('#CHROM'):
			filtered_vcf = open('mutect_filtered.vcf', 'a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif line.startswith('#CHROM'):
			headers = line.strip().split('\t')
			Tmut = headers.index(sample1)
			Nmut = headers.index(sample2)
			filtered_vcf = open('mutect_filtered.vcf','a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif not line.startswith('#'):
			columns = line.strip().split('\t')
			Filter = columns[6]
			if re.search('PASS', Filter):
				normal_coverage = int(line.split('\t')[Nmut].split(':')[1].split(',')[0]) + int(line.split('\t')[Nmut].split(':')[1].split(',')[1])
				tumor_coverage = int(line.split('\t')[Tmut].split(':')[1].split(',')[0]) + int(line.split('\t')[Tmut].split(':')[1].split(',')[1])
				tumor_var_depth = int(line.split('\t')[Tmut].split(':')[1].split(',')[1])
				tumor_var_freq = float(float(line.split('\t')[Tmut].split(':')[4])*100)
				normal_var_freq = float(float(line.split('\t')[Nmut].split(':')[4])*100)
				if normal_var_freq != 0:
					t2n_ratio = float(tumor_var_freq) / float(normal_var_freq)
				else:
					t2n_ratio = 5
				if normal_coverage >=10 and tumor_coverage >= 10 and tumor_var_freq >= 2.5 and tumor_var_depth >= 3 and t2n_ratio >= 5:
					filtered_vcf=open('mutect_filtered.vcf','a')
					filtered_vcf.write(line)
					filtered_vcf.close()
	vcf.close()

	############################################################
	##Filtering Strelka snv calls and adding GT information in##
	############################################################
	print("Filtering Strelka SNV")
	filtered_vcf = open('strelka_filtered.vcf','w')
	filtered_vcf.close()
	vcf = open('Strelka_output/results/passed.somatic.snvs.vcf')
	for line in vcf:
		if line.startswith('#') and not re.search('##FORMAT=<ID=DP,', line) and not line.startswith('#CHROM'):
			filtered_vcf = open('strelka_filtered.vcf', 'a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif line.startswith('#') and re.search('##FORMAT=<ID=DP,', line):
			filtered_vcf = open('strelka_filtered.vcf', 'a')
			filtered_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' + line)
			filtered_vcf.close()
		elif line.startswith('#CHROM'):
			headers=line.strip().split('\t')
			Tst=headers.index('TUMOR')
			Nst=headers.index('NORMAL')
			filtered_vcf = open('strelka_filtered.vcf', 'a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif not line.startswith('#'):
			columns = line.strip().split('\t')
			ref = columns[3]
			alt = columns[4]
			Filter = columns[6]
			INFO = columns[7]
			Format = columns[8]
			Normal = columns[9]
			Tumor = columns[10]
			if re.search('PASS', Filter):
				if alt == 'A':
					normal_variant_depth = int(line.split('\t')[Nst].split(':')[4].split(',')[0])
					tumor_variant_depth = int(line.split('\t')[Tst].split(':')[4].split(',')[0])
				elif alt =='C':
					normal_variant_depth = int(line.split('\t')[Nst].split(':')[5].split(',')[0])
					tumor_variant_depth = int(line.split('\t')[Tst].split(':')[5].split(',')[0])
				elif alt == 'G':
					normal_variant_depth = int(line.split('\t')[Nst].split(':')[6].split(',')[0])
					tumor_variant_depth = int(line.split('\t')[Tst].split(':')[6].split(',')[0])
				elif alt == 'T':
					normal_variant_depth = int(line.split('\t')[Nst].split(':')[7].split(',')[0])
					tumor_variant_depth = int(line.split('\t')[Tst].split(':')[7].split(',')[0])
				n_cov = int(line.split('\t')[Nst].split(':')[0])
				t_cov = int(line.split('\t')[Tst].split(':')[0])
				T_freq = float((tumor_variant_depth/t_cov)*100)
				if normal_variant_depth !=0:
					N_freq = float(normal_variant_depth/n_cov)
				else:
					N_freq = 0
				if N_freq == 0:
					t2n_ratio = 5
				else:
					t2n_ratio = T_freq / N_freq
				if n_cov >= 10 and t_cov>=10 and tumor_variant_depth >=3 and T_freq >=5 and t2n_ratio >=5:
					Format = 'GT:' + Format
					INFO_split = INFO.split(';')
					SGT = INFO_split[3].replace('SGT=','').split('->')
					Normal_GT = ''
					Tumor_GT = ''
					i = 1
					for x in SGT[0]:
						if x == ref and i==1:
							Normal_GT = Normal_GT + '0'
							i = i + 1
						elif x == ref and i==2:
							Normal_GT = Normal_GT + '/0'
							i = i + 1
						elif x!= ref and i==1:
							Normal_GT = Normal_GT + '1'
							i = i + 1
						elif x != ref and i==2:
							Normal_GT = Normal_GT + '/1'
							i = i + 1
					i = 1
					for x in SGT[1]:
						if x == ref and i==1:
							Tumor_GT = Tumor_GT + '0'
							i = i + 1
						elif x == ref and i==2:
							Tumor_GT = Tumor_GT + '/0'
							i = i + 1
						elif x != ref and i == 1:
							Tumor_GT = Tumor_GT + '1'
							i = i + 1
						elif x != ref and i==2:
							Tumor_GT = Tumor_GT + '/1'
							i = i + 1
					filtered_vcf = open('strelka_filtered.vcf', 'a')
					filtered_vcf.write(str('\t'.join(columns[0:8])) + '\t' + Format + '\t' + Normal_GT + ':' + Normal + '\t' + Tumor_GT + ':' + Tumor + '\n')
					filtered_vcf.close()
	vcf.close()

	##############################
	##Filter SomaticSniper calls##
	##############################
	print("Filtering SomaticSniper SNV")
	filtered_vcf = open('somaticsniper_filtered.vcf','w')
	filtered_vcf.close()
	vcf = open('SS.vcf')
	for line in vcf:
		if line.startswith('#') and not line.startswith('#CHROM'):
			filtered_vcf = open('somaticsniper_filtered.vcf','a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif line.startswith('#CHROM'):
			headers = line.strip().split('\t')
			Tss = headers.index('TUMOR')
			Nss = headers.index('NORMAL')
			filtered_vcf = open('somaticsniper_filtered.vcf','a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif not line.startswith('#'):
			normal_coverage = int(line.split('\t')[Nss].split(':')[2])
			tumor_coverage = int(line.split('\t')[Tss].split(':')[2])
			somatic_status = int(line.split('\t')[Tss].split(':')[11])
			variant_count = int(line.split('\t')[Tss].split(':')[3].split(',')[2]) + int(line.split('\t')[Tss].split(':')[3].split(',')[3])
			normal_variant_count = int(line.split('\t')[Nss].split(':')[3].split(',')[2]) + int(line.split('\t')[Nss].split(':')[3].split(',')[3])
			Tumor_variant_freq = float((variant_count / tumor_coverage) * 100)
			if normal_variant_count !=0:
				normal_freq = float((normal_variant_count / normal_coverage) * 100)
				t2nratio = float(Tumor_variant_freq / normal_freq)
			else:
				normal_freq = 0
				t2nratio = 5
			if normal_coverage >= 10 and tumor_coverage >=10 and somatic_status ==2 and variant_count >= 3 and Tumor_variant_freq >= 5 and t2nratio >=5:
				filtered_vcf=open('somaticsniper_filtered.vcf','a')
				filtered_vcf.write(line)
				filtered_vcf.close()
	vcf.close()

	#########################################################
	##Edit varscan snp vcf and filter for somatic mutations##
	#########################################################'
	print("Filering Varscan SNP")
	filtered_vcf = open('varscan_filtered.vcf','w')
	filtered_vcf.close()
	vcf = open('snp.vcf')
	for line in vcf:
		if line.startswith('#') and re.search(r'DP4', line):
			new_DP4 = line.replace('ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev',
								   'ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases"')
			filtered_vcf = open('varscan_filtered.vcf','a')
			filtered_vcf.write(new_DP4)
			filtered_vcf.close()
		elif line.startswith('#') and not re.search(r'DP4', line) and not line.startswith('#CHROM'):
			filtered_vcf = open('varscan_filtered.vcf','a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif line.startswith('#CHROM'):
			headers = line.strip().split('\t')
			Tvs=headers.index('TUMOR')
			Nvs=headers.index('NORMAL')
			filtered_vcf = open('varscan_filtered.vcf','a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif re.search(r'SOMATIC', line, re.IGNORECASE):
			normal_coverage = int(line.split('\t')[Nvs].split(':')[2])
			tumor_coverage = int(line.split('\t')[Tvs].split(':')[2])
			tumor_var_depth = int(line.split('\t')[Tvs].split(':')[4])
			tumor_var_freq = float(line.split('\t')[Tvs].split(':')[5].replace('%',''))
			normal_var_freq = float(line.split('\t')[Nvs].split(':')[5].replace('%',''))
			p_value = float(line.split('\t')[7].split(';')[-1].replace('SPV=',''))
			if float(line.split('\t')[Nvs].split(':')[5].replace('%','')) !=0:
				t2n_ratio = float(line.split('\t')[Tvs].split(':')[5].replace('%','')) / float(line.split('\t')[Nvs].split(':')[5].replace('%',''))
			else:
				t2n_ratio = 5
			if normal_coverage >=10 and tumor_coverage >= 10 and tumor_var_depth >= 3 and t2n_ratio >= 5 and tumor_var_freq >=2.5:
				filtered_vcf = open('varscan_filtered.vcf','a')
				filtered_vcf.write(line)
				filtered_vcf.close()
	vcf.close()

	################################################################
	## Use GATK to combine all of the variants from various callers##
	################################################################
	print('Combining variants')
	cmd_GATK = GATK + ' -T CombineVariants -R ' + genome + ' -V:varscan varscan_filtered.vcf -V:mutect mutect_filtered.vcf '\
			   '-V:strelka strelka_filtered.vcf -V:somaticsniper somaticsniper_filtered.vcf -o combined_calls.vcf -genotypeMergeOptions UNIQUIFY'
	p_GATK = subprocess.Popen(cmd_GATK, shell=True)
	p_GATK.wait()

	####################################
	## Run annovar to annotate variants##
	####################################
	print('Running annovar')
	cmd1 = 'convert2annovar.pl -format vcf4old combined_calls.vcf --withzyg -outfile snp.av'
	cmd2 = 'table_annovar.pl snp.av ' + annovar_db + ' -out snp.sum' + ' -remove -protocol ' + annovar_anno
	p1 = subprocess.Popen(cmd1, shell=True)
	p1.wait()
	p2 = subprocess.Popen(cmd2, shell=True)
	p2.wait()

	######################################
	#Extract coverage info from vcf file##
	#to add to annotation data          ##
	######################################
	print("Extracting coverage from combined VCF")
	vcf = open('combined_calls.vcf')
	vcf_cov_dict = {}
	for line in vcf:
		if line.startswith('#CHROM'):
			headers = line.strip().split('\t')
			varscanT = headers.index('TUMOR.varscan')
			varscanN = headers.index('NORMAL.varscan')
			strelkaT = headers.index('TUMOR.strelka')
			strelkaN = headers.index('NORMAL.strelka')
			mutectT = headers.index(sample1+'.mutect')
			mutectN = headers.index(sample2+'.mutect')
			sniperT = headers.index('TUMOR.somaticsniper')
			sniperN = headers.index('NORMAL.somaticsniper')
		if not line.startswith('#'):
			columns = line.split('\t')
			chrm = columns[0]
			pos = columns[1]
			ref = columns[3]
			alt = columns[4]
			info = columns[7]
			DictID = chrm + ':' + pos
			callers = info.strip().split(';')[-1].replace('set=', '')
			if re.search('Intersection', callers) or re.search('varscan', callers):
				if callers == 'Intersection':
					caller_count = 4
					callers = 'varscan-strelka-mutect-somaticsniper'
				else:
					caller_count = callers.count('-') + 1
				if re.search('SPV=',info):
					info_split = info.split(';')
					for x in info_split:
						if re.search('SPV=',x):
							p_val = x.replace('SPV=','')
				else:
					p_val = '.'
				form = columns[8].split(':')
				DP4 = form.index('DP4')
				Freq = form.index('FREQ')
				trfor = int(columns[varscanT].split(':')[DP4].split(',')[0])
				trrev = int(columns[varscanT].split(':')[DP4].split(',')[1])
				tvfor = int(columns[varscanT].split(':')[DP4].split(',')[2])
				tvrev = int(columns[varscanT].split(':')[DP4].split(',')[3])
				tumor_read1 = trfor + trrev
				tumor_read2 = tvfor + tvrev
				tcov = trfor + trrev + tvfor + tvrev
				tfreq = columns[varscanT].split(':')[Freq]
				nrfor = int(columns[varscanN].split(':')[DP4].split(',')[0])
				nrrev = int(columns[varscanN].split(':')[DP4].split(',')[1])
				nvfor = int(columns[varscanN].split(':')[DP4].split(',')[2])
				nvrev = int(columns[varscanN].split(':')[DP4].split(',')[3])
				normal_read1 = nrfor + nrrev
				normal_read2 = nvfor + nvrev
				ncov = nrfor + nrrev + nvfor + nvrev
				nfreq = columns[varscanN].split(':')[Freq]
				vcf_cov_dict[chrm + ':' + pos] = {}
				vcf_cov_dict[chrm + ':' + pos]['pval'] = p_val
				vcf_cov_dict[chrm + ':' + pos]['Note'] = str(caller_count)+':'+callers
				vcf_cov_dict[chrm + ':' + pos]['trfor'] = trfor
				vcf_cov_dict[chrm + ':' + pos]['trrev'] = trrev
				vcf_cov_dict[chrm + ':' + pos]['tvfor'] = tvfor
				vcf_cov_dict[chrm + ':' + pos]['tvrev'] = tvrev
				vcf_cov_dict[chrm + ':' + pos]['nrfor'] = nrfor
				vcf_cov_dict[chrm + ':' + pos]['nrrev'] = nrrev
				vcf_cov_dict[chrm + ':' + pos]['nvfor']=nvfor
				vcf_cov_dict[chrm + ':' + pos]['nvrev']=nvrev
				vcf_cov_dict[chrm + ':' + pos]['tumor_read1'] = tumor_read1
				vcf_cov_dict[chrm + ':' + pos]['tumor_read2'] = tumor_read2
				vcf_cov_dict[chrm + ':' + pos]['normal_read1'] = normal_read1
				vcf_cov_dict[chrm + ':' + pos]['normal_read2'] = normal_read2
				vcf_cov_dict[chrm + ':' + pos]['tumor_coverage'] = tcov
				vcf_cov_dict[chrm + ':' + pos]['tumor_freq'] = tfreq
				vcf_cov_dict[chrm + ':' + pos]['normal_coverage'] = ncov
				vcf_cov_dict[chrm + ':' + pos]['normal_freq'] = nfreq
			elif re.search('somaticsniper', callers):
				caller_count=callers.count('-') + 1
				form = columns[8].split(':')
				DP4 = form.index('DP4')
				trfor = int(columns[sniperT].split(':')[DP4].split(',')[0])
				trrev = int(columns[sniperT].split(':')[DP4].split(',')[1])
				tvfor = int(columns[sniperT].split(':')[DP4].split(',')[2])
				tvrev = int(columns[sniperT].split(':')[DP4].split(',')[3])
				tumor_read1 = trfor+trrev
				tumor_read2 = tvfor+tvrev
				tcov = trfor + trrev + tvfor + tvrev
				nrfor = int(columns[sniperN].split(':')[DP4].split(',')[0])
				nrrev = int(columns[sniperN].split(':')[DP4].split(',')[1])
				nvfor = int(columns[sniperN].split(':')[DP4].split(',')[2])
				nvrev = int(columns[sniperN].split(':')[DP4].split(',')[3])
				normal_read1 = nrfor+nrrev
				normal_read2 = nvfor+nvrev
				ncov = nrfor + nrrev + nvfor + nvrev
				tfreq = str(round((tumor_read2/tcov)*100,2))+'%'
				nfreq = str(round((normal_read2/ncov)*100,2))+'%'
				vcf_cov_dict[chrm + ':' + pos] = {}
				vcf_cov_dict[chrm + ':' + pos]['pval'] = '.'
				vcf_cov_dict[chrm + ':' + pos]['Note'] = str(caller_count) + ':' + callers
				vcf_cov_dict[chrm + ':' + pos]['trfor'] = trfor
				vcf_cov_dict[chrm + ':' + pos]['trrev'] = trrev
				vcf_cov_dict[chrm + ':' + pos]['tvfor'] = tvfor
				vcf_cov_dict[chrm + ':' + pos]['tvrev'] = tvrev
				vcf_cov_dict[chrm + ':' + pos]['nrfor'] = nrfor
				vcf_cov_dict[chrm + ':' + pos]['nrrev'] = nrrev
				vcf_cov_dict[chrm + ':' + pos]['nvfor'] = nvfor
				vcf_cov_dict[chrm + ':' + pos]['nvrev'] = nvrev
				vcf_cov_dict[chrm + ':' + pos]['tumor_read1'] = tumor_read1
				vcf_cov_dict[chrm + ':' + pos]['tumor_read2'] = tumor_read2
				vcf_cov_dict[chrm + ':' + pos]['normal_read1'] = normal_read1
				vcf_cov_dict[chrm + ':' + pos]['normal_read2'] = normal_read2
				vcf_cov_dict[chrm + ':' + pos]['tumor_coverage'] = tcov
				vcf_cov_dict[chrm + ':' + pos]['tumor_freq'] = tfreq
				vcf_cov_dict[chrm + ':' + pos]['normal_coverage'] = ncov
				vcf_cov_dict[chrm + ':' + pos]['normal_freq'] = nfreq
			elif re.search('strelka', callers):
				caller_count = callers.count('-')+1
				form = columns[8].split(':')
				AU = form.index('AU')
				CU = form.index('CU')
				GU = form.index('GU')
				TU = form.index('TU')
				trfor = '.'
				trrev = '.'
				tvfor = '.'
				tvrev = '.'
				nrfor = '.'
				nrrev = '.'
				nvfor = '.'
				nvrev = '.'
				p_val = '.'
				if ref == 'A':
					tumor_read1 = int(columns[strelkaT].split(':')[AU].split(',')[0])
					normal_read1 = int(columns[strelkaN].split(':')[AU].split(',')[0])
				elif ref == 'C':
					tumor_read1 = int(columns[strelkaT].split(':')[CU].split(',')[0])
					normal_read1 = int(columns[strelkaN].split(':')[CU].split(',')[0])
				elif ref =='G':
					tumor_read1 = int(columns[strelkaT].split(':')[GU].split(',')[0])
					normal_read1 = int(columns[strelkaN].split(':')[GU].split(',')[0])
				elif ref =='T':
					tumor_read1 = int(columns[strelkaT].split(':')[TU].split(',')[0])
					normal_read1 = int(columns[strelkaN].split(':')[TU].split(',')[0])
				else:
					tumor_read1 = '.'
					normal_read1 = '.'
				if alt == 'A':
					tumor_read2 = int(columns[strelkaT].split(':')[AU].split(',')[0])
					normal_read2 = int(columns[strelkaN].split(':')[AU].split(',')[0])
				elif alt == 'C':
					tumor_read2 = int(columns[strelkaT].split(':')[CU].split(',')[0])
					normal_read2 = int(columns[strelkaN].split(':')[CU].split(',')[0])
				elif alt =='G':
					tumor_read2 = int(columns[strelkaT].split(':')[GU].split(',')[0])
					normal_read2 = int(columns[strelkaN].split(':')[GU].split(',')[0])
				elif alt =='T':
					tumor_read2 = int(columns[strelkaT].split(':')[TU].split(',')[0])
					normal_read2 = int(columns[strelkaN].split(':')[TU].split(',')[0])
				else:
					tumor_read2 = '.'
					normal_read2 = '.'
				if tumor_read2 != '.':
					tcov = tumor_read1 + tumor_read2
					ncov = normal_read1+normal_read2
					tfreq = str(round((tumor_read2/tcov)*100,2))+'%'
					nfreq = str(round((normal_read2/ncov)*100,2))+'%'
				else:
					tcov = '.'
					ncov = '.'
					tfreq = '.'
					nfreq = '.'
				vcf_cov_dict[chrm + ':' + pos] = {}
				vcf_cov_dict[chrm + ':' + pos]['Note'] = str(caller_count) + ':' + callers
				vcf_cov_dict[chrm + ':' + pos]['pval'] = '.'
				vcf_cov_dict[chrm + ':' + pos]['trfor'] = trfor
				vcf_cov_dict[chrm + ':' + pos]['trrev'] = trrev
				vcf_cov_dict[chrm + ':' + pos]['tvfor'] = tvfor
				vcf_cov_dict[chrm + ':' + pos]['tvrev'] = tvrev
				vcf_cov_dict[chrm + ':' + pos]['nrfor'] = nrfor
				vcf_cov_dict[chrm + ':' + pos]['nrrev'] = nrrev
				vcf_cov_dict[chrm + ':' + pos]['nvfor'] = nvfor
				vcf_cov_dict[chrm + ':' + pos]['nvrev'] = nvrev
				vcf_cov_dict[chrm + ':' + pos]['tumor_read1'] = tumor_read1
				vcf_cov_dict[chrm + ':' + pos]['tumor_read2'] = tumor_read2
				vcf_cov_dict[chrm + ':' + pos]['normal_read1'] = normal_read1
				vcf_cov_dict[chrm + ':' + pos]['normal_read2'] = normal_read2
				vcf_cov_dict[chrm + ':' + pos]['tumor_coverage'] = tcov
				vcf_cov_dict[chrm + ':' + pos]['tumor_freq'] = tfreq
				vcf_cov_dict[chrm + ':' + pos]['normal_coverage'] = ncov
				vcf_cov_dict[chrm + ':' + pos]['normal_freq'] = nfreq
			elif re.search('mutect', callers):
				caller_count = callers.count('-') + 1
				form = columns[8].split(':')
				if not re.search(',', alt):
					AD = form.index('AD')
				DP = form.index('DP')
				FA = form.index('FA')
				trfor = '.'
				trrev = '.'
				tvfor = '.'
				tvrev = '.'
				nrfor = '.'
				nrrev = '.'
				nvfor = '.'
				nvrev = '.'
				p_val = '.'
				if not re.search(',', alt):
					tumor_read1 = int(columns[mutectT].split(':')[AD].split(',')[0])
					tumor_read2 = int(columns[mutectT].split(':')[AD].split(',')[1])
					normal_read1 = int(columns[mutectN].split(':')[AD].split(',')[0])
					normal_read2 = int(columns[mutectN].split(':')[AD].split(',')[1])
					tcov = tumor_read1 + tumor_read2
					ncov = normal_read1 + normal_read2
					tfreq = str(round((tumor_read2/tcov)*100,2))+'%'
					nfreq = str(round((normal_read2/ncov)*100,2))+'%'
				else:
					tumor_read1 = 0
					tumor_read2 = 0
					normal_read1 = 0
					normal_read2 = 0
					tcov = 0
					ncov = 0
					tfreq = str(0) + '%'
					nfreq = str(0) + '%'
				vcf_cov_dict[chrm + ':' + pos]={}
				vcf_cov_dict[chrm + ':' + pos]['Note'] = str(caller_count) + ':' + callers
				vcf_cov_dict[chrm + ':' + pos]['pval'] = '.'
				vcf_cov_dict[chrm + ':' + pos]['trfor'] = trfor
				vcf_cov_dict[chrm + ':' + pos]['trrev'] = trrev
				vcf_cov_dict[chrm + ':' + pos]['tvfor'] = tvfor
				vcf_cov_dict[chrm + ':' + pos]['tvrev'] = tvrev
				vcf_cov_dict[chrm + ':' + pos]['nrfor'] = nrfor
				vcf_cov_dict[chrm + ':' + pos]['nrrev'] = nrrev
				vcf_cov_dict[chrm + ':' + pos]['nvfor'] = nvfor
				vcf_cov_dict[chrm + ':' + pos]['nvrev'] = nvrev
				vcf_cov_dict[chrm + ':' + pos]['tumor_read1'] = tumor_read1
				vcf_cov_dict[chrm + ':' + pos]['tumor_read2'] = tumor_read2
				vcf_cov_dict[chrm + ':' + pos]['normal_read1'] = normal_read1
				vcf_cov_dict[chrm + ':' + pos]['normal_read2'] = normal_read2
				vcf_cov_dict[chrm + ':' + pos]['tumor_coverage'] = tcov
				vcf_cov_dict[chrm + ':' + pos]['tumor_freq'] = tfreq
				vcf_cov_dict[chrm + ':' + pos]['normal_coverage'] = ncov
				vcf_cov_dict[chrm + ':' + pos]['normal_freq'] = nfreq
	vcf.close()

	########################################
	##Now that I have a dict I need to    ##
	##produce a final sheet with coverage ##
	##that Will go into our database      ##
	########################################
	print("Formatting coverage")
	nonsyn_snv = open(SNP_ANNO)
	nonsyn_file = open('nonsyn_SQL_insert.txt', 'w')
	nonsyn_file.close()
	date = datetime.datetime.now().replace(microsecond=0)
	date_trunc = datetime.datetime.now().date()
	first_line = True
	for line in nonsyn_snv:
		if first_line:
			first_line= False
			continue
		columns = line.strip().split('\t')
		if (re.search(r'nonsynonymous', columns[8])or re.search(r'frame', columns[8])or re.search(r'stop', columns[8])\
			or re.search(r'nonsynonymous', columns[13])or re.search(r'frame', columns[13]) or re.search(r'stop', columns[13])\
			or re.search(r'nonsynonymous', columns[18]) or re.search(r'frame', columns[18])or re.search(r'stop', columns[18])):
			Chr = columns[0]
			start = columns[1]
			ID = Chr + ':' + start
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
			COSMIC = columns[26]
			mrn = mrn
			seq_center = seq_center
			sampleID = sampleID
			sample_gDNA = sampleID + ' chr' + Chr + ':' + start
			gDNA = 'chr' + Chr + ':' + start
			tumor_type = tumor_type
			source = source
			sample_note = sample_note
			sample_center = sample_note + ' ' + source + ' ' + seq_center
			variant_key = str(Chr) + ':' + str(start) + '-' + str(end) + ' ' + str(ref) + '>' + str(alt)
			if ref_gene_detail != 'NA':
				columns[9] = columns[7]
			if known_gene_detail != 'NA':
				columns[14] = columns[12]
			if ens_gene_detail != 'NA':
				columns[19] = columns[17]
			p_val = vcf_cov_dict[ID]['pval']
			Note = vcf_cov_dict[ID]['Note']
			trfor = vcf_cov_dict[ID]['trfor']
			trrev = vcf_cov_dict[ID]['trrev']
			tvfor = vcf_cov_dict[ID]['tvfor']
			tvrev = vcf_cov_dict[ID]['tvrev']
			nrfor = vcf_cov_dict[ID]['nrfor']
			nrrev = vcf_cov_dict[ID]['nrrev']
			nvfor = vcf_cov_dict[ID]['nvfor']
			nvrev = vcf_cov_dict[ID]['nvrev']
			tumor_read1 = vcf_cov_dict[ID]['tumor_read1']
			tumor_read2 = vcf_cov_dict[ID]['tumor_read2']
			normal_read1 = vcf_cov_dict[ID]['normal_read1']
			normal_read2 = vcf_cov_dict[ID]['normal_read2']
			tcov = vcf_cov_dict[ID]['tumor_coverage']
			tfreq = vcf_cov_dict[ID]['tumor_freq']
			ncov = vcf_cov_dict[ID]['normal_coverage']
			nfreq = vcf_cov_dict[ID]['normal_freq']
			nonsyn_file = open('nonsyn_SQL_insert.txt', 'a')
			nonsyn_file.write(str(mrn)+"\t"+str(seq_center)+"\t"+str(sampleID)+"\t"+str('\t'.join(columns[0:5]))+'\t-\t-\t'+str(columns[20])+'\t'+str(tumor_read1)+'\t'+str(tumor_read2)+'\t'\
							+str('\t'.join(columns[5:7]))+'\t'+str('\t'.join(columns[8:12]))+'\t'+str('\t'.join(columns[13:17]))+'\t'+str('\t'.join(columns[18:20]))+'\t'\
							+str('\t'.join(columns[21:26]))+'\t'+str(normal_read1)+'\t'+str(normal_read2)+'\t'+str(trfor)+'\t'+str(trrev)+'\t'+str(tvfor)+'\t'+str(tvrev)+'\t'\
							+str(nrfor)+'\t'+str(nrrev)+'\t'+str(nvfor)+'\t'+str(nfreq)+'\t'+str(nvrev)+'\t'+str(tfreq)+'\tSomatic\t'+str(p_val)+'\t'+str(sampleID)+' chr'+str(Chr)+':'+str(start)+'\t'+str(Note)+'\t'\
							+str(gDNA)+'\t'+str(tumor_type)+'\t'+str(source)+'\t'+str(sample_note)+' '+str(source)+' '+str(seq_center)+'\t'+str(tcov)+'\t'+str(ncov)+'\t'+str(sample_note)+'\t'+str(variant_key)+'\t'+str(COSMIC)+'\t'+str(date)+'\t'+resection_date+'\t'+run_date+'\t'+sequencer+'\t'+kit+'\t'+note+'\t'+index+'\n')
			nonsyn_file.close()
		else:
			Chr = columns[0]
			start = columns[1]
			ID = Chr + ':' + start
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
			COSMIC = columns[26]
			mrn = mrn
			seq_center = seq_center
			sampleID = sampleID
			sample_gDNA = sampleID + ' chr' + Chr + ':' + start
			gDNA = 'chr' + Chr + ':' + start
			tumor_type = tumor_type
			source = source
			sample_note = sample_note
			sample_center = sample_note + ' ' + source + ' ' + seq_center
			variant_key = str(Chr) + ':' + str(start) + '-' + str(end) + ' ' + str(ref) + '>' + str(alt)
			if ref_gene_detail != 'NA':
				columns[9] = columns[7]
			if known_gene_detail != 'NA':
				columns[14] = columns[12]
			if ens_gene_detail != 'NA':
				columns[19] = columns[17]
			p_val = vcf_cov_dict[ID]['pval']
			Note = vcf_cov_dict[ID]['Note']
			trfor = vcf_cov_dict[ID]['trfor']
			trrev = vcf_cov_dict[ID]['trrev']
			tvfor = vcf_cov_dict[ID]['tvfor']
			tvrev = vcf_cov_dict[ID]['tvrev']
			nrfor = vcf_cov_dict[ID]['nrfor']
			nrrev = vcf_cov_dict[ID]['nrrev']
			nvfor = vcf_cov_dict[ID]['nvfor']
			nvrev = vcf_cov_dict[ID]['nvrev']
			tumor_read1 = vcf_cov_dict[ID]['tumor_read1']
			tumor_read2 = vcf_cov_dict[ID]['tumor_read2']
			normal_read1 = vcf_cov_dict[ID]['normal_read1']
			normal_read2 = vcf_cov_dict[ID]['normal_read2']
			tcov = vcf_cov_dict[ID]['tumor_coverage']
			tfreq = vcf_cov_dict[ID]['tumor_freq']
			ncov = vcf_cov_dict[ID]['normal_coverage']
			nfreq = vcf_cov_dict[ID]['normal_freq']
			all_file = open('all_other_mutations.txt', 'a')
			all_file.write(str(mrn)+"\t"+str(seq_center)+"\t"+str(sampleID)+"\t"+str('\t'.join(columns[0:5]))+'\t-\t-\t'+str(columns[20])+'\t'+str(tumor_read1)+'\t'+str(tumor_read2)+'\t'\
							+str('\t'.join(columns[5:7]))+'\t'+str('\t'.join(columns[8:12]))+'\t'+str('\t'.join(columns[13:17]))+'\t'+str('\t'.join(columns[18:20]))+'\t'\
							+str('\t'.join(columns[21:26]))+'\t'+str(normal_read1)+'\t'+str(normal_read2)+'\t'+str(trfor)+'\t'+str(trrev)+'\t'+str(tvfor)+'\t'+str(tvrev)+'\t'\
							+str(nrfor)+'\t'+str(nrrev)+'\t'+str(nvfor)+'\t'+str(nfreq)+'\t'+str(nvrev)+'\t'+str(tfreq)+'\tSomatic\t'+str(p_val)+'\t'+str(sampleID)+' chr'+str(Chr)+':'+str(start)+'\t'+str(Note)+'\t'\
							+str(gDNA)+'\t'+str(tumor_type)+'\t'+str(source)+'\t'+str(sample_note)+' '+str(source)+' '+str(seq_center)+'\t'+str(tcov)+'\t'+str(ncov)+'\t'+str(sample_note)+'\t'+str(variant_key)+'\t'+str(COSMIC)+'\t'+str(date)+'\t'+resection_date+'\t'+run_date+'\t'+sequencer+'\t'+kit+'\t'+note+'\t'+index+'\n')
			all_file.close()

	########################################
	##create a final vcf file to be stored##
	########################################
	print("Formatting final VCF")
	vcf = open('combined_calls.vcf')
	vcf_final = open('final.vcf', 'w')
	vcf_final.close()
	for line in vcf:
		if line.startswith('##file') or line.startswith('##FILTER') or line.startswith('##FORMAT=<ID=GT')\
		or line.startswith('##FORMAT=<ID=AD') or line.startswith('##FORMAT=<ID=BQ') or line.startswith('##FORMAT=<ID=DP,')\
		or line.startswith('##FORMAT=<ID=DP') or line.startswith('##FORMAT=<ID=FA') or line.startswith('##INFO') or line.startswith('##contig'):
			vcf_final = open('final.vcf', 'a')
			vcf_final.write(line)
			vcf_final.close()
		if line.startswith('#CHROM'):
			headers = line.strip().split('\t')
			varscanT = headers.index('TUMOR.varscan')
			varscanN = headers.index('NORMAL.varscan')
			strelkaT = headers.index('TUMOR.strelka')
			strelkaN = headers.index('NORMAL.strelka')
			mutectT = headers.index(sample1+'.mutect')
			mutectN = headers.index(sample2+'.mutect')
			sniperT = headers.index('TUMOR.somaticsniper')
			sniperN = headers.index('NORMAL.somaticsniper')
			vcf_final = open('final.vcf', 'a')
			vcf_final.write(str('\t'.join(headers[0:9])) + '\tTUMOR\tNORMAL\n')
			vcf_final.close()
		elif not line.startswith('#'):
			columns = line.split('\t')
			chrm = columns[0]
			pos = columns[1]
			ref = columns[3]
			alt = columns[4]
			info = columns[7]
			callers = info.strip().split(';')[-1].replace('set=', '')
			if re.search('Intersection', callers) or re.search('varscan', callers):
				if not re.search(',', alt):
					form = columns[8].split(':')
					GT = form.index('GT')
					AD = form.index('AD')
					DP = form.index('DP')
					FREQ = form.index('FREQ')
					vcf_final = open('final.vcf', 'a')
					vcf_final.write(str('\t'.join(columns[0:8]))+'\tGT:AD:DP:FREQ\t'+columns[varscanT].split(':')[GT]+':'+ columns[varscanT].split(':')[AD]+':'+columns[varscanT].split(':')[DP]+':'+columns[varscanT].split(':')[FREQ]+'\t'+\
									columns[varscanN].split(':')[GT]+':'+ columns[varscanN].split(':')[AD]+':'+columns[varscanN].split(':')[DP]+':'+columns[varscanN].split(':')[FREQ]+'\n')
					vcf_final.close()
				else:
					vcf_final = open('final.vcf', 'a')
					vcf_final.write(str('\t'.join(columns[0:8]))+'\tGT:AD:DP:FREQ\t.:.:.:.\t.:.:.:.\n')
					vcf_final.close()
			elif re.search('somaticsniper', callers):
				if not re.search(',', alt):
					form = columns[8].split(':')
					DP4 = form.index('DP4')
					GT = form.index('GT')
					trfor = int(columns[sniperT].split(':')[DP4].split(',')[0])
					trrev = int(columns[sniperT].split(':')[DP4].split(',')[1])
					tvfor = int(columns[sniperT].split(':')[DP4].split(',')[2])
					tvrev = int(columns[sniperT].split(':')[DP4].split(',')[3])
					tumor_read1 = trfor + trrev
					tumor_read2 = tvfor + tvrev
					tcov = trfor + trrev + tvfor + tvrev
					nrfor = int(columns[sniperN].split(':')[DP4].split(',')[0])
					nrrev = int(columns[sniperN].split(':')[DP4].split(',')[1])
					nvfor = int(columns[sniperN].split(':')[DP4].split(',')[2])
					nvrev = int(columns[sniperN].split(':')[DP4].split(',')[3])
					normal_read1 = nrfor + nrrev
					normal_read2 = nvfor + nvrev
					ncov = nrfor + nrrev + nvfor + nvrev
					tfreq = str(round((tumor_read2/tcov)*100,2))+'%'
					nfreq = str(round((normal_read2/ncov)*100,2))+'%'
					vcf_final = open('final.vcf', 'a')
					vcf_final.write(str('\t'.join(columns[0:8]))+'\tGT:AD:DP:FREQ\t'+str(columns[sniperT].split(':')[GT])+':'+ str(tumor_read2)+':'+str(tcov)+':'+str(tfreq)+'\t'+\
									str(columns[sniperN].split(':')[GT])+':'+ str(normal_read2)+':'+str(ncov)+':'+str(nfreq)+'\n')
					vcf_final.close()
				else:
					vcf_final = open('final.vcf', 'a')
					vcf_final.write(str('\t'.join(columns[0:8]))+'\tGT:AD:DP:FREQ\t.:.:.:.\t.:.:.:.\n')
					vcf_final.close()
			elif re.search('mutect', callers):
				if not re.search(',', alt):
					form = columns[8].split(':')
					GT = form.index('GT')
					AD = form.index('AD')
					DP = form.index('DP')
					FA = form.index('FA')
					vcf_final = open('final.vcf', 'a')
					vcf_final.write(str('\t'.join(columns[0:8]))+'\tGT:AD:DP:FREQ\t'+columns[mutectT].split(':')[GT]+':'+ columns[mutectT].split(':')[AD].split(',')[1]+':'+columns[mutectT].split(':')[DP]+':'+str(float(columns[mutectT].split(':')[FA])*100)+'%\t'+\
							columns[mutectN].split(':')[GT]+':'+ columns[mutectN].split(':')[AD].split(',')[1]+':'+columns[mutectN].split(':')[DP]+':'+str(float(columns[mutectN].split(':')[FA])*100)+'%\n')
					vcf_final.close()
				else:
					vcf_final = open('final.vcf', 'a')
					vcf_final.write(str('\t'.join(columns[0:8]))+'\tGT:AD:DP:FREQ\t.:.:.:.\t.:.:.:.\n')
					vcf_final.close()
			elif re.search('strelka', callers):
				if not re.search(',', alt):
					form = columns[8].split(':')
					GT = form.index('GT')
					AU = form.index('AU')
					CU = form.index('CU')
					GU = form.index('GU')
					TU = form.index('TU')
					if ref == 'A':
						tumor_read1 = int(columns[strelkaT].split(':')[AU].split(',')[0])
						normal_read1 = int(columns[strelkaN].split(':')[AU].split(',')[0])
					elif ref == 'C':
						tumor_read1 = int(columns[strelkaT].split(':')[CU].split(',')[0])
						normal_read1 = int(columns[strelkaN].split(':')[CU].split(',')[0])
					elif ref =='G':
						tumor_read1 = int(columns[strelkaT].split(':')[GU].split(',')[0])
						normal_read1 = int(columns[strelkaN].split(':')[GU].split(',')[0])
					elif ref =='T':
						tumor_read1 = int(columns[strelkaT].split(':')[TU].split(',')[0])
						normal_read1 = int(columns[strelkaN].split(':')[TU].split(',')[0])
					else:
						tumor_read1 = '.'
						normal_read1 = '.'
					if alt == 'A':
						tumor_read2 = int(columns[strelkaT].split(':')[AU].split(',')[0])
						normal_read2 = int(columns[strelkaN].split(':')[AU].split(',')[0])
					elif alt == 'C':
						tumor_read2 = int(columns[strelkaT].split(':')[CU].split(',')[0])
						normal_read2 = int(columns[strelkaN].split(':')[CU].split(',')[0])
					elif alt =='G':
						tumor_read2 = int(columns[strelkaT].split(':')[GU].split(',')[0])
						normal_read2 = int(columns[strelkaN].split(':')[GU].split(',')[0])
					elif alt =='T':
						tumor_read2 = int(columns[strelkaT].split(':')[TU].split(',')[0])
						normal_read2 = int(columns[strelkaN].split(':')[TU].split(',')[0])
					else:
						tumor_read2 = '.'
						normal_read2 = '.'
					if tumor_read2 != '.':
						tcov = tumor_read1 + tumor_read2
						ncov = normal_read1 + normal_read2
						tfreq = str(round((tumor_read2/tcov)*100,2))+'%'
						nfreq = str(round((normal_read2/ncov)*100,2))+'%'
					vcf_final = open('final.vcf', 'a')
					vcf_final.write(str('\t'.join(columns[0:8]))+'\tGT:AD:DP:FREQ\t'+str(columns[strelkaT].split(':')[GT])+':'+ str(tumor_read2)+':'+str(tcov)+':'+str(tfreq)+'\t'+\
									str(columns[strelkaN].split(':')[GT])+':'+ str(normal_read2)+':'+str(ncov)+':'+str(nfreq)+'\n')
					vcf_final.close()
				else:
					vcf_final=open('final.vcf', 'a')
					vcf_final.write(str('\t'.join(columns[0:8]))+'\tGT:AD:DP:FREQ\t.:.:.:.\t.:.:.:.\n')
					vcf_final.close()
	vcf.close()

	###########################################
	##need to repeat everything done for snvs##
	##for indels.  only two of the callers   ##
	##used call indels                       ##
	###########################################
	print("Filtering Strelka indels")
	filtered_vcf = open('strelka_indel_filtered.vcf','w')
	filtered_vcf.close()
	vcf = open('Strelka_output/results/passed.somatic.indels.vcf')
	for line in vcf:
		if line.startswith('#') and not re.search('##FORMAT=<ID=DP,', line) and not line.startswith('#CHROM'):
			filtered_vcf = open('strelka_indel_filtered.vcf', 'a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif line.startswith('#') and re.search('##FORMAT=<ID=DP,', line):
			filtered_vcf = open('strelka_indel_filtered.vcf', 'a')
			filtered_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
			filtered_vcf.close()
		elif line.startswith('#CHROM'):
			headers = line.strip().split('\t')
			Tst=headers.index('TUMOR')
			Nst=headers.index('NORMAL')
			filtered_vcf = open('strelka_indel_filtered.vcf', 'a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif not line.startswith('#'):
			columns = line.strip().split('\t')
			ref = columns[3]
			alt = columns[4]
			Filter = columns[6]
			INFO = columns[7]
			Format = columns[8]
			Normal = columns[9]
			Tumor = columns[10]
			if re.search('PASS', Filter):
				normal_variant_depth = int(line.split('\t')[Nst].split(':')[3].split(',')[1])
				tumor_variant_depth = int(line.split('\t')[Tst].split(':')[3].split(',')[1])
				n_cov = int(line.split('\t')[Nst].split(':')[0])
				t_cov = int(line.split('\t')[Tst].split(':')[0])
				T_freq = float((tumor_variant_depth/t_cov)*100)
				if normal_variant_depth !=0:
					N_freq = float(normal_variant_depth/n_cov)
				else:
					N_freq = 0
				if N_freq == 0:
					t2n_ratio = 5
				else:
					t2n_ratio = T_freq/N_freq
				if n_cov >= 10 and t_cov>=10 and tumor_variant_depth >=3 and T_freq >=5 and t2n_ratio >=5:
					Format = 'GT:' + Format
					INFO_split = INFO.split(';')
					SGT_index = index_column_substring(INFO_split, 'SGT')
					SGT = INFO_split[SGT_index].replace('SGT=','').split('->')
					Normal_GT = '0/0'
					if SGT[1] == 'het':
						Tumor_GT = '0/1'
					else:
						Tumor_GT = '1/1'
					filtered_vcf=open('strelka_indel_filtered.vcf', 'a')
					filtered_vcf.write(str('\t'.join(columns[0:8]))+'\t'+Format+'\t'+Normal_GT+':'+Normal+'\t'+Tumor_GT+':'+Tumor+'\n')
					filtered_vcf.close()
	vcf.close()

	#####################################
	##Edit varscan indel vcf and filter##
	##for somatic mutations            ##
	#####################################
	print("Filtering Varscan indels")
	filtered_vcf = open('varscan_filtered_indel.vcf','w')
	filtered_vcf.close()
	vcf = open('indel.vcf')
	for line in vcf:
		if line.startswith('#') and re.search(r'DP4', line):
			new_DP4=line.replace('ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev',
								 'ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases"')
			filtered_vcf = open('varscan_filtered_indel.vcf','a')
			filtered_vcf.write(new_DP4)
			filtered_vcf.close()
		elif line.startswith('#') and not re.search(r'DP4', line) and not line.startswith('#CHROM'):
			filtered_vcf = open('varscan_filtered_indel.vcf','a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif line.startswith('#CHROM'):
			headers = line.strip().split('\t')
			Tvs = headers.index('TUMOR')
			Nvs = headers.index('NORMAL')
			filtered_vcf = open('varscan_filtered_indel.vcf','a')
			filtered_vcf.write(line)
			filtered_vcf.close()
		elif re.search(r'SOMATIC', line, re.IGNORECASE):
			normal_coverage = int(line.split('\t')[Nvs].split(':')[2])
			tumor_coverage = int(line.split('\t')[Tvs].split(':')[2])
			tumor_var_depth = int(line.split('\t')[Tvs].split(':')[4])
			tumor_var_freq = float(line.split('\t')[Tvs].split(':')[5].replace('%',''))
			normal_var_freq = float(line.split('\t')[Nvs].split(':')[5].replace('%',''))
			p_value = float(line.split('\t')[7].split(';')[-1].replace('SPV=',''))
			if float(line.split('\t')[Nvs].split(':')[5].replace('%','')) !=0:
				t2n_ratio = float(line.split('\t')[Tvs].split(':')[5].replace('%','')) / float(line.split('\t')[Nvs].split(':')[5].replace('%',''))
			else:
				t2n_ratio = 5
			if normal_coverage >=10 and tumor_coverage >= 10 and tumor_var_depth >= 3 and t2n_ratio >= 5 and tumor_var_freq >=2.5:
				filtered_vcf = open('varscan_filtered_indel.vcf','a')
				filtered_vcf.write(line)
				filtered_vcf.close()
	vcf.close()

	##################################
	##Use GATK to combine all of the##
	##indel variants from callers   ##
	##################################
	print('Combining indels variants')
	cmd_GATK = GATK + ' -T CombineVariants -R ' + genome + ' -V:varscan varscan_filtered_indel.vcf ' +\
			   '-V:strelka strelka_indel_filtered.vcf -o combined_indel_calls.vcf -genotypeMergeOptions UNIQUIFY '
	p_GATK = subprocess.Popen(cmd_GATK, shell=True)
	p_GATK.wait()

	############################
	##Run annovar to annotate ##
	## indel variants         ##
	############################
	print('Annotating combined indels with annovar')
	cmd1 = 'convert2annovar.pl -format vcf4old combined_indel_calls.vcf --withzyg -outfile indel.av'
	cmd2 = 'table_annovar.pl indel.av ' + annovar_db + ' -out indel.sum -remove -protocol ' + annovar_anno
	p1 = subprocess.Popen(cmd1, shell=True)
	p1.wait()
	p2 = subprocess.Popen(cmd2, shell=True)
	p2.wait()

	######################################
	##extract coverage info from vcf file#
	######################################
	print("Extracting coverage info")
	vcf = open('combined_indel_calls.vcf')
	vcf_cov_dict = {}
	for line in vcf:
		if line.startswith('#CHROM'):
			headers = line.strip().split('\t')
			varscanT = headers.index('TUMOR.varscan')
			varscanN = headers.index('NORMAL.varscan')
			strelkaT = headers.index('TUMOR.strelka')
			strelkaN = headers.index('NORMAL.strelka')
		if not line.startswith('#'):
			columns = line.split('\t')
			chrm = columns[0]
			pos = columns[1]
			ref = columns[3]
			alt = columns[4]
			info = columns[7]
			DictID = chrm + ':' + pos
			callers = info.strip().split(';')[-1].replace('set=', '')
			if re.search('Intersection', callers) or re.search('varscan', callers):
				if callers == 'Intersection':
					caller_count = 2
					callers = 'varscan-strelka'
				else:
					caller_count = callers.count('-')+1
				if re.search('SPV=',info):
					info_split = info.split(';')
					for x in info_split:
						if re.search('SPV=',x):
							p_val = x.replace('SPV=','')
				else:
					p_val = '.'
				form = columns[8].split(':')
				DP4 = form.index('DP4')
				Freq = form.index('FREQ')
				trfor = int(columns[varscanT].split(':')[DP4].split(',')[0])
				trrev = int(columns[varscanT].split(':')[DP4].split(',')[1])
				tvfor = int(columns[varscanT].split(':')[DP4].split(',')[2])
				tvrev = int(columns[varscanT].split(':')[DP4].split(',')[3])
				tumor_read1 = trfor+trrev
				tumor_read2 = tvfor+tvrev
				tcov = trfor+trrev+tvfor+tvrev
				tfreq = columns[varscanT].split(':')[Freq]
				nrfor = int(columns[varscanN].split(':')[DP4].split(',')[0])
				nrrev = int(columns[varscanN].split(':')[DP4].split(',')[1])
				nvfor = int(columns[varscanN].split(':')[DP4].split(',')[2])
				nvrev = int(columns[varscanN].split(':')[DP4].split(',')[3])
				normal_read1 = nrfor+nrrev
				normal_read2 = nvfor+nvrev
				ncov = nrfor+nrrev+nvfor+nvrev
				nfreq = columns[varscanN].split(':')[Freq]
				vcf_cov_dict[chrm + ':' + pos] = {}
				vcf_cov_dict[chrm + ':' + pos]['pval'] = p_val
				vcf_cov_dict[chrm + ':' + pos]['Note'] = str(caller_count) + ':' + callers
				vcf_cov_dict[chrm + ':' + pos]['trfor'] = trfor
				vcf_cov_dict[chrm + ':' + pos]['trrev'] = trrev
				vcf_cov_dict[chrm + ':' + pos]['tvfor'] = tvfor
				vcf_cov_dict[chrm + ':' + pos]['tvrev'] = tvrev
				vcf_cov_dict[chrm + ':' + pos]['nrfor'] = nrfor
				vcf_cov_dict[chrm + ':' + pos]['nrrev'] = nrrev
				vcf_cov_dict[chrm + ':' + pos]['nvfor'] = nvfor
				vcf_cov_dict[chrm + ':' + pos]['nvrev'] = nvrev
				vcf_cov_dict[chrm + ':' + pos]['tumor_read1'] = tumor_read1
				vcf_cov_dict[chrm + ':' + pos]['tumor_read2'] = tumor_read2
				vcf_cov_dict[chrm + ':' + pos]['normal_read1'] = normal_read1
				vcf_cov_dict[chrm + ':' + pos]['normal_read2'] = normal_read2
				vcf_cov_dict[chrm + ':' + pos]['tumor_coverage'] = tcov
				vcf_cov_dict[chrm + ':' + pos]['tumor_freq'] = tfreq
				vcf_cov_dict[chrm + ':' + pos]['normal_coverage'] = ncov
				vcf_cov_dict[chrm + ':' + pos]['normal_freq'] = nfreq
			elif re.search('strelka', callers):
				caller_count = callers.count('-')+1
				form = columns[8].split(':')
				DP = form.index('DP')
				TIR = form.index('TIR')
				trfor = '.'
				trrev = '.'
				tvfor = '.'
				tvrev = '.'
				nrfor = '.'
				nrrev = '.'
				nvfor = '.'
				nvrev = '.'
				p_val = '.'
				normal_read2 = int(line.split('\t')[strelkaN].split(':')[TIR].split(',')[1])
				tumor_read2 = int(line.split('\t')[strelkaT].split(':')[TIR].split(',')[1])
				n_cov = int(line.split('\t')[strelkaN].split(':')[DP])
				t_cov = int(line.split('\t')[strelkaT].split(':')[DP])
				T_freq = float((tumor_read2/t_cov)*100)
				if normal_read2 != 0:
					N_freq = float(normal_read2/n_cov)
				else:
					N_freq = 0
				tumor_read1 = t_cov-tumor_read2
				normal_read1 = n_cov-normal_read2
				vcf_cov_dict[chrm + ':' + pos] = {}
				vcf_cov_dict[chrm + ':' + pos]['Note'] = str(caller_count) + ':' + callers
				vcf_cov_dict[chrm +':' + pos]['pval'] = '.'
				vcf_cov_dict[chrm + ':' + pos]['trfor'] = trfor
				vcf_cov_dict[chrm + ':' + pos]['trrev'] = trrev
				vcf_cov_dict[chrm + ':' + pos]['tvfor'] = tvfor
				vcf_cov_dict[chrm + ':' + pos]['tvrev'] = tvrev
				vcf_cov_dict[chrm + ':' + pos]['nrfor'] = nrfor
				vcf_cov_dict[chrm + ':' + pos]['nrrev'] = nrrev
				vcf_cov_dict[chrm + ':' + pos]['nvfor'] = nvfor
				vcf_cov_dict[chrm + ':' + pos]['nvrev'] = nvrev
				vcf_cov_dict[chrm + ':' + pos]['tumor_read1'] = tumor_read1
				vcf_cov_dict[chrm + ':' + pos]['tumor_read2'] = tumor_read2
				vcf_cov_dict[chrm + ':' + pos]['normal_read1'] = normal_read1
				vcf_cov_dict[chrm + ':' + pos]['normal_read2'] = normal_read2
				vcf_cov_dict[chrm + ':' + pos]['tumor_coverage'] = tcov
				vcf_cov_dict[chrm + ':' + pos]['tumor_freq'] = tfreq
				vcf_cov_dict[chrm + ':' + pos]['normal_coverage'] = ncov
				vcf_cov_dict[chrm + ':' + pos]['normal_freq'] = nfreq
	vcf.close()

	############################################
	##Now that I have a dict I need to produce##
	##a final sheet with coverage that Will go##
	## into our database                      ##
	############################################
	nonsyn_snv = open(INDEL_ANNO)
	nonsyn_file = open('nonsyn_SQL_insert.txt', 'a')
	nonsyn_file.close()
	date_trunc = datetime.datetime.now().date()
	first_line = True
	for line in nonsyn_snv:
		if first_line:
			first_line= False
			continue
		columns = line.strip().split('\t')
		if  (re.search(r'nonsynonymous', columns[8])or re.search(r'frame', columns[8])or re.search(r'stop', columns[8])\
			or re.search(r'nonsynonymous', columns[13])or re.search(r'frame', columns[13]) or re.search(r'stop', columns[13])\
			or re.search(r'nonsynonymous', columns[18]) or re.search(r'frame', columns[18])or re.search(r'stop', columns[18])):
			Chr = columns[0]
			start = columns[1]
			end = columns[2]
			ref = columns[3]
			alt = columns[4]
			if re.search('-', alt):
				ID = Chr + ':' + str(int(start)-1)
			else:
				ID = Chr + ':' + start
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
			COSMIC = columns[26]
			mrn = ""
			seq_center = "VHIO"
			sample_gDNA = sampleID + ' chr' + Chr + ':' + start
			gDNA = 'chr' + Chr + ':' + start
			tumor_type = tumor_type
			source = ""
			sample_note = ""
			sample_center = sample_note + ' ' + source + ' ' + seq_center
			variant_key = str(Chr) + ':' + str(start) + '-' + str(end) + ' ' + str(ref) + '>' + str(alt)
			if ref_gene_detail != 'NA':
				columns[9] = columns[7]
			if known_gene_detail != 'NA':
				columns[14] = columns[12]
			if ens_gene_detail != 'NA':
				columns[19] = columns[17]
			p_val = vcf_cov_dict[ID]['pval']
			Note = vcf_cov_dict[ID]['Note']
			trfor = vcf_cov_dict[ID]['trfor']
			trrev = vcf_cov_dict[ID]['trrev']
			tvfor = vcf_cov_dict[ID]['tvfor']
			tvrev = vcf_cov_dict[ID]['tvrev']
			nrfor = vcf_cov_dict[ID]['nrfor']
			nrrev = vcf_cov_dict[ID]['nrrev']
			nvfor = vcf_cov_dict[ID]['nvfor']
			nvrev = vcf_cov_dict[ID]['nvrev']
			tumor_read1 = vcf_cov_dict[ID]['tumor_read1']
			tumor_read2 = vcf_cov_dict[ID]['tumor_read2']
			normal_read1 = vcf_cov_dict[ID]['normal_read1']
			normal_read2 = vcf_cov_dict[ID]['normal_read2']
			tcov = vcf_cov_dict[ID]['tumor_coverage']
			tfreq = vcf_cov_dict[ID]['tumor_freq']
			ncov = vcf_cov_dict[ID]['normal_coverage']
			nfreq = vcf_cov_dict[ID]['normal_freq']
			nonsyn_file = open('nonsyn_SQL_insert.txt', 'a')
			nonsyn_file.write(str(mrn)+"\t"+str(seq_center)+"\t"+str(sampleID)+"\t"+str('\t'.join(columns[0:5]))+'\t-\t-\t'+str(columns[20])+'\t'+str(tumor_read1)+'\t'+str(tumor_read2)+'\t'\
							+str('\t'.join(columns[5:7]))+'\t'+str('\t'.join(columns[8:12]))+'\t'+str('\t'.join(columns[13:17]))+'\t'+str('\t'.join(columns[18:20]))+'\t'\
							+str('\t'.join(columns[21:26]))+'\t'+str(normal_read1)+'\t'+str(normal_read2)+'\t'+str(trfor)+'\t'+str(trrev)+'\t'+str(tvfor)+'\t'+str(tvrev)+'\t'\
							+str(nrfor)+'\t'+str(nrrev)+'\t'+str(nvfor)+'\t'+str(nfreq)+'\t'+str(nvrev)+'\t'+str(tfreq)+'\tSomatic\t'+str(p_val)+'\t'+str(sampleID)+' chr'+str(Chr)+':'+str(start)+'\t'+str(Note)+'\t'\
							+str(gDNA)+'\t'+str(tumor_type)+'\t'+str(source)+'\t'+str(sample_note)+' '+str(source)+' '+str(seq_center)+'\t'+str(tcov)+'\t'+str(ncov)+'\t'+str(sample_note)+'\t'+str(variant_key)+'\t'+str(COSMIC)+'\t'+str(date)+'\t'+resection_date+'\t'+run_date+'\t'+sequencer+'\t'+kit+'\t'+note+'\t'+index+'\n')
			nonsyn_file.close()
		else:
			Chr = columns[0]
			start = columns[1]
			end = columns[2]
			ref = columns[3]
			alt = columns[4]
			if re.search('-', alt):
				ID = Chr + ':' + str(int(start)-1)
			else:
				ID = Chr + ':' + start
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
			COSMIC = columns[26]
			mrn = " "
			seq_center = "VHIO"
			sample_gDNA = sampleID + ' chr' + Chr + ':' + start
			gDNA = 'chr' + Chr + ':' + start
			tumor_type = tumor_type
			source = " "
			sample_note = " "
			sample_center = sample_note + ' ' + source + ' ' + seq_center
			variant_key = str(Chr) + ':' + str(start) + '-' + str(end) + ' ' + str(ref) + '>' + str(alt)
			if ref_gene_detail != 'NA':
				columns[9] = columns[7]
			if known_gene_detail != 'NA':
				columns[14] = columns[12]
			if ens_gene_detail != 'NA':
				columns[19] = columns[17]
			p_val = vcf_cov_dict[ID]['pval']
			Note = vcf_cov_dict[ID]['Note']
			trfor = vcf_cov_dict[ID]['trfor']
			trrev = vcf_cov_dict[ID]['trrev']
			tvfor = vcf_cov_dict[ID]['tvfor']
			tvrev = vcf_cov_dict[ID]['tvrev']
			nrfor = vcf_cov_dict[ID]['nrfor']
			nrrev = vcf_cov_dict[ID]['nrrev']
			nvfor = vcf_cov_dict[ID]['nvfor']
			nvrev = vcf_cov_dict[ID]['nvrev']
			tumor_read1 = vcf_cov_dict[ID]['tumor_read1']
			tumor_read2 = vcf_cov_dict[ID]['tumor_read2']
			normal_read1 = vcf_cov_dict[ID]['normal_read1']
			normal_read2 = vcf_cov_dict[ID]['normal_read2']
			tcov = vcf_cov_dict[ID]['tumor_coverage']
			tfreq = vcf_cov_dict[ID]['tumor_freq']
			ncov = vcf_cov_dict[ID]['normal_coverage']
			nfreq = vcf_cov_dict[ID]['normal_freq']
			all_file = open('all_other_mutations.txt', 'a')
			all_file.write(str(mrn)+"\t"+str(seq_center)+"\t"+str(sampleID)+"\t"+str('\t'.join(columns[0:5]))+'\t-\t-\t'+str(columns[20])+'\t'+str(tumor_read1)+'\t'+str(tumor_read2)+'\t'\
							+str('\t'.join(columns[5:7]))+'\t'+str('\t'.join(columns[8:12]))+'\t'+str('\t'.join(columns[13:17]))+'\t'+str('\t'.join(columns[18:20]))+'\t'\
							+str('\t'.join(columns[21:26]))+'\t'+str(normal_read1)+'\t'+str(normal_read2)+'\t'+str(trfor)+'\t'+str(trrev)+'\t'+str(tvfor)+'\t'+str(tvrev)+'\t'\
							+str(nrfor)+'\t'+str(nrrev)+'\t'+str(nvfor)+'\t'+str(nfreq)+'\t'+str(nvrev)+'\t'+str(tfreq)+'\tSomatic\t'+str(p_val)+'\t'+str(sampleID)+' chr'+str(Chr)+':'+str(start)+'\t'+str(Note)+'\t'\
							+str(gDNA)+'\t'+str(tumor_type)+'\t'+str(source)+'\t'+str(sample_note)+' '+str(source)+' '+str(seq_center)+'\t'+str(tcov)+'\t'+str(ncov)+'\t'+str(sample_note)+'\t'+str(variant_key)+'\t'+str(COSMIC)+'\t'+str(date)+'\t'+resection_date+'\t'+run_date+'\t'+sequencer+'\t'+kit+'\t'+note+'\t'+index+'\n')
			all_file.close()
	nonsyn_snv.close()


	# Extract peptides
	print("Extracting pepdides")
	sample_sheet = open('nonsyn_SQL_insert.txt')
	epitope_file = open('Formatted_epitope_variant.txt','a')
	for line in sample_sheet:
		columns = line.strip().split('\t')
		mrn = columns[0]
		seq_center = columns[1]
		sampleID = columns[2]
		source = columns[48]
		tumor_type = columns[47]
		sample_note = columns[52]
		sample_gDNA = columns[44]
		gDNA = columns[46]
		sample_center = columns[49]
		variant_key = columns[53]
		chrom = columns[3]
		start = columns[4]
		stop = columns[5]
		ref = columns[6]
		alt = columns[7]
		func_ref_gene = columns[13]
		gene_ref_gene = columns[14]
		exonic_func_ref = columns[15]
		AA_change_refGene = columns[16].split(',')
		func_UCSC_gene = columns[17]
		gene_UCSC_gene = columns[18]
		exonic_func_UCSC = columns[19]
		AA_change_UCSCGene = columns[20].split(',')
		func_ens_gene = columns[21]
		gene_ens_gene = columns[22]
		exonic_func_ens = columns[23]
		AA_change_ensGene = columns[24].split(',')
		if re.search(r'nonsynonymous', exonic_func_ref) or re.search(r'frame', exonic_func_ref):
			for entry in AA_change_refGene:
				epitope_file.write(mrn+'\t'+seq_center+'\t'+sampleID +'\t'+source+'\t'+tumor_type+'\t'+sample_note+'\t'+\
								   sample_gDNA+'\t'+gDNA+'\t'+sample_center+'\t'+variant_key+'\t'+chrom+'\t'+start+'\t'+\
								stop+'\t'+ref+'\t'+alt+'\t'+func_ref_gene+'\t'+exonic_func_ref+'\t'+\
								   str(sub(':','\t',(entry)))+'\n')
		if re.search(r'nonsynonymous', exonic_func_UCSC) or re.search(r'frame', exonic_func_UCSC):
			for entry in AA_change_UCSCGene:
				epitope_file.write(mrn+'\t'+seq_center+'\t'+sampleID +'\t'+source+'\t'+tumor_type+'\t'+sample_note+'\t'+\
								   sample_gDNA+'\t'+gDNA+'\t'+sample_center+'\t'+variant_key+'\t'+chrom+'\t'+start+'\t'+\
								   stop+'\t'+ref+'\t'+alt+'\t'+func_UCSC_gene+'\t'+exonic_func_UCSC+'\t'+\
								   str(sub(':','\t',(entry)))+'\n')
		if re.search(r'nonsynonymous', exonic_func_ens) or re.search(r'frame', exonic_func_ens):
			for entry in AA_change_ensGene:
				epitope_file.write(mrn+'\t'+seq_center+'\t'+sampleID +'\t'+source+'\t'+tumor_type+'\t'+sample_note+'\t'+\
								   sample_gDNA+'\t'+gDNA+'\t'+sample_center+'\t'+variant_key+'\t'+chrom+'\t'+start+'\t'+\
								   stop+'\t'+ref+'\t'+alt+'\t'+func_ens_gene+'\t'+exonic_func_ens+'\t'+\
								   str(sub(':','\t',(entry)))+'\n')
	epitope_file.write('\n')
	epitope_file.close()
	sample_sheet.close()
	print('Formatted file created')

	################################################################
	## Now I need to run through this list of transcripts and mut ##
	## and create the epitope files.                              ##
	################################################################
	################################################################
	## create two dictionaries of transcript IDS and their AA seq ##
	## or thier cDNA seq                                          ##
	################################################################
	print('Creating epitopes')
	AA_seq = {}
	dict1 = open(FASTA_AA_DICT)
	for line in dict1:
		entry = line.rstrip("\n").split(":")
		key, values = entry[0], entry[1]
		AA_seq[key] = values
	cDNA_seq = {}
	dict2 = open(FASTA_cDNA)
	for line in dict2:
		entry = line.rstrip("\n").split(":")
		key, values = entry[0], entry[1]
		cDNA_seq[key] = values
	epitope_file = open('SQL_Epitopes.txt', 'a')
	input_file = open('Formatted_epitope_variant.txt')
	for line in input_file:
		columns = line.rstrip('\n').split('\t')
		if len(columns) > 21:
			#############################
			## Load all pertinent info ##
			#############################
			mrn = columns[0].strip()
			seq_center = columns[1].strip()
			sampleID = columns[2].strip()
			source = columns[3].strip()
			tumor_type = columns[4].strip()
			sample_note = columns[5].strip()
			sample_gDNA = columns[6].strip()
			gDNA = columns[7].strip()
			sample_center = columns[8].strip()
			variant_key = columns[9].strip()
			chrom = columns[10].strip()
			spot = columns[11].strip()
			stop = columns[12].strip()
			ref = columns[13].strip()
			alt = columns[14].strip()
			func_gene = columns[15].strip()
			exonic_func = columns[16].strip()
			gene_name = columns[17].strip()
			transcriptID = columns[18].strip()
			exon = columns[19].strip()
			cDNA_raw = columns[20].strip()
			errors = ''
			protein_raw = columns[21].replace(' ','')
			################################################
			## Nonsynonymous point mutations to 25 mers   ##
			################################################
			if exonic_func == 'nonsynonymous SNV' and re.search(r'^p\.', protein_raw):
				protein_strip = protein_raw.strip()
				##################################
				## extract the AA change info   ##
				##################################
				position = int(protein_strip[(protein_strip.find('.')+2):len(protein_strip)-1])
				ref_AA = protein_strip[(protein_strip.find('.')+1)]
				var_AA = protein_strip[len(protein_strip)-1]
				##################################
				## gather AA seq for transcript ##
				##################################
				protein_seq = AA_seq.get(transcriptID, 'AA_seq not present for this transcript')
				################################################
				## I need a quick dirty way to check and make ##
				## sure that the annotation is correct.  The  ##
				## simplest way for me to do that is to make  ##
				## sure that the AA is what it is supposed to ##
				## be at that position.                       ##
				################################################
				FASTA_AA = protein_seq[position-1:position]
				if FASTA_AA == ref_AA:
					#######################################################
					## Create wt 25 mer by extracting the  positions     ##
					## depending on certian conditions i.e. mutations    ##
					## that fall close to the begining or end of the     ##
					## Amino acid.                                       ##
					##                                                   ##
					## Create the nonsynonymous epitope by extracting 12 ##
					## Amino acids before the mutation, the mutation     ##
					## and then 12 amino acids after the mutation        ##
					#######################################################
					if position >= 13:
						WT_25mer = protein_seq[(position-13):position+12]
						Mut_25mer = protein_seq[(position-13):position-1] + var_AA + protein_seq[(position):(position+12)]
					elif position < 13:
						WT_25mer = protein_seq[0 : position+12]
						Mut_25mer = protein_seq[0:position-1] + var_AA + protein_seq[(position):(position+12)]
				#####################################################
				## in event the AA doesn't match the reference at  ##
				## that position I don't want to make a false      ##
				## epitope and I want to flag the transcript so I  ##
				## can try to figure out what happened             ##
				#####################################################
				elif FASTA_AA != ref_AA and errors.startswith('AA_seq not'):
					WT_25mer = ''
					Mut_25mer = ''
				elif FASTA_AA != ref_AA and not errors.startswith('AA_seq not'):
					WT_25mer = ''
					Mut_25mer = ''
				#########################################################
				## Open and write the epitopes to a final epitope file ##
				#########################################################
				epitope_file=open('SQL_Epitopes.txt', 'a')
				epitope_file.write(str('\t'.join(columns[0:]))+'\t'+str(position)+'\t'+ errors+'\t'+WT_25mer+'\t'+Mut_25mer+'\t'+str(variant_key)+' '+transcriptID+'\n')
				epitope_file.close()
			#############################################################################################
			## The above section handled all of the nonsynonymous point mutations.  The remanider will ##
			## tackle the problem of insertion and deletion mutations.  The strategy taken is to       ##
			## take the cDNA for the appropriate transcript, make the noted change, then translate the ##
			## sequence, after that I need to either 1) extract 12 AA before the indel and everything  ##
			## after until 1st stop codon or end of transcript or 2) for non frameshifts just take the ##
			## 12 amino acids before or after the indel no matter the size of the indel.               ##
			#############################################################################################
			##############################
			## 1st frameshift deletions ##
			##############################
			elif exonic_func == 'frameshift deletion' and re.search(r'^p\.', protein_raw):
				protein_strip = protein_raw.strip()
				#########################################################
				##Extract the cDNA seq to the corresponding transcript ##
				#########################################################
				ref_cDNA_seq = cDNA_seq.get(transcriptID, 'cDNA not present for this transcript')
				##############################
				##count of bases deleted so ##
				##I can make the change in  ##
				##the cDNA.                 ##
				##############################
				len_del = len(ref)
				cDNA_strip = cDNA_raw.strip()
				###############################
				## Single base deletions are ##
				## annotated differently so  ##
				## i need to extract the     ##
				## position from the cDNA    ##
				## change column differently ##
				###############################
				if len_del >1 and protein_strip.startswith('p.') and not protein_strip.startswith('p.X'):
					cDNA_pos = cDNA_strip[(cDNA_strip.find('.')+1):cDNA_strip.find('_')]
					position = int(protein_strip[(protein_strip.find('.')+2):protein_strip.find('f')])
				elif len_del == 1 and protein_strip.startswith('p.')and not protein_strip.startswith('p.X'):
					cDNA_pos = cDNA_strip[(cDNA_strip.find('.')+1):cDNA_strip.find('del')]
					position = int(protein_strip[(protein_strip.find('.')+2):protein_strip.find('f')])
				if not protein_strip.startswith('p.'):
					position = 0
				elif protein_strip.startswith('p.X'):
					position = 0
				#######################################################
				##Uses that position info obtained above to extract  ##
				##cDNA regions to the right and left and stitch them ##
				##together to create the new mutant cDNA             ##
				#######################################################
				mut_cDNA_left = ref_cDNA_seq[0 : int(cDNA_pos)-1]
				mut_cDNA_right = ref_cDNA_seq[int(cDNA_pos)+int(len_del)-1:]
				mut_cDNA_seq = mut_cDNA_left+mut_cDNA_right
				###############################################################
				##Translate the ref and mut_cDNA and then extract the regions##
				###############################################################
				ref_FASTA = translate_dna.translate_dna(ref_cDNA_seq)
				mut_FASTA = translate_dna.translate_dna(mut_cDNA_seq)
				mut_stop = int(mut_FASTA.find('X'))
				##########################################################
				##Extract the regions depending upon where the position ##
				## is same as was done with the nonsynonymous mutations ##
				##########################################################
				if position >= 13  and mut_stop > 0:
					WT_25mer_temp = ref_FASTA[(position-13):position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer = mut_FASTA[(position-13):mut_stop]
				elif position < 13 and position > 0 and mut_stop > 0:
					WT_25mer_temp = ref_FASTA[0 : position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer = mut_FASTA[0:mut_stop]
				elif position >= 13  and mut_stop < 0:
					WT_25mer_temp = ref_FASTA[(position-13):position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer = mut_FASTA[(position-13):]
				elif position < 13 and position > 0 and mut_stop < 0:
					WT_25mer_temp = ref_FASTA[0 : position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer = mut_FASTA [0:]
				elif position == 0 or ref_cDNA_seq.startswith('cDNA'):
					errors = errors+' can not code for this mutated AA_position '
					WT_25mer = ''
					Mut_25mer = ''
				epitope_file = open('SQL_Epitopes.txt', 'a')
				epitope_file.write(str('\t'.join(columns[0:]))+'\t'+str(position)+'\t'+ errors+'\t'+WT_25mer+'\t'+Mut_25mer+'\t'+str(variant_key)+' '+transcriptID+'\n')
				epitope_file.close()
			###############################
			## 2nd frameshift insertions ##
			###############################
			elif exonic_func == 'frameshift insertion' and re.search(r'^p\.', protein_raw):
				protein_strip = protein_raw.strip()
				#########################################################
				##Extract the cDNA seq to the corresponding transcript ##
				#########################################################
				ref_cDNA_seq = cDNA_seq.get(transcriptID, 'cDNA not present for this transcript')
				#################################
				##count of bases inserteded so ##
				##I can make the change in     ##
				##the cDNA.                    ##
				#################################
				len_ins = len(alt)
				cDNA_strip = cDNA_raw.strip()
				#####################################
				## frameshift insertions are either##
				## annotated as duplications or    ##
				## large insertions of a different ##
				## base and this requires them to  ##
				## be stripped differently to get  ##
				## to the position of the insertion##
				#####################################
				if  re.search(r'dup', cDNA_strip):
					cDNA_pos = cDNA_strip[int(cDNA_strip.find('.'))+1:int(cDNA_strip.find('dup'))]
					ins = cDNA_strip[int(cDNA_strip.find('dup'))+3:]
				elif re.search(r'_',cDNA_strip):
					cDNA_pos = cDNA_strip[int(cDNA_strip.find('.'))+1:int(cDNA_strip.find('_'))]
					ins = cDNA_strip[int(cDNA_strip.find('ins'))+3:]
				######################################
				## after I have the insertion point ##
				## I can peice together the cDNA    ##
				## like done for the deletion above ##
				######################################
				mut_cDNA_left = ref_cDNA_seq[0 : int(cDNA_pos)]
				mut_cDNA_right = ref_cDNA_seq[int(cDNA_pos):]
				mut_cDNA_seq = mut_cDNA_left+ins+mut_cDNA_right
				##################################
				## need to get the AA position  ##
				##################################
				if protein_strip.startswith('p.') and re.search(r'fs', protein_strip) and not protein_strip.startswith('p.X'):
					position = int(protein_strip[(protein_strip.find('.')+2):protein_strip.find('f')])
				elif protein_strip.startswith('p.') and re.search(r'delins', protein_strip) and not protein_strip.startswith('p.X'): ##different versions of annovar annotated this differently this has corrected it thus far.
					position = int(protein_strip[(protein_strip.find('.')+2):protein_strip.find('delins')])
				elif protein_strip.startswith('p.X'):
					position = 0
				elif not protein_strip.startswith('p.'):
					postion = 0
				###############################################
				##Translate new cDNA and extract new epitopes##
				###############################################
				ref_FASTA = translate_dna.translate_dna(ref_cDNA_seq)
				mut_FASTA = translate_dna.translate_dna(mut_cDNA_seq)
				mut_stop = int(mut_FASTA.find('X'))
				if position >= 13 and mut_stop > 0:
					WT_25mer_temp = ref_FASTA[(position-13):position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer = mut_FASTA[(position-13):mut_stop]
				elif position < 13 and position > 0 and mut_stop > 0:
					WT_25mer_temp = ref_FASTA[0 : position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer = mut_FASTA [0:mut_stop]
				elif position >= 13  and mut_stop < 0:
					WT_25mer_temp = ref_FASTA[(position-13):position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer = mut_FASTA[(position-13):]
				elif position < 13 and position > 0 and mut_stop < 0:
					WT_25mer_temp = ref_FASTA[0 : position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer = mut_FASTA [0:]
				elif position == 0 or ref_cDNA_seq.startswith('cDNA'):
					WT_25mer = ''
					Mut_25mer = ''
				epitope_file = open('SQL_Epitopes.txt', 'a')
				epitope_file.write(str('\t'.join(columns[0:]))+'\t'+str(position)+'\t'+ errors+'\t'+WT_25mer+'\t'+Mut_25mer+'\t'+str(variant_key)+' '+transcriptID+'\n')
				epitope_file.close()
			##########################################
			##3rd nonframeshift deletions to 25mers ##
			##########################################
			elif exonic_func == 'nonframeshift deletion' and re.search(r'^p\.', protein_raw):
				protein_strip = protein_raw.strip()
				#########################################################
				##Extract the cDNA seq to the corresponding transcript ##
				#########################################################
				ref_cDNA_seq = cDNA_seq.get(transcriptID, 'cDNA not present for this transcript')
				##############################
				##count of bases deleted so ##
				##I can make the change in  ##
				##the cDNA.                 ##
				##############################
				len_del = len(ref)
				cDNA_strip = cDNA_raw.strip()
				###############################
				## Single base deletions are ##
				## annotated differently so  ##
				## i need to extract the     ##
				## position from the cDNA    ##
				## change column differently ##
				###############################
				if re.search(r'_', cDNA_strip):
					cDNA_pos = cDNA_strip[(cDNA_strip.find('.')+1):cDNA_strip.find('_')]
				else:
					cDNA_pos = cDNA_strip[(cDNA_strip.find('.')+1):cDNA_strip.find('del')]
				if protein_strip.startswith('p.') and not protein_strip.startswith('p.X'):
					position = int(protein_strip[(protein_strip.find('.')+1):protein_strip.find('_')])
				elif not protein_strip.startswith('p.'):
					position = 0
				elif protein_strip.startswith('p.X'):
					position = 0
				#######################################################
				##Uses that position info obtained above to extract  ##
				##cDNA regions to the right and left and stitch them ##
				##together to create the new mutant cDNA             ##
				#######################################################
				mut_cDNA_left = ref_cDNA_seq[0 : int(cDNA_pos)-1]
				mut_cDNA_right = ref_cDNA_seq[int(cDNA_pos)+int(len_del)-1:]
				mut_cDNA_seq = mut_cDNA_left+mut_cDNA_right
				###############################################################
				##Translate the ref and mut_cDNA and then extract the regions##
				###############################################################
				ref_FASTA = translate_dna.translate_dna(ref_cDNA_seq)
				mut_FASTA = translate_dna.translate_dna(mut_cDNA_seq)
				mut_stop = int(mut_FASTA.find('X'))
				if position >= 13 :
					WT_25mer_temp = ref_FASTA[(position-13):position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer_temp = mut_FASTA[(position-13):position+12]
					Mut_25mer = Mut_25mer_temp.replace('X','')
				elif position < 13 and position > 0:
					WT_25mer_temp = ref_FASTA[0 : position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer_temp = mut_FASTA [0 : position+12]
					Mut_25mer = Mut_25mer_temp.replace('X','')
				elif position == 0 or ref_cDNA_seq.startswith('cDNA'):
					WT_25mer = ''
					Mut_25mer = ''
				epitope_file=open('SQL_Epitopes.txt', 'a')
				epitope_file.write(str('\t'.join(columns[0:]))+'\t'+str(position)+'\t'+ errors+'\t'+WT_25mer+'\t'+Mut_25mer+'\t'+str(variant_key)+' '+transcriptID+'\n')
				epitope_file.close()
			###########################################
			##4th nonframeshift insertions to 25mers ##
			###########################################
			elif exonic_func == 'nonframeshift insertion' and re.search(r'^p\.', protein_raw):
				protein_strip = protein_raw.strip()
				#########################################################
				##Extract the cDNA seq to the corresponding transcript ##
				#########################################################
				ref_cDNA_seq = cDNA_seq.get(transcriptID, 'cDNA not present for this transcript')
				#################################
				##count of bases inserteded so ##
				##I can make the change in     ##
				##the cDNA.                    ##
				#################################
				len_ins = len(alt)
				cDNA_strip = cDNA_raw.strip()
				cDNA_pos = cDNA_strip[int(cDNA_strip.find('.'))+1:int(cDNA_strip.find('_'))]
				ins = cDNA_strip[int(cDNA_strip.find('ins'))+3:]
				#######################################################
				##Uses that position info obtained above to extract  ##
				##cDNA regions to the right and left and stitch them ##
				##together to create the new mutant cDNA             ##
				#######################################################
				mut_cDNA_left = ref_cDNA_seq[0 : int(cDNA_pos)]
				mut_cDNA_right = ref_cDNA_seq[int(cDNA_pos):]
				mut_cDNA_seq = mut_cDNA_left+ins+mut_cDNA_right
				if protein_strip.startswith('p.') and not protein_strip.startswith('p.X'):
					position = int(protein_strip[(protein_strip.find('.')+2):protein_strip.find('del')])
				elif not protein_strip.startswith('p.'):
					position = 0
				elif protein_strip.startswith('p.X'):
					position = 0
				###############################################################
				##Translate the ref and mut_cDNA and then extract the regions##
				###############################################################
				ref_FASTA = translate_dna.translate_dna(ref_cDNA_seq)
				mut_FASTA = translate_dna.translate_dna(mut_cDNA_seq)
				mut_stop = int(mut_FASTA.find('X'))
				if position >= 13 :
					WT_25mer_temp = ref_FASTA[(position-13):position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer_temp = mut_FASTA[(position-13):position+12]
					Mut_25mer = Mut_25mer_temp.replace('X','')
				elif position < 13 and position > 0 :
					WT_25mer_temp = ref_FASTA[0 : position+12]
					WT_25mer = WT_25mer_temp.replace('X','')
					Mut_25mer_temp = mut_FASTA [0 : position+12]
					Mut_25mer = Mut_25mer_temp.replace('X','')
				elif position == 0 or ref_cDNA_seq.startswith('cDNA'):
					WT_25mer = 'can not code for this mutated AA postion'
					Mut_25mer = ''
				epitope_file = open('SQL_Epitopes.txt', 'a')
				epitope_file.write(str('\t'.join(columns[0:]))+'\t'+str(position)+'\t'+ errors+'\t'+WT_25mer+'\t'+Mut_25mer+'\t'+str(variant_key)+' '+transcriptID+'\n')
				epitope_file.close()
	print('Epitopes have been created...')
	input_file.close()
