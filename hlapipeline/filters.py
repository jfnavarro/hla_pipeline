"""
@author: jfnavarro
"""
import gzip
import re
from hlapipeline.common import index_column_substring

def mutect2_filter(input, output, sample1_ID, sample2_ID):
    filtered_vcf = open(output, 'w')
    vcf = open(input)
    for line in vcf:
        if line.startswith('#') and not line.startswith('#CHROM'):
            filtered_vcf.write(line)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            Tmut = headers.index(sample1_ID)
            Nmut = headers.index(sample2_ID)
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            Filter = columns[6]
            if re.search('PASS', Filter):
                # AD variable
                n_split = columns[Nmut].split(':')[1].split(',')
                t_split = columns[Tmut].split(':')[1].split(',')
                normal_coverage = int(n_split[0]) + int(n_split[1])
                tumor_coverage = int(t_split[0]) + int(t_split[1])
                tumor_var_depth = int(t_split[1])
                # AF variable
                tumor_var_freq = float(float(columns[Tmut].split(':')[2]) * 100)
                normal_var_freq = float(float(columns[Nmut].split(':')[2]) * 100)
                if normal_var_freq != 0:
                    t2n_ratio = tumor_var_freq / normal_var_freq
                else:
                    t2n_ratio = 5
                # NOTE this filter seems to too strict with Mutect2 (where variants are already filtered)
                if normal_coverage >= 10 and tumor_coverage >= 10 and tumor_var_freq >= 2.5 and tumor_var_depth >= 3 and t2n_ratio >= 5:
                    filtered_vcf.write(line)
    filtered_vcf.close()
    vcf.close()

def strelka2_filter(input, output):
    alt_index = {'A': 4, 'C': 5, 'G': 6, 'T': 7}
    filtered_vcf = open(output, 'w')
    vcf = gzip.open(input, 'rt')
    for line in vcf:
        if line.startswith('#') and not '##FORMAT=<ID=DP,' in line and not line.startswith('#CHROM'):
            filtered_vcf.write(line)
        elif line.startswith('#') and '##FORMAT=<ID=DP,' in line:
            filtered_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' + line)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            Tst = headers.index('TUMOR')
            Nst = headers.index('NORMAL')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            ref = columns[3]
            alt = columns[4]
            Filter = columns[6]
            INFO = columns[7]
            Format = columns[8]
            Normal = columns[Nst]
            Tumor = columns[Tst]
            if 'PASS' in Filter:
                n_split = Normal.split(':')
                t_split = Tumor.split(':')
                normal_variant_depth = int(n_split[alt_index[alt]].split(',')[0])
                tumor_variant_depth = int(t_split[alt_index[alt]].split(',')[0])
                n_cov = int(n_split[0])
                t_cov = int(t_split[0])
                T_freq = float((tumor_variant_depth / t_cov) * 100)
                if normal_variant_depth != 0:
                    N_freq = float(normal_variant_depth / n_cov)
                    t2n_ratio = T_freq / N_freq
                else:
                    t2n_ratio = 5
                if n_cov >= 10 and t_cov >= 10 and tumor_variant_depth >= 3 and T_freq >= 5 and t2n_ratio >= 5:
                    Format = 'GT:' + Format
                    INFO_split = INFO.split(';')
                    SGT = INFO_split[6].replace('SGT=', '').split('->')
                    Normal_GT = ''
                    Tumor_GT = ''
                    i = 1
                    for x in SGT[0]:
                        if x == ref and i == 1:
                            Normal_GT = Normal_GT + '0'
                            i = i + 1
                        elif x == ref and i == 2:
                            Normal_GT = Normal_GT + '/0'
                            i = i + 1
                        elif x != ref and i == 1:
                            Normal_GT = Normal_GT + '1'
                            i = i + 1
                        elif x != ref and i == 2:
                            Normal_GT = Normal_GT + '/1'
                            i = i + 1
                    i = 1
                    for x in SGT[1]:
                        if x == ref and i == 1:
                            Tumor_GT = Tumor_GT + '0'
                            i = i + 1
                        elif x == ref and i == 2:
                            Tumor_GT = Tumor_GT + '/0'
                            i = i + 1
                        elif x != ref and i == 1:
                            Tumor_GT = Tumor_GT + '1'
                            i = i + 1
                        elif x != ref and i == 2:
                            Tumor_GT = Tumor_GT + '/1'
                            i = i + 1
                    filtered_vcf.write('{}\t{}\t{}:{}\t{}:{}\n'.format('\t'.join(columns[0:8]), Format, Normal_GT, Normal, Tumor_GT, Tumor))
    vcf.close()
    filtered_vcf.close()

def somaticSniper_filter(input, output):
    filtered_vcf = open(output, 'w')
    vcf = open(input)
    for line in vcf:
        if line.startswith('#') and not line.startswith('#CHROM'):
            filtered_vcf.write(line)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            Tss = headers.index('TUMOR')
            Nss = headers.index('NORMAL')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            n_split = columns[Nss].split(':')
            t_split = columns[Tss].split(':')
            normal_coverage = int(n_split[2])
            tumor_coverage = int(t_split[2])
            somatic_status = int(t_split[11])
            variant_count = int(t_split[3].split(',')[2]) + int(t_split[3].split(',')[3])
            normal_variant_count = int(n_split[3].split(',')[2]) + int(n_split[3].split(',')[3])
            Tumor_variant_freq = float((variant_count / tumor_coverage) * 100)
            if normal_variant_count != 0:
                normal_freq = float((normal_variant_count / normal_coverage) * 100)
                t2nratio = float(Tumor_variant_freq / normal_freq)
            else:
                t2nratio = 5
            if normal_coverage >= 10 and tumor_coverage >= 10 and somatic_status == 2 \
                    and variant_count >= 3 and Tumor_variant_freq >= 5 and t2nratio >= 5:
                filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()

def strelka2_filter_indels(input, output):
    filtered_vcf = open(output, 'w')
    vcf = gzip.open(input, 'rt')
    for line in vcf:
        if line.startswith('#') and not '##FORMAT=<ID=DP,' in line and not line.startswith('#CHROM'):
            filtered_vcf.write(line)
        elif line.startswith('#') and '##FORMAT=<ID=DP,' in line:
            filtered_vcf.write(
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            Tst = headers.index('TUMOR')
            Nst = headers.index('NORMAL')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            Filter = columns[6]
            INFO = columns[7]
            Format = columns[8]
            Normal = columns[Nst]
            Tumor = columns[Tst]
            if 'PASS' in Filter:
                n_split = columns[Nst].split(':')
                t_split = columns[Tst].split(':')
                normal_variant_depth = int(n_split[3].split(',')[1])
                tumor_variant_depth = int(t_split[3].split(',')[1])
                n_cov = int(n_split[0])
                t_cov = int(t_split[0])
                T_freq = float((tumor_variant_depth / t_cov) * 100)
                if normal_variant_depth != 0:
                    N_freq = float(normal_variant_depth / n_cov)
                    t2n_ratio = T_freq / N_freq
                else:
                    t2n_ratio = 5
                if n_cov >= 10 and t_cov >= 10 and tumor_variant_depth >= 3 and T_freq >= 5 and t2n_ratio >= 5:
                    Format = 'GT:' + Format
                    INFO_split = INFO.split(';')
                    SGT_index = index_column_substring(INFO_split, 'SGT')
                    SGT = INFO_split[SGT_index].replace('SGT=', '').split('->')
                    Normal_GT = '0/0'
                    if SGT[1] == 'het':
                        Tumor_GT = '0/1'
                    else:
                        Tumor_GT = '1/1'
                    filtered_vcf.write('{}\t{}\t{}:{}\t{}:{}\n'.format('\t'.join(columns[0:8]), Format,
                                                                       Normal_GT, Normal, Tumor_GT, Tumor))
    vcf.close()
    filtered_vcf.close()

def varscan_filter(input, output):
    filtered_vcf = open(output, 'w')
    vcf = open(input)
    for line in vcf:
        if line.startswith('#') and 'DP4' in line:
            # Ugly hack so CombineVariants works
            new_DP4 = line.replace(
                'ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev"',
                'ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases"')
            filtered_vcf.write(new_DP4)
        elif line.startswith('#') and not line.startswith('#CHROM'):
            filtered_vcf.write(line)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            Tvs = headers.index('TUMOR')
            Nvs = headers.index('NORMAL')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            Filter = columns[6]
            INFO = columns[7]
            if 'SOMATIC' in INFO and 'PASS' in Filter:
                n_split = columns[Nvs].split(':')
                t_split = columns[Tvs].split(':')
                normal_coverage = int(n_split[2])
                tumor_coverage = int(t_split[2])
                tumor_var_depth = int(t_split[4])
                tumor_var_freq = float(t_split[5].replace('%', ''))
                normal_var_freq = float(n_split[5].replace('%', ''))
                if normal_var_freq != 0:
                    t2n_ratio = tumor_var_freq / normal_var_freq
                else:
                    t2n_ratio = 5
                if normal_coverage >= 10 and tumor_coverage >= 10 and tumor_var_depth >= 3\
                        and t2n_ratio >= 5 and tumor_var_freq >= 2.5:
                    filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()