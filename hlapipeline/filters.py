import gzip

def index_column_substring(your_list, substring):
    for i, s in enumerate(your_list):
        if substring in s:
            return i
    return -1

def mutect2_filter(input, output, sample1_ID, sample2_ID):
    """
    Filters a Mutect2 VCF to only keep variants that PASS
    :param input:
    :param output:
    :param sample1_ID:
    :param sample2_ID:
    :return:
    """
    filtered_vcf = open(output, 'w')
    vcf = open(input)
    for line in vcf:
        if line.startswith('#') and not line.startswith('#CHROM'):
            line.replace(sample1_ID, "TUMOR")
            line.replace(sample2_ID, "NORMAL")
            filtered_vcf.write(line)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if 'PASS' in columns[headers.index('FILTER')]:
                filtered_vcf.write(line)
        else:
            filtered_vcf.write(line)
    filtered_vcf.close()
    vcf.close()

def strelka2_filter(input, output):
    """
    Filters a Strelka2 SNVs VCF to add the genotype and keep variants that PASS
    :param input:
    :param output:
    :return:
    """
    filtered_vcf = open(output, 'w')
    vcf = gzip.open(input, 'rt')
    for line in vcf:
        if '##FORMAT=<ID=DP,' in line:
            filtered_vcf.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n' + line)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            ref = columns[headers.index('REF')]
            if 'PASS' not in columns[headers.index('FILTER')]:
                continue
            INFO = columns[headers.index('INFO')]
            Format = columns[headers.index('FORMAT')]
            Tumor = columns[headers.index('TUMOR')]
            Normal = columns[headers.index('NORMAL')]
            Format = 'GT:' + Format
            INFO_split = INFO.split(';')
            SGT = INFO_split[6].replace('SGT=', '').split('->')
            Normal_GT = ''
            Tumor_GT = ''
            i = 1
            for x in SGT[0]:
                if x == ref and i == 1:
                    Normal_GT += '0'
                    i = i + 1
                elif x == ref and i == 2:
                    Normal_GT += '/0'
                    i = i + 1
                elif x != ref and i == 1:
                    Normal_GT += '1'
                    i = i + 1
                elif x != ref and i == 2:
                    Normal_GT += '/1'
                    i = i + 1
            i = 1
            for x in SGT[1]:
                if x == ref and i == 1:
                    Tumor_GT += '0'
                    i = i + 1
                elif x == ref and i == 2:
                    Tumor_GT += '/0'
                    i = i + 1
                elif x != ref and i == 1:
                    Tumor_GT += '1'
                    i = i + 1
                elif x != ref and i == 2:
                    Tumor_GT += '/1'
                    i = i + 1
            filtered_vcf.write('{}\t{}\t{}:{}\t{}:{}\n'.format('\t'.join(columns[0:8]),
                                                               Format,
                                                               Normal_GT,
                                                               Normal,
                                                               Tumor_GT,
                                                               Tumor))
        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()

def somaticSniper_filter(input, output):
    """
    Filters a SomaticSniper VCF to keep only somatic variants
    :param input:
    :param output:
    :return:
    """
    filtered_vcf = open(output, 'w')
    vcf = open(input)
    for line in vcf:
        if line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            format = columns[headers.index('FORMAT')]
            tumor = columns[headers.index('TUMOR')]
            somatic_status_index = format.split(':').index('SS')
            somatic_status = int(columns[tumor].split(':')[somatic_status_index])
            if somatic_status == 2:
                filtered_vcf.write(line)
        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()

def strelka2_filter_indels(input, output):
    """
    Filters a Strelka2 indels VCF to add the genotype and keep variants that PASS
    :param input:
    :param output:
    :return:
    """
    filtered_vcf = open(output, 'w')
    vcf = gzip.open(input, 'rt')
    for line in vcf:
        if '##FORMAT=<ID=DP,' in line:
            filtered_vcf.write(
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            if 'PASS' not in columns[headers.index('FILTER')]:
                continue
            INFO = columns[headers.index('INFO')]
            Format = columns[headers.index('FORMAT')]
            Tumor = columns[headers.index('TUMOR')]
            Normal = columns[headers.index('NORMAL')]
            Format = 'GT:' + Format
            INFO_split = INFO.split(';')
            SGT_index = index_column_substring(INFO_split, 'SGT')
            if SGT_index != -1:
                SGT = INFO_split[SGT_index].replace('SGT=', '').split('->')
                Normal_GT = '0/0'
                Tumor_GT = '0/1' if SGT[1] == 'het' else '1/1'
                filtered_vcf.write('{}\t{}\t{}:{}\t{}:{}\n'.format('\t'.join(columns[0:8]),
                                                                   Format,
                                                                   Normal_GT,
                                                                   Normal,
                                                                   Tumor_GT,
                                                                   Tumor))
        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()

def varscan_filter(input, output):
    """
    Filters a Varscan VCF reformat the DP4 field and keep only SOMATIC variants that PASS
    :param input:
    :param output:
    :return:
    """
    filtered_vcf = open(output, 'w')
    vcf = open(input)
    for line in vcf:
        if line.startswith('#') and 'DP4' in line:
            new_DP4 = line.replace(
                'ID=DP4,Number=1,Type=String,Description="Strand read counts: ref/fwd, ref/rev, var/fwd, var/rev"',
                'ID=DP4,Number=4,Type=Integer,Description="# high-quality ref-forward bases, ref-reverse, alt-forward and alt-reverse bases"')
            filtered_vcf.write(new_DP4)
        elif line.startswith('#CHROM'):
            headers = line.strip().split('\t')
            filtered_vcf.write(line)
        elif not line.startswith('#'):
            columns = line.strip().split('\t')
            Filter = columns[headers.index('FILTER')]
            INFO = columns[headers.index('INFO')]
            if 'SOMATIC' in INFO and 'PASS' in Filter:
                filtered_vcf.write(line)
        else:
            filtered_vcf.write(line)
    vcf.close()
    filtered_vcf.close()