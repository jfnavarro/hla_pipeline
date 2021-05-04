"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com
"""
from hlapipeline.epitopes import create_epitope_varcode
from collections import namedtuple
import numpy as np
import vcfpy
from pyensembl import EnsemblRelease

#  A convenience namedtuple to store the information of an epitope
Epitope = namedtuple('Epitope', 'transcript gene func dnamut aamut flags wtseq mutseq')

#  A convenience namedtuple to store the information of an annotated record (VEP)
Record_INFO = namedtuple('INFO', 'Allele Consequence SYMBOL Gene Feature_type Feature BIOTYPE \
    EXON INTRON HGVSc HGVSp cDNA_position CDS_position Protein_position Existing_variation FLAGS gnomAD_AF')


class Variant:
    def __init__(self):
        self.chrom = None
        self.start = None
        self.ref = None
        self.alt = None
        self.epitopes = None
        self.callers = None
        self.num_callers = None
        self.status = None
        self.type = None
        self.dbsnp = None
        self.gnomad = None
        self.cosmic = None
        self.gene = None

    @property
    def key(self):
        return "{}:{} {}>{}".format(self.chrom, self.start, self.ref, self.alt)

    def __str__(self):
        return '{}:{} {}>{} {} {}'.format(self.chrom, self.start, self.ref, self.alt, self.type, self.status)


def epitopes(record, info, ens_data):
    """
    This function computes the epitopes (mutated and wt peptides) of
    a VEP annotated variant (record from vcfpy) using the effects and
    their isoforms from Ensembl.
    The function only considers nonsynonymous and framshift effects.
    :param record: A vcfpy record containing the variant information from VEP
    :param info: 
    :param ens_data: Ensembl's database cache version used to annotate the VCF
    :return:
        A list of unique epitopes detected in the variant
        Epitope (transcript gene func dnamut aamut flags wtseq mutseq)
    """

    funcensGene = info.Consequence
    epitopes = list()
    if 'missense' in funcensGene or 'frame' in funcensGene:
        gene = info.SYMBOL
        transcript = info.Feature
        # sequence = ens_data.transcript_by_id(info.Feature)
        mut_dna = info.HGVSc.split(':')[1] if len(info.HGVSc.split(':')) > 1 else ''
        mut_aa = info.HGVSp.split(':')[1] if len(info.HGVSp.split(':')) > 1 else ''

        # TODO this should return a list 
        pos, flags, wtmer, mutmer = create_epitope_varcode(record.CHROM,
                                                           record.POS,
                                                           record.REF,
                                                           info.Allele,
                                                           ens_data,
                                                           transcript)
        epitopes.append(Epitope(transcript, gene, funcensGene, mut_dna, mut_aa, flags, wtmer, mutmer))
    return epitopes


def filter_variants_rna(file, tumor_coverage, tumor_var_depth,
                        tumor_var_freq, num_callers, ensembl_version):
    """
    This function processes a list of annotated RNA variants from VEP (VCF).
    It then applies some filters to the variants and computes the epitopes of each of
    the variants nonsynonymous and frameshift effects.
    The input is expected to contain HaplotypeCaller and Varscan RNA variants.
    It returns a list of Variant() objects.
    :param file: the VEP annotated RNA variants
    :param tumor_coverage: filter value for the number of total reads (DP)
    :param tumor_var_depth: filter value for the number of allelic reads (AD)
    :param tumor_var_freq: filter value for the Variant Allele Frequency (VAF)
    :param num_callers: filter value for the number of variant callers
    :param ensembl_version: ensembl version to use for the generation of epitoes
    :return:
        A list of Variant() objects
    """
    ens_data = EnsemblRelease(int(ensembl_version))
    variants = list()
    reader = vcfpy.Reader.from_path(file)
    for record in reader:
        for info in record.INFO['CSQ']:
            record_INFO = Record_INFO(*info.split('|'))
            funcensGene = record_INFO.Consequence
            has_func_ens = 'missense' in funcensGene or 'frame' in funcensGene
            avsnp150 = record_INFO.Existing_variation.split('&')[0] if 'rs' in record_INFO.Existing_variation else 'NA'
            gnomad_AF = record_INFO.gnomAD_AF if record_INFO.gnomAD_AF != '' else 'NA'
            cosm_count = record_INFO.Existing_variation.count('COSV')
            cosmic70 = ';'.join(record_INFO.Existing_variation.split('&')[-cosm_count::]) if cosm_count > 0 else 'NA'
            gene = record_INFO.SYMBOL

            if has_func_ens:
                called = {x.sample: x.data for x in record.calls if x.called}
                filtered = dict()
                pass_variants = 0
                try:
                    if 'HaplotypeCaller' in called and 'PASS' in record.FILTER:
                        tumor_DP = int(called['HaplotypeCaller']['DP'])
                        tumor_AD = int(called['HaplotypeCaller']['AD'][1])
                        tumor_VAF = np.around(tumor_AD / float(tumor_DP) * 100, 3) if tumor_DP > 0.0 else 0.0
                        if tumor_DP >= tumor_coverage and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth:
                            pass_variants += 1
                        filtered['HaplotypeCaller'] = '{};{};{}'.format(tumor_DP, tumor_AD, tumor_VAF)
                    if 'varscan' in called and 'PASS' in record.FILTER:
                        tumor_DP = int(called['varscan']['DP'])
                        tumor_AD = int(called['varscan']['AD'][0])
                        tumor_VAF = float(called['varscan']['FREQ'][0].replace('%', '')) if tumor_DP > 0.0 else 0.0
                        if tumor_DP >= tumor_coverage and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth:
                            pass_variants += 1
                        filtered['varscan'] = '{};{};{}'.format(tumor_DP, tumor_AD, tumor_VAF)
                except KeyError:
                    continue

                variant_epitopes = epitopes(record, record_INFO, ens_data)
                variant = Variant()
                variant.chrom = record.CHROM
                variant.start = record.POS
                variant.ref = record.REF
                variant.alt = record_INFO.Allele
                variant.callers = '|'.join(['{}:{}'.format(key, value) for key, value in filtered.items()])
                variant.num_callers = len(filtered)
                variant.status = pass_variants >= num_callers
                variant.epitopes = variant_epitopes
                variant.dbsnp = avsnp150
                variant.gnomad = gnomad_AF
                variant.cosmic = cosmic70
                variant.type = 'rna'
                variant.gene = gene
                variants.append(variant)

    return variants


def filter_variants_dna(file, normal_coverage, tumor_coverage, tumor_var_depth,
                        tumor_var_freq, normal_var_freq, t2n_ratio, num_callers,
                        num_callers_indel, ensembl_version):
    """
    This function processes a list of annotated DNA variants from VEP (VCF).
    It then applies some filters to the variants and computes the epitopes of each of
    the variants nonsynonymous and frameshift effects.
    The input is expected to contain Mutect2, Strelka, SomaticSniper and Varscan DNA variants.
    It returns a list of Variant() objects.
    :param file: the VEP annotated somatic variants
    :param normal_coverage: filter value for the number of normal total reads (DP)
    :param tumor_coverage: filter value for the number of tumor total reads (DP)
    :param tumor_var_depth: filter value for the number tumor alleic reads (AD)
    :param tumor_var_freq: filter value for the tumor Variant Allele Frequency (VAF)
    :param normal_var_freq: filter value for the normal Variant Allele Frequency (VAF)
    :param t2n_ratio: filter value for the ratio between tumor and normal VAFs
    :param num_callers: filter value for the number of variant callers
    :param num_callers_indel: filter value for the number of variant callers (indels)
    :param ensembl_version: ensembl version to use for epitope obtention
    :return:
        A list of Variant() objects
    """

    ens_data = EnsemblRelease(int(ensembl_version))
    variants = list()
    reader = vcfpy.Reader.from_path(file)
    for record in reader:
        for info in record.INFO['CSQ']:
            record_INFO = Record_INFO(*info.split('|'))
            funcensGene = record_INFO.Consequence
            has_func_ens = 'missense' in funcensGene or 'frame' in funcensGene
            avsnp150 = record_INFO.Existing_variation.split('&')[0] if 'rs' in record_INFO.Existing_variation else 'NA'
            gnomad_AF = record_INFO.gnomAD_AF if record_INFO.gnomAD_AF != '' else 'NA'
            cosm_count = record_INFO.Existing_variation.count('COSV')
            cosmic70 = ';'.join(record_INFO.Existing_variation.split('&')[-cosm_count::]) if cosm_count > 0 else 'NA'
            gene = record_INFO.SYMBOL

            if has_func_ens:
                called = {x.sample: x.data for x in record.calls if x.called}
                filtered = dict()
                pass_snp = 0
                pass_indel = 0
                try:
                    if 'NORMAL.mutect' in called and 'TUMOR.mutect' in called and 'PASS' in record.FILTER:
                        normal_DP = int(called['NORMAL.mutect']['DP'])
                        normal_AD = int(called['NORMAL.mutect']['AD'][1])
                        normal_VAF = np.around(float(called['NORMAL.mutect']['AF'][0]) * 100,
                                               3) if normal_DP > 0.0 else 0.0
                        tumor_DP = int(called['TUMOR.mutect']['DP'])
                        tumor_AD = int(called['TUMOR.mutect']['AD'][1])
                        tumor_VAF = np.around(float(called['TUMOR.mutect']['AF'][0]) * 100, 3)
                        tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                        if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                                and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio:
                            pass_snp += 1
                        filtered['mutect'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                        normal_AD,
                                                                        normal_VAF,
                                                                        tumor_DP,
                                                                        tumor_AD,
                                                                        tumor_VAF)
                    if 'NORMAL.somaticsniper' in called and 'TUMOR.somaticsniper' in called:
                        normal_DP = int(called['NORMAL.somaticsniper']['DP'])
                        normal_AD = sum(called['NORMAL.somaticsniper']['DP4'][2:])
                        normal_VAF = np.around((normal_AD / float(normal_DP)) * 100, 3) if normal_DP > 0.0 else 0.0
                        tumor_DP = int(called['TUMOR.somaticsniper']['DP'])
                        tumor_AD = sum(called['TUMOR.somaticsniper']['DP4'][2:])
                        tumor_VAF = np.around((tumor_AD / float(tumor_DP)) * 100, 3)
                        tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                        is_somatic = int(called['TUMOR.somaticsniper']['SS']) == 2
                        if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                                and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio and is_somatic:
                            pass_snp += 1
                        if is_somatic:
                            filtered['somaticsniper'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                                   normal_AD,
                                                                                   normal_VAF,
                                                                                   tumor_DP,
                                                                                   tumor_AD,
                                                                                   tumor_VAF)

                    if ('NORMAL.varscan' in called and 'TUMOR.varscan' in called) \
                            or ('NORMAL.varscan_indel' in called and 'TUMOR.varscan_indel' in called) \
                            and 'PASS' in record.FILTER and 'SOMATIC' in record.INFO:
                        label_index = 'varscan' if 'NORMAL.varscan' in called else 'varscan_indel'
                        normal_DP = int(called['NORMAL.{}'.format(label_index)]['DP'])
                        normal_AD = sum(called['NORMAL.{}'.format(label_index)]['DP4'][2:])
                        normal_VAF = float(called['NORMAL.{}'.format(label_index)]['FREQ'][0].replace('%', ''))
                        tumor_DP = int(called['TUMOR.{}'.format(label_index)]['DP'])
                        tumor_AD = sum(called['TUMOR.{}'.format(label_index)]['DP4'][2:])
                        tumor_VAF = float(called['TUMOR.{}'.format(label_index)]['FREQ'][0].replace('%', ''))
                        tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                        if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                                and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio:
                            if 'indel' in label_index:
                                pass_indel += 1
                            else:
                                pass_snp += 1
                        filtered[label_index] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                           normal_AD,
                                                                           normal_VAF,
                                                                           tumor_DP,
                                                                           tumor_AD,
                                                                           tumor_VAF)
                    if 'NORMAL.strelka' in called and 'TUMOR.strelka' in called and 'PASS' in record.FILTER:
                        ref_index = record.REF + 'U'
                        alt_index = str(record.ALT[0].serialize()) + 'U'
                        # normal_DP = int(called['NORMAL.strelka']['DP'])
                        normal_AD1 = int(called['NORMAL.strelka'][ref_index][0])
                        normal_AD2 = int(called['NORMAL.strelka'][alt_index][0])
                        normal_DP = normal_AD1 + normal_AD2
                        normal_VAF = np.around((normal_AD2 / float(normal_DP)) * 100, 3) if normal_DP > 0.0 else 0.0
                        # tumor_DP = int(called['TUMOR.strelka']['DP'])
                        tumor_AD1 = int(called['TUMOR.strelka'][ref_index][0])
                        tumor_AD2 = int(called['TUMOR.strelka'][alt_index][0])
                        tumor_DP = tumor_AD1 + tumor_AD2
                        tumor_VAF = np.around((tumor_AD2 / float(tumor_DP)) * 100, 3)
                        tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                        if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                and tumor_VAF >= tumor_var_freq and tumor_AD2 >= tumor_var_depth \
                                and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio:
                            pass_snp += 1
                        filtered['strelka'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                         normal_AD2,
                                                                         normal_VAF,
                                                                         tumor_DP,
                                                                         tumor_AD2,
                                                                         tumor_VAF)
                    if 'NORMAL.strelka_indel' in called and 'TUMOR.strelka_indel' in called and 'PASS' in record.FILTER:
                        # normal_DP = int(called['NORMAL.strelka_indel']['DP'])
                        normal_AD1 = int(called['NORMAL.strelka_indel']['TAR'][0])
                        normal_AD2 = int(called['NORMAL.strelka_indel']['TIR'][0])
                        normal_DP = normal_AD1 + normal_AD2
                        normal_VAF = np.around((normal_AD2 / float(normal_DP)) * 100, 3) if normal_DP > 0.0 else 0.0
                        # tumor_DP = int(called['TUMOR.strelka_indel']['DP'])
                        tumor_AD1 = int(called['TUMOR.strelka_indel']['TAR'][0])
                        tumor_AD2 = int(called['TUMOR.strelka_indel']['TIR'][0])
                        tumor_DP = tumor_AD1 + tumor_AD2
                        tumor_VAF = np.around((tumor_AD2 / float(tumor_DP)) * 100, 3)
                        tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                        if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                                and tumor_VAF >= tumor_var_freq and tumor_AD2 >= tumor_var_depth \
                                and normal_VAF <= normal_var_freq and tumor_normal_ratio >= t2n_ratio:
                            pass_indel += 1
                        filtered['strelka_indel'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                               normal_AD2,
                                                                               normal_VAF,
                                                                               tumor_DP,
                                                                               tumor_AD2,
                                                                               tumor_VAF)
                except KeyError:
                    continue

                variant_epitopes = epitopes(record, record_INFO, ens_data)
                variant = Variant()
                variant.chrom = record.CHROM
                variant.start = record.POS
                variant.ref = record.REF
                variant.alt = record.ALT[0].serialize()
                variant.callers = '|'.join(['{}:{}'.format(key, value) for key, value in filtered.items()])
                variant.num_callers = len(filtered)
                variant.status = pass_snp >= num_callers or pass_indel >= num_callers_indel
                variant.epitopes = variant_epitopes
                variant.dbsnp = avsnp150
                variant.gnomad = gnomad_AF
                variant.cosmic = cosmic70
                variant.type = 'dna'
                variant.gene = gene
                variants.append(variant)

    return variants
