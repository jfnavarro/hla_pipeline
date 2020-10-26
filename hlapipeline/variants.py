"""
@author: jfnavarro
"""
from hlapipeline.epitopes import create_epitope
from collections import defaultdict, namedtuple
import numpy as np
import vcfpy

Epitope = namedtuple('Epitope', 'transcript dnamut aamut flags wtseq mutseq')

class Variant:
    def __init__(self):
        self.chrom = None
        self.start = None
        self.ref = None
        self.alt = None
        self.effects = None
        self.epitopes = None
        self.callers = None
        self.num_callers = None
        self.status = None
        self.type = None

    @property
    def key(self):
        return "{}:{} {}>{}".format(self.chrom, self.start, self.ref, self.alt)

    @property
    def ensemblGene(self):
        gene = None
        for effect in self.effects:
            if 'ensGene' in effect:
                gene = effect.split("_")[1]
                break
        return gene

    def __str__(self):
        return '{}:{} {}>{} {} {}'.format(self.chrom, self.start, self.ref, self.alt, self.type, self.status)

def effects(record, cDNA_seq_dict, AA_seq_dict):
    """
    This function computes the effects of a variant (record from vcfpy) using
    3 databases (Ensembl, NCBI and UCSC). The function will also compute the
    mutated peptide of the variant for each effect/transcript.
    The function only considers nonsynonymous and framshift effects.
    :param record: A vcfpy record containing the variant infomation from Annovar
    :param cDNA_seq_dic: a dictionary of cDNA sequences of the transcripts
    :param AA_seq: a dictionary of AA sequences of the transcripts
    :return:
        A list of unique epitopes detected in the variant (dna_mut, aa_mut, flags, wt_peptide, mut_peptide)
        A dict of unique effects detected in each of the variants unique mutated peptide
    """
    funcensGene = ''.join(record.INFO['ExonicFunc.ensGene'])
    has_func_ens = 'nonsynonymous' in funcensGene or 'frame' in funcensGene
    funcknownGene = ''.join(record.INFO['ExonicFunc.knownGene'])
    has_func_known = 'nonsynonymous' in funcknownGene or 'frame' in funcknownGene
    funcRefGene = ''.join(record.INFO['ExonicFunc.refGene'])
    has_func_ref = 'nonsynonymous' in funcRefGene or 'frame' in funcRefGene
    already_processed = set()
    epitopes = list()
    effects_seen_in = defaultdict(set)
    if has_func_ens:
        for mutation in record.INFO['AAChange.ensGene']:
            if len(mutation.split(':')) == 5:
                gene, transcript, exon, mut_dna, mut_aa = mutation.split(':')
                cDNA_seq = cDNA_seq_dict.get(transcript, 'None').strip()
                AA_seq = AA_seq_dict.get(transcript, 'None').strip()
                pos, flags, wtmer, mutmer = create_epitope(record.REF, funcensGene, mut_dna, mut_aa, cDNA_seq, AA_seq)
                if mutmer not in already_processed and mutmer != "-":
                    already_processed.add(mutmer)
                    epitopes.append(Epitope(transcript, mut_dna, mut_aa, flags, wtmer, mutmer))
                if mutmer != "-":
                    effects_seen_in[mutmer].add('ensGene_{}_{}'.format(gene, funcRefGene))
    if has_func_known:
        for mutation in record.INFO['AAChange.knownGene']:
            if len(mutation.split(':')) == 5:
                gene, transcript, exon ,mut_dna, mut_aa = mutation.split(':')
                cDNA_seq = cDNA_seq_dict.get(transcript, 'None').strip()
                AA_seq = AA_seq_dict.get(transcript, 'None').strip()
                pos, flags, wtmer, mutmer = create_epitope(record.REF, funcknownGene, mut_dna, mut_aa, cDNA_seq, AA_seq)
                if mutmer not in already_processed and mutmer != "-":
                    already_processed.add(mutmer)
                    epitopes.append(Epitope(transcript, mut_dna, mut_aa, flags, wtmer, mutmer))
                if mutmer != "-":
                    effects_seen_in[mutmer].add('knownGene_{}_{}'.format(gene, funcRefGene))
    if has_func_ref:
        for mutation in record.INFO['AAChange.refGene']:
            if len(mutation.split(':')) == 5:
                gene, transcript, exon, mut_dna, mut_aa = mutation.split(':')
                cDNA_seq = cDNA_seq_dict.get(transcript, 'None').strip()
                AA_seq = AA_seq_dict.get(transcript, 'None').strip()
                pos, flags, wtmer, mutmer = create_epitope(record.REF, funcRefGene, mut_dna, mut_aa, cDNA_seq, AA_seq)
                if mutmer not in already_processed and mutmer != "-":
                    already_processed.add(mutmer)
                    epitopes.append(Epitope(transcript, mut_dna, mut_aa, flags, wtmer, mutmer))
                if mutmer != "-":
                    effects_seen_in[mutmer].add('refGene_{}_{}'.format(gene, funcRefGene))

    return epitopes, effects_seen_in


def filter_variants_germline(file, tumor_coverage, tumor_var_depth,
                             tumor_var_freq, num_callers, cDNA_seq_dict, AA_seq_dict):
    """
    This function parses a list of annotated germline variants from Annovar.
    It then apply some filters to the variants and computes the eptiopes of the variants effects.
    It returns a list of Variant() objects.
    :param file: the annovar annotated germline variants
    :param tumor_coverage: filter value for the number of total reads (DP)
    :param tumor_var_depth: filter value for the number of allelic reads (AD)
    :param tumor_var_freq: filter value for the Variant Allele Frequency (VAF)
    :param num_callers: filter value for the number of variant callers
    :param cDNA_seq_dict: dictionary of transcripts to cDNA sequences
    :param AA_seq_dict: dictionary of AA to cDNA sequences
    :return:
        A list of Variant() objects
    """
    variants = list()
    reader = vcfpy.Reader.from_path(file)
    for record in reader:

        funcensGene = ''.join(record.INFO['ExonicFunc.ensGene'])
        has_func_ens = 'nonsynonymous' in funcensGene or 'frame' in funcensGene
        funcknownGene = ''.join(record.INFO['ExonicFunc.knownGene'])
        has_func_known = 'nonsynonymous' in funcknownGene or 'frame' in funcknownGene
        funcRefGene = ''.join(record.INFO['ExonicFunc.refGene'])
        has_func_ref = 'nonsynonymous' in funcRefGene or 'frame' in funcRefGene

        if has_func_ens or has_func_known or has_func_ref:
            called = {x.sample: x.data for x in record.calls if x.called}
            filtered = dict()
            try:
                if 'HaplotypeCaller' in called and 'PASS' in record.FILTER:
                    tumor_DP = int(called['HaplotypeCaller']['DP'])
                    tumor_AD = int(called['HaplotypeCaller']['AD'][1])
                    tumor_VAF = np.around(tumor_AD / float(tumor_DP) * 100, 3)
                    if tumor_DP >= tumor_coverage and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth:
                        filtered['HaplotypeCaller'] = '{};{};{}'.format(tumor_DP, tumor_AD, tumor_VAF)
                if 'varscan' in called and 'PASS' in record.FILTER:
                    tumor_DP = int(called['varscan']['DP'])
                    tumor_AD = int(called['varscan']['AD'][0])
                    tumor_VAF = float(called['varscan']['FREQ'][0].replace('%', ''))
                    if tumor_DP >= tumor_coverage and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth:
                        filtered['varscan'] = '{};{};{}'.format(tumor_DP, tumor_AD, tumor_VAF)
            except KeyError:
                continue

            is_valid = len(filtered.keys()) >= num_callers
            variant_effects, dbs_seen_in = effects(record, cDNA_seq_dict, AA_seq_dict)
            variant = Variant()
            variant.chrom = record.CHROM
            variant.start = record.POS
            variant.ref = record.REF
            variant.alt = record.ALT[0].serialize()
            variant.callers = '|'.join(['{}:{}'.format(key, value) for key, value in filtered.items()])
            variant.num_callers = len(filtered)
            variant.status = is_valid
            variant.effects = [';'.join(dbs_seen_in[r[-1]]) for r in variant_effects]
            variant.epitopes = variant_effects
            variant.type = 'germline'
            variants.append(variant)

    return variants

def filter_variants_somatic(file, normal_coverage, tumor_coverage, tumor_var_depth,
                            tumor_var_freq, t2n_ratio, num_callers, num_callers_indel,
                            cDNA_seq, AA_seq):
    """
    This function parses a list of annotated somatic variants from Annovar.
    It then apply some filters to the variants and computes the eptiopes of the variants effects.
    It returns a list of Variant() objects.
    :param file: the annovar annotated somatic variants
    :param normal_coverage: filter value for the number of normal total reads (DP)
    :param normal_var_depth: filter value for the number normal allelic reads (AD)
    :param normal_var_freq: filter value for the normal Variant Allele Frequency (VAF)
    :param tumor_coverage: filter value for the number of tumor total reads (DP)
    :param tumor_var_depth: filter value for the number tumor alleic reads (AD)
    :param tumor_var_freq: filter value for the tumor Variant Allele Frequency (VAF)
    :param t2n_ratio: filter value for the ratio between tumor and normal VAFs
    :param num_callers: filter value for the number of variant callers
    :param num_callers_indel: filter value for the number of variant callers (indels)
    :param cDNA_seq_dict: dictionary of transcripts to cDNA sequences
    :param AA_seq_dict: dictionary of AA to cDNA sequences
    :return:
        A list of Variant() objects
    """
    variants = list()
    reader = vcfpy.Reader.from_path(file)
    for record in reader:
        funcensGene = ''.join(record.INFO['ExonicFunc.ensGene'])
        has_func_ens = 'nonsynonymous' in funcensGene or 'frame' in funcensGene
        funcknownGene = ''.join(record.INFO['ExonicFunc.knownGene'])
        has_func_known = 'nonsynonymous' in funcknownGene or 'frame' in funcknownGene
        funcRefGene = ''.join(record.INFO['ExonicFunc.refGene'])
        has_func_ref = 'nonsynonymous' in funcRefGene or 'frame' in funcRefGene
        if has_func_ens or has_func_known or has_func_ref:
            called = {x.sample: x.data for x in record.calls if x.called}
            filtered = dict()
            try:
                if 'NORMAL.mutect' in called and 'TUMOR.mutect' in called and 'PASS' in record.FILTER:
                    normal_DP = int(called['NORMAL.mutect']['DP'])
                    normal_AD = int(called['NORMAL.mutect']['AD'][1])
                    normal_VAF = np.around(float(called['NORMAL.mutect']['AF'][0]) * 100, 3)
                    tumor_DP = int(called['TUMOR.mutect']['DP'])
                    tumor_AD = int(called['TUMOR.mutect']['AD'][1])
                    tumor_VAF = np.around(float(called['TUMOR.mutect']['AF'][0]) * 100, 3)
                    tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                    if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                            and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                            and tumor_normal_ratio >= t2n_ratio:
                        filtered['mutect'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                        normal_AD,
                                                                        normal_VAF,
                                                                        tumor_DP,
                                                                        tumor_AD,
                                                                        tumor_VAF)
                if 'NORMAL.somaticsniper' in called and 'TUMOR.somaticsniper' in called:
                    normal_DP = int(called['NORMAL.somaticsniper']['DP'])
                    normal_AD = sum(called['NORMAL.somaticsniper']['DP4'][2:])
                    normal_VAF = np.around((normal_AD / float(normal_DP)) * 100, 3)
                    tumor_DP = int(called['TUMOR.somaticsniper']['DP'])
                    tumor_AD = sum(called['TUMOR.somaticsniper']['DP4'][2:])
                    tumor_VAF = np.around((tumor_AD / float(tumor_DP)) * 100, 3)
                    tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                    is_somatic = int(called['TUMOR.somaticsniper']['SS']) == 2
                    if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                            and tumor_VAF >= tumor_var_freq and tumor_AD >= tumor_var_depth \
                            and tumor_normal_ratio >= t2n_ratio and is_somatic:
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
                            and tumor_normal_ratio >= t2n_ratio:
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
                    normal_VAF = np.around((normal_AD2 / float(normal_DP)) * 100, 3)
                    # tumor_DP = int(called['TUMOR.strelka']['DP'])
                    tumor_AD1 = int(called['TUMOR.strelka'][ref_index][0])
                    tumor_AD2 = int(called['TUMOR.strelka'][alt_index][0])
                    tumor_DP = tumor_AD1 + tumor_AD2
                    tumor_VAF = np.around((tumor_AD2 / float(tumor_DP)) * 100, 3)
                    tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                    if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                            and tumor_VAF >= tumor_var_freq and tumor_AD2 >= tumor_var_depth \
                            and tumor_normal_ratio >= t2n_ratio:
                        filtered['strelka'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                         normal_AD,
                                                                         normal_VAF,
                                                                         tumor_DP,
                                                                         tumor_AD,
                                                                         tumor_VAF)
                if 'NORMAL.strelka_indel' in called and 'TUMOR.strelka_indel' in called and 'PASS' in record.FILTER:
                    # normal_DP = int(called['NORMAL.strelka_indel']['DP'])
                    normal_AD1 = int(called['NORMAL.strelka_indel']['TAR'][0])
                    normal_AD2 = int(called['NORMAL.strelka_indel']['TIR'][0])
                    normal_DP = normal_AD1 + normal_AD2
                    normal_VAF = np.around((normal_AD2 / float(normal_DP)) * 100, 3)
                    # tumor_DP = int(called['TUMOR.strelka_indel']['DP'])
                    tumor_AD1 = int(called['TUMOR.strelka_indel']['TAR'][0])
                    tumor_AD2 = int(called['TUMOR.strelka_indel']['TIR'][0])
                    tumor_DP = tumor_AD1 + tumor_AD2
                    tumor_VAF = np.around((tumor_AD2 / float(tumor_DP)) * 100, 3)
                    tumor_normal_ratio = tumor_VAF / normal_VAF if normal_VAF != 0 else t2n_ratio
                    if normal_DP >= normal_coverage and tumor_DP >= tumor_coverage \
                            and tumor_VAF >= tumor_var_freq and tumor_AD2 >= tumor_var_depth \
                            and tumor_normal_ratio >= t2n_ratio:
                        filtered['strelka'] = '{};{};{};{};{};{}'.format(normal_DP,
                                                                         normal_AD,
                                                                         normal_VAF,
                                                                         tumor_DP,
                                                                         tumor_AD,
                                                                         tumor_VAF)
            except KeyError:
                continue

            is_valid = len([x for x in filtered.keys() if 'indel' not in x]) >= num_callers \
                       or len([x for x in filtered.keys() if 'indel' in x]) >= num_callers_indel

            variant_effects, dbs_seen_in = effects(record, cDNA_seq, AA_seq)
            variant = Variant()
            variant.chrom = record.CHROM
            variant.start = record.POS
            variant.ref = record.REF
            variant.alt = record.ALT[0].serialize()
            variant.callers = '|'.join(['{}:{}'.format(key, value) for key, value in filtered.items()])
            variant.num_callers = len(filtered)
            variant.status = is_valid
            variant.effects = [';'.join(dbs_seen_in[r[-1]]) for r in variant_effects]
            variant.epitopes = variant_effects
            variant.type = 'somatic'
            variants.append(variant)

    return variants