"""
@author: jfnavarro
"""
import re
from Bio.Seq import translate

def translate_dna(seq):
    return translate(seq, to_stop=True)

def create_epitope(ref, exonic_func, cDNA_mut, protein_mut, cDNA_seq, protein_seq):
    """
    This function computes an epitope (mutated peptide) from a given mutation (variant)
    :param ref: REF (reference) of the variant
    :param exonic_func: the Annovar exonic function of the variant
    :param cDNA_mut: the cDNA mutation of the variant
    :param protein_mut: the AA mutation of the variant
    :param cDNA_seq: the cDNA sequence of the transcript where the variant is located
    :param protein_seq: the protein sequence of the transcript where the variant is located
    :return:
        The AA position
        Error creating flags if any
        The WT peptide
        The mutated peptide (+12 -12 trimmed)
    """
    errors = 'Flags:'
    WT_25mer = '-'
    Mut_25mer = '-'
    position = int(re.findall(r'\d+', protein_mut)[0])
    cDNA_pos = int(re.findall(r'\d+', cDNA_mut)[0])
    # Nonsynonymous point mutations
    if exonic_func == 'nonsynonymous_SNV':
        # extract the AA change info
        ref_AA = protein_mut[protein_mut.find('.') + 1]
        var_AA = protein_mut[len(protein_mut) - 1]
        # gather AA seq for transcript
        if protein_seq == 'None':
            errors += ' AA_seq not present for this transcript'
        else:
            # Check annotation is correct and
            # create mut and wt epitopes by extracting 12
            # AAs before the mutation, the mutation
            # and then 12 AAs after the mutation
            FASTA_AA = protein_seq[position - 1:position]
            if FASTA_AA == ref_AA:
                if position == 0:
                    errors += ' can not code for this mutated AA_position'
                else:
                    end = position + 12 if position + 12 < len(protein_seq) else None
                    start = position - 13 if position > 13 else 0
                    WT_25mer = protein_seq[start:end]
                    Mut_25mer = protein_seq[start:position - 1] + var_AA + protein_seq[position:end]
                if position == 1:
                    errors += ' mutation occurs in start codon'
            else:
                errors += ' Ref in AA_seq does not match file Ref {} {}'.format(FASTA_AA, ref_AA)
    # frameshift/non-frameshift insertions/deletions/substitutions
    elif 'frameshift' in exonic_func:
        if cDNA_seq == 'None':
            errors += ' cDNA not present for this transcript'
        else:
            if exonic_func in ['frameshift_deletion', 'nonframeshift_deletion']:
                len_del = len(ref)
                mut_cDNA_left = cDNA_seq[0:cDNA_pos - 1]
                mut_cDNA_right = cDNA_seq[cDNA_pos + len_del - 1:]
                mut_cDNA_seq = mut_cDNA_left + mut_cDNA_right
            elif exonic_func in ['frameshift_insertion', 'nonframeshift_insertion']:
                if 'dup' in cDNA_mut:
                    ins = cDNA_mut[int(cDNA_mut.find('dup')) + 3:]
                elif 'ins' in cDNA_mut:
                    ins = cDNA_mut[int(cDNA_mut.find('ins')) + 3:]
                else:
                    errors += ' could not find mutation in cDNA'
                    ins = ''
                mut_cDNA_left = cDNA_seq[0:cDNA_pos]
                mut_cDNA_right = cDNA_seq[cDNA_pos:]
                mut_cDNA_seq = mut_cDNA_left + ins + mut_cDNA_right
            elif exonic_func in ['frameshift_substitution', 'nonframeshift_substitution']:
                if 'delins' in cDNA_mut:
                    subs = cDNA_mut[int(cDNA_mut.find('delins')) + 6:]
                else:
                    errors += ' could not find mutation in cDNA'
                    subs = ''
                mut_cDNA_left = cDNA_seq[0:cDNA_pos]
                mut_cDNA_right = cDNA_seq[cDNA_pos:]
                mut_cDNA_seq = mut_cDNA_left + subs + mut_cDNA_right
            else:
                mut_cDNA_seq = cDNA_seq
                errors += ' unknown exonic function {}'.format(exonic_func)
            if not protein_mut.startswith('p.'):
                position = 0
            elif protein_mut.startswith('p.X'):
                position = 0
                errors += ' mutation occurs in stop codon'
            # Create mut and wt epitopes by extracting 12
            # AAs before and after the mutation
            if position == 0:
                errors += ' can not code for this mutated AA_position'
            else:
                ref_FASTA = str(translate_dna(cDNA_seq.replace(' ', '')))
                mut_FASTA = str(translate_dna(mut_cDNA_seq.replace(' ', '')))
                end = position + 12 if position + 12 < len(ref_FASTA) else None
                start = position - 13 if position > 13 else 0
                WT_25mer = ref_FASTA[start:end]
                Mut_25mer = mut_FASTA[start:]
            if not cDNA_seq.startswith('ATG'):
                errors += ' no ATG start codon for this transcript cDNA'
            if position == 1:
                errors += ' mutation occurs in start codon'
    elif 'stop' in exonic_func:
        errors += ' stop mutation'
    else:
        errors += ' unknown exonic function {}'.format(exonic_func)
    return position, errors, WT_25mer, Mut_25mer
