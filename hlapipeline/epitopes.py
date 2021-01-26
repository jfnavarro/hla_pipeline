"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com
"""
import re
from Bio.Seq import translate
from varcode import Variant


def translate_dna(seq):
    return translate(seq, to_stop=True)


def create_epitope_varcode(chrm, start, ref, alt, transcript, db='GRCh37'):
    """
    This function computes an epitope (mutated peptide) from a given mutation (variant)
    It uses the package Varcode to obtain the WT and mutated sequences as well as the
    mutation position.
    :param db: the Ensembl database to use (GRCh37 or GRCh38)
    :param chrm: the chromosome of the variant
    :param start: the position of the variant
    :param transcript: the transcript ID of the variant
    :param ref: REF (reference allele) of the variant
    :param alt: ALT (alternative allele) of the variant
    :return:
        The AA position
        Error creating flags if any
        The WT peptide (+12 -12)
        The mutated peptide (+12 -12)
    """
    errors = list()
    wt_mer = '-'
    mut_mer = '-'
    pos = -1

    # Retrieve variant info
    vinfo = Variant(contig=chrm, start=start, ref=ref, alt=alt, ensembl=db)
    effects = vinfo.effects()
    effect = None
    for e in effects:
        if e is not None and e.transcript_id == transcript:
            effect = e
    if effect is None:
        errors.append('Could not find effects for this transcript (using top effect)')
        effect = effects.top_priority_effect()

    # Retrieve effect type
    protein_mut = effect.short_description
    if protein_mut is None:
        errors.append('Could not retrieve AA mutation for this effect')
    elif not protein_mut.startswith('p.'):
        errors.append('Invalid AA mutation: {}'.format(protein_mut))
    elif protein_mut.startswith('p.X'):
        errors.append('AA mutation occurs in stop codon')
    else:
        # Retrieve AA mut pos (it is already ajudted to 0-based)
        pos = effect.aa_mutation_start_offset
        if pos is None:
            errors.append('Could not find the position for this mutation')
        elif pos < 0:
            errors.append('Can not code for this mutated position')
        elif pos == 0:
            errors.append('Mutation occurs in start codon')
        else:
            if effect.mutant_protein_sequence is None or effect.original_protein_sequence is None:
                errors.append('Could not retrieve protein sequences')
            else:
                # Type of effect
                effect_type = type(effect).__name__
                if 'Stop' in effect_type:
                    errors.append('Stop mutation')
                elif 'FrameShift' in effect_type or 'Substitution' in effect_type \
                        or 'Insertion' in effect_type or 'Deletion' in effect_type:
                    end_wt = pos + 13
                    if end_wt > len(effect.original_protein_sequence):
                        errors.append('End of sequence is shorter than 12aa from mutation (WT)')
                        end_wt = len(effect.original_protein_sequence)
                    start = pos - 12
                    if start < 0:
                        errors.append('Start of sequence is shorter than 12aa from mutation')
                        start = 0
                    wt_mer = effect.original_protein_sequence[start:end_wt]
                    if 'FrameShift' in effect_type:
                        mut_mer = effect.mutant_protein_sequence[start:]
                    else:
                        end_mut = pos + 13
                        if end_mut > len(effect.mutant_protein_sequence):
                            errors.append('End of sequence is shorter than 12aa from mutation (MUT)')
                            end_mut = len(effect.mutant_protein_sequence)
                        mut_mer = effect.mutant_protein_sequence[start:end_mut]
                else:
                    errors.append('Unknown exonic function {}'.format(effect_type))
    return pos, ';'.join(errors), wt_mer, mut_mer


def create_epitope(ref, alt, exonic_func, cDNA_mut, protein_mut, cDNA_seq, protein_seq):
    """
    This function computes an epitope (mutated peptide) from a given mutation (variant)
    :param ref: REF (reference allele) of the variant
    :param alt: ALT (alternative allele) of the variant
    :param exonic_func: the Annovar exonic function of the variant
    :param cDNA_mut: the cDNA mutation of the variant
    :param protein_mut: the AA mutation of the variant
    :param cDNA_seq: the cDNA sequence of the transcript where the variant is located
    :param protein_seq: the protein sequence of the transcript where the variant is located
    :return:
        The AA position
        Error creating flags if any
        The WT peptide (+12 position -12)
        The mutated peptide (+12 position -12)
    """
    errors = list()
    WT_25mer = '-'
    Mut_25mer = '-'
    # position is given in a starting by 1 coordinate system
    cDNA_pos = int(re.findall(r'\d+', cDNA_mut)[0]) - 1
    position = int(re.findall(r'\d+', protein_mut)[0]) - 1
    # nonsynonymous point mutations
    if exonic_func == 'nonsynonymous_SNV':
        # extract the AA change info
        ref_AA = protein_mut[protein_mut.find('.') + 1]
        var_AA = protein_mut[len(protein_mut) - 1]
        # gather AA seq for transcript
        if protein_seq == 'None':
            errors.append('AA seq not present for this transcript')
        else:
            # create mut and wt epitopes by extracting 12
            # AAs before the mutation, the mutation
            # and then 12 AAs after the mutation
            if position < 0:
                errors.append('Can not code for this mutated position')
            elif position == 0:
                errors.append('Mutation occurs in start codon')
            else:
                end = position + 13
                if end > len(protein_seq):
                    errors.append('End of sequence is shorter than 12aa from mutation')
                    end = len(protein_seq)
                start = position - 12
                if start < 0:
                    errors.append('Start of sequence is shorter than 12aa from mutation')
                    start = 0
                WT_25mer = protein_seq[start:end]
                Mut_25mer = protein_seq[start:position] + var_AA + protein_seq[position + 1:end]
    # frameshift/non-frameshift insertions/deletions/substitutions
    elif 'frameshift' in exonic_func:
        if cDNA_seq == 'None':
            errors.append('cDNA seq not present for this transcript')
        else:
            if exonic_func in ['frameshift_deletion', 'nonframeshift_deletion']:
                len_del = len(ref) - len(alt)
                mut_cDNA_seq = cDNA_seq[:cDNA_pos] + cDNA_seq[cDNA_pos + len_del:]
            elif exonic_func in ['frameshift_insertion', 'nonframeshift_insertion']:
                key = 'dup' if 'dup' in cDNA_mut else 'ins'
                ins = cDNA_mut[int(cDNA_mut.find(key)) + 3:]
                mut_cDNA_seq = cDNA_seq[:cDNA_pos] + ins + cDNA_seq[cDNA_pos + 1:]
            elif exonic_func in ['frameshift_substitution', 'nonframeshift_substitution']:
                subs = cDNA_mut[int(cDNA_mut.find('delins')) + 6:]
                len_subs = len(subs)
                mut_cDNA_left = cDNA_seq[:cDNA_pos]
                mut_cDNA_right = cDNA_seq[cDNA_pos + len_subs:]
                mut_cDNA_seq = mut_cDNA_left + subs + mut_cDNA_right
            else:
                # this should never happen
                position = -1
            if not protein_mut.startswith('p.'):
                errors.append('Invalid mutation: {}'.format(protein_mut))
            elif protein_mut.startswith('p.X'):
                errors.append('Mutation occurs in stop codon')
            elif position < 0:
                errors.append('Can not code for this mutated position')
            elif position == 0:
                errors.append('Mutation occurs in start codon')
            else:
                # create mut and wt epitopes by extracting 12
                # AAs before the mutation, the mutation
                # and then 12 AAs after the mutation
                # in the case of the mut epitope end end is
                # extended until the first stop codon
                ref_FASTA = str(translate_dna(cDNA_seq.replace(' ', '')))
                mut_FASTA = str(translate_dna(mut_cDNA_seq.replace(' ', '')))
                end_wt = position + 13
                if end_wt > len(ref_FASTA):
                    errors.append('End of sequence is shorter than 12aa from mutation')
                    end_wt = len(ref_FASTA)
                start = position - 12
                if start < 0:
                    errors.append('Start of sequence is shorter than 12aa from mutation')
                    start = 0
                WT_25mer = ref_FASTA[start:end_wt]
                if 'nonframeshift' in exonic_func:
                    end_mut = position + 13
                    if end_mut > len(mut_FASTA):
                        errors.append('End of sequence is shorter than 12aa from mutation (MUT)')
                        end_mut = len(mut_FASTA)
                    Mut_25mer = mut_FASTA[start:end_mut]
                else:
                    Mut_25mer = mut_FASTA[start:]
                if not cDNA_seq.startswith('ATG'):
                    errors.append('No ATG start codon for this transcript cDNA')
    elif 'stop' in exonic_func:
        errors.append('Stop mutation')
    else:
        errors.append('Unknown exonic function {}'.format(exonic_func))
    return position, ';'.join(errors), WT_25mer, Mut_25mer
