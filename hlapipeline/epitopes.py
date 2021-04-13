"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com
"""
import re
from Bio.Seq import translate
from Bio.Data.IUPACData import protein_letters_3to1
from varcode import Variant


def translate_dna(seq):
    return translate(seq, to_stop=True)


def create_epitope_varcode(chrm, start, ref, alt, db, transcript):
    # Retrieve variant info
    vinfo = Variant(contig=chrm, start=start, ref=ref, alt=alt, ensembl=db)
    effect = [effect for effect in vinfo.effects() if effect.transcript_id == transcript][0]
    errors = "Flags:"
    wt_mer = '-'
    mut_mer = '-'
    pos = -1
    if effect is None:
        errors += ' could not infer the effect'
    else:
        # Retrieve effect type
        protein_mut = effect.short_description
        if protein_mut is None:
            errors += ' could not retrieve AA mutation'
        elif not protein_mut.startswith('p.'):
            errors += ' invalid mutation {}'.format(protein_mut)
        elif protein_mut.startswith('p.X'):
            errors += ' mutation occurs in stop codon'
        else:
            # Retrieve pos
            pos = effect.aa_mutation_start_offset
            if pos is None:
                errors += ' could not find the position for this mutation'
            elif pos == 0:
                errors += ' can not code for this mutated position'
            elif pos == 1:
                errors += ' mutation occurs in start codon'
            else:
                if effect.mutant_protein_sequence is None or effect.original_protein_sequence is None:
                    errors += ' could not retrieve protein sequence'
                else:
                    # Type of effect
                    effect_type = type(effect).__name__
                    if 'Stop' in effect_type:
                        errors += ' stop mutation'
                    elif 'FrameShift' in effect_type:
                        wt_mer = effect.original_protein_sequence[pos-12:pos+13]
                        mut_mer = effect.mutant_protein_sequence[pos-12:]
                    elif 'Substitution' in effect_type \
                            or 'Deletion' in effect_type:
                        wt_mer = effect.original_protein_sequence[pos-12:pos+13]
                        mut_mer = effect.mutant_protein_sequence[pos-12:pos+13]
                    elif 'Insertion' in effect_type:
                        size = int(abs(len(ref) - len(alt)) / 3)
                        wt_mer = effect.original_protein_sequence[pos-12:pos+13+size]
                        mut_mer = effect.mutant_protein_sequence[pos-12:pos+13+size]
                    else:
                        errors += ' unknown exonic function {}'.format(effect_type)
    return pos, errors, wt_mer, mut_mer