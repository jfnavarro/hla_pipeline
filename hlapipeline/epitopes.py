# -*- coding: utf-8 -*-

import re
from re import sub
from Bio import SeqIO
from hlapipeline.common import translate_dna

def extract_peptides(input, output, sampleID):
    # Extract peptides
    print("Extracting peptides")
    snv = open(input)
    snv_lines = snv.readlines()
    header = snv_lines.pop(0).strip().split('\t')
    epitope_file = open(output, 'w')
    for line in snv_lines:
        columns = line.strip().split('\t')
        variant_key = columns[header.index('VARIANT-KEY')]
        chrom = columns[header.index('CHR')]
        start = columns[header.index('START')]
        stop = columns[header.index('END')]
        ref = columns[header.index('REF')]
        alt = columns[header.index('ALT')]
        func_ref_gene = columns[header.index('Func.refGene')]
        exonic_func_ref = columns[header.index('ExonicFunc.refGene')]
        AA_change_refGene = columns[header.index('AAChange.refGene')].split(',')
        func_UCSC_gene = columns[header.index('Func.knownGene')]
        exonic_func_UCSC = columns[header.index('ExonicFunc.knownGene')]
        AA_change_UCSCGene = columns[header.index('AAChange.knownGene')].split(',')
        func_ens_gene = columns[header.index('Func.ensGene')]
        exonic_func_ens = columns[header.index('ExonicFunc.ensGene')]
        AA_change_ensGene = columns[header.index('AAChange.ensGene')].split(',')
        if re.search(r'nonsynonymous', exonic_func_ref) or re.search(r'frame', exonic_func_ref):
            # RefSeq Annotation
            for entry in AA_change_refGene:
                epitope_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sampleID,
                                                                                     variant_key,
                                                                                     chrom,
                                                                                     start,
                                                                                     stop,
                                                                                     ref,
                                                                                     alt,
                                                                                     func_ref_gene,
                                                                                     exonic_func_ref,
                                                                                     sub(':', '\t', entry)))
            # UCSC annotation
        if re.search(r'nonsynonymous', exonic_func_UCSC) or re.search(r'frame', exonic_func_UCSC):
            for entry in AA_change_UCSCGene:
                epitope_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sampleID,
                                                                                     variant_key,
                                                                                     chrom,
                                                                                     start,
                                                                                     stop,
                                                                                     ref,
                                                                                     alt,
                                                                                     func_UCSC_gene,
                                                                                     exonic_func_UCSC,
                                                                                     sub(':', '\t', entry)))
            # Ensembl annotation
        if re.search(r'nonsynonymous', exonic_func_ens) or re.search(r'frame', exonic_func_ens):
            for entry in AA_change_ensGene:
                epitope_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(sampleID,
                                                                                     variant_key,
                                                                                     chrom,
                                                                                     start,
                                                                                     stop,
                                                                                     ref,
                                                                                     alt,
                                                                                     func_ens_gene,
                                                                                     exonic_func_ens,
                                                                                     sub(':', '\t', entry)))
    epitope_file.close()
    snv.close()

def create_epitopes(input, output, FASTA_AA_DICT, FASTA_cDNA_DICT):
    # Create list of AA and cDNA sequences
    print('Creating epitopes')
    AA_seq = dict()
    with open(FASTA_AA_DICT, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            AA_seq[record.id.split("|")[1]] = str(record.seq)
    cDNA_seq = dict()
    with open(FASTA_cDNA_DICT, "rU") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            cDNA_seq[record.id.split("|")[0]] = str(record.seq)
    epitope_file = open('SQL_Epitopes.txt', 'w')
    input_file = open('Formatted_epitope_variant.txt')
    header = 'SAMPLE_ID\tVARIANT-KEY\tCHR\tSTART\tSTOP\tREF\tALT\tfunc_ref_gene\texonic_func_ref\tGene\t' \
             'Transcript_ID\tExon_Numb\tNT_CHANGE\tAA_CHANGE\tPOSITION\tERRORS\tWT25MER\tMUT25MER\n'
    epitope_file.write(header)
    for line in input_file:
        columns = line.rstrip('\n').split('\t')
        if len(columns) < 14:
            print("Formatted epitote with wrong number of columns {}".format(','.join(columns)))
            continue
        try:
            ref = columns[5].strip()
            exonic_func = columns[8].strip()
            transcriptID = columns[10].strip()
            cDNA_strip = columns[12].strip()
            protein_strip = columns[13].strip()
            errors = 'Flags:'
            WT_25mer = '-'
            Mut_25mer = '-'
            position = int(re.findall(r'\d+', protein_strip)[0])
            cDNA_pos = int(re.findall(r'\d+', cDNA_strip)[0])
            protein_seq = AA_seq.get(transcriptID, 'AA_seq not present for this transcript').strip()
            ref_cDNA_seq = cDNA_seq.get(transcriptID, 'cDNA not present for this transcript').strip()
            # Nonsynonymous point mutations to 25 mers
            if exonic_func == 'nonsynonymous SNV':
                # extract the AA change info
                ref_AA = protein_strip[protein_strip.find('.') + 1]
                var_AA = protein_strip[len(protein_strip) - 1]
                # gather AA seq for transcript
                if protein_seq == 'AA_seq not present for this transcript':
                    errors += ' AA_seq not present for this transcript'
                else:
                    # Check annotation is correct
                    # Create mut and wt epitopes by extracting 12
                    # AAs before the mutation, the mutation
                    # and then 12 AAs after the mutation
                    # TODO we may be extending over protein_seq
                    # TODO we should check for that
                    FASTA_AA = protein_seq[position - 1:position]
                    if FASTA_AA == ref_AA:
                        if position >= 13:
                            WT_25mer = protein_seq[position - 13:position + 12]
                            Mut_25mer = protein_seq[position - 13:position - 1] + var_AA + protein_seq[position:position + 12]
                        elif position < 13 and position > 0:
                            WT_25mer = protein_seq[0:position + 12]
                            Mut_25mer = protein_seq[0:position - 1] + var_AA + protein_seq[position:position + 12]
                        elif position == 0:
                            errors += ' can not code for this mutated AA_position'
                        if position == 1:
                            errors += ' mutation occurs in start codon'
                    else:
                        errors += ' Ref in AA_seq does not match file Ref'
            # frameshift/non-frameshift insertions/deletions/substitutions
            elif 'frameshift' in exonic_func:
                if ref_cDNA_seq == 'cDNA not present for this transcript':
                    errors += ' cDNA not present for this transcript'
                else:
                    if exonic_func in ['frameshift deletion', 'nonframeshift deletion']:
                        len_del = len(ref)
                        mut_cDNA_left = ref_cDNA_seq[0:cDNA_pos - 1]
                        mut_cDNA_right = ref_cDNA_seq[cDNA_pos + len_del - 1:]
                        mut_cDNA_seq = mut_cDNA_left + mut_cDNA_right
                    elif exonic_func in ['frameshift insertion', 'nonframeshift insertion']:
                        if 'dup' in cDNA_strip:
                            ins = cDNA_strip[int(cDNA_strip.find('dup')) + 3:]
                        elif 'ins' in cDNA_strip:
                            ins = cDNA_strip[int(cDNA_strip.find('ins')) + 3:]
                        else:
                            errors += ' could not find mutation in cDNA'
                            ins = ''
                        mut_cDNA_left = ref_cDNA_seq[0:cDNA_pos]
                        mut_cDNA_right = ref_cDNA_seq[cDNA_pos:]
                        mut_cDNA_seq = mut_cDNA_left + ins + mut_cDNA_right
                    elif exonic_func in ['frameshift substitution', 'nonframeshift substitution']:
                        if 'delins' in cDNA_strip:
                            subs = cDNA_strip[int(cDNA_strip.find('delins')) + 6:]
                        else:
                            errors += ' could not find mutation in cDNA'
                            subs = ''
                        mut_cDNA_left = ref_cDNA_seq[0:cDNA_pos]
                        mut_cDNA_right = ref_cDNA_seq[cDNA_pos:]
                        mut_cDNA_seq = mut_cDNA_left + subs + mut_cDNA_right
                    else:
                        mut_cDNA_seq = ref_cDNA_seq
                        errors += ' unknown exonic function'
                    if not protein_strip.startswith('p.'):
                        position = 0
                    elif protein_strip.startswith('p.X'):
                        position = 0
                        errors += ' mutation occurs in stop codon'
                    # Obtain protein sequences from DNA
                    ref_FASTA = str(translate_dna(ref_cDNA_seq.replace(' ', '')))
                    mut_FASTA = str(translate_dna(mut_cDNA_seq.replace(' ', '')))
                    # TODO we may be extending over ref_FASTA and mut_FASTA
                    # TODO we should check for that
                    # Create mut and wt epitopes by extracting 12
                    # AAs before and after the mutation
                    if position >= 13:
                        WT_25mer = ref_FASTA[position - 13:position + 12]
                        Mut_25mer = mut_FASTA[position - 13:]
                    elif position < 13 and position > 0:
                        WT_25mer = ref_FASTA[0:position + 12]
                        Mut_25mer = mut_FASTA[0:]
                    elif position == 0:
                        errors += ' can not code for this mutated AA_position'
                    if not ref_cDNA_seq.startswith('ATG'):
                        errors += ' no ATG start codon for this transcript cDNA'
                    if position == 1:
                        errors += ' mutation occurs in start codon'
            elif re.search(r'^stop', exonic_func):
                position = ''.join([s for s in protein_strip if s.isdigit()])
                errors += ' stop mutation'
            else:
                errors += ' unknown exonic function'
            # NOTE just to make sure that they are not empty
            WT_25mer = "-" if not WT_25mer else WT_25mer
            mut_FASTA = "-" if not mut_FASTA else mut_FASTA
            epitope_file.write('{}\t{}\t{}\t{}\t{}\n'.format('\t'.join(columns[0:]),
                                                             position,
                                                             errors,
                                                             WT_25mer,
                                                             Mut_25mer))
        except Exception as err:
            print('An unknown error {}\t happened parsing epitope {}'.format(err, ','.join(columns)))
    input_file.close()
    epitope_file.close()