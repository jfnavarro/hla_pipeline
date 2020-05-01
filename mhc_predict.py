from collections import Counter
from _collections import  defaultdict
import json
import argparse
from common import exec_command

def compute_MHC(hla_exome_cancer, hla_exome_normal, hla_rna, overlap_final):
    HLA_dict = defaultdict(list)

    # First parse hla_exome_cancer and normal (HLA-LA format)
    # MHCflurry format (HLA-A*02:01)
    print('Loading DNA tumor HLAs..')
    with open(hla_exome_cancer) as f:
        for line in f.readlines():
            columns = line.strip().split('\t')
            hla = columns[1].split("_")[-1]
            alleles = columns[1:]
            HLA_dict[hla].extend(alleles)
    print('Loading DNA normal HLAs..')
    with open(hla_exome_normal) as f:
        for line in f.readlines():
            columns = line.strip().split('\t')
            hla = columns[1].split("_")[-1]
            alleles = columns[1:]
            HLA_dict[hla].extend(alleles)

    # Parse RNA hlas (arcasHLA JSON format)
    # MHCflurry format (HLA-A*02:01)
    print('Loading RNA HLAs..')
    with open(hla_rna) as f:
        HLA_dict.update(json.load(f))

    # Filter HLAs by occurrences
    filtered_hla = []
    for alleles in HLA_dict.values():
        f = Counter(alleles)
        items = f.most_common()
        cutoff = items[0][1] if len(items) == 1 else items[1][1]
        filtered_hla += [y[0] for y in items if y[1] >= cutoff]

    # Create protein FASTA file
    print('Creating protein sequencess..')
    with open('protein_sequences.fasta', 'w') as fwrite:
        with open(overlap_final, 'r') as fread:
            lines = fread.readlines()
            header = lines.pop(0).strip().split('\t')
            for line in lines:
                columns = line.strip().split('\t')
                if int(columns[header.index('Number of Exomes samples (passing)')]) > 0\
                        or columns[header.index('Cosmic Info')] != '-':
                    protein_name = '{} {} {} {}'.format(columns[header.index('Variant key')],
                                                        columns[header.index('transcript ID')],
                                                        columns[header.index('cDNA change')],
                                                        columns[header.index('AA change')])
                    protein_seq = columns[header.index('Mut Epitope')].strip()
                    if protein_seq != '-':
                        fwrite.write('>{}\n{}\n'.format(protein_name, protein_seq))


    # Run prediction
    print('Predicting MHCs..')
    cmd = 'mhcflurry-predict-scan protein_sequences.fasta --alleles {} ' \
          '--peptide-lengths 8 9 10 11 12'.format(','.join(filtered_hla))
    exec_command(cmd)
    print('Completed')

parser = argparse.ArgumentParser(description='Script to predict MHC using MHCflurry and data from from Jareds pipeline '
                                             '(adjusted by Jose Fernandez) <jc.fernandes.navarro@gmail.com>',
                                 prog='mhc_predict.py',
                                 usage='mhc_predict.py [options] '
                                       '--hla-dna-normal [file with HLA predictions from DNA (Normal)] '
                                       '--hla-dna-tumor [file with HLA predictions from DNA (Tumor)] '
                                       '--hla-rna [file with HLA predictions from RNA]'
                                       '--final-table [file with the final variants generated with merge_results.py]')

parser.add_argument('--hla-dna-normal', default=None, required=True,
                    help='A file containing predicted HLAs from normal DNA (table format)')
parser.add_argument('--hla-dna-tumor', default=None, required=True,
                    help='A file containing predicted HLAs from tumor DNA (table format)')
parser.add_argument('--hla-rna', default=None, required=True,
                    help='A file containing predicted HLAs from RNA (JSON format)')
parser.add_argument('--final-table', default=None, required=True,
                    help='A file with the final variants generated with merge_results.py (table format)')

args = parser.parse_args()
compute_MHC(args.hla_dna_tumor, args.hla_dna_normal, args.hla_rna, args.final_table)