#! /usr/bin/env python
"""
This tools uses the output of merge_results.py and HLAs predicted with either
the dna_pipeline.py and/or the rna_pipeline.py to predict neoantigens.
The tool extracts the WT and MUT peptides and makes affinity binding
predictions for the HLAs (class I). The tools uses MHC-flurry for the predictions.

@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
from hlapipeline.common import exec_command
from collections import Counter
from _collections import defaultdict
import json
from argparse import ArgumentParser, RawDescriptionHelpFormatter
import sys

def main(hla_dna, hla_rna, overlap_final, alleles_file, mode, results, results_filter):

    if not hla_dna and not hla_rna:
        sys.stderr.write("Error, need HLAs as input.\n")
        sys.exit(1)

    HLA_dict = defaultdict(list)

    # First parse hla_exome_cancer and normal (HLA-LA format)
    if hla_dna:
        print('Loading DNA HLAs..')
        for file in hla_dna:
            with open(file) as f:
                for line in f.readlines():
                    columns = line.strip().split('\t')
                    hla = columns[1].split("_")[-1]
                    alleles = columns[2:]
                    HLA_dict[hla].extend(alleles)

    # Parse RNA hlas (arcasHLA JSON format)
    if hla_rna:
        print('Loading RNA HLAs..')
        for file in hla_rna:
            with open(file) as f:
                local_dict = json.load(f)
                for hla, alleles in local_dict.items():
                    HLA_dict[hla].extend(alleles)

    # Filter HLAs by occurrences
    filtered_hla = []
    for alleles in HLA_dict.values():
        f = Counter(alleles)
        items = f.most_common()
        cutoff = items[0][1] if len(items) == 1 else items[1][1]
        filtered_hla += [y[0] for y in items if y[1] >= cutoff]
    # MHCflurry format (HLA-A*02:01)
    filtered_hla = ['HLA-{}'.format(':'.join(x.split(':')[0:2])) for x in filtered_hla]
    print('Alleles found: {}'.format(','.join(filtered_hla)))

    # Filter HLAs by allowed alleles in MHCflurry
    allowed_alleles = set()
    with open(alleles_file) as f:
        for x in f.readlines():
            allowed_alleles.add(x.strip())
    filtered_hla = [x for x in filtered_hla if x in allowed_alleles]
    print('Alleles allowed in MHCflurry: {}'.format(','.join(filtered_hla)))

    if len(filtered_hla) == 0:
        sys.stderr.write("Error, list of HLAs is empty.\n")
        sys.exit(1)

    # Create protein FASTA file
    print('Creating protein sequencess..')
    with open('protein_sequences_mu.fasta', 'w') as fwrite_mu:
        with open('protein_sequences_wt.fasta', 'w') as fwrite_wt:
            with open(overlap_final) as fread:
                lines = fread.readlines()
                header = lines.pop(0).strip().split('\t')
                added_proteins_mu = set()
                added_proteins_wt = set()
                for line in lines:
                    columns = line.strip().split('\t')
                    pass_somatic = int(columns[header.index('Number of DNA samples (passing)')]) > 0
                    pass_germline = int(columns[header.index('Number of RNA samples (passing)')]) > 0
                    if (mode == "both" and pass_somatic and pass_germline) or (mode == "dna" and pass_somatic)\
                        or (mode == "rna" and pass_germline) or (mode == "either" and (pass_somatic or pass_germline)):
                        protein_name = '{}_{}_{}'.format(''.join(columns[header.index('Variant key')].split()),
                                                         columns[header.index('cDNA change')],
                                                         columns[header.index('AA change')])
                        protein_seq_mu = columns[header.index('Mut Epitope')].strip()
                        protein_seq_wt = columns[header.index('Wt Epitope')].strip()
                        # Should probably make sure that all the letters in the protein seq are alpha (isalpha())
                        if protein_seq_mu != '-' and '.' not in protein_seq_mu and protein_seq_mu not in added_proteins_mu:
                            fwrite_mu.write('>{}\n{}\n'.format(protein_name, protein_seq_mu))
                            added_proteins_mu.add(protein_seq_mu)
                        # Should probably make sure that all the letters in the protein seq are alpha (isalpha())
                        if protein_seq_wt != '-' and '.' not in protein_seq_wt and protein_seq_wt not in added_proteins_wt:
                            fwrite_wt.write('>{}\n{}\n'.format(protein_name, protein_seq_wt))
                            added_proteins_wt.add(protein_seq_wt)
                del added_proteins_mu
                del added_proteins_wt

    results_filter = '' if results == 'all' else results_filter

    # Run predictions
    print('Predicting MHCs with MUT peptides..')
    cmd = 'mhcflurry-predict-scan protein_sequences_mu.fasta --alleles {} ' \
          '--no-throw --results-{} {} --out predictions_mut.csv --peptide-lengths 8-12'.format(' '.join(filtered_hla),
                                                                                               results, results_filter)
    exec_command(cmd)

    print('Predicting MHCs with WT peptides..')
    cmd = 'mhcflurry-predict-scan protein_sequences_wt.fasta --alleles {} ' \
          '--no-throw --results-{} {} --out predictions_wt.csv --peptide-lengths 8-12'.format(' '.join(filtered_hla),
                                                                                              results, results_filter)
    exec_command(cmd)

    print('Completed')

if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__, formatter_class=RawDescriptionHelpFormatter)
    parser.add_argument('--hla-dna', nargs='+', default=None, required=False,
                        help='A file or files containing predicted HLAs from normal DNA (HLA-LA table format)')
    parser.add_argument('--hla-rna', nargs='+', default=None, required=False,
                        help='A file or files containing predicted HLAs from RNA (arcasHLA JSON format)')
    parser.add_argument('--variants', default=None, required=True,
                        help='A file with the final variants generated with merge_results.py (table format)')
    parser.add_argument('--alleles', default=None, required=True,
                        help='A file containing the allowed alleles in MHCflurry')
    parser.add_argument('--mode', default='either',
                        help='Mode to use to extract peptides from the variants (both, dna (default), rna, either)',
                        choices=['both', 'dna', 'rna', 'either'])
    parser.add_argument('--results', default='best',
                        help='Whether to include all results for each peptide or only the best one (default=best)',
                        choices=['all', 'best'])
    parser.add_argument('--results-filter', default='affinity',
                        help='What filtering criteria to use when using --results best (default=affinity)',
                        choices=['presentation_score', 'processing_score', 'affinity', 'affinity_percentile'])
    args = parser.parse_args()
    main(args.hla_dna, args.hla_rna, args.variants, args.alleles, args.mode, args.results, args.results_filter)
