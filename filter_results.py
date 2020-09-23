#! /usr/bin/env python
"""
@author: Jose Fernandez Navarro <jc.fernandez.navarro@gmail.com>
"""
import argparse
import numpy as np
from _collections import defaultdict
import sys


def main(somatic, epitopes, germline, rna_counts):
    if not somatic and not germline:
        sys.stderr.write("Error, no variants given as input (somatic or germline).\n")
        sys.exit(1)

    variant_dict = {}


parser = argparse.ArgumentParser(
    description='Script that filters variants created with merge_results.py\n'
                'Created by Jose Fernandez <jc.fernandes.navarro@gmail.com>',
    prog='filter_results.py',
    usage='filter_results.py [options] results_file.txt\n'
          '--somatic [somatic variants results files]\n'
          '--epitope [epitopes results files]\n'
          '--germine [germine variants results files]\n'
          '--counts [rna gene counts results]')
parser.add_argument('--tcov', nargs='+', default=None, required=False,
                    help='List of files with the variants obtained with the somatic pipeline')
parser.add_argument('--treads', nargs='+', default=None, required=True,
                    help='List of files with the the epitotes (somatic and/or germline)')
parser.add_argument('--vaf', nargs='+', default=None, required=False,
                    help='List of files with the variants obtained with the germline pipeline')
parser.add_argument('--callers', nargs='+', default=None, required=True,
                    help='List of files with the gene counts results of the samples')

args = parser.parse_args()
main(args.somatic, args.epitope, args.rna, args.germine)
