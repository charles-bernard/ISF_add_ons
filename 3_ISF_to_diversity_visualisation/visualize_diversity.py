#!/usr/local/bin/python3.5


import argparse
import os
from glob import glob
from ete3 import *

ncbi = NCBITaxa()

parser = argparse.ArgumentParser(
    description='This script runs rpsblast with the ISF output gene family in fasta '
                'against COG db in order to assign functions to all the representative '
                'sequences of the family')
parser.add_argument('-i', '--input_dir', dest='input_dir', type=str,
                    help='specify that path to a directory that must store the '
                         'Multitwin_edges files of gene families (*_Multitwin_edges.csv)')
parser.add_argument('-t', '--reference_tree', dest='reference_tree', type=str,
                    help='specify that path to a reference phylogenetic tree in Newick format')
args = parser.parse_args()

# LOAD EDGES FILES
list_edges_files = glob(os.path.join(args.input_dir, '*MultiTwin_edges.csv'))
list_family_names = list()
dict_species = dict()
dict_taxons = dict()

for i in range(0, len(list_edges_files)):
    list_family_names.append(os.path.basename(list_edges_files[i]).split("_MultiTwin_edges.csv")[0])
    list_species = list()
    list_taxons = list()
    with open(list_edges_files[i]) as f:
        next(f)
        for row in f:
            curr_species = row.split("\t")[0]
            try:
                int(curr_species)
                list_taxons.append(curr_species)
                list_species.appendlist(ncbi.get_taxid_translator([curr_species]).values())[0]
            except:
                list_species.append(curr_species)
                translated_name = list(ncbi.get_name_translator([curr_species]).values())
                if translated_name:
                    list_taxons.append(translated_name[0][0])
                else:
                    list_taxons.append('')
    f.close()
    dict_species[list_family_names[i]] = list_species
    dict_taxons[list_family_names[i]] = list_taxons
    print(dict_species['CDPS'])
    print(dict_taxons['CDPS'])

# LOAD REFERENCE TREE


