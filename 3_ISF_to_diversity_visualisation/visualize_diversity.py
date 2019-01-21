#!/usr/local/bin/python3.5


import argparse
import os
from glob import glob
from ete3 import *

ncbi = NCBITaxa()


def assign_lineage_to_leaf(node):
    if node.is_leaf():
        if node.name != 'NA':
            try:
                node_lineage = ncbi.get_lineage(int(node.name))
            except:
                node_lineage = ['NA']
        else:
            node_lineage = ['NA']
        node.lineage = node_lineage
    else:
        node_children = node.children
        for i in range(0, len(node_children)):
            assign_lineage_to_node(node_children[i])


parser = argparse.ArgumentParser(
    description='This script runs rpsblast with the ISF output gene family in fasta '
                'against COG db in order to assign functions to all the representative '
                'sequences of the family')
parser.add_argument('-i', '--input_dir', dest='input_dir', type=str,
                    help='specify that path to a directory that must store the '
                         'Multitwin_edges files of gene families (*_Multitwin_edges.csv)')
parser.add_argument('-t', '--reference_tree', dest='tree', type=str,
                    help='specify that path to a reference phylogenetic tree in Newick format'
                         '(the name of each leaf must correspond to a NCBI taxonomic id)')
args = parser.parse_args()

# LOAD EDGES FILES
list_edges_files = glob(os.path.join(args.input_dir, '*MultiTwin_edges.csv'))
list_family_names = list()
dict_species = dict()
dict_taxons = dict()
dict_lineage = dict()

for i in range(0, len(list_edges_files)):
    list_family_names.append(os.path.basename(list_edges_files[i]).split("_MultiTwin_edges.csv")[0])
    list_species = list()
    list_taxons = list()
    list_lineages = list()
    with open(list_edges_files[i]) as f:
        next(f)
        for row in f:
            curr_species = row.split("\t")[0]
            try:
                int(curr_species)
                list_taxons.append(curr_species)
                list_species.appendlist(ncbi.get_taxid_translator([curr_species]).values())[0]
                list_lineages.append(ncbi.get_lineage(list_taxons[-1]))
            except:
                list_species.append(curr_species)
                translated_name = list(ncbi.get_name_translator([curr_species]).values())
                if translated_name:
                    list_taxons.append(translated_name[0][0])
                    list_lineages.append(ncbi.get_lineage(list_taxons[-1]))
                else:
                    list_taxons.append('NA')
                    list_lineages.append('[NA]')
    f.close()
    dict_species[list_family_names[i]] = list_species
    dict_taxons[list_family_names[i]] = list_taxons
    print(dict_species['CDPS'])
    print(dict_taxons['CDPS'])

# LOAD REFERENCE TREE
tree = Tree(args.tree, format=1)
assign_lineage_to_leaf(tree)



