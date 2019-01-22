#!/usr/local/bin/python3.5


import argparse
import os
from glob import glob
from ete3 import *

ncbi = NCBITaxa()


def load_species_entries_of_family(edges_file):
    list_species = list()
    with open(edges_file) as f:
        # Skip header
        next(f)
        for row in f:
            # Species entry is given by the 1st column
            # of the edge file
            list_species.append(row.split("\t")[0])
    f.close()
    return list_species


def species_entries_to_names_and_taxid(list_species):
    list_names = list()
    list_taxons = list()
    for i in range(0, len(list_species)):
        species_entry = list_species[i]
        try:
            # if species entry is an integer, then it might
            # be a taxonomic id
            int(species_entry)
            try:
                # if species entry is indeed a taxonomic id
                # use it to retrieve the species name
                list_names.append(list(ncbi.get_taxid_translator([species_entry]).values())[0])
                list_taxons.append(species_entry)
            except:
                # else keep the entry as species name
                # and assign 'NA' as taxon
                list_names.append(species_entry)
                list_taxons.append('NA')
        except:
            # if species entry is a string, then it might
            # be a species name
            try:
                translated_name = list(ncbi.get_name_translator([species_entry]).values())
                list_taxons.append(translated_name[0][0])
                list_names.append(species_entry)
            except:
                list_names.append(species_entry)
                list_taxons.append('NA')
    return [list_names, list_taxons]


def species_taxons_to_tax_info(list_taxons):
    list_lineages = list()
    list_ranks = list()
    for i in range(0, len(list_taxons)):
        # Use taxonomic id to extract lineages and ranks from
        # NCBI taxonomy db
        if list_taxons[i] != 'NA':
            list_lineages.append(ncbi.get_lineage(list_taxons[i]))
            list_ranks.append(ncbi.get_rank(list_lineages[-1]))
        else:
            list_lineages.append('[NA]')
            list_ranks.append([{'NA': 'NA'}])
    return [list_lineages, list_ranks]


def load_tax_info_of_families(list_edges_files):
    super_dict = dict()
    super_dict['name'] = dict()
    super_dict['taxon'] = dict()
    super_dict['lineage'] = dict()
    super_dict['ranks'] = dict()
    for i in range(0, len(list_edges_files)):
        curr_family_name = os.path.basename(list_edges_files[i]).split("_MultiTwin_edges.csv")[0]
        species_entries = load_species_entries_of_family(list_edges_files[i])
        super_dict['name'][curr_family_name], super_dict['taxon'][curr_family_name] = \
            species_entries_to_names_and_taxid(species_entries)
        super_dict['lineage'][curr_family_name], super_dict['ranks'][curr_family_name] = \
            species_taxons_to_tax_info(super_dict['taxon'][curr_family_name])
    return super_dict


def load_leaves_entries_of_tree(tree):
    list_leaves = list()
    list_names = list()
    for node in tree.traverse():
        if node.is_leaf():
            list_leaves.append(node)
            list_names.append(node.name)
    return [list_leaves, list_names]



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
family_dict = load_tax_info_of_families(list_edges_files)
family_names = family_dict['taxon'].keys()

# LOAD REFERENCE TREE
tree = Tree(args.tree, format=1)
tree_leaves, tree_names = load_leaves_entries_of_tree(tree)
tree_dict = dict()
tree_dict['name'], tree_dict['taxon'] = \
    species_entries_to_names_and_taxid(tree_names)
tree_dict['lineage'], tree_dict['ranks'] = \
    species_taxons_to_tax_info(tree_dict['taxon'])

