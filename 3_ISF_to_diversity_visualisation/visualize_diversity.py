#!/usr/local/bin/python3.5


import argparse
import os
from glob import glob
from ete3 import *
import numpy as np
import pandas as pd

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
    # for each family of genes, get the corresponding list of species names,
    # taxons, lineages and ranks
    for i in range(0, len(list_edges_files)):
        # Family name given by the prefix of the edge file
        curr_family_name = os.path.basename(list_edges_files[i]).split("_MultiTwin_edges.csv")[0]
        species_entries = load_species_entries_of_family(list_edges_files[i])
        super_dict['name'][curr_family_name], super_dict['taxon'][curr_family_name] = \
            species_entries_to_names_and_taxid(species_entries)
        super_dict['lineage'][curr_family_name], super_dict['ranks'][curr_family_name] = \
            species_taxons_to_tax_info(super_dict['taxon'][curr_family_name])
    return super_dict


def load_leaves_entries_of_tree(tree):
    # get address and name of each leaf
    list_leaves = list()
    list_names = list()
    for node in tree.traverse():
        if node.is_leaf():
            list_leaves.append(node)
            list_names.append(node.name)
    return [list_leaves, list_names]


def taxon_to_leaves(list_leaves, list_lineage):
    # create a dictionary that tell for any rank
    # of a lineage the leaves that are concerned
    taxon_to_leaf = dict()
    for i in range(0, len(list_leaves)):
        curr_lineage = list_lineage[i]
        for j in range(0, len(curr_lineage)):
            if curr_lineage[j] in taxon_to_leaf:
                taxon_to_leaf[curr_lineage[j]].append(list_leaves[i])
            else:
                taxon_to_leaf[curr_lineage[j]] = [list_leaves[i]]
    return taxon_to_leaf


def find_best_matches(family_dict, taxon_to_leaf, hist_df, lin_thresh = 2):
    family_names = family_dict['taxon'].keys()
    matched_leaves = dict()
    unanchored_taxons = dict()
    for i in range(0, family_names):
        # Consider the current family of genes and access the
        # lineages of all its representative species
        matched_leaves[family_names[i]] = list()
        unanchored_taxons[family_names[i]] = list()
        fam_lineages = family_dict['lineage'][family_names[i]]
        for j in range(0, len(fam_lineages)):
            # Consider the lineage of the current species in the family
            curr_lineage = fam_lineages[j]
            for k in range(len(curr_lineage)-1, lin_thresh, 1):
                # Visit the lineage from descendants to ascendants
                # until reaching the threshold (a threshold is set because
                # we don't want a species to be projected in all bacteria for instance)
                if curr_lineage[k] in taxon_to_leaf:
                    # if descendant is found in one or several leaves,
                    # then retrieves the name of these leaves and stop
                    # exploring the lineage
                    curr_matched_leaves = taxon_to_leaf[curr_lineage[k]]
                    matched_leaves[family_names[i]].append(curr_matched_leaves)
            if not matched_leaves[family_names[i]]:
                # If no matched leaves for the current lineage, then add the species
                # to the list of unanchored species of the current family
                unanchored_taxons.append(family_dict['taxon'][family_names[i]])
                hist_df[family_names[i]]['UNANCHORED'] += 1
            else:
                n = len(curr_matched_leaves)
                for l in range(0, n):
                    hist_df[family_names[i]][curr_matched_leaves[l]] += 1 / float(n)
    return [hist_df, matched_leaves, unanchored_taxons]


##################################################################
# MAIN ###########################################################
##################################################################

# PARSE ARGUMENTS
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
taxon_to_leaves_dict = taxon_to_leaves(tree_leaves, tree_dict['lineage'])


# FIND BEST LEAVES TO ANCHOR EACH SPECIES OF A
# FAMILY OF GENES IN THE REFERENCE TREE
# Initialize matrix of count
df = pd.DataFrame(np.zeros((len(tree_leaves)+1, len(family_names))),
                  columns=family_names, index=tree_names + ['UNANCHORED'])



