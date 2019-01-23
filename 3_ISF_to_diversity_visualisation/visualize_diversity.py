#!/usr/local/bin/python3.5


import argparse
import os
from glob import glob
from ete3 import *
import numpy as np
import pandas as pd

ncbi = NCBITaxa()


######################################################
# FUNCTIONS TO EXTRACT, CONVERT TAXONOMIC ATTRIBUTES #
######################################################
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
            list_lineages.append('NA')
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
    list_taxons = list()
    list_names = list()
    for node in tree.traverse():
        if node.is_leaf():
            list_leaves.append(node)
            list_names.append(node.name)
            list_taxons.append(node.taxon)
    return [list_leaves, list_names, list_taxons]


def taxon_to_leaves(list_leaves, list_lineage):
    # create a dictionary that tell for any rank
    # of a lineage the leaves that are concerned
    taxon_to_leaf = dict()
    for i in range(0, len(list_leaves)):
        curr_lineage = list_lineage[i]
        if curr_lineage[0] != 'NA':
            for j in range(0, len(curr_lineage)):
                if curr_lineage[j] in taxon_to_leaf:
                    taxon_to_leaf[curr_lineage[j]] += [list_leaves[i]]
                else:
                    taxon_to_leaf[curr_lineage[j]] = [list_leaves[i]]
    return taxon_to_leaf


######################################################
# FUNCTION TO ANCHOR A SPECIES ON TREE NODES #########
######################################################
def find_best_matches(family_dict, taxon_to_leaves, hist_df, lin_thresh = 2):
    family_names = list(family_dict['taxon'])
    matched_leaves = dict()
    unanchored_taxons = dict()
    for i in range(0, len(family_names)):
        # Consider the current family of genes and access the
        # lineages of all its representative species
        matched_leaves[family_names[i]] = dict()
        matched_leaves[family_names[i]]['NA'] = list()
        unanchored_taxons[family_names[i]] = list()
        fam_lineages = family_dict['lineage'][family_names[i]]
        for j in range(0, len(fam_lineages)):
            # Consider the lineage of the current species in the family
            curr_lineage = fam_lineages[j]
            if curr_lineage != 'NA':
                matched_leaves[family_names[i]][curr_lineage[-1]] = list()
                for k in range(len(curr_lineage)-1, lin_thresh, -1):
                    # Visit the lineage from descendants to ascendants until reaching the
                    # threshold rank (a threshold is set because we don't want a species
                    # to be projected in all bacteria leaves for instance)
                    if curr_lineage[k] in taxon_to_leaves:
                        # if descendant is found in one or several leaves,
                        # then retrieves the name of these leaves and stop
                        # exploring the lineage
                        curr_matched_leaves = taxon_to_leaves[curr_lineage[k]]
                        matched_leaves[family_names[i]][curr_lineage[-1]] += curr_matched_leaves
                        break
                if not matched_leaves[family_names[i]][curr_lineage[-1]]:
                    # If no matched leaves for the current lineage, then add the species
                    # to the list of unanchored species of the current family
                    unanchored_taxons[family_names[i]] += [curr_lineage[-1]]
                    hist_df[family_names[i]]['UNANCHORED'] += 1
                else:
                    # Update the occurrence matrix of each leaf for the current
                    # family of gene. If one species of the family is projected in several
                    # leaves, divide the occurrence of these leaves by the number of matches
                    n = len(curr_matched_leaves)
                    for l in range(0, n):
                        hist_df[family_names[i]][curr_matched_leaves[l].name] += 1 / float(n)
            else:
                matched_leaves[family_names[i]]['NA'] += ['NA']
                unanchored_taxons[family_names[i]] += ['NA']
                hist_df[family_names[i]]['UNANCHORED'] += 1
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
parser.add_argument('-o', '--output_dir', dest='output_dir', type=str,
                    help='specify that path to the output directory')
parser.add_argument('-t', '--reference_tree', dest='tree', type=str,
                    help='specify that path to a reference phylogenetic tree in Newick format'
                         '(each leaf must have a "taxon" attribute corresponding to its NCBI taxonomic id)')
args = parser.parse_args()
if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

# LOAD EDGES FILES
list_edges_files = glob(os.path.join(args.input_dir, '*MultiTwin_edges.csv'))
family_dict = load_tax_info_of_families(list_edges_files)
family_names = list(family_dict['taxon'])


# LOAD REFERENCE TREE
tree = Tree(args.tree, format=1)
tree_leaves, tree_names, tree_taxons = load_leaves_entries_of_tree(tree)
tree_dict = dict()
tree_dict['name'], tree_dict['taxon'] = \
    species_entries_to_names_and_taxid(tree_taxons)
tree_dict['lineage'], tree_dict['ranks'] = \
    species_taxons_to_tax_info(tree_dict['taxon'])
taxon_to_leaves_dict = taxon_to_leaves(tree_leaves, tree_dict['lineage'])


# FIND BEST LEAVES TO ANCHOR EACH SPECIES OF A
# FAMILY OF GENES IN THE REFERENCE TREE
# Initialize matrix of count
hist_df = pd.DataFrame(np.zeros((len(tree_leaves)+1, len(family_names))),
                       columns=family_names, index=tree_names + ['UNANCHORED'])
hist_df, matched_leaves, unanchored_taxons = \
    find_best_matches(family_dict, taxon_to_leaves_dict, hist_df)

# OUTPUTS
hist_df = np.floor(hist_df)
hist_df.to_csv(path_or_buf=os.path.join(args.output_dir, "hist_mat.csv"),
               sep=",", header=True, index=True)
