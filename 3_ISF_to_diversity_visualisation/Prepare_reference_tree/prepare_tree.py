#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
# import pandas
# import re
from ete3 import *
import subprocess
import argparse
ncbi = NCBITaxa()


def create_output_dir(tree_path):
    # Creates an output directory one directory ahead of tree_path
    parent_dir = os.path.join(os.path.dirname(tree_path), '..')
    output_dir = os.path.join(parent_dir, "Prepared_trees")
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    return output_dir


def convert_itol_to_nhx(tree_path):
    # An Itol tree includes bootstrap values inside brackets but this cannot
    # be read by ete3. A conversion of format is therefore required to make ete3
    # consider these values as node annotations
    nhex_tree = os.path.join(os.path.dirname(tree_path),
                             os.path.basename(tree_path).split(".")[0] + "_NHX.txt")
    output_file = open(nhex_tree, "w")
    subprocess.call(['sed', 's/\[/\[\&\&NHX:bootstrap=/g', tree_path],
                    stdout=output_file)
    output_file.close()
    return nhex_tree


def name_internal_nodes(node, k):
    if not node.is_leaf():
        node.name = "Node_" + str(k)
        k += 1
        children = node.get_children()
        for i in range(0, len(children)):
            k = name_internal_nodes(children[i], k)
    return k


def read_node_dict(dict_path):
    # Associate to each node name, its corresponding taxonomic attributes
    tree_dict = dict()
    tree_dict['species'] = dict()
    tree_dict['taxid'] = dict()
    tree_dict['lineage'] = dict()
    tree_dict['rank'] = dict()
    with open(dict_path, "r") as input_file:
        # Skip header
        next(input_file)
        for row in input_file:
            row_fields = row.split("\t")
            node_name = row_fields[0]
            taxon = row_fields[1].strip()
            try:
                # if taxon is an integer, then it might be a taxonomic id
                int(taxon)
                try:
                    # if species entry is indeed a taxonomic id
                    # use it to retrieve the species name
                    tree_dict['species'][node_name] = list(ncbi.get_taxid_translator([taxon]).values())[0]
                    tree_dict['taxid'][node_name] = taxon
                except:
                    tree_dict['species'][node_name] = 'NA'
                    tree_dict['taxid'][node_name] = 'NA'
            except:
                # if taxon is a string, then it might be a species name
                try:
                    translated_name = list(ncbi.get_name_translator([taxon]).values())
                    tree_dict['species'][node_name] = taxon
                    tree_dict['taxid'][node_name] = (translated_name[0][0])
                except:
                    tree_dict['species'][node_name] = 'NA'
                    tree_dict['taxid'][node_name] = 'NA'
            # Use tax id to retrieve lineage and ranks
            if tree_dict['taxid'][node_name] != 'NA':
                tree_dict['lineage'][node_name] = ncbi.get_lineage(tree_dict['taxid'][node_name])
                tree_dict['rank'][node_name] = ncbi.get_rank(tree_dict['lineage'][node_name])
            else:
                tree_dict['lineage'][node_name] = ['NA']
                tree_dict['rank'][node_name] = [{'NA': 'NA'}]
    input_file.close()
    return tree_dict


def assign_br_avg_length(node, nb_attached_leaves, parent_dist):
    sum_distance_to_leaves = 0
    if node.is_leaf():
        nb_attached_leaves += 1
        sum_distance_to_leaves += node.dist + parent_dist
        node.br_avg_length = 0
    else:
        children = node.get_children()
        parent_dist += node.dist
        for i in range(0, len(children)):
            [curr_nb_attached_leaves, curr_sum_distance_to_leaves] = \
                assign_br_avg_length(children[i], 0, parent_dist)
            nb_attached_leaves += curr_nb_attached_leaves
            sum_distance_to_leaves += curr_sum_distance_to_leaves
        node.br_avg_length = (sum_distance_to_leaves - nb_attached_leaves * parent_dist) / nb_attached_leaves
    return [nb_attached_leaves, sum_distance_to_leaves]


def toto(tree):
    a = 0
    k = 0
    for node in tree.traverse():
        if node.is_leaf():
            a += node.get_distance(tree)
            k +=1
    print(a)
    print(k)


def collapse_by_br_length(node, cutoff, nodes_to_collapse):
    if not node.is_leaf():
        if node.br_avg_length < cutoff:
            nodes_to_collapse.append(node.name)
        else:
            children = node.get_children()
            for i in range(0, len(children)):
                collapse_by_br_length(children[i], cutoff, nodes_to_collapse)


##################################################################
# MAIN ###########################################################
##################################################################
# PARSE ARGUMENTS
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--tree', dest='input_tree', type=str,
                    help='specify that path to a reference tree')
parser.add_argument('-d', '--tree_dict', dest='tree_dict', type=str,
                    help='specify that path to a tabular dictionary '
                         'that associates a taxonomic id to each tree leaf name')
parser.add_argument('--collapse_method', dest='collapse', type=str,
                    help='specify the method to collapse branches '
                         '(can be either "rank deduplication" or "branch avg length")')
parser.add_argument('--length_cutoff', dest="length_cutoff", type=float,
                    help='specify the cutoff length if --collapse_method has been set '
                         'to "branch avg length" (0.65 by default)')
parser.add_argument('--rank_to_deduplicate', dest="rank_cutoff", type=float, default = "phylum",
                    help='specify the taxonomic rank to consider to collapse nodes if '
                         '--collapse_method has been set to "rank deduplication" ("phylum" by default)')
args = parser.parse_args()


# SETUP
output_dir = create_output_dir(os.path.abspath(args.tree))
nhx_tree = convert_itol_to_nhx(args.tree)
tree = Tree(nhx_tree, format=1)
name_internal_nodes(tree, 0)
tree.write(features=['name'], format=1, outfile=os.path.join(output_dir, "with_internal_nodes.txt"))
tree_dict = read_node_dict(args.tree_dict)
method_recognized = False

# COLLAPSE BRANCHES BY AVG LENGTH
if args.collapse == "branch avg length":
    method_recognized = True
    if args.length_cutoff:
        tree.br_avg_length = 0
        assign_br_avg_length(tree, 0, -tree.dist)
        nodes_to_collapse = list()
        collapse_by_br_length(tree, args.length_cutoff, nodes_to_collapse)
    else:
        sys.exit()

with open(os.path.join(output_dir, 'collapse.txt'), "w") as f:
    for i in range(0, len(nodes_to_collapse)):
        f.write("%s\n" % nodes_to_collapse[i])
    f.close()


# ############################################################
# # INPUT FILE 1: taxonomic attributes of nodes ##############
# ############################################################
# # node_name_to_attribute = sys.arv[1]
# node_attributes_file = "/home/charles/Downloads/infer_age_dev/leaves_id/final_dictionary.csv"
#
# # read all the taxonomic information related to a node name
# node_attributes = pandas.read_csv(node_attributes_file, sep="\t")
# node_attributes.columns = ['node', 'ncbi_accession_number', 'jgi_submission_id',
#                            'taxon', 'phylum', 'class', 'order', 'family', 'genus', 'species']
#
#
# ############################################################
# # INPUT FILE 2: phylo. tree in Newick ######################
# ############################################################
# # tree = sys.argv[2]
# # tree_file = "/home/charles/Downloads/infer_age_dev/input_trees/input_tree_based_on_ribosomal_proteins.txt"
# tree_file = "/home/charles/Downloads/collapsed_newick.txt"
# tree = Tree(tree_file, format=1)
#
# ############################################################
# # Add new taxonomic attributes to nodes ####################
# ############################################################
# for node in tree.traverse():
#     if node.is_leaf():
#         leaf = node
#         leaf_index = (node_attributes['node'] == leaf.name)
#         tax_id = node_attributes.loc[leaf_index, 'taxon']
#         if not tax_id.empty:
#             try:
#                 leaf.taxon = int(tax_id.values[0])
#             except:
#                 leaf.taxon = 'NA'
#         phylum = node_attributes.loc[leaf_index, 'phylum']
#         if not phylum.empty and phylum.values[0] != "<not present>":
#             phylum_name = list(ncbi.get_taxid_translator([phylum.values[0]]).values())[0]
#         else:
#             phylum_name = re.search(r'([A-Za-z]+_)+', leaf.name).group(0)
#         leaf.phylum = phylum_name
#
#
# ############################################################
# # Print attribute phylum along the tree ####################
# ############################################################
# with open("/home/charles/Downloads/printed_tree1.txt", 'w') as f:
#     f.write(tree.get_ascii(attributes=['name']))
#     f.close()
#
# # ############################################################
# # # Compute average branch length for each leaf ##############
# # ############################################################
# # nodes_to_retain = list()
# # for node in tree.traverse(strategy = 'postorder'):
# #     if node.is_leaf():
# #         parent = node.up
# #         parent_children = parent.children
# #         dist_sum = node.dist
# #         n_children_leaves = 1
# #         for i in range(0, len(parent_children)):
# #             if parent_children[i] != node and parent_children[i].is_leaf():
# #                 dist_sum += parent_children[i].dist
# #                 n_children_leaves += 1
# #         branch_avg_length = dist_sum / float(n_children_leaves)
# #         if branch_avg_length < 0.65:
# #             nodes_to_retain.append(node)
# #
# #
# # tree.prune(nodes_to_retain)
# #
# # # ############################################################
# # # # Print attribute phylum along the tree ####################
# # # ############################################################
# # with open("/home/charles/Documents/UPMC/infer_age_dev/leaves_id/printed_tree2.txt", 'w') as f:
# #     f.write(tree.get_ascii(attributes=['phylum']))
# #     f.close()
#
# tree.write(features=['name','taxon'], format=1,
#            outfile="/home/charles/Downloads/infer_age_dev/input_trees/reference_tree.txt")
#
#
#
# ############################################################
# # Collapse brother nodes into one if from the same phylum ##
# ############################################################
# def who_to_kill(node):
#     # If node is not a leaf, then visit its children
#     if not node.is_leaf():
#         children = node.children
#         n_children = len(children)
#         for i in range(0, n_children):
#             child = children[i]
#             who_to_kill(child)
#     # else visit its brothers
#     else:
#         parent = node.up
#         brothers = parent.children
#         n_brothers = len(brothers)
#         for i in range(0, n_brothers):
#             if brothers[i] != node:
#                 brother = brothers[i]
#                 # list the descendants of the visited brother
#                 if not brother.is_leaf():
#                     brother_desc = brother.get_descendants()
#                     brother_n_desc = len(brother_desc)
#                 else:
#                     brother_desc = [brother]
#                     brother_n_desc = 1
#                 # and find the major phylum among them
#                 brother_desc_phyla = list()
#                 for j in range(0, brother_n_desc):
#                     if brother_desc[j].is_leaf():
#                         brother_desc_phyla.append(brother_desc[j].phylum)
#                 brother_desc_major_phylum = max(set(brother_desc_phyla),
#                                                 key=brother_desc_phyla.count)
#                 # if major phylum in brother descendance is the same as visited node
#                 # then delete the leaves of the brother descendance corresponding to
#                 # the phylum
#                 if node.phylum == brother_desc_major_phylum:
#                     for j in range(0, brother_n_desc):
#                         if brother_desc[j].is_leaf():
#                             if brother_desc[j].phylum == brother_desc_major_phylum:
#                                 if brother_n_desc == 1:
#                                     try:
#                                         node.do_kill
#                                     except:
#                                         brother_desc[j].do_kill = True
#                                 else:
#                                     brother_desc[j].do_kill = True
#
#
#
#
#
#
# def collapse_twins(tree):
#     n = 1
#     while n > 0:
#         who_to_kill(tree)
#         nodes_to_kill = tree.search_nodes(do_kill=True)
#         n = len(nodes_to_kill)
#         for i in range(0, n):
#             nodes_to_kill[i].delete()
#
#
# collapse_twins(tree)
#
# ############################################################
# # Print attribute phylum along the tree ####################
# ############################################################
# with open("/home/charles/Documents/UPMC/infer_age_dev/leaves_id/printed_tree3.txt", 'w') as f:
#     f.write(tree.get_ascii(attributes=['phylum']))
#     f.close()
#
# tree.render("toto.pdf")