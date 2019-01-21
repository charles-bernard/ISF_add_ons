#!/usr/bin/env python

import sys
import argparse

parser = argparse.ArgumentParser(
    description='This script runs rpsblast with the ISF output gene family in fasta '
                'against COG db in order to assign functions to all the representative '
                'sequences of the family')
parser.add_argument('--input_dir', dest='input_dir', type=str,
					help = 'specify that path to a directory that must store the '
						   'Multitwin_Edges files of gene families (*_Multitwin_Edges.csv)')








# Arguments

#familyFile, fasta format file with the family ID and sequence IDs 
# e.g. :
# >F1
# 33
# 56
#

# familyFile = sys.argv[1]
#
# #annotationFile, TAB delimited file with the node ID and Taxonony information (need to be exactly as in the reference tree in newick format)
# annotationFile = sys.argv[2]

#treeFile, reference tree in newick format
# treeFile = sys.argv[1]
#
# # Read Tree
#
# from ete3 import *
#
# tree = Tree(treeFile, format=0, quoted_node_names=True)
#
# # Add internal nodes
#
# n = 0
#
# for node in tree.traverse():
# 	if not node.is_leaf():
# 		n+=1
# 		node.name = "C" + str(n)
#
#
# # Print the tree
# print tree.get_ascii(show_internal=True)
# print tree
#
# # tree.render('tree.pdf')
#
# treeNodes = []
# for node in tree.traverse("levelorder"):
# 	treeNodes.append(node.name)
#


# Annotation dico
# Annotation={}
#
# annotInfo = open(annotationFile,"r")
#
# for annot in annotInfo:
# 	annot=annot.rstrip()
# 	seqID,seqAnnot = annot.split("\t")
# 	Annotation[seqID]=seqAnnot
#
# annotInfo.close()
#
# # Start family annotation ( place on tree )
#
# familyInfo = open(familyFile,"r")
#
# output = open(familyFile+".treeAnnot",'w')
#
# members = set()
# NbFam = 0
#
# for info in familyInfo:
# 	info = info.rstrip()
# 	if ">" in info:
# 		if NbFam != 0 :
# 			#print members
# 			if len(members) == 1:
# 				output.write(famID+"\t"+list(members)[0]+"\n")
# 			else:
# 				try:
# 					posTree = tree.get_common_ancestor(members)
# 				except:
# 					posTree = "Unknown"
# 				try:
# 					output.write(famID+"\t"+posTree.name+"\n")
# 				except:
# 					output.write(famID+"\t"+posTree+"\n")
# 				#print famID+"\t"+posTree
# 			#print members
# 		famID = info[1:]
# 		members = set()
# 		NbFam+=1
# 	else:
# 		if Annotation[info] in treeNodes:
# 			members.add(Annotation[info])
#
# # Print last family
# if len(members) == 1:
# 	output.write(famID+"\t"+list(members)[0]+"\n")
# else:
# 	try:
# 		posTree = tree.get_common_ancestor(members)
# 	except:
# 		posTree = "Unknown"
# 	try:
# 		output.write(famID+"\t"+posTree.name+"\n")
# 	except:
# 		output.write(famID+"\t"+posTree+"\n")
#
#
# output.close()
# familyInfo.close()
#END
