#!/usr/local/bin/python3.5

# Note: this script parses the header of each fasta with the aim
# to retrieve the corresponding species-sequence relationship(s).
# This parser doesn't handle yet all the peculiar syntaxes that may
# occur in a fasta header but should do it eventually. To this end,
# if this script fails at producing the expected output, I would be
# thankful if you could report the incident for further improvements.

import argparse
import os
import re

parser = argparse.ArgumentParser(
    description='This script retrieves the list of the unique species which possess ' 
                'at least one gene/protein of the family of genes extended by ISF and '
                'summarizes this information in a tabular file. '
                'Several of such tabular files, thus representative of different family of genes, '
                'may eventually be concatenated together to serve as input for '
                'https://github.com/TeamAIRE/MultiTwin to construct bipartite graphs. '
                'Additionally, this script summarizes the species-sequence relationship(s) '
                'of each and every sequence returned by ISF as a tabular file of the same '
                'structure as above. This might be useful if ISF was run against nr, because '
                'one sequence might be identical in many different species.')
parser.add_argument('-d', dest='isf_dir', type=str,
                    help='specify the path to the ISF output directory')
parser.add_argument('--family_name', dest="family_name", type=str,
                    help='specify the name of the gene family whose retrieved sequences belong to')
parser.add_argument('--consider_strains', dest="consider_strains", action='store_true', default=False,
                    help='tell the program that different strains of one species have to be '
                         'considered as as many different species (optional)')
parser.add_argument('--dictionary', dest='dictionary', type=str,
                    help='If you submitted to ISF a fasta file for which the header of the sequences '
                         'exhibits only a sequence id ; and if you have beside a tabular dictionary with a '
                         'description associated to each of these sequence ids, then you can specify its path '
                         'in order for the program to retrieve the corresponding sequence-species '
                         'relationship(s) (optional)')
# parser.add_argument('--tx_db', dest="taxonomy_db", type=str,
#                     help='specify the path to the \'fullnamelineage.dmp\' of the NCBI taxonomy database '
#                     'in case you would like to substitute the name of the species by their corresponding '
#                     'taxonomic id (optional)')
args = parser.parse_args()

######################################################################
# Read the comprehensive list of fasta ###############################
######################################################################
seq_list = os.path.join(args.isf_dir, 'sequence_found_and_bases.faa')
n_records = 0
id = list()
species = list()

with open(seq_list, mode='r') as infile:
    for line in infile:
        # retrieve fasta header
        header = re.search(r'^>.*$', line)

        if header:
            header = header.group(0)

            # get sequence id
            # (identified as the string that comes straight after the ">" char)
            tmp_id = re.search(r'^> *([a-zA-Z]|[0-9]|[-_.])+', header)
            tmp_id = re.split(r' *>', tmp_id.group(0))
            curr_id = tmp_id[1]

            # if the species is not mentioned directly in the fasta header but inside a
            # dictionary
            if args.dictionary:
                with open(args.dictionary, mode='r') as dic:
                    for dic_line in dic:
                        header = re.match(str(curr_id) + '.*$', dic_line)
                        if header:
                            break
                dic.close()
                header = header.group(0)


            # retrieve all unique species listed in the header
            # (identified as strings enclosed in brackets)
            tmp_species = re.findall(r'\[[A-Z][a-z] ?.*?]', header)

            if not tmp_species:
                curr_unique_species = list([''])
            else:
                if not args.consider_strains:
                    for i in range(0, len(tmp_species)):
                        # get rid of characters related to strains
                        # while adding exception for sp. cases
                        if not re.match(r'\[[A-Z][^A-Z0-9.-]*?\]', tmp_species[i]):
                            curr_string = re.match(r'^\[[A-Z][a-z]+( |_)?(sp\. ?[^\]]*)?([A-Z]?[a-z]+( |_)?)?', tmp_species[i])
                            curr_string = curr_string.group(0)
                            tmp_species[i] = curr_string + ']'
                for i in range(0, len(tmp_species)):
                    curr_tmp_species = tmp_species[i]
                    tmp_species[i] = curr_tmp_species[1:-1].rstrip()
                curr_unique_species = list(set(tmp_species))

            species = species + curr_unique_species
            # repeat sequence id as many times as there are related species
            for i in range(0, len(curr_unique_species)):
                id.append(curr_id)
infile.close()
unique_species = sorted(set(species))

######################################################################
# Write retrieved lists of species and ids into a tabular file #######
######################################################################
multi_twin_infile = os.path.join(args.isf_dir, args.family_name + '_MultiTwin_edges.csv')
with open(multitwin_infile, mode='w') as multi_twin_outfile:
    multi_twin_outfile.write('species\tgene_family\n')
    for i in range(0, len(unique_species)):
        if unique_species[i]:
            multi_twin_outfile.write('%s\t%s\n' % (unique_species[i], args.family_name))
multi_twin_outfile.close()

spe_seq_file = os.path.join(args.isf_dir, args.family_name + 'species_seq_relationships.csv')
with open(spe_seq_file, mode='w') as spe_seq_outfile:
    spe_seq_outfile.write('species\tsequence_id\n')
    for i in range(0, len(id)):
        spe_seq_outfile.write('%s\t%s\n' % (species[i], id[i]))
spe_seq_outfile.close()

######################################################################
# Turn species name into tx id #######################################
######################################################################
# if args.taxonomy_db:
#     print('You have chosen to convert species name to taxonomic id\n'
#           'This may take a little while ...')
#     multitwin_infile = os.path.join(args.isf_dir, 'MultiTwin_edges_with_txid.csv')
#     with open(multitwin_infile, mode='w') as outfile:
#         outfile.write('species\tgene_family\n')
#         for i in range(0, len(unique_species)):
#             with open(args.taxonomy_db, mode='r') as tx_file:
#                 for line in tx_file:
#                     fields = re.split('\t|\t', line)
#                     curr_species = fields[2]
#                     if curr_species == species[i]:
#                         unique_species[i] = fields[0]
#                         break
#             tx_file.close()
#             outfile.write('%s\t%s\n' % (unique_species[i], args.family_name))
#             print('retrieve species %d / %d\r' % (i+1, len(unique_species)),
#                   sep='', end='', flush=True)
#     outfile.close()
