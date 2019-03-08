#!/usr/local/bin/python3.5

# Note: this script parses the header of each fasta with the aim
# to retrieve the corresponding species-sequence relationship(s).
# This parser doesn't handle yet all the peculiar syntaxes that may
# occur in a fasta header but should do it eventually. To this end,
# if this script fails at producing the expected output, I would be
# thankful if you could report the incident for further improvements.

import argparse
import os
import logging
import re
import sys
import tempfile
import subprocess
from ete3 import *
ncbi = NCBITaxa()

def redirect_msg(msg, level = "info"):
    print(msg)
    if level == "info":
        logging.info(msg)


def retrieve_lineages(species):
    phylum = lineage = ""
    try:
        taxid = list(ncbi.get_name_translator([species]).values())[0][0]
    except:
        return ["", ""]

    lineage_list = ncbi.get_lineage(taxid)
    lineage = list(ncbi.get_taxid_translator([lineage_list[0]]).values())[0]
    lineage_ranks_dict = ncbi.get_rank(lineage_list)
    for i in range(1, len(lineage_list)):
        curr_taxid = lineage_list[i]
        curr_name = list(ncbi.get_taxid_translator([curr_taxid]).values())[0]
        curr_rank = lineage_ranks_dict[curr_taxid]
        if curr_rank == "phylum":
            phylum = curr_name
        lineage = lineage + '; ' + curr_name
    if lineage == "{}":
        return["", ""]
    else:
        return [lineage, phylum]


parser = argparse.ArgumentParser(
    description='For an ensemble of families of genes, this script associates to each family name a list of'
                ' representatives species.'
                ' This information for all families is eventually summarized as a tabular file that can serve as'
                ' input for https://github.com/TeamAIRE/MultiTwin to construct bipartite graphs.'
                'Additionally, this script summarizes the species-sequence relationship(s) '
                'of each and every sequence returned by ISF as a tabular file of the same '
                'structure as above. This might be useful if ISF was run against nr, because '
                'one sequence might be identical in many different species.')
parser.add_argument('-i', '--input_dir', dest='fa_dir', type=str,
                    help='specify the path to the directory '
                         'that contains several gene families in fasta (the filename of '
                         'each fasta will be used to name the corresponding gene family)')
parser.add_argument('-o', '--output_dir', dest="output_dir", type=str,
                    help='specify the path to the output directory')
parser.add_argument('-l', dest="get_lineages", action="store_true", default=False,
                    help='whether to output lineage of the species as well')
parser.add_argument('--consider_strains', dest="consider_strains", action='store_true', default=True,
                    help='tell the program that different strains of one species have to be '
                         'considered as as many different species (optional)')
parser.add_argument('--dictionary', dest='dictionary', type=str,
                    help='If you submitted to ISF a fasta file for which the header of the sequences '
                         'exhibits only a sequence id ; and if you have beside a tabular dictionary with a '
                         'description associated to each of these sequence ids, then you can specify its path '
                         'in order for the program to retrieve the corresponding sequence-species '
                         'relationship(s) (optional)')
args = parser.parse_args()

##########################################
# SETUP ##################################
##########################################
# get the directory of the executing script
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

# save the abspath of the output directory
args.output_dir = os.path.abspath(args.output_dir)

# initialize logfile
log_file = os.path.join(args.output_dir, 'logfile.log')
if os.path.isfile(log_file):
    os.remove(log_file)
logging.basicConfig(filename=os.path.join(args.output_dir, 'logfile.log'), level=logging.DEBUG)

# get family names and fill the list of paths to fasta
fasta_files = list()
family_names = list()
dir_files = os.listdir(args.fa_dir)
dir_files.sort()
for dir_file in dir_files:
    if re.search(r'\.fa(a)?(sta)?$', dir_file):
        fasta_files.append(os.path.join(args.fa_dir, dir_file))
        family_names.append(os.path.basename(dir_file).split(".")[0])

# define output
multi_twin_in = os.path.join(args.output_dir, "MultiTwin_input_file.tsv")
comprehensive_output_dir = os.path.join(args.output_dir, "species_sequences_relationships")
if not os.path.exists(comprehensive_output_dir):
    os.makedirs(comprehensive_output_dir)

# define dictionary of families.
multi_twin_dict = dict()
multi_twin_dict["id"] = dict()
multi_twin_dict["species"] = dict()
if args.get_lineages:
    multi_twin_dict["lineage"] = dict()
    multi_twin_dict["phylum"] = dict()

##########################################
# JOB ####################################
##########################################
for i in range(0, len(fasta_files)):
    curr_fa = fasta_files[i]
    curr_family = family_names[i]
    multi_twin_dict["species"][curr_family] = list()
    multi_twin_dict["id"][curr_family] = list()
    redirect_msg("  * processing family %d/%d: %s" % (i+1, len(fasta_files), curr_family))

    with open(fasta_files[i], mode="r") as f:

        for line in f:

            is_header = re.search(r'^>.*$', line)
            if is_header:
                header = is_header.group(0)

                # (id is identified as the string that comes straight after the ">" char)
                id = list()
                tmp_id = re.search(r'^> *([a-zA-Z]|[0-9]|[-._])+', header)
                tmp_id = re.split(r' *>', tmp_id.group(0))
                curr_id = tmp_id[1]

                redirect_msg("    * processing id %s" % (curr_id))

                curr_species = re.findall(r'\[.+\]', header)
                if args.get_lineages:
                    curr_lineages = list()
                    curr_phylum = list()

                if not curr_species:
                    curr_species_name = subprocess.check_output(
                        "efetch -db protein -id %s -format gpc \
                         | egrep 'INSDSeq_organism' \
                         | awk -F \"</?INSDSeq_organism>\" '{ printf($2); }'" % curr_id,
                        shell=True).decode('utf-8')
                    curr_species = curr_species_name

                j = 0
                while j < len(curr_species):
                    if not args.consider_strains:
                        # get rid of characters related to strains while adding exception for sp. cases
                        if not re.match(r'\[[A-Z][^A-Z0-9.-]*?\]', curr_species[j]):
                            curr_string = re.match(r'^\[[A-Z][a-z]+( |_)?(sp\. ?[^\]]*)?([A-Z]?[a-z]+( |_)?)?',
                                                   curr_species[j])
                            try:
                                curr_species[j] = curr_string.group(0)
                            except:
                                curr_species[j] = curr_species[j]
                    curr_species[j] = curr_species[j][1:-1].rstrip()
                    id.append(curr_id)
                    if args.get_lineages:
                        try:
                            tmp = multi_twin_dict["lineage"][curr_species[j]]
                        except:
                            lineage_out = retrieve_lineages(curr_species[j])

                            if not lineage_out[0]:
                                curr_species_name = subprocess.check_output(
                                    "efetch -db protein -id %s -format gpc \
                                     | egrep 'INSDSeq_organism' \
                                     | awk -F \"</?INSDSeq_organism>\" '{ printf($2); }'" % curr_id,
                                    shell=True).decode('utf-8')
                                curr_species = [curr_species_name]
                                j = 0
                                lineage_out = retrieve_lineages(curr_species[0])
                            multi_twin_dict["lineage"][curr_species[j]] = dict()
                            multi_twin_dict["phylum"][curr_species[j]] = dict()
                            multi_twin_dict["lineage"][curr_species[j]] = lineage_out[0]
                            multi_twin_dict["phylum"][curr_species[j]] = lineage_out[1]
                    j = j+1

                redirect_msg("      found [%s]" % curr_species[0])

                multi_twin_dict["species"][curr_family] = multi_twin_dict["species"][curr_family] + curr_species
                multi_twin_dict["id"][curr_family] = multi_twin_dict["id"][curr_family] + id

    f.close()


#############################################
# OUTPUT ####################################
#############################################
with open(multi_twin_in, mode='w') as f1:
    if not args.get_lineages:
        f1.write('#family_name\tspecies\n')
    else:
        f1.write('#family_name\tspecies\tphylum\tlineage\n')

    for i in range(0, len(fasta_files)):
        curr_fa = fasta_files[i]
        curr_family = family_names[i]
        list_species = multi_twin_dict["species"][curr_family]
        unique_species = sorted(set(list_species))
        list_ids = multi_twin_dict["id"][curr_family]

        for j in range(0, len(unique_species)):
            if unique_species[j]:
                if not args.get_lineages:
                    f1.write('%s\t%s\n' % (curr_family, unique_species[j]))
                else:
                    f1.write('%s\t%s\t%s\t%s\n' % (curr_family, unique_species[j],
                             multi_twin_dict["phylum"][unique_species[j]],
                             multi_twin_dict["lineage"][unique_species[j]]))

        curr_file = os.path.join(comprehensive_output_dir, curr_family + 'species_sequences_relationships.tsv')
        with open(curr_file, mode='w') as f2:
            if not args.get_lineages:
                f2.write('#id\tspecies\n')
            else:
                f2.write('#family_name\tspecies\tphylum\tlineage\n')
            for j in range(0, len(list_ids)):
                if not args.get_lineages:
                    f2.write('%s\t%s\n' % (list_ids[j], list_species[j]))
                else:
                    f2.write('%s\t%s\t%s\t%s\n' % (list_ids[j], list_species[j],
                                                   multi_twin_dict["phylum"][list_species[j]],
                                                   multi_twin_dict["lineage"][list_species[j]]))
        f2.close()

f1.close()

