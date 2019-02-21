#!/usr/local/bin/python3.5

import sys
import argparse
import re
import os
import subprocess
import tempfile


def get_fasta_ids(fasta):
    ids = list()
    with open(fasta, mode="r") as f:
        for line in f:
            is_header = re.search(r'^>.*$', line)
            if is_header:
                ids.append(is_header.group(0)[1:])
    f.close()
    return ids


def write_ids(ids, output_file):
    with open(output_file, mode="w") as f:
        for id in ids:
            f.write("%s\n" % id)
    f.close()


def get_remote_seqs_from_iteration1(ids, edge_file, output_file):
    f = open(output_file, mode="w")
    for id in ids:
        subprocess.call(['egrep', '^' + id, edge_file], stdout=f)
    f.close()


def filter_edges(edges_it1_file, query_ids_file, output_file, thresh):
    f = open(output_file, mode="w")
    subprocess.call(['awk', "-f", 'filter_ISF_edges.awk', "-v",
                     "thresh=" + thresh, query_ids_file,
                     edges_it1_file], stdout=f)
    f.close()


parser = argparse.ArgumentParser()
parser.add_argument('-e', dest="edge_file", type=str,
                    help='specify the path to the edge_file')
parser.add_argument('-f', dest="family_file", type=str,
                    help='specify the path to the initial family (fasta)')
parser.add_argument('-o', dest='output_dir', type=str,
                    help='specify the path to the output directory')
parser.add_argument('-t', dest='pident_thr', type=str,
                    help='pourcentage_identity_threshold')
args = parser.parse_args()

if not os.path.exists(args.output_dir):
    os.makedirs(args.output_dir)

ids = get_fasta_ids(args.family_file)
query_ids_file = os.path.join(args.output_dir, "query_ids.txt")
write_ids(ids, query_ids_file)

edges_it1_file = os.path.join(args.output_dir, "edges_first_iteration.tsv")
get_remote_seqs_from_iteration1(ids, args.edge_file, edges_it1_file)

filt_edges_file = os.path.join(args.output_dir, "filtered_edges_first_iteration.tsv")
filter_edges(edges_it1_file, query_ids_file, filt_edges_file, args.pident_thr)
