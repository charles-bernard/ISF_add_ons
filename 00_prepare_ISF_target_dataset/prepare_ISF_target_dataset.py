#!/usr/local/bin/python3.5

import argparse
import os
import re
import subprocess


def get_assemblies(assembly_directory):
    print(' * Retrieving name and number of proteins of assemblies...')
    assembly_dictionary = dict()
    files = os.listdir(assembly_directory)
    for file in files:
        if re.search(r'\.fa(a)?(sta)?$', file):
            assembly_name = re.match(r'^.*?_.*?\.[0-9]', file).group(0)
            assembly_dictionary[assembly_name] = dict()
            assembly_dictionary[assembly_name]['fasta'] = os.path.join(assembly_directory, file)
            nb_proteins = subprocess.check_output(['grep', '-c', '>', assembly_dictionary[assembly_name]['fasta']])
            assembly_dictionary[assembly_name]['nb_proteins'] = int(nb_proteins.strip())     
    return(assembly_dictionary)
    

def assembly_name_to_taxid(assembly_dictionary):
    print(' * Retrieving the taxid of each assembly...')
    assemblies_names = list(assembly_dictionary.keys())
    for name in assemblies_names:
        assembly_link = 'https://www.ncbi.nlm.nih.gov/assembly/' + name + '/?report=xml'
        taxid = subprocess.check_output([
        'curl %s | egrep \'SpeciesTaxid\' | egrep -o \'[0-9]+\'' % assembly_link ], shell = True)
        taxid = int(taxid.strip())
        assembly_dictionary[name]['taxid'] = taxid
    return(assembly_dictionary)
    

def taxid_to_lineages(assembly_dictionary):
    
    
    
def write_dictionary(assembly_dictionary, output_directory):
    output_file = os.path.join(output_di)
        
    


parser = argparse.ArgumentParser(description='This script prepares the target dataset of ISF')
parser.add_argument('--assembly_directory', '-a', dest='assembly_directory', type=str,
                    help='specify the path to the directory where the assemblies in fasta are located')
parser.add_argument('--output_directory', '-o', dest='output_directory', type=str,
                    help='specify the path to the directory where the target dataset of ISF will be located')
parser.add_argument('--dataset_name', '-n', dest='dataset_name', type=str,
                    help='specify the name of the target dataset')
args = parser.parse_args()

assembly_dictionary = get_assemblies(args.assembly_directory)
assembly_dictionary = assembly_name_to_taxid(assembly_dictionary)
