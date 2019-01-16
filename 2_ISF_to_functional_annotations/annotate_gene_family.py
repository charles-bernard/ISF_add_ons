#!/usr/local/bin/python3.5

import argparse
import subprocess
import logging
import os
import sys
import time
import re


def run_rpsblast(args, fa, fa_basename):
    cmd = 'rpsblast+ -query \"%s\" -db %s '\
          '-out %s -evalue %f -outfmt 6' %\
          (fa, os.path.basename(os.path.normpath(args.db)).split("_")[0],
           os.path.join(args.outdir, fa_basename + '_rpsblast.out'), args.e_value)

    logging.info('Executing rpsblast+ as: %s' % cmd)
    print('Executing rpsblast+ as: %s' % cmd)

    start_time = time.time()
    rpsblast_result = subprocess.Popen(args=cmd, shell=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    elapsed_time = time.time() - start_time

    rpsblast_out, rpsblast_err = rpsblast_result.communicate()
    if rpsblast_err:
        err = "* rpsblast+ generated the following error:\n%s" % rpsblast_err
        logging.critical(err)
        print(err)
    else:
        msg = '* rpsblast+ successfully ran in %s' % \
              time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        logging.info(msg)
        print(msg)


def concat_seq_with_no_matches_to_rspblast_out(args, fa, fa_basename, script_dir):
    cmd = 'awk -f \"%s\" \"%s\" \"%s\" > \"%s\"' %\
          (os.path.join(script_dir, 'Subscripts', 'left_join_fa_cdd.awk'),
           fa, os.path.join(args.outdir, fa_basename + '_rpsblast.out'),
           os.path.join(args.outdir, fa_basename + '_comprehensive_rpsblast.out'))

    logging.info('Concatenate sequences that have no cdd matches with rpsblast output : %s' % cmd)
    print('Concatenate sequences that have no cdd matches with rpsblast output : %s' % cmd)

    start_time = time.time()
    concat_result = subprocess.Popen(args=cmd, shell=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elapsed_time = time.time() - start_time

    concat_out, concat_err = concat_result.communicate()
    if concat_err:
        err = "* concatenation generated the following error:\n%s" % concat_err
        logging.critical(err)
        print(err)
    else:
        msg = '* concatenation successfully ran in %s' % \
              time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        logging.info(msg)
        print(msg)


def convert_cddid_to_cogid(args, fa_basename, script_dir):
    cmd = 'awk -v comprehensive_out=\"%s\" -v summary_out=\"%s\" '\
          '-f \"%s\" \"%s\" \"%s\" \"%s\" \"%s\"' %\
          (os.path.join(args.outdir, 'COG_comprehensive_output.txt'),
           os.path.join(args.outdir, 'COG_enrichment.txt'),
           os.path.join(script_dir, 'Subscripts', 'cdd_2_cog.awk'),
           os.path.join(script_dir, '../required_files', 'fun.txt'),
           os.path.join(args.outdir, fa_basename + '_comprehensive_rpsblast.out'),
           os.path.join(script_dir, '../required_files', 'cddid.tbl'),
           os.path.join(script_dir, '../required_files', 'whog'))

    logging.info('Executing cdd 2 cog as: %s' % cmd)
    print('Executing cdd 2 cog as: %s' % cmd)

    start_time = time.time()
    cdd2cog_result = subprocess.Popen(args=cmd, shell=True,
                                       stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elapsed_time = time.time() - start_time

    cdd2cog_out, cdd2cog_err = cdd2cog_result.communicate()
    if cdd2cog_err:
        err = "* cdd2cog+ generated the following error:\n%s" % cdd2cog_err
        logging.critical(err)
        print(err)
    else:
        msg = '* cdd2cog successfully ran in %s' % \
              time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
        logging.info(msg)
        print(msg)


parser = argparse.ArgumentParser(
    description='This script runs rpsblast with the ISF output gene family in fasta '
                'against COG db in order to assign functions to all the representative '
                'sequences of the family')
parser.add_argument('-o', '--output_dir', dest="outdir", type=str,
                    help='specify the path to the output directory')
parser.add_argument('--input_file', dest='fa_file', type=str,
                    help='specify the path to the fasta file corresponding to the gene family '
                         '(the filename will be used to name the gene family)')
parser.add_argument('--input_dir', dest='fa_dir', type=str,
                    help='(alternative to --input_file), specify the path to the directory '
                         'that contains several gene families in fasta (the filename of '
                         'each fasta will be used to name the corresponding gene family)')
parser.add_argument('--db', dest='db', type=str,
                    help='specify the path to the target database. '
                         'You may want to download the Conserved Domain Database '
                         '(wget -r ftp://ftp.ncbi.nlm.nih.gov/pub/mmdb/cdd/little_endian) '
                         'and specify the path to the directory of one of its child databases '
                         '(COG, Pfam, Tigr etc...) e.g: -db \"~/cdd/little_endian/COG_LE/\"')
parser.add_argument('--evalue', dest='e_value', type=float,
                    help='Expectation value (E) threshold for saving hits', default=0.1)
parser.add_argument('--cog_enrichment', dest='cog_enrich', action='store_true', default=False,
                    help='weither you want to compute the enrichment in each COG functional category '
                         'of the gene family (only relevant if --db was set to */COG_LE/)')
args = parser.parse_args()


##########################################
# SETUP ##################################
##########################################
# get the directory of the executing script
script_dir = os.path.dirname(os.path.realpath(sys.argv[0]))

# get working directory
working_directory = os.getcwd()

# save the abspath of the output directory
args.outdir = os.path.abspath(args.outdir)

# initialize logfile
logging.basicConfig(filename=os.path.join(args.outdir, 'logfile.log'), level=logging.DEBUG)

# It is important to change to the db directory because rpsblast+
# doesn't handle absolute path to db as arg
logging.info('Changing into %s' % args.db)
print('Changing into %s' % args.db)
os.chdir(args.db)

# Build the list of paths to fasta
fasta_files = list()
if args.fa_file:
    fasta_files.append(args.fa_file)
else:
    dir_files = os.listdir(args.fa_dir)
    for dir_file in dir_files:
        if re.search(r'\.fa(a)?(sta)?$', dir_file):
            fasta_files.append(os.path.abspath(dir_file))

##########################################
# RUN RPSBLAST+ & CONVERT CDD TO COG #####
##########################################
for fa in fasta_files:
    fa_basename = os.path.basename(fa).split(".")[0]
    run_rpsblast(args, fa, fa_basename)
    concat_seq_with_no_matches_to_rspblast_out(args, fa, fa_basename, script_dir)
    if args.cog_enrich:
        convert_cddid_to_cogid(args, fa_basename, script_dir)


