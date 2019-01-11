#!/usr/local/bin/python3.5


import argparse
import logging
import os
import subprocess
import time
from glob import glob
import sys


def read_list_of_paths(path_list):
    k = 0;
    in_files = list()
    with open(path_list, mode='r') as f:
        for line in f:
            in_files.append(line.strip())
            k += 1
    f.close()
    return in_files


def create_logger(name, filename):
    try:
        os.remove(filename)
    except OSError:
        pass
    # Initialize the logger
    logger = logging.getLogger(name)
    logger.setLevel(logging.DEBUG)
    # Create a file handler
    handler = logging.FileHandler(filename, mode="a+")
    handler.setLevel(logging.DEBUG)
    # Create a logging format
    formatter = logging.Formatter('%(levelname)s:%(name)s:%(message)s')
    handler.setFormatter(formatter)
    # Add the handlers to the logger
    logger.addHandler(handler)
    return logger


def create_command(args):
    cmd = 'python3.5 \"%s\" -query \"%s\" -db_fa \"%s\" '\
          '-th %s -run %s -pident_thr %s -cov_thr %s -eval_thr %s '\
          '-min_size %s -max_size %s -test_all_chain %s '\
          '-mkdb_ %s -faa_split %s -wd %s' % \
          (args.isf_path, curr_fa, curr_db,
           args.th, args.run, args.pident_thr, args.cov_thr, args.eval_thr,
           args.min_size, args.max_size, args.test_all_chain,
           args.mkdb_, args.faa_split, curr_outdir)
    if args.nr_db:
        cmd = "%s -nr_db %s" % (cmd, args.nr_db)
    if args.diamond:
        cmd = "%s -diamond %s" % (cmd, args.diamond)
    else:
        cmd = "%s -blast_ %s" % (cmd, args.blast_)
    return cmd


parser = argparse.ArgumentParser(description='This script runs several instances of ISF')
# ISF_batch_run_args
parser.add_argument('--list_of_fasta_paths', dest='list_fa', type=str,
                    help='specify the path to the list of paths to input fasta files')
parser.add_argument('--list_of_db_paths', dest='list_db', type=str,
                    help='specify the path to the list of paths to target databases' 
                         '(the list can be only one path')
parser.add_argument('--output_dir', dest='output_dir', type=str,
                    help='specify the path to directory that will store all the ISF outputs')
parser.add_argument('--ISF_path', dest='isf_path', type=str,
                    help='specify the path to isf_main.py')
# ISF_args
parser.add_argument('-th', help='number of thread to use', type=int, required=True)
parser.add_argument('-run', help='number of run to do', type=int, default=10**6)
parser.add_argument('-pident_thr', help='threshold limit for identity default = 30.0, value > 0.0',
                    type=float, default=30.0)
parser.add_argument('-cov_thr', help='threshold limit fo coverage default = 80.0, value > 0.0',
                    type=float, default=80.0)
parser.add_argument('-eval_thr', help='threshold limit fo evalue default = 0.00001, value > 0.0',
                    type=float, default=0.00001)
parser.add_argument('-min_size', help=' %% of the shortest query size use to limit by size blasted sequence'
                                      'default 80', default=80, type=float)
parser.add_argument('-max_size', help=' %% of the longest query size use to limit by size blasted sequence'
                                      'default 120', default=120, type=float)
parser.add_argument('-test_all_chain', help='does we test conservation of alignment position, default previous',
                    action='store', choices=["no", "previous", "all"], default="previous")
parser.add_argument('-diamond', help="path to diamond binary if set the software will use diamond instead of blast.")
parser.add_argument('-blast_', help=" where is blast ? default use: blastp if not in PATH, indicate the path to exe",
                    default='blastp')
parser.add_argument('-mkdb_', help='where is makeblastdb ? default use: makeblastdb if not in PATH'
                                   ', indicate the path to exe',
                    default='makeblastdb')
parser.add_argument('-faa_split', help='where is exonerate fastasplit, default use: fastasplit',
                    default='fastasplit')
parser.add_argument('-nr_db', help='path to nr database, if set when the iteration finish \
results are blast against nr', type=str)
args = parser.parse_args()

#####################################################
# SETUP #############################################
#####################################################
# Read list of paths to input fasta and db
fa_paths = read_list_of_paths(args.list_fa)
db_paths = read_list_of_paths(args.list_db)
# Initialize logfile
general_log = create_logger('general_log', os.path.join(args.output_dir, "ISF_batch.log"))
# Create a summary directory that will eventually store all the gene families
summary_dir = os.path.join(args.output_dir, "00_summary")
os.makedirs(summary_dir, exist_ok=True)

#####################################################
# BATCH RUN #########################################
#####################################################
for i in range(0, len(fa_paths)):
    # define target database
    if len(db_paths) == 1:
        curr_db = db_paths[0]
    else:
        curr_db = db_paths[i]

    # retrieve the fasta file to be processed
    curr_fa = fa_paths[i]
    fa_name = os.path.splitext(os.path.basename(curr_fa))[0]

    # log name of the current record
    record_msg = "Executing ISF for \"%s\" (seq n°%s)" % (fa_name, i+1)
    general_log.info(record_msg)
    print(record_msg)

    # Create an output subdirectory for the current instance of ISF
    curr_outdir = os.path.join(args.output_dir, fa_name)
    os.makedirs(curr_outdir, exist_ok=True)

    # initialize logfile for the current instance of ISF
    current_log = create_logger('current_log', os.path.join(curr_outdir, "ISF.log"))

    # define ISF command line
    cmd = create_command(args)
    current_log.info(cmd)

    # run isf, compute time elapsed, and log everything !
    start_time = time.time()

    ISF_result = subprocess.Popen(args=cmd, shell=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in ISF_result.stdout:
        print(line)
        current_log.info(line)

    elapsed_time = time.time() - start_time

    ISF_out, ISF_err = ISF_result.communicate()
    if ISF_err:
        err = "* ISF generated the following error:\n%s\n" % ISF_err
        print(err)
        current_log.critical("The following error has been produced:\n%s" % ISF_err)
        general_log.critical("The following error has been produced for this query:\n%s" % ISF_err)
    else:
        general_log.info("Gene family successfully created in %s" %
                         time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
        # create a symbolic link of the comprehensive gene family into the summary directory
        summary_files = glob(os.path.join(curr_outdir, "sequence_found_and_bases.faa*"))
        for f in summary_files:
            os.symlink(f, os.path.join(summary_dir, fa_name + os.path.splitext(f)[1]))

    sys.stdout.flush()
