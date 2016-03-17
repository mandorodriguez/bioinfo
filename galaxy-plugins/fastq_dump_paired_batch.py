#!/bin/env python

import os
import sys
import argparse
import pdb
import subprocess
import shutil


"""
Program takes a file with a list of SRA identifiers and performs
a fastqdump on each and stores them in a directory.

"""

#-------------------------------------------------------------------------------

def fastq_dump(accession_id, directory):
    """
    This function will perform a fastq-dump on the given accession number and put the
    downloaded files in the given directory.

    """
    print "Doing a fastq-dump on accession %s:" % accession_id


    out = subprocess.Popen(["fastq-dump","--log-level","fatal","--split-3", "--accession", accession_id,
                            "--ncbi_error_report","never"],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.STDOUT)

    if not out.stdout is None:
        for line in iter(out.stdout.readline,b''):
            print line.rstrip()

        print "\n"
        
    # move to the given directory
    if directory != ".":

        if os.path.isfile("./%s_1.fastq" % accession_id):
            
            shutil("./%s_1.fastq" % accession_id, "%s/%s_R1.fastq" % (directory,accession_id))
            shutil("./%s_2.fastq" % accession_id, "%s/%s_R2.fastq" % (directory,accession_id))
            
        else:
            
            raise Exception("fastq-dump failed for %s. No files downloaded" % accession_id)
        
#-------------------------------------------------------------------------------


parser = argparse.ArgumentParser()

parser.add_argument("-a","--accession-file", type=str,
                    help="File with a list of accession identifiers")
parser.add_argument("-d","--directory", type=str,default=".",
                    help="")
parser.add_argument('-s', '--accession', nargs='*',
                    help="A list of space separated accession numbers")

args = parser.parse_args()


#-------------------------------------------------------------------------------


accession_ids = []
directory = args.directory

# We load any ids we get from the command line (easiest)

if not args.accession is None:
    accession_ids += args.accession 


# Now we read in ids from a file.
if not args.accession_file is None:
    try:

        with open(args.accession_file,"rU") as input_handle:

            for line in input_handle:

                accession_ids.append(line.rstrip())

    except Exception,e:
        print "Error reading accession ids from file %s" % args.accession_file

accession_ids = list(set(accession_ids))

print "Attempting to download from %d accession ids\n" % len(accession_ids)

# Now we just loop through and use the fastq function on each

for aid in accession_ids:

    try:
        fastq_dump(aid, directory)
    except Exception,e:
        print "Error downloading %s: %s" % (aid, e) 


print "Done!"

