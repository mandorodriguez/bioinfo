#!/usr/bin/env python

import argparse, os, shutil, sys
import pdb

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

"""

Just opens a fasta file and checks to see if a sequence isn't a multiple
of 3. If not it prints the ID and description to stdout.
"""


#-------------------------------------------------------------------------------

def __main__():

    #Parse Command Line
    parser = argparse.ArgumentParser(description="""Will parse a fasta file and print out the id
    and description of any fasta records that are not a multiple of 3.

    """)

    parser.add_argument("fasta_file", type=str,
                        help="Multifasta file input",metavar="my_file.fasta")
    

    args = parser.parse_args()

    fasta_handle = open(args.fasta_file,"rU")

    for fasta_record in SeqIO.parse(fasta_handle, "fasta"):

        if (len(fasta_record.seq) % 3) != 0:

            print "%s\t%s" % (fasta_record.id, fasta_record.description)

    fasta_handle.close()

            
if __name__=="__main__": __main__()
