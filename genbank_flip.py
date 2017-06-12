#!/usr/bin/env python

__version__="1.0"

"""
Name: genbank_flip.py

Description : Changes the non-coding regions of the genbank to genes

Author: Mando Rodriguez

"""
import argparse, os, sys
import pdb

from Bio import SeqIO
from Bio.Seq import Seq

#-------------------------------------------------------------------------------
def parse_genbank_file(genbank_file):
    
    genbank_recs = []

    genbank_input_handle = open(genbank_file, "rU")
    for gbrec in  SeqIO.parse(genbank_input_handle, "genbank"):

        genbank_recs.append(gbrec)

    genbank_input_handle.close()

    print "Read in %i genbank records from file %s" % (len(genbank_recs), genbank_file)

    return genbank_recs[0]

#-------------------------------------------------------------------------------
# Main function call
def __main__():

    parser = argparse.ArgumentParser()
    parser.add_argument("genbank_file", type=str,
                        help="The genbank file to flip", metavar="file.gb")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="flipped.gb",metavar="flipped_genbank.gb")

    args = parser.parse_args()

    genbank_file = args.genbank_file

    if args.outfile == "flipped.gb":
        output_file = "flipped_%s" % genbank_file
    
    genbank_rec = parse_genbank_file(genbank_file)

    

    pdb.set_trace()


#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
