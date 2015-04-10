#!/usr/bin/env python

import argparse, os, sys
import pdb


import operator

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIStandalone
from Bio import SearchIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqFeature import FeatureLocation





#-------------------------------------------------------------------------------
#
# full example for parsing blast files here:
#   http://bugzilla.open-bio.org/attachment.cgi?id=293&action=view
#
def _parse_blast_list(blast_file):

    b_list = []
    for qresult in SearchIO.parse(blast_file, 'blast-text'):

        b_list.append(qresult)
            
    return b_list

    
#-------------------------------------------------------------------------------

def _parse_blast_dict(blast_file):
    """
    Loads a blast file into a dictionary structure.
    """
    b_dict = {}
    for qresult in SearchIO.parse(blast_file, 'blast-text'):

        if b_dict.has_key(qresult.id):
            print "Blast read collision on id '%s' in file '%s'" % (qresult.id, blast_file)

        b_dict[qresult.id] = qresult
            
    return b_dict

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Main function call
def __main__():

    parser = argparse.ArgumentParser()
    #parser.add_argument("snp_table", type=str,
    #                    help="The snp table to input")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="blast_fetch.txt")
    parser.add_argument("-b", "--blast", type=str,
                        help="A raw Blast text file")


    args = parser.parse_args()

    output_file = args.outfile
    blast_file =args.blast

    blast_data = _parse_blast_dict(blast_file)

    pdb.set_trace()


#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
