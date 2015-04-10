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


"""

Gets the true position from the blast hit via the snp position.


-----

filename:molecule filename:pos ... filename:molecule filename:pos

"""


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

class QueryContig:

    #---------------------------------------------------------------

    def __init__(self, qfile, blast_data):


        self.qfile = qfile
        
        self.names = self._parse_query_names(qfile)

        self.hits = []
        
        for name in self.names:
            
            self.hits.append( self.get_blast_hits(name, blast_data) )

        pdb.set_trace()
    #---------------------------------------------------------------
   
    def _parse_query_names(self, query_file):
    
        fasta_seqgen = SeqIO.parse(open(query_file), 'fasta')

        query_names = []
    
        for f in fasta_seqgen:

            query_names.append(f.name)

        return query_names

    #---------------------------------------------------------------

    def get_blast_hits(self, name, blast_data):
        """
        Gets all blast hits for this query contig name
        """

        blast_hits = []
        
        for query_result in blast_data:
        
            for hit in [h for h in query_result.hits if h.id == name]:

                for hsp in hit.hsps:

                    if hsp.hit_id == name:

                        blast_hits.append(hsp)

        qname = None
        if len(blast_hits) > 0:
            qname = blast_hits[0].query_id
            
        return (name, qname, blast_hits)

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
    parser.add_argument('-q', '--query', nargs='*', help="list of query files separated by spaces.")


    args = parser.parse_args()

    output_file = args.outfile
    blast_file =args.blast
    query_files = args.query

    blast_data = _parse_blast_list(blast_file)

    query_contigs = [QueryContig(qf, blast_data) for qf in query_files]
    
    pdb.set_trace()


#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
