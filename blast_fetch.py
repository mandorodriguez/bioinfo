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

def _parse_snp_positions_dict(snp_file):

    num_snps = 0

    snp_positions = {}
    
    with open(snp_file, "rU") as input_handle:

        for line in input_handle:

            if not line.isspace():

                snp_txt = line.split()

                if len(snp_txt) != 2:
                    print "Error parsing entry \"%s\", not in form 'locus position'" % line
                    continue

                else:

                    if snp_positions.get(snp_txt[0]) is None:

                        try:
                            snp_positions[snp_txt[0]] = [int(snp_txt[1])]
                        except Exception,e:
                            raise Exception("Error with snp: %s" % e)

                    else:

                        # if this refpos is in the positions for this locus, when we ignore it,
                        # otherwise It needs to be added.
                        try:
                            ref_pos_value = int(snp_txt[1])
                        except Exception, e:
                            print "Can't add '%s' to locus %s: %s" % (snp_txt[1], snp_txt[0], e)
                            continue
                                    
                        if not ref_pos_value in snp_positions[snp_txt[0]]:

                               
                            snp_positions[snp_txt[0]].append(ref_pos_value)

                        else:

                            print "SNP position %s is already in %s, skipping" % (snp_txt[1], snp_txt[0])
                    
                num_snps += 1

            

        print "Parsed %d snp positions from file '%s'" % (num_snps,snp_file)
        
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

    parser.add_argument("-b", "--blast", type=str,
                        help="A raw Blast text file")
    parser.add_argument("-s", "--snp_panel", type=str,
                        help="A SNP panel file")
    parser.add_argument('-q', '--query', nargs='*',
                        help="list of query files separated by spaces.")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="blast_fetch.txt")


    args = parser.parse_args()

    output_file = args.outfile
    blast_file =args.blast
    query_files = args.query

    blast_data = _parse_blast_list(blast_file)

    query_contigs = [QueryContig(qf, blast_data) for qf in query_files]
    
    pdb.set_trace()


#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
