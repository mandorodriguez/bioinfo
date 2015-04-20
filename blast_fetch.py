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
    """
    Parses a blast file into a list
    """
    b_list = []
    for qresult in SearchIO.parse(blast_file, 'blast-text'):

        b_list.append(qresult)
            
    return b_list

    

#-------------------------------------------------------------------------------

def _parse_snp_positions_dict(snp_file):
    """
    This function just loads a snp panel file into a dictionary
    that hashes the list of positions to the molecule.
    """
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

        return snp_positions
    
#---------------------------------------------------------------
   
def _parse_query_names(query_file):
    
    fasta_seqgen = SeqIO.parse(open(query_file), 'fasta')

    query_names = []
    
    for f in fasta_seqgen:
            
        query_names.append(f.name)

    return query_names
    
#-------------------------------------------------------------------------------

class SNP:
    """
    This SNP class just contains the molecule, position, and each contig file
    to process so that it can all be output on one line. The contigs are processed
    in the order of the list query_contig_files.
    """
    def __init__(self, molecule, pos, query_contigs, blast_data):


        self.molecule = molecule
        self.pos = pos
        self.query_id = "%s_%d_SUBSEQ" % (molecule, pos)

        self.query_contigs = [QueryContig(qc, self.query_id, blast_data) for qc in query_contigs]

    #---------------------------------------------------------------

    def to_string(self):

        return "\t".join( [self.molecule, str(self.pos)] + [qc.get_names() for qc in self.query_contigs] )
        
    
#-------------------------------------------------------------------------------


class QueryContig:

    #---------------------------------------------------------------

    def __init__(self, query_tuple, query_id, blast_data):


        self.qfile = query_tuple[0]
        
        self.names = query_tuple[1]

        self.query_id = query_id
        
        self.hits = []
        
        for name in self.names:

            bhit_tuple = self.get_blast_hits(name, query_id, blast_data)
            
            # if we get any blast hits that match on the tuple then we
            # append
            
            if len(bhit_tuple[1]) > 0:
                
                self.hits.append( bhit_tuple )

    #---------------------------------------------------------------
    
    def get_names(self):

        return ",".join([ht[0] for ht in self.hits])

    #---------------------------------------------------------------

    def get_positions(self):

        hit_starts = []

        for htuple in self.hits:
            
            for h in htuple[1]:

                hit_starts.append( "%s:%s" % (h.hit_id, h.hit_start) )

        return ",".join( hit_starts )
        
    #---------------------------------------------------------------

    def get_blast_hits(self, name, query_id, blast_data):
        """
        Gets all blast hits for this query contig name
        """

        blast_hits = []
        
        for query_result in blast_data:
        
            for hit in [h for h in query_result.hits if h.id == name]:

                for hsp in hit.hsps:

                    if hsp.hit_id == name and hsp.query_id == query_id:

                        blast_hits.append(hsp)


        return (name, blast_hits)

#-------------------------------------------------------------------------------

def get_header(query_contigs):

    return "\t".join(["molecule", "pos"] + [qc[0] for qc in query_contigs])

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
    snp_positions = _parse_snp_positions_dict(args.snp_panel)

    # loads the contigs with the filename into a list of tuples
    # so they preserve the order given. Only saving the base filename.
    query_contigs = [(os.path.splitext(os.path.basename(qf))[0], _parse_query_names(qf)) for qf in query_files]

    snps_w_data = []

    # Iterate over the molecule and each position in each molecule.
    for mol,positions in snp_positions.items():

        for pos in positions:
            
            snps_w_data.append( SNP(mol, pos, query_contigs, blast_data) )
        
    pdb.set_trace()


#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
