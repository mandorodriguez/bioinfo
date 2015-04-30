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

from multiprocessing.dummy import Pool as ThreadPool

"""

Gets the true position from the blast hit via the snp position.


"""


#************************************************************************************
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

    

#*************************************************************************************

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

    
#********************************************************************************

class SNP:
    """
    This SNP class just contains the molecule, position, and each contig file
    to process so that it can all be output on one line. The contigs are processed
    in the order of the list query_contig_files.
    """
    def __init__(self, molecule, pos, query_files, blast_data, flanking_bases=20):


        self.molecule = molecule
        self.pos = pos
        self.query_id = "%s_%d_SUBSEQ" % (molecule, pos)

        self.query_contigs = [QueryContig(qf, self.query_id, blast_data, flanking_bases) for qf in query_files]

    #---------------------------------------------------------------

    def to_string(self):

        return "\t".join( [self.molecule, str(self.pos)] + [qc.get_ids_w_start() for qc in self.query_contigs] )
        

#-------------------------------------------------------------------------------

def compute_snp(st):
    """
    This will compute a snp for the blast fetch in threaded mode.
    Nothing fancy. Work is all in the SNP class.
    """
    try: 
        return SNP(st[0], st[1], st[2], st[3]) 
    except Exception,e:
        err_str= "Unable to compute snp entry for %s position %s: %s" % (st[0], st[1], e)
        print err_str
        raise Exception(err_str)

#***********************************************************************************


class QueryContig:

    #---------------------------------------------------------------

    def __init__(self, query_file, query_id, blast_data, flanking_bases=20):


        self.qfile = query_file.metafilename
        
        self.names = query_file.names

        self.query_id = query_id

        self.flanking_bases = flanking_bases
        
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

    def get_ids_w_start(self):

        hit_ids = []
        hit_starts = []
        true_snp_pos = []
        
        for htuple in self.hits:
            
            for h in htuple[1]:

                hit_ids.append ( h.hit_id )
                hit_starts.append( str(h.hit_start) )
                true_snp_pos.append( self.snp_offset(h) )

        return "%s\t%s" % (",".join( hit_ids ), ",".join( true_snp_pos ) )
        
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
    def snp_offset(self, hit):
        """
        This calculates the offset to get the true snp position in the
        blast hit.
        """


        offset = self.flanking_bases - hit.query_start + 1

        # here we slide over if we have any indels on the left hand side of the snp
        if str(hit.query.seq).count('-') > 0:


            indel_indexes = [i for i,k in enumerate(str(hit.query.seq)) if k == '-']

            for i in indel_indexes:
                if i < offset:
                    offset += 1
                            
        if offset > hit.hit_end or offset < 0:
            return "Error[Q(%s,%s), H(%s,%s)]" % (hit.query_start, hit.query_end, hit.hit_start, hit.hit_end)
            #offset = -1

        # we add one to the offset just before returning so we have the position from start 1
        return str(hit.hit_start + offset+1)


#*********************************************************************************

class QueryFile:
    """
    Just a class to keep track of the query contig files and the names
    that they load.
    """
    #--------------------------------------------------------------------------

    def __init__(self, query_filename):


        self.filename = ""
        self.metafilename = ""
        
        lineparts = query_filename.split(':')

        if len(lineparts) == 2:

            self.filename = lineparts[0]
            self.metafilename = lineparts[1].replace(' ','_')

        elif len(lineparts) == 1:

            self.filename = lineparts[0]

            try:
                self.metafilename = os.path.splitext(os.path.basename(lineparts[0]))[0]
            except Exception:
                self.metafilename = lineparts[0]

        else:
                    
            raise Exception("Can't parse '%s', not a valid query contig FASTA filename" % line)


        self.names = self._parse_query_names(self.filename)

    #--------------------------------------------------------------------------

   
    def _parse_query_names(self,query_file):
    
        fasta_seqgen = SeqIO.parse(open(query_file), 'fasta')

        query_names = []
    
        for f in fasta_seqgen:
            
            query_names.append(f.name)

        return query_names
    
#*******************************************************************************

def get_header(query_files):

    return "\t".join(["molecule", "pos"] + ["hit_id:%s\thit_pos:%s" % (qf.metafilename,qf.metafilename) for qf in query_files])

#-------------------------------------------------------------------------------
# Main function call
def __main__():

    parser = argparse.ArgumentParser()

    parser.add_argument("-b", "--blast", type=str,
                        help="A raw Blast text file")
    parser.add_argument("-s", "--snp-panel", type=str,
                        help="A SNP panel file")
    parser.add_argument('-q', '--query', nargs='*',
                        help="list of query files separated by spaces.")
    parser.add_argument("-n", "--num-threads", type=int,
                        help="Number of threads to use", default=0)
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="blast_fetch.txt")


    args = parser.parse_args()

    output_file = args.outfile
    blast_file =args.blast
    query_files = args.query
    num_threads = args.num_threads

    blast_data = _parse_blast_list(blast_file)
    snp_positions = _parse_snp_positions_dict(args.snp_panel)

    # loads the contig files into a list of objects. The order
    # of these will be the processing and output table order.
    query_files = [QueryFile(qf) for qf in query_files]

    snps_w_data = []

    # here we run things in threaded mode to speed things up if the option is given
    if num_threads > 0:

        # iterate over each mol,position and create input tuples
        inputs = []
        
        for mol,positions in snp_positions.items():

            print "Setting up %d SNP jobs for molecule '%s'" % (len(positions),mol)

            for pos in positions:
            
                inputs.append( (mol, pos, query_files, blast_data) )


        print "Running %d individual input jobs for processing on %d threads..." % (len(inputs), num_threads)
            
        pool = ThreadPool(num_threads)
        pool_output = pool.map(compute_snp, inputs)
        pool.close()
        pool.join()
            
        print "Collecting results..."

        inputs = None # not needed anymore

        snps_w_data = [r for r in pool_output]

    else:
        
        # Iterate over the molecule and each position in each molecule.
        for mol,positions in snp_positions.items():

            for pos in positions:
            
                snps_w_data.append( SNP(mol, pos, query_files, blast_data) )


    with open(output_file, 'w') as output_handle:
        
        output_handle.write( get_header(query_files) + "\n")
        
        for s in snps_w_data:
            output_handle.write(s.to_string() + "\n")


    print "Done writing table to '%s'" % output_file

#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
