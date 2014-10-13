#!/usr/bin/env python



import argparse, os, sys
import pdb
import operator

#- just stuff I imported into other stuff
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Blast import NCBIXML
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIStandalone
from Bio import SearchIO
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA
from Bio.SeqFeature import FeatureLocation




#-------------------------------------------------------------------------------
# Function gets the maximum bitscore contained within a Hit object and its fragments
#-------------------------------------------------------------------------------
def get_max_id_and_bitscore(hit):
    hmax = max(hit.hsps, key=operator.attrgetter('bitscore'))

    
    return (hmax.hit_id, hmax.bitscore)

#-------------------------------------------------------------------------------

def filter_max_bitscore(qr):

    hit_maxes = [get_max_id_and_bitscore(h) for h in qr.hits]

    # get the max bitscore from the list of hit maxes
    has_max_bitscore = max(hit_maxes,key=operator.itemgetter(1))

    max_filter = lambda hsp: hsp.hit_id == has_max_bitscore[0] and hsp.bitscore == has_max_bitscore[1]

    filtered_qr = qr.hsp_filter(max_filter)

    return filtered_qr
    
#-------------------------------------------------------------------------------
# Going to do all of the work in this __main__ function
# 
#-------------------------------------------------------------------------------
def __main__():


    #Parse Command Line
    parser = argparse.ArgumentParser()


    #- a default argument for the blast file
    parser.add_argument('blast_file', type=str, help="A blast file to process")
    parser.add_argument('-o', type=str, help="A name for the output file")
    args = parser.parse_args()



    blast_data = []
    
    # check for a valid file
    if not os.path.isfile( args.blast_file ):
        sys.exit("\n%s is not a valid file!\n" % args.blast_file)

    else:

        # this sets the output file
        output_file = ""
        if args.o is None:
            output_file = "filtered_blast.txt"
        else:
            output_file = args.o
        
        #
        # Class object for what is being read in is:
        #   http://biopython.org/DIST/docs/api/Bio.SearchIO._model.query.QueryResult-class.html
        #
        num_total_qresults = 0
        num_filtered_qresults = 0
        
        for qresult in SearchIO.parse(args.blast_file, 'blast-text'):

            # if no hits we just keep going
            if len(qresult.hits) == 0:
                continue
            
            # this counts the total number of hits loaded
            num_total_qresults += sum([len(h) for h in qresult.hits])

            filtered_qresult = filter_max_bitscore(qresult)
                
            num_filtered_qresults += sum([len(h) for h in filtered_qresult.hits])

            blast_data.append(filtered_qresult)

            

        print "Filted out %d hits from %d total hits" % (num_total_qresults - num_filtered_qresults, num_total_qresults)

        output_handle = open(output_file, "w")
        
        count = SearchIO.write(blast_data, output_handle, 'blast-tab') # not working just yet
        print "Wrote out %s blast records to '%s'" % (count,output_file)

        output_handle.close()




#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()

