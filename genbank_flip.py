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


old_string = """
     gene            %s..%s
                     /locus_tag="%s"
     CDS             %s..%s
                     /locus_tag="%s"
                     /product="Non-coding region"
"""
#-------------------------------------------------------------------------------

gene_string = """
      gene            %s..%s
      /locus_tag="%s"
"""



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

    output_file = args.outfile
    
    genbank_rec = parse_genbank_file(genbank_file)

    
    source_start = genbank_rec.features[0].location.start
    source_end = genbank_rec.features[0].location.end

    gene_coords = []
    
    next_start = 0

    for feature in genbank_rec.features[1:]:

        if feature.type == "gene":
            
            this_start = int(feature.location.start)
            this_end = int(feature.location.end)

            if this_start > next_start:

                gene_coords.append( (next_start,this_start) ) 

                next_start = this_end-1
        
                
    gene_coords.append( (next_start, int(source_end)) )
        

    # now we read in the genbank file and prepare to reoutput it.
    genbank_lines = []
    
    with open(genbank_file, "rU") as input_handle:

        for line in input_handle:
            genbank_lines.append(line)

    
    # now we output the genbank
    with open(output_file,"w") as output_handle:

        for line in genbank_lines:
            output_handle.write(line+"\n")

            if "gene" in line:
                break


        for coords in gene_coords:
            output_handle.write( gene_string % (coords[0],coords[1]) )

        start_printing = False
        for line in genbank_lines:

            if "ORIGIN" in line:
                start_printing = True

            if start_printing:
                output_handle.write(line+"\n")
                
    
    pdb.set_trace()


#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
