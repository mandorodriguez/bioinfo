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
    parser.add_argument("-p", "--prefix", type=str, 
                        help="Prefix for the converted genes", default="NC", metavar="NC")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="flipped.gb",metavar="flipped_genbank.gb")
    parser.add_argument("-m", "--min", type=str,
                        help="Minimum Length to allow for a gene", default=0,metavar="MIN_LENGTH")
    parser.add_argument("genbank_file", type=str,
                        help="The genbank file to flip", metavar="file.gb")

    args = parser.parse_args()

    genbank_file = args.genbank_file

    output_file = args.outfile
    
    genbank_rec = parse_genbank_file(genbank_file)

    gene_prefix = args.prefix

    min_length = args.min
    
    source = genbank_rec.features[0]
    source_start = genbank_rec.features[0].location.start
    source_end = genbank_rec.features[0].location.end

    gene_coords = []
    
    next_start = 0

    # collect what will be the new gene coords
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
    genbank_text = ""
    with open(genbank_file, "rU") as input_handle:

        for line in input_handle:
            genbank_lines.append(line)

    genbank_text = "".join(genbank_lines)

    genbank_parts = genbank_text.split("FEATURES")

    genbank_prefix = genbank_parts[0]

    genbank_source_seq = genbank_parts[1].split("ORIGIN")[1]

    nc_gene_num = 0


    # now we output the genbank

    with open(output_file,"w") as output_handle:

        output_handle.write(genbank_prefix)

        output_handle.write("FEATURES             Location/Qualifiers\n")

        output_handle.write(format(' ',"5") + "source" + format(' ',"11") + "%s..%s" % (source_start+1,source_end+1) + "\n")

        
        for key,val in source.qualifiers.iteritems():
            output_handle.write(format(' ',"22")+"/%s=\"%s\"\n" % (key, ",".join(val)))

        for coords in gene_coords:
            
            if coords[1]-coords[0] > min_length:
                output_handle.write( gene_string % (coords[0]+1,coords[1]+1,"%s_%s" % (gene_prefix,format(nc_gene_num,"06"))) )
                nc_gene_num += 1

        output_handle.write("ORIGIN\n")
        output_handle.write(genbank_source_seq)
                
    

    print "Flipped genbank is written to %s" % output_file
    
#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
