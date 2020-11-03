#!/usr/bin/env python3

__version__="1.3"

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

    genbank_input_handle = open(genbank_file, "r")
    for gbrec in  SeqIO.parse(genbank_input_handle, "genbank"):

        genbank_recs.append(gbrec)

    genbank_input_handle.close()

    print("Read in %i genbank records from file %s" % (len(genbank_recs), genbank_file))

    return genbank_recs[0]

#-------------------------------------------------------------------------------

def in_gene(n,coords,verbose=False):


    for c in coords:

        # if n[0] >= c[0] and n[0] <= c[1]:

        #     print("first coordinate in {} < {} < {}".format(c[0],n[0],c[1]))
        #     return True


        # if n[1] >= c[0] and n[1] <= c[1]:

        #     print("second coordinate in {} < {} < {}".format(c[0],n[1],c[1]))
        #     return True


        if (n[0] >= c[0] and n[0] <= c[1]) or (n[1] >= c[0] and n[1] <= c[1]):

            if verbose:
                print_coord(n,c)

            return True

    
    return False

#-------------------------------------------------------------------------------

def print_coord(n,c):

    if n[0] == c[0]:
        print("First position {},{} is equal to starting coord {},{}".format(n[0],n[1],c[0],c[1]))
    elif n[1] == c[0]:
        print("Second position {},{} is equal to starting coord {},{}".format(n[0],n[1],c[0],c[1]))
    elif n[0] == c[1]:
        print("First position {},{} is equal to ending coord {},{}".format(n[0],n[1],c[0],c[1]))
    elif n[1] == c[1]:
        print("Second position {},{} is equal to ending coord {},{}".format(n[0],n[1],c[0],c[1]))
    elif (n[0] > c[0] and n[0] < c[1]) and (n[1] > c[0] and n[1] < c[1]):
        print("Both positions {},{} are between coords {},{}".format(n[0],n[1],c[0],c[1]))
    elif n[0] > c[0] and n[0] < c[1]:
        print("First position {},{} is between coords {},{}".format(n[0],n[1],c[0],c[1]))
    elif n[1] > c[0] and n[1] < c[1]:
        print("Second position {},{} is between coords {},{}".format(n[0],n[1],c[0],c[1]))
    else:
        print("Coords aren't intergenic {} - {},{} - {}".format(n[0],c[0],c[1],n[1]))


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
    parser.add_argument("-v", "--verbose", action="store_true", 
                        help="Provides more output to the screen")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="flipped.gb",metavar="flipped_genbank.gb")
    parser.add_argument("-m", "--min", type=int,
                        help="Minimum Length to allow for a gene", default=0,metavar="MIN_LENGTH")
    parser.add_argument("genbank_file", type=str,
                        help="The genbank file to flip", metavar="file.gb")

    args = parser.parse_args()

    genbank_file = args.genbank_file

    output_file = args.outfile
    
    genbank_rec = parse_genbank_file(genbank_file)

    gene_prefix = args.prefix

    verbose = args.verbose

    min_length = int(args.min)
    
    source = genbank_rec.features[0]
    source_start = genbank_rec.features[0].location.start
    source_end = genbank_rec.features[0].location.end

    noncoding_coords = []
    gene_coords = []
    new_coords = []

    next_start = 0

    # collect all the current gene coord positions
    for feature in genbank_rec.features[1:]:

        if feature.type == "gene":

            gene_coords.append( (int(feature.location.start), int(feature.location.end)) )

    # sort the coords by first coord just in case they're out of order (normally not the case)
    gene_coords.sort(key=lambda tup: tup[0])


    # Loop through the coords and collect the 'gaps' in between them
    for gc in gene_coords:

        this_start = gc[0]
        this_end = gc[1]

        if this_start > next_start:

            noncoding_coords.append( (next_start,this_start-1) ) 

            next_start = this_end+1
        

    noncoding_coords.append( (next_start, int(source_end)) )
        

    # filter out any overlaps here

    for ncc in noncoding_coords:

        if not in_gene(ncc,gene_coords,verbose):
            new_coords.append(ncc)

    # now we read in the genbank file and prepare to reoutput it.
    genbank_lines = []
    genbank_text = ""
    with open(genbank_file, "r") as input_handle:

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
        
        for key,val in source.qualifiers.items():
            output_handle.write(format(' ',"22")+"/%s=\"%s\"\n" % (key, ",".join(val)))

        for coords in new_coords:
            
            if coords[1]-coords[0] >= min_length:

                this_coord_1 = coords[0]+1
                this_coord_2 = coords[1]+1
                
                output_handle.write( gene_string % (this_coord_1,this_coord_2,"%s_%s" % (gene_prefix,format(nc_gene_num,"06"))) )
                nc_gene_num += 1

        output_handle.write("ORIGIN\n")
        output_handle.write(genbank_source_seq)
                
    

    print("Flipped genbank is written to %s" % output_file)
    
#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
