#!/usr/bin/env python

import argparse
import pdb

#*** SNP class object ******************************************
# Bulding functionality into this class object
# to save a bunch of work.
#***************************************************************

class SNP:
    
    #--------------------------------------------
    def __init__(self, qindexes, line, ref_base_index=3):

        self.snp_data = line.split('\t')

        self.qindexes = qindexes

        self.ref_base = self.snp_data[ref_base_index]
        
        self.qbases = [self.snp_data[qi] for qi in self.qindexes]

#*** end SNP class *********************************************


#-------------------------------------------------------------------------------

def load_table(table_file):

    infile = open(table_file, 'r')

    table_data = infile.readlines()

    infile.close()

    # The header should always be the first line
    header = table_data[0]

    # collect the qindexes from the header with this loop within a list
    qindexes = [indx for indx,colname in enumerate(header.split('\t')) if "qbase:" in colname]


    # put each line representing a SNP into a SNP object which was declared earlier
    snp_objects = []
    
    for snp_line in table_data[1:]:

        snp_objects.append( SNP(qindexes, snp_line) )


    # returning everything in a tuple
    
    return (header, qindexes, snp_objects)

#-------------------------------------------------------------------------------
# Main function call
def __main__():

    parser = argparse.ArgumentParser()
    parser.add_argument("snp_table", type=str,
                        help="The snp table to input")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="snp_genotype_out.txt")
    args = parser.parse_args()

    # open the snp table
    
    header, qindexes, snp_objects = load_table(args.snp_table)
    
    pdb.set_trace()

#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
