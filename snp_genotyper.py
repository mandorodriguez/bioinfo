#!/usr/bin/env python

import argparse
import pdb

#*******************************************************************************
# Based on the original code here:
#   https://github.com/Henshina/Genotyper/blob/master/Geno1.py
#*******************************************************************************




#*** SNP class object ******************************************
# Bulding functionality into this class object
# to save a bunch of work.
#***************************************************************
class SNP:
    
    #--------------------------------------------
    def __init__(self, qindexes, line, ref_base_index=3):

        self.snp_data = line.split('\t')

        self.qindexes = qindexes

        # just the ref base
        self.ref_base = self.snp_data[ref_base_index]

        # collect the qbases from the SNP
        self.qbases = [self.snp_data[qi] for qi in self.qindexes]

        # get the ref counts from the SNP
        self.ref_counts = self.get_ref_counts()

        # Convert our ref counts to a string as an ID for hashing.
        self.count_id = ''.join([str(b) for b in self.ref_counts])
        
    #--------------------------------------------

    def ref_matches(self):
        """
        This method just returns the number of times the refbase
        matches all qbases.
        """
        return self.qbases.count( self.ref_base )

    #--------------------------------------------

    def get_ref_counts(self):
        """
        performs the walk along the qbases and does the count if
        refbase bases match the qbase and resets if it doesn't.
        Only need to call this once on init since the SNP is immutable.
        """
        base_count = []
        count = 0
        for qb in self.qbases:

            if qb == self.ref_base:
                count += 1
            else:
                count = 0

            base_count.append(count)

        return base_count

#*** end SNP class *********************************************


#***** Genotyper class object ***********************************

class Genotyper:
    
    #--------------------------------------------
    def __init__(self):
        pass
    
#****************************************************************




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
    print "Loaded %d SNPs from file '%s' with %d qbases per SNP.\n" % (len(snp_objects), table_file, len(qindexes))
    
    return (header, qindexes, snp_objects)


#************* Genotyper class ******************************************


#-------------------------------------------------------------------------------
# Main function call
def __main__():

    parser = argparse.ArgumentParser()
    parser.add_argument("snp_table", type=str,
                        help="The snp table to input")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="snp_genotype_out.txt")
    args = parser.parse_args()

    # open the snp table and load it into some data types.
    
    header, qindexes, snp_objects = load_table(args.snp_table)
    
p    pdb.set_trace()

#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
