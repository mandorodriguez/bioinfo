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

        # molecule is at index 0
        self.molecule = self.snp_data[0]
        
        self.qindexes = qindexes

        # just the ref base
        self.ref_base = self.snp_data[ref_base_index]

        # collect the qbases from the SNP
        self.qbases = [self.snp_data[qi] for qi in self.qindexes]

        # get the ref counts from the SNP
        self.ref_counts = self.get_ref_counts()

        # Convert our ref counts to a string as an ID for hashing.
        self.count_id = ''.join([str(b) for b in self.ref_counts])

        # A class variable for the assigned encoding to set later
        self.code = None

    
    #--------------------------------------------

    def first_half(self):
        return self.snp_data[0:self.qindexes[0]]

    #--------------------------------------------

    def second_half(self):
        return self.snp_data[self.qindexes[0]:]
        
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

    #--------------------------------------------

    def get_num_refs(self, num):
        """
        Returns a list of the bases where all counts match 'num'
        """
        
        refs = []
        
        for indx, count in enumerate(self.ref_counts):

            if count == num:
                refs.append( self.qbases[indx] )

        return sorted(refs)

    #--------------------------------------------

    def get_count_cols(self, num):
        """
        Returns a list with all bases that match up for
        a column count.

        """

        count_columns = []
        
        for i in range(1,num+1):

            nrefs = self.get_num_refs(i)

            if len(nrefs) > 0:
                
                count_columns.append( ''.join(nrefs) )

            else:

                count_columns.append('--')

        return count_columns
    
#*** end SNP class *********************************************




#***** AlphaIncrementer class object ***********************************
class AlphaIncrementer:

    #--------------------------------------------
    def __init__(self, start=chr(ord('A') - 1)):

        self.current_char = 0


    #--------------------------------------------

    def get_next(self):
        """
        Increments the current character to the next one
        and returns the new current.
        """
        
        #nchar = chr(ord(self.current_char)+1)

        #self.current_char = nchar

        self.current_char += 1

        return str(self.current_char)
    
#****************************************************************



#****************************************************************
#
# A class that contains the dictionary for hashing against the molecule
# name.
class MoleculeDict:

    def __init__(self):
        self.molecule_dict = {}

    #--------------------------------------------
    def add(self, snp):

        if not self.molecule_dict.has_key(snp.count_id):

            self.molecule_dict[snp.count_id] = []
            

        self.molecule_dict[snp.count_id].append(snp.molecule)


    #--------------------------------------------
    def get(self, snp):

        if self.molecule_dict.has_key(snp.count_id):

            return self.molecule_dict[snp.count_id]

        else:

            return []
    #--------------------------------------------
    def keys(self):
        return self.molecule_dict.keys()
    #-----------------------------------------------
    def get_string(self, snp, unique=False):

        if len(self.get(snp)) > 0:
            
            if unique:
                return ",".join(list(set(self.get(snp))))
            else:
                return ",".join(self.get(snp))
            
        else:
            return "--"
            
#****************************************************************


#****************************************************************
# 
#
#****************************************************************
class RefposDict:
    def __init__(self):
        self.refpos_dict = {}

    #--------------------------------------------
    def add(self, snp):

        if not self.refpos_dict.has_key(snp.count_id):

            self.refpos_dict[snp.count_id] = []
            

        self.refpos_dict[snp.count_id].append(snp.ref_base)


    #--------------------------------------------
    def get(self, snp):

        if self.refpos_dict.has_key(snp.count_id):

            return self.refpos_dict[snp.count_id]

        else:

            return []

    #-----------------------------------------------
    def get_string(self, snp, unique=False):

        if len(self.get(snp)) > 0:
            
            if unique:
                return ",".join(list(set(self.get(snp))))
            else:
                return ",".join(self.get(snp))
            
        else:
            return "--"
        
#****************************************************************



        
#-------------------------------------------------------------------------------

def load_table(table_file):
    """
    Loads the table file into a tuple that returns the header, qindexes
    and a list of all snp objects.
    """
    
    infile = open(table_file, 'r')

    table_data = infile.readlines()

    infile.close()

    # The header should always be the first line
    header = table_data[0].split('\t')

    # collect the qindexes from the header with this loop within a list
    qindexes = [indx for indx,colname in enumerate(header) if "qbase:" in colname]


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
    parser.add_argument("-s", "--sort", action="store_true",
                        help="Sort snps by code")
    parser.add_argument("-c", "--counts", type=int,
                        help="Number of count columns to include in the table", default=3)
    parser.add_argument("-m", "--min_table", action="store_true",
                        help="Output a minimum table with only the counts and codes.")

    args = parser.parse_args()

    output_file = args.outfile

    # open the snp table and load it into some data types.
    
    header, qindexes, snp_objects = load_table(args.snp_table)

    # Here I'm going to encode the snp objects in a hash
    # using the alpha incrementer.
    alphainc = AlphaIncrementer()
    encoding_dict = {}
    molecule_dict = MoleculeDict()
    refpos_dict = RefposDict()

    # loop through snp objects and set each one
    for snp in snp_objects:

        # if this is not given a key, then we add it and set it.
        if not encoding_dict.has_key(snp.count_id):

            encoding_dict[snp.count_id] = alphainc.get_next()

        snp.code = encoding_dict[snp.count_id]


        # add this snp info to the molecule and refpos dicts
        molecule_dict.add(snp)
        refpos_dict.add(snp)
        
    #------- Finished loop for encoding ---------------------------------------



    #
    # From here we just output with the new columns, same loop order as before.
    #
    with open(output_file, 'w') as of:

        newcols_start = qindexes[0] 
        newcols_end = qindexes[0] + 1

        # First write the header.
        of.write( "\t".join(header[0:newcols_start] + ["Group", "Count", "Reference Positions", "Molecule Number"] + header[newcols_end:]) )

        
        for snp in snp_objects:

            line = "\t".join( snp.first_half() + [ snp.code, "".join([ str(r) for r in snp.get_ref_counts()]), refpos_dict.get_string(snp), molecule_dict.get_string(snp)] + snp.second_half() ) 

            of.write(line)
            
    print "Output %d SNPs to file %s" % (len(snp_objects), output_file)

#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
