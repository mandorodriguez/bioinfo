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
    def __init__(self, qindexes, line, query_genes, ref_pos_index=1, ref_base_index=3):

        self.snp_data = line.split('\t')

        # molecule is at index 0
        self.molecule = self.snp_data[0]
        
        self.qindexes = qindexes

        # the ref pos
        self.ref_pos = self.snp_data[ref_pos_index]

        # just the ref base
        self.ref_base = self.snp_data[ref_base_index]

        # collect the qbases from the SNP
        self.qbases = [self.snp_data[qi] for qi in self.qindexes]

        self.query_genes = query_genes

        # get the ref counts from the SNP
        self.pattern_list = self.get_pattern()

        # Convert our ref counts to a string as an ID for hashing.
        self.pattern = ''.join([str(b) for b in self.pattern_list])


        # the occurance of different numbers in the pattern
        self.zeroes = self.pattern_list.count(0)
        self.ones = self.pattern_list.count(1)
        self.twos = self.pattern_list.count(2)
        self.threes = self.pattern_list.count(3)

        self.info = self.set_informative()


        # lists of the genes at the positions with the number in the pattern
        self.genes_w_one = self.list_to_string(self.get_genes(1))
        self.genes_w_two = self.list_to_string(self.get_genes(2))
        self.genes_w_three = self.list_to_string(self.get_genes(3))

        # A class variable for the assigned encoding to set later
        self.group = None

    #--------------------------------------------

    def list_to_string(self, l):

        if len(l) == 0:
            return "--"
        else:
            return ",".join(l)
    
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

    def set_informative(self):
        """
        Just going to use the global ones, twos and threes variables for the check.

        Simple check to find out if it's informative or not.
        """

        if self.ones == 1 or self.zeroes == 1:

            if self.twos == 0 or (self.twos == 1 and self.threes <= 1):

                return "NI"


        return "PI"

    #--------------------------------------------

    def get_genes(self, num):
        """
        retreives all genes in the self.query_genes list that match up
        with the postion of the number in the pattern_list.
        """
        indexes = [i for i,val in enumerate(self.pattern_list) if val == num]

        return [self.query_genes[i] for i in indexes]
        
    #--------------------------------------------

    def get_pattern(self):
        """
        performs the walk along the qbases and does the count if
        refbase bases match the qbase and resets if it doesn't.
        Only need to call this once on init since the SNP is immutable.
        """
        encoding = []
        seen = {}
        count = 0

        # set the first for refbase
        seen[self.ref_base] = 0
                
        for qb in self.qbases:

            if not seen.has_key(qb):
                count += 1
                seen[qb] = count


            encoding.append(seen[qb])


        return encoding


    #--------------------------------------------

    def get_num_refs(self, num):
        """
        Returns a list of the bases where all counts match 'num'
        """
        
        refs = []
        
        for indx, count in enumerate(self.pattern_list):

            if count == num:
                refs.append( self.qbases[indx] )

        return sorted(refs)

    
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
# name and saves the reference position.
class MoleculeDict:

    def __init__(self):
        self.molecule_dict = {}

    #--------------------------------------------
    def add(self, snp):

        if not self.molecule_dict.has_key(snp.molecule):

            self.molecule_dict[snp.pattern] = []
            

        self.molecule_dict[snp.pattern].append(str(snp.ref_pos))


    #--------------------------------------------
    def get(self, snp):

        if self.molecule_dict.has_key(snp.molecule):

            return self.molecule_dict[snp.molecule]

        else:

            return []
    #--------------------------------------------
    def keys(self):
        return self.molecule_dict.keys()
            
#****************************************************************


#****************************************************************
# 
#
#****************************************************************
class GroupDict:
    def __init__(self):
        self._dict = {}

    #--------------------------------------------
    def add(self, snp):

        if not self._dict.has_key(snp.pattern):

            self._dict[snp.pattern] = {}

        if not self._dict[snp.pattern].has_key(snp.molecule):

            self._dict[snp.pattern][snp.molecule] = []

            
        #self._dict[snp.pattern].append(snp)
        self._dict[snp.pattern][snp.molecule].append(snp)


    #--------------------------------------------
    def get(self, snp):

        reflist = []
        
        if self._dict.has_key(snp.pattern):

            pattern = self._dict[snp.pattern]

            
            
            for k in pattern.keys():

                refs = ",".join([str(s.ref_pos) for s in pattern[k]])

                molpos = "%s: %s" % (k,refs)

                reflist.append(molpos)

        

        return reflist


    #-----------------------------------------------
    def get_string(self, snp, unique=False):
        
        if len(self.get(snp)) > 0:

            items = []
            if unique:
                
                items = list(set(self.get(snp)))

            else:

                items = self.get(snp)

            return " | ".join(items)
            
        else:
            return "--"
        
#****************************************************************

# class GroupDict:
#     def __init__(self):
#         self._dict = {}

#     #--------------------------------------------
#     def add(self, snp):

#         if not self._dict.has_key(snp.pattern):

#             self._dict[snp.pattern] = []

            
#         self._dict[snp.pattern].append(snp)


#     #--------------------------------------------
#     def get(self, snp):

#         if self._dict.has_key(snp.pattern):

#             return self._dict[snp.pattern]

#         else:

#             return []

#     #-----------------------------------------------
#     def get_string(self, snp, unique=False):
        
#         if len(self.get(snp)) > 0:

#             items = []
#             if unique:
                
#                 items = list(set(self.get(snp)))

#             else:

#                 items = self.get(snp)

#             return ",".join([s.ref_pos for s in items if snp.molecule == s.molecule])
            
#         else:
#             return "--"
        



#****************************************************************
# 
#
#****************************************************************
class TypeDict:
    def __init__(self):
        self._dict = {}



    #--------------------------------------------
    def add(self, snp):

        keys = list(set(snp.pattern_list))

        for k in keys:

            self.add_to_dict(snp,k)
            
    #--------------------------------------------
    def add_to_dict(self, snp, key):

        if not self._dict.has_key(key):

            self._dict[key] = []
            

        self._dict[key].append(snp)


    #--------------------------------------------
    def get(self, type):
        """
        Returns all snps collected with the number matching 'type'
        in its code.
        """

        if self._dict.has_key(type):

            return self._dict[type]

        else:

            return []


        
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

    query_genes = [header[qi].replace("qbase:","").rstrip() for qi in qindexes]


    # put each line representing a SNP into a SNP object which was declared earlier
    snp_objects = []
    
    for snp_line in table_data[1:]:

        snp_objects.append( SNP(qindexes, snp_line, query_genes) )


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
#    parser.add_argument("-m", "--min_table", action="store_true",
#                        help="Output a minimum table with only the patterns and codes.")

    args = parser.parse_args()

    output_file = args.outfile

    # open the snp table and load it into some data types.
    
    header, qindexes, snp_objects = load_table(args.snp_table)

    # Here I'm going to encode the snp objects in a hash
    # using the alpha incrementer.
    alphainc = AlphaIncrementer()
    encoding_dict = {}
    group_dict = GroupDict()
    type_dict = TypeDict()

    # loop through snp objects and set each one
    for snp in snp_objects:

        # if this is not given a key, then we add it and set it.
        if not encoding_dict.has_key(snp.pattern):

            encoding_dict[snp.pattern] = alphainc.get_next()

        snp.group = encoding_dict[snp.pattern]


        # add this snp info to the molecule and refpos dicts
        group_dict.add(snp)
        type_dict.add(snp)
        
    #------- Finished loop for encoding ---------------------------------------


    #
    # From here we just output with the new columns, same loop order as before.
    #
    with open(output_file, 'w') as of:

        newcols_start = qindexes[0] 
        newcols_end = qindexes[0] 

        # First write the header.
        of.write( "\t".join(header[0:newcols_start] + ["Pattern", "Group", "Informative", "Genomes_w_1", "Genomes_w_2", "Genomes_w_3", "Reference Positions"] + header[newcols_end:]) )

        for snp in snp_objects:

            line = "\t".join( snp.first_half() + [ snp.pattern, snp.group, snp.info, snp.genes_w_one, snp.genes_w_two, snp.genes_w_three, group_dict.get_string(snp)] + snp.second_half() ) 

            of.write(line)
            
    print "Output %d SNPs to file %s" % (len(snp_objects), output_file)

#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
