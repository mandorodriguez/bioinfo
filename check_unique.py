#!/usr/bin/env python

import argparse
import pdb

#-------------------------------------------------------------------------------

def load_table(csv_file):
    """
    Loads the table file into a tuple that returns the header, qindexes
    and a list of all snp objects.
    """
    
    infile = open(csv_file, 'r')

    csv_data = infile.readlines()

    infile.close()

    return csv_data

#-------------------------------------------------------------------------------

def __main__():

    parser = argparse.ArgumentParser()
    parser.add_argument("csv_file", type=str,
                        help="The csv to input")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="check_unique_out.txt")

    args = parser.parse_args()

    output_file = args.outfile

    csv_data = load_table(args.csv_file)
    
    pdb.set_trace()
#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()



