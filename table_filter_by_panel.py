#!/usr/bin/env python

import argparse
import pdb


#-------------------------------------------------------------------------------
# Main function call
def __main__():

    parser = argparse.ArgumentParser()
    parser.add_argument("panel", type=str,
                        help="Panel file with molecue and positions",metavar="panel_file")
    parser.add_argument("snp_table", type=str,
                        help="The snp table to input",metavar="snp_table")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="snp_table_out.txt")

    args = parser.parse_args()


    snp_table_file = args.snp_table
    panel_file = args.panel
    out_file = args.outfile

    panel_dict = {}

    # read in all of the panel file
    with open(panel_file,'r') as pf:

        for line in pf:

            (molecule, pos) = line.split()

            if panel_dict.has_key(molecule):

                panel_dict[molecule].append(pos)

            else:

                panel_dict[molecule] = [pos]

    #----- end loading panel ----------------------

    out_file_descriptor = open(out_file,'w')

    # loop through each line in the table and check if the mol and pos is present
    # if it is then we write to the output file, otherwise we ignore it.
    with open(snp_table_file,'r') as stf:

        for line in stf:

            parts = line.split()

            this_mol = parts[0]
            this_pos = parts[1]

            if panel_dict.has_key(this_mol):
                if this_pos in panel_dict[this_mol]:
                    out_file_descriptor.write(line)


    # Done
    out_file_descriptor.close()
    print "Finished writing table '%s'" % out_file
    
#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
