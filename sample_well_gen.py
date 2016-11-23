#!/usr/bin/env python
"""


"""


import argparse
import pdb
import csv

#-------------------------------------------------------------------------------

def read_file(csv_file):
    with open(csv_file, 'r') as f:
        data = [row for row in csv.reader(f.read().splitlines())]

    return data

    
#-------------------------------------------------------------------------------

class Indexes:

    def __init__(self, index_file):

        data = read_file(index_file)

        self.headers = data[0]

        self.alpha = [(l[0],l[1]) for l in data[1:] if l[0] != '' and l[1] != '']
        self.num = [(l[2],l[3]) for l in data[1:] if l[2] != '' and l[3] != '']

        self.indexes = []
        self.lines = []

        curr_char = 'A'
        curr_num = 1
        
        for a in self.alpha:
            for n in self.num:

                sample_well = "%s%s" % (curr_char,str(curr_num))
                
                self.indexes.append( ( sample_well, a, n) )
                self.lines.append( [sample_well, a[0], a[1], n[0], n[1]] )

                curr_num += 1

            curr_char = chr(ord(curr_char) + 1)
            curr_num = 1
            
                    

#-------------------------------------------------------------------------------

def __main__():

    parser = argparse.ArgumentParser(description="""
Program loads two hiseq tables from the machine and transfers all of the
sample well data from the first into the second and outputs a new csv file into
a file called new_sheet.csv
    """)
    parser.add_argument("index_file", type=str,
                        help="The csv file with the two sets of indexes to use to fill in the sample wells ")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file, by defualt new_sheet.csv", default="new_sheet.csv")

    args = parser.parse_args()

    output_file = args.outfile

    
    indx = Indexes(args.index_file)


    first_line = ["Sample_Well"]+indx.headers
    
    with open(output_file,'w') as outf:

        outf.write(",".join(first_line)+"\n")

        for l in indx.lines:


            outf.write(",".join(l)+"\n")
            

    print "CSV file has been written to %s" % output_file
    
#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()



