#!/usr/bin/env python
"""

Program loads two hiseq tables from the machine and transfers all of the
sample well data from the first into the second. 

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

class NewSampleWellTable:

    def __init__(self, csv_file, sw):

        self.data = read_file(csv_file)

        self.sample_wells = {}

        self.header_indexes = {}

        sample_well_index = None

        found_data = False

        for line in self.data:

            if "Sample_ID" in line:
                
                for i,item in enumerate(line):
                        
                        self.header_indexes[item] = i
 
                sample_well_index = self.header_indexes["Sample_Well"]
                
                found_data = True


            if found_data:

                sample_well_id = line[sample_well_index]

                if not sw.sample_wells.has_key(sample_well_id):

                    print "Sample Well %s for '%s' not found" % (sample_well_id,",".join(line))
                    continue
                
                else:

                    
                    for i,j in enumerate(self.header_indexes):


                        if "ndex" in j:
                            this_index = self.header_indexes[j]

                            value = sw.get_data(sample_well_id,j)

                            line[this_index] = value


                    
#-------------------------------------------------------------------------------

class SampleWells:

    def __init__(self, csv_file):

        self.data = read_file(csv_file)

        self.sample_wells = {}

        self.header_indexes = {}

        self.sample_wells = {}

        self.sample_well_index = None

        found_data = False

        # now go through the data to get the sample wells
        for line in self.data:

            if "Sample_ID" in line:
                
                for i,item in enumerate(line):
                    
                    self.header_indexes[item] = i

                self.sample_well_index = self.header_indexes["Sample_Well"]
                
                found_data = True


            if found_data and not "Sample_ID" in line:

                sample_well_id = line[self.sample_well_index]

                if self.sample_wells.has_key(sample_well_id):

                    print "Sample well %s is not unique" % sample_well_id

                else:

                    self.sample_wells[sample_well_id] = line

        print "Read in %d Sample Wells" % len(self.sample_wells.keys())

    #---------------------------------------------------------------------------

    def get_data(self, sw_id, key):

        line = self.sample_wells[sw_id]

        this_index = self.header_indexes[key]

        return line[this_index]


#-------------------------------------------------------------------------------

def __main__():

    parser = argparse.ArgumentParser(description="""
Program loads two hiseq tables from the machine and transfers all of the
sample well data from the first into the second and outputs a new csv file into
a file called new_sheet.csv
    """)
    parser.add_argument("sample_well_file", type=str,
                        help="The csv file with the sample well data")
    parser.add_argument("empty_well_file", type=str,
                        help="The csv file with parts to fill in from the sample well file ")
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file, by defualt new_sheet.csv", default="new_sheet.csv")

    args = parser.parse_args()

    output_file = args.outfile

    sw = SampleWells(args.sample_well_file)
    nsw = NewSampleWellTable(args.empty_well_file,sw)


    with open(output_file,'w') as outf:

        for line in nsw.data:

            outf.write(",".join(line)+"\n")
            

    print "CSV file has been written to %s" % output_file
    
#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()



