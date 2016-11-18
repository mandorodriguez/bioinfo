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

    if len(csv_data) == 1 and not '\n' in csv_data and not '\r' in csv_data:

        #print "Condensing list"
        
        csv_data = csv_data[0]

    if '\r' in csv_data:
        
        #print "Splitting data by return carriages"

        csv_data = csv_data.split('\r')

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

    found_data = False

    header_indexes = []
    unique_indexes = {}
    indexes = []

    for line in csv_data:

        if "Sample_ID" in line:
            for indx,tag in enumerate(line.split(',')):
                if "ndex" in tag:
                    indexes.append(indx)
                    header_indexes.append(tag)

            found_data = True

        if found_data:

            line_parts = line.split(',')

            index_id = ",".join([line_parts[i] for i in indexes])

            if unique_indexes.has_key(index_id):

                unique_indexes[index_id].append(line_parts[0])
                
            else:

                unique_indexes[index_id] = [line_parts[0]]
        
    outfile = open(output_file, 'w')

    outfile.write("Sample_ID\tUnique\t"+'\t'.join( header_indexes )+'\n')
       
    for index_id,sample_names in unique_indexes.iteritems():

        is_or_not = "unique" if len(sample_names) == 1 else "not unique"

        index_tags = index_id.split(',')
        
        for s in sample_names:

            print "%s\t%s" % (s,is_or_not)

            outfile.write("%s\t%s\t%s\n" % (s, is_or_not,'\t'.join(index_tags)) )

            

    outfile.close()

    print "CSV file has been written to %s" % output_file
    
#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()



