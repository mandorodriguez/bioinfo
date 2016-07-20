#!/usr/bin/env python

import argparse
import os
import pdb

#*******************************************************************************
# Program: Spades batch
#
# Author: Mando Rodriguez
#
#
# Program that creates a batch of scripts to submit to the sun grid engine for
# running spades assemblies.
#
#*******************************************************************************


#-------------------------------------------------------------------------------

def write_job_script(filename):

    with open(filename, 'w') as js:

        js.write("""
#$ -q bigmem.q
#$ -l bigmem
#$ -N spades
#$ -o spades_ISU-1002_S22_$JOB_ID.log
#$ -j y
#$ -cwd


spades.py --pe1-1 ~/MRSA_assemblies/sequences/ISU-1002_S22_L001_R1_001.fastq --pe1-2 ~/MRSA_assemblies/sequences/ISU-1002_S22_L001_R2_001.fastq --careful -k 21,33,55,77,99,127 -t 32 -o spades_results_ISU-1002_S22 
exit 0
""")
        
#-------------------------------------------------------------------------------
        
def fasta_files(directory):

    fasta_files = []
    
    for dirpath,_,filenames in os.walk(directory):
        
        for f in filenames:

            fullname = os.path.abspath(os.path.join(dirpath, f))

            if fullname.endswith(".fasta"):

                fasta_files.append(fullname)
           
#-------------------------------------------------------------------------------
# Main function call
def __main__():

    parser = argparse.ArgumentParser()
    parser.add_argument("-o", "--outfile", type=str,
                        help="The output file", default="snp_genotype_out.txt")
    parser.add_argument("-t", "--translation_table", type=int, default=1,
                        help="The translation table to use for amino translations")
    parser.add_argument("-g", "--groupfile", type=str,
                        help="Output file indexed by pattern groups")


    args = parser.parse_args()

    output_file = args.outfile

    pdb.set_trace()


#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
