#!/usr/bin/env python

import argparse
import os
import subprocess
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

def write_job_script(outdir, seqid, r1, r2):


    job_filename = "/".join([outdir, "%s.job" % seqid])

    with open(job_filename, 'w') as js:

        js.write("""
#$ -q bigmem.q
#$ -l bigmem
#$ -N spades
#$ -o %s_$JOB_ID.log
#$ -j y
#$ -cwd


spades.py --pe1-1 %s --pe1-2 %s --careful -k 21,33,55,77,99,127 -t 32 -o spades_results_%s 
exit 0
""" % (seqid, r1, r2, seqid))

    return job_filename

#-------------------------------------------------------------------------------
        
def fasta_files(directory, extension=".fasta.gz"):
    """
    returns a dict of all the files with the R1 and R2 file in a tuple
    """
    fasta_files = []
    
    for dirpath,_,filenames in os.walk(directory):
        
        for f in filenames:

            fullname = os.path.abspath(os.path.join(dirpath, f))

            if True or fullname.endswith(extension):

                fasta_files.append(fullname)

    
    return fasta_files
           
#-------------------------------------------------------------------------------
# Main function call
def __main__():

    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--fasta-dir", type=str,
                        help="Directory with all the fastq sequences")
    parser.add_argument("-o", "--out", type=str,
                        help="Result output directory")
    parser.add_argument("-s", "--submit",  action='store_true',default=False,
                        help="Submit the job scripts after creating them")


    args = parser.parse_args()

    fasta_dir = args.fasta_dir
    result_dir = args.out
    submit = args.submit

    all_files = fasta_files(fasta_dir)

    sequences_dict = {}

    for f in all_files:

        just_filename = f.split('/')[-1]

        seqid = ""
        
        if "_R1_" in just_filename:
            seqid = just_filename.split("_R1_")[0]

            if not sequences_dict.has_key(seqid):
                sequences_dict[seqid] = {}
                
            sequences_dict[seqid]["R1"] = f
            
        elif "_R2_" in just_filename:
            seqid = just_filename.split("_R2_")[0]
            
            if not sequences_dict.has_key(seqid):
                sequences_dict[seqid] = {}

            sequences_dict[seqid]["R2"] = f

                    
    #----- now we see about writing the job files ----------------

    if not os.path.exists(result_dir):
        os.mkdir(result_dir)
        print "Created result directory %s" % result_dir
    else:
        print "Result directory %s is ready" % result_dir

    job_files = []

    for k in sequences_dict.keys():

        job_file = write_job_script(result_dir, k, sequences_dict[k]["R1"], sequences_dict[k]["R2"])

        job_files.append(job_file)


    #----- Submit jobs? -------------------------------------------

    if submit:

        print "Submitting %d jobs" % len(job_files)        

        for j in job_files:
            subprocess.call(["qsub", j])


#-------------------------------------------------------------------------------
if __name__=="__main__": __main__()
