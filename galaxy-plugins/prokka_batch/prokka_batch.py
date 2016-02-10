#!/bin/env python

import os
import sys
import argparse
import pdb
import subprocess 


"""
Program is a wrapper for prokka to allow the passing of a
file with arguments to send to prokka in a system call
so that it can be used in a batch mode with lists in galaxy.

"""



#-------------------------------------------------------------------------------


parser = argparse.ArgumentParser()

parser.add_argument("arg_table", type=str,
                    help="The argument table to pass to prokka")

parser.add_argument("--contigs", type=str)
"""
parser.add_argument("--outdir", type=str)
parser.add_argument("--prefix", type=str)
parser.add_argument("--cpus", type=str)
parser.add_argument("--increment", type=str)
parser.add_argument("--gffver", type=str)
parser.add_argument("--addgenes", type=str)
parser.add_argument("--mincontig", type=str)
parser.add_argument("--compliant", type=str)
parser.add_argument("--centre", type=str)
parser.add_argument("--genus", type=str)
parser.add_argument("--species", type=str)
parser.add_argument("--strain", type=str)
parser.add_argument("--plasmid", type=str)
parser.add_argument("--kingdom", type=str)
parser.add_argument("--gcode", type=str)
parser.add_argument("--usegenus", type=str)
parser.add_argument("--proteins", type=str)
parser.add_argument("--metagenome", type=str)
parser.add_argument("--fast", type=str)
parser.add_argument("--evalue", type=str)
parser.add_argument("--rfam", type=str)
parser.add_argument("--norrna", type=str)
parser.add_argument("--notra", type=str)
"""
args = parser.parse_args()

metaname = args.contigs.rstrip(".fasta")

#-------------------------------------------------------------------------------
# Here open up the tabular file and put it into a dict
#-------------------------------------------------------------------------------
prokka_args = {}

with open(args.arg_table, 'r') as at:

    table_data = at.readlines()

    for line in table_data:
        
        line_parts = [i.rstrip() for i in line.split('\t')]

        if not prokka_args.has_key(line_parts[0]):

            prokka_args[line_parts[0]] = line_parts

        else:

            print "Data for '%s' already loaded in table" % line_parts[0]


#------------------- End of load table ------------------------------------------

pdb.set_trace()
# collect our arguments and run prokka

cmd = ["prokka"]

cmd += ["--cpus", args.cpus]
cmd += ["--quiet"]
cmd += ["--outdir","outdir"]
cmd += ["--prefix", args.prefix]
cmd += ["--increment", args.increment]
cmd += ["--gffver", args.gffver]
cmd += ["--addgenes", args.addgenes]
cmd += ["--mincontig", args.mincontig]
cmd += ["--complaint",args.compliant]

#cmd += ["--gcode", args.gcode]
cmd += ["--usegenus",args.usegenus]
cmd += ["--proteins",args.proteins]
cmd += ["--metagenome",args.metagenome]
cmd += ["--fast", args.fast]
cmd += ["--evalue", args.evalue]
cmd += ["--rfam", args.rfam]
cmd += ["--norrna", args.norrna]
cmd += ["--notrna", args.notrna]

# these we fill in from the argument table. If there's no match the key error
# should let the script quit with an error
run_args = prokka_args[metaname]

cmd += ["--strain", run_args[0]]
cmd += ["--locustag", run_args[1]]
cmd += ["--centre", run_args[2]]
cmd += ["--genus", run_args[3]]
cmd += ["--species", run_args[4]]

if not run_args[5]=="":
    cmd += ["--plasmid", run_args[5]]

if not run_args[6]=="":
    cmd += ["--gcode", run_args[6]]


print "running prokka:\n  %s\n" % " ".join(cmd)

out = subprocess.Popen(cmd,stdout=subprocess.PIPE,
                     stderr=subprocess.STDOUT, shell=True)

if not out.stdout is None:
    for line in iter(out.stdout.readline,b''):
        print line.rstrip()
        
#--------------------------------------------------------------------------------

    





