#!/bin/env python

import os
import sys
import argparse
import pdb
from subprocess import call


"""
Program is a wrapper for prokka to allow the passing of a
file with arguments to send to prokka in a system call
so that it can be used in a batch mode with lists in galaxy.

"""

#-------------------------------------------------------------------------------
def run_prokka(arg_file, cpus, outdir, prefix, locustag, increment, gffver, addgenes,
               mincontig, compliant, centre, genus, species,
               strain, plasmid, kingdom, gcode, usegenus, proteins, metagenome,
               fast, evalue, rfam, norrna, notrna):
    out = call(["ls","-l"])

#-------------------------------------------------------------------------------


parser = argparse.ArgumentParser()

parser.add_argument("arg_table", type=str,
                    help="The argument table to pass to prokka")
parser.add_argument("--outdir", type=str)
parser.add_argument("--prefix", type=str)
parser.add_argument("--locustag", type=str)
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

args = parser.parse_args()

pdb.set_trace()





