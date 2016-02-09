#!/bin/env python

import os
import sys
import argparse
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
parser.add_argument("-o", "--outfile", type=str,
                    help="The output file", default="snp_genotype_out.txt")
parser.add_argument("-t", "--translation_table", type=int, default=1,
                    help="The translation table to use for amino translations")
parser.add_argument("-g", "--groupfile", type=str,
                    help="Output file indexed by pattern groups")






