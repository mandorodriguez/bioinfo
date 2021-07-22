#!/usr/bin/env python3

import os, sys
import pdb
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('file', type=str, help="A file with the numbers.")

args = parser.parse_args()

num_file = args.file
tmp_nums = []

with open(num_file,'r') as nf:

    tmp_nums = nf.readlines()

nums = list(set([int(n) for n in tmp_nums]))

nums.sort()

groups = []
this_group = []
num_len = len(nums)

for i,n in enumerate(nums):

    try:
        this_group.append(n)

        if nums[i]+1 != nums[i+1]:
            if len(this_group) > 1:
                groups.append( this_group )
            this_group = []

    except IndexError:
        if len(this_group) > 1:
            groups.append(this_group)
        break


    
with open("grouped-numbers.txt","w") as of:
    
    for g in groups:
        of.write("{}-{}\n".format( g[0], g[-1] ))


