#!/usr/bin/python3

################################################################################
# summarize_repeats.py
# Author: Linnéa Smeds
# Date: 27 September 2023

# Summarise repeat information from output from instersectBed -wo (regions of
# interests overlapped with repeats). Makes sure each base is only counted once,
# even if there are overlaps in the repeat annotation. 

################################################################################
##### Import libraries
import argparse
import re


################################################################################
##### Parse commandline
parser = argparse.ArgumentParser(description = "This script takes tab-separated \
output from intersectBed -wo with genomic regions of interest first, and the \
repeats they overlap with second. The regions are summarized per base, so that \
any overlap is calculated only once")
parser.add_argument('-f', '--file', help = 'Overlap from intersectBed', type = str, required = True)
parser.add_argument('-o', '--output', help = 'Output file name', type = str, required = True)
args = parser.parse_args()


################################################################################
##### Main code

# Go through all bases of interest and assign them as "NA" (non repetetive),
# then go through the corresponding repeat and reset the values that overlap
bp_dict={}
with open(args.file, 'r') as infile:
    for line in infile:
        line = line.rstrip()
        tabs = line.split("\t")
        #Add chromosome if it doesn't exists
        if tabs[0] not in bp_dict:
            bp_dict[tabs[0]]={}
        # Go through all positions in the variant
        start=int(tabs[1])+1
        end=int(tabs[2])+1
        for i in range(start, end, 1):
            if i not in bp_dict[tabs[0]]:
                bp_dict[tabs[0]][i]="NA"
        # Go through overlap, only add bps that are in the dict
        rstart=int(tabs[5])+1
        rend=int(tabs[6])+1
        for j in range(rstart, rend, 1):
            if j in bp_dict[tabs[0]]:
                bp_dict[tabs[0]][j]=tabs[7]
                if bp_dict[tabs[0]][j] != tabs[7]:
                    print("WARNING: Overlapping repeats! "+bp_dict[tabs[0]][j]+" is replaced with "+tabs[7])


# Print outfile (one line per base pair)
with open(args.output, 'w') as outfile:
    for chr in sorted(bp_dict.keys()):
        for pos in sorted(bp_dict[chr].keys()):
            start=pos-1
            print(chr+"\t"+str(start)+"\t"+str(pos)+"\t"+bp_dict[chr][pos], file=outfile)
