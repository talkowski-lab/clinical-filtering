#!/usr/bin/env python
# -*- coding: utf-8 -*-


#
# This is a simple script to create an upset plot from a list of gene lists in the same directory.
# Inputs: Path to file of geneLists.txt, this should be a file with <name of gene list> and <path to genelist>
# Output: Upset plot for gene lists
#

#For input, need gene list name as column names and gene names as row names. 

# Imports
import sys
import re
import UpSetPlot
import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--list", help="List of Gene Names")
parser.add_argument("-o", "--output", help="Upset Plot")

args = parser.parse_args()

geneLists = open(args.list, mode = 'r')
geneLines = geneLists.readlines()
geneLists.close()

for line in GeneLines:
