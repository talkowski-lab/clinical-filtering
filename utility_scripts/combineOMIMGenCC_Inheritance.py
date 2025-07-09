#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# This takes a genCC and OMIM Inheritance File
# and creates a file with unique gene names and inheritance codes
# 12 March 2025
#
# Updates: 8 July 2025
# Updated flag text

# Imports
import sys
import re
import argparse

def separate_digits(number):
    """Separates a number into its individual digits.

    Args:
        number: The number to separate.

    Returns:
        A list of digits.
    """

    digits = []
    while number > 0:
        digit = number % 10
        digits.insert(0, digit)  # Insert at the beginning to maintain order
        number //= 10
    return digits

parser = argparse.ArgumentParser()

#-i INTERVAL_LIST -g GENEFILE -o OUTPUT
parser.add_argument("-i", "--genCC_curated", help="genCC Curated Input (GenCC Processed)")
parser.add_argument("-g", "--OMIM_curated", help="OMIM Curated Input (OMIM_Step2)")
parser.add_argument("--a", dest="gen_CC_complete", help="genCC complete gene list as .txt")
parser.add_argument("-o", "--output", help="Output File")
parser.add_argument("--genCC", help="Flag: Add OMIM inheritance ONLY IF genCC gene is 1) non-autosomal OR 2) gene not in genCC complete gene list. --a required if this flag is ticked.", action='store_true')
parser.add_argument("--limit5", help="Add flag if only want to report 5 if only output", action='store_true')

parser.set_defaults(genCC=False)
parser.set_defaults(limit5=False)

args = parser.parse_args()

genCC = args.genCC
limit5 = args.limit5

#Make list of input gene names
OMIM_File = open(args.OMIM_curated, mode = 'r')
OMIMLines = OMIM_File.readlines()
OMIM_File.close()

genCC_File = open(args.genCC_curated, mode = 'r')
genCCLines = genCC_File.readlines()
genCC_File.close()

genCCDict={}
OMIMDict={}

genCC_List = []

#Make list of genCC genes
if genCC==True: 
    genCCCompFile = open(args.gen_CC_complete, mode = 'r')
    genCC_Comp_Lines = genCCCompFile.readlines()
    counter = 0
    for line in genCCLines:
        if counter == 0:
            counter = counter+1
            continue
        genCC_List.append(line.strip('\n'))
    genCCCompFile.close()

#Output file
output_File = open(args.output, mode = 'w')

output_File.write("approvedGeneSymbol"+"\t"+"inheritance_code" + "\n")

#Read in OMIM and genCC curated gene lists

counter = 0
for line in genCCLines:
    if counter == 0:
        counter = counter+1
        continue
    line = line.strip('\n')
    valueList = line.split('\t')

    geneSymbol = valueList[0]
    inheritanceValue = valueList[1]

    genCCDict[geneSymbol]=inheritanceValue

counter = 0
for line in OMIMLines:
    if counter == 0:
        counter = counter+1
        continue
    line = line.strip('\n')
    valueList = line.split('\t')

    geneSymbol = valueList[0]
    inheritanceValue = valueList[1]

#only include OMIM genes if (not in genCC AND in OMIM strong) OR (if in genCC AND in OMIM strong AND have inheritance 3 or 4)
    if genCC == True:
        if (geneSymbol not in genCC_List) or (geneSymbol in genCC_List and (3 in separate_digits(inheritanceValue) or 4 in separate_digits(inheritanceValue))):
            OMIMDict[geneSymbol]=inheritanceValue
    if genCC == False: 
        OMIMDict[geneSymbol]=inheritanceValue

#Sorted master gene list
master_keys = sorted(list(set(genCCDict.keys()).union(set(OMIMDict.keys()))))


for key in master_keys:
    if key in OMIMDict.keys():
        OMIM_inh = set(separate_digits(int(OMIMDict[key])))
    else:
        OMIM_inh = set()

    if key in genCCDict.keys():
        genCC_inh = set(separate_digits(int(genCCDict[key])))
    else:
        genCC_inh = set()

    joint_inh_list = sorted(list(OMIM_inh.union(genCC_inh)))
  
    string_list = [str(item) for item in joint_inh_list]

    if limit5==True and len(string_list)>1 and ("5" in string_list):
        string_list.remove("5")

    string_inh = "".join(string_list)
    output_File.write(key + "\t" + string_inh + "\n")

output_File.close()
exit()