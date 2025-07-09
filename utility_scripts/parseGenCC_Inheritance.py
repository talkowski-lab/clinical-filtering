#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# This takes a genCC TSV File
# and creates a file with unique gene names and inheritance codes
# 12 March 2025
#
# Updates: 8 July 2025
# Changed text surround flags. 
# Added a high-priority flag to include only strong/definitive entries, unless no strong/definitive entries exist. 
# Added a flag to output list of genCC genes, and strong/definitive Gen CC genes

# Imports
import sys
import re
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--genCC", help="genCC tsv Input file")
parser.add_argument("-o", "--output", help="Parsed genCC Output")
parser.add_argument("--high", help="Flag: Return genes/inheritance for strong/definitive entries only", action='store_true')
parser.add_argument("--limit5", help="Flag: Only output a 5 inheritance code if no other inheritance code", action='store_true')
parser.add_argument("--high_priority", help="Flag: Include only strong/definitive entries, unless no strong/definitive entries exist.", action='store_true')

parser.add_argument("--output_strong_list", help="Flag: Output list of genes with strong/definitive annotations.", action='store_true')
parser.add_argument("--output_all", help="Flag: Output all genCC genes", action='store_true')

parser.set_defaults(high=False)
parser.set_defaults(limit5=False)
parser.set_defaults(high_priority=False)

parser.set_defaults(output_strong_list=False)
parser.set_defaults(output_all=False)

args = parser.parse_args()

genCCFile = open(args.genCC, mode = 'r')
genCCLines = genCCFile.readlines()
genCCFile.close()

high = args.high
limit5 = args.limit5
high_priority = args.high_priority
output_strong = args.output_strong_list
output_all = args.output_all

#open output files
outputFile = open(args.output, mode = 'w')

if output_strong == True: 
    strongFile = open("genCC_strongDef_genes.txt", mode = 'w')

if output_all == True: 
    allFile = open("genCC_all_genes.txt", mode = 'w')

#Dictionary with genenames as keys (so no duplicates), and separate for strong/definitive entries only
geneDict = {}
geneDict_high = {}

evidenceList = ["Definitive", "Strong"]

#Possible inheritances

#NB X-linked ONLY is lumped with other
inheritanceDict = {"Autosomal dominant":1, "Autosomal dominant inheritance":1, "Autosomal dominant inheritance with maternal imprinting":1, "Autosomal dominant inheritance with maternal imprinting HP:0012275":1, "Autosomal dominant inheritance with paternal imprinting":1, "AD":1, "Autosomal recessive":2, "Autosomal recessive inheritance":2, "AR":2, "X-linked dominant inheritance":3, "X-linked dominant":3, "X-linked recessive inheritance":4, "X-linked recessive":4, "X-linked inheritance":5, "Unknown":5, "Somatic Somatic mosaicismism":5, "Digenic inheritance":5, "Somatic mosaicism":5, "Unknown inheritance":5, "Mitochondrial inheritance":5, "X-linked":5, "Semi-dominant mode of inheritance":5, "(Literally mode of inheritance but we map to) Unknown inheritance":5, "Y-linked inheritance":5, "Semidominant mode of inheritance":5, "Semidominant":5, "Mitochondrial":5, "Digenic inheritance HP:0010984":5, "":5}

lineCounter = 0

outputFile.write("approvedGeneSymbol"+"\t"+"inheritance_code" + "\n")

for line in genCCLines:
    # Skip header line and count lines
    if lineCounter == 0:
        lineCounter = lineCounter+1
        continue
    else:
        lineCounter = lineCounter+1

    # Strip trailing new line
    line = line.strip('\n')

    # Get the values
    valueList = line.split('\t')

    if len(valueList)<11:
        continue

    # Get the fields
    geneSymbol = valueList[2].strip().strip("\"")
    evidence = valueList[8].strip().strip("\"")
    inheritance = valueList[10].strip().strip("\"")

    #Code inheritance
    #1 AD, 2 AR, 3 XLD, 4 XLR, 5 other
    #Looks to see if already has inheritance type in dictionary, if not, adds.

#Gene dict high will have genes and inheritances with high or definitive evidence. Gene Dict will have all genCC genes/inheritances.

    if evidence in evidenceList: 
        if geneSymbol in geneDict_high.keys():
            geneDict_high[geneSymbol].add(inheritanceDict[inheritance])
        else: 
            geneDict_high[geneSymbol] = {inheritanceDict[inheritance]}
        
    if geneSymbol in geneDict.keys():
        geneDict[geneSymbol].add(inheritanceDict[inheritance])
    else: 
        geneDict[geneSymbol] = {inheritanceDict[inheritance]}


    #Keys: 
    #1 AD only
    #2 AR only
    #3 XLD only
    #4 XLR only
    #5 other only
    #12 AD, AD
    #15 AD, other
    #25 AR, other
    #125 AD, AR, other
    #34 XLD, XLR
    #35 XLD, other
    #45 XLR, other
    #345 XLD, XLR, other

high_keys = sorted(geneDict_high.keys())

if high_priority == True and high == False: 
    for key in sorted(geneDict.keys()):
        if key in high_keys: 
            string_list = [str(item) for item in sorted(list(geneDict_high[key]))]
        else:
            string_list = [str(item) for item in sorted(list(geneDict[key]))]

        if limit5==True and len(string_list)>1 and ("5" in string_list):
            string_list.remove("5")
        
        string_inh = "".join(string_list)
    
        outputFile.write(key + "\t" + string_inh + "\n")

elif (high_priority == True and high == True) or (high_priority == False and high == True):
        for key in sorted(geneDict.keys()):
            if key in high_keys: 
                string_list = [str(item) for item in sorted(list(geneDict_high[key]))]

            if limit5==True and len(string_list)>1 and ("5" in string_list):
                string_list.remove("5")
        
        string_inh = "".join(string_list)
    
        outputFile.write(key + "\t" + string_inh + "\n")

else: 
    for key in sorted(geneDict.keys()):
        string_list = [str(item) for item in sorted(list(geneDict[key]))]

        if limit5==True and len(string_list)>1 and ("5" in string_list):
           string_list.remove("5")
        
        string_inh = "".join(string_list)
    
        outputFile.write(key + "\t" + string_inh + "\n")

if output_strong == True: 
    for key in sorted(geneDict_high.keys()):
        strongFile.write(key + "\n")

if output_all == True: 
    for key in sorted(geneDict.keys()):
        allFile.write(key + "\n")

outputFile.close()
exit()