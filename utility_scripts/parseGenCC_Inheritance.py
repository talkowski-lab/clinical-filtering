#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# This takes a genCC TSV File
# and creates a file with unique gene names and inheritance codes
# 12 March 2025
#

# Imports
import sys
import re
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("-g", "--genCC", help="genCC tsv Input file")
parser.add_argument("-o", "--output", help="Parsed genCC Output")
parser.add_argument("--high", help="Add flag if only want to return genes/inheritance for strong/definitive entries", action='store_true')
parser.set_defaults(high=False)

args = parser.parse_args()

genCCFile = open(args.genCC, mode = 'r')
genCCLines = genCCFile.readlines()
genCCFile.close()

high - args.high

outputFile = open(args.output, mode = 'w')

#Dictionary with genenames as keys (so no duplicates)
geneDict = {}

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
    #Looks to see if already has inheritance type in dictionary, if not, adds

    #Adds individual inheritance code to a set associated with the gene symbol from each. 
    #Does not add if not sufficient evidence AND high flag is marked. 

    if high == True and evidence in evidenceList:
        if geneSymbol in geneDict.keys():
            geneDict[geneSymbol].add(inheritanceDict[inheritance])
        else: 
            geneDict[geneSymbol] = {inheritanceDict[inheritance]}

    if high == False: 
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

for key in sorted(geneDict.keys()):
    string_list = [str(item) for item in sorted(list(geneDict[key]))]
    string_inh = "".join(string_list)
    outputFile.write(key + "\t" + string_inh + "\n")

outputFile.close()
exit()