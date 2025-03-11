#!/usr/bin/env python
# -*- coding: utf-8 -*-


#
# This takes the output of parseGeneMap2
# and creates a file with unique gene names and inheritance codes
# 22 April 2024
#


# Imports
import sys
import re


# Read from stdin

print("approvedGeneSymbol"+"inheritance_code")

#Dictionary with genenames as keys (so no duplicates)
phenoDict = {}

#Possible inheritances
#NB: X-linked is X-linked recessive

inheritanceDict = {"Autosomal dominant":1, "?Autosomal dominant":1, "Autosomal recessive":2, "?Autosomal recessive":2, "X-linked dominant":3, "?X-linked dominant":3, "X-linked recessive":4, "?X-linked recessive":4, "Digenic dominant":5, "Digenic recessive":5, "Inherited chromosomal imbalance":5, "Isolated cases":5, "Mitochondrial":5, "Multifactorial":5, "Pseudoautosomal dominant":5, "Pseudoautosomal recessive":5, "Somatic mosaicism":5, "Somatic mutation":5, "X-linked":4, "Y-linked":5}

lineCounter = 0

for line in sys.stdin:

    # Skip comments
    if line.startswith('#'):
        continue

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

    lenList = len(valueList)

    #skip any gene without inheritance information
    if lenList<13:
        continue

    # Get the fields
    chromosome = valueList[0]
    genomicPositionStart = valueList[1]
    genomicPositionEnd = valueList[2]
    mimNumber = valueList[3]
    geneSymbols = valueList[4]
    geneName = valueList[5]
    approvedGeneSymbol = valueList[6]
    entrezGeneID = valueList[7]
    ensemblGeneID = valueList[8]
    phenotype = valueList[9]
    phenotypeMimNumber = valueList[10]
    phenotypeMappingKey = valueList[11]

    #Skip any gene without an approved gene symbol
    if approvedGeneSymbol=="":
        continue

    ###CHANGE IF NEEDED
    #Skip any provisional or suceptibility genes
    if phenotype[0]=="{" or phenotype[0]=="|" or phenotype[0]=="?" or phenotype[0]=="[":
        continue

    #Code inheritance
    #1 AD, 2 AR, 3 XLD, 4 XLR (XL), 5 other
    #Looks to see if already has inheritance type in dictionary, if not, adds

    #Empty inheritance set
    inheritance=set()

    if lenList==13:
        inheritance1 = valueList[12]

        if inheritance1 in inheritanceDict.keys():
            inheritance1=inheritanceDict[inheritance1]
        else: 
            inheritance1=5

        inheritance.add(inheritance1)

    elif lenList==14:
        inheritance1 = valueList[12]
        inheritance2 = valueList[13]

        if inheritance1 in inheritanceDict.keys():
            inheritance1=inheritanceDict[inheritance1]
        else: 
            inheritance1=5

        if inheritance2 in inheritanceDict.keys():
            inheritance2=inheritanceDict[inheritance2]
        else: 
            inheritance2=5

        inheritance.add(inheritance1)
        inheritance.add(inheritance2)

    elif lenList==15:
        inheritance1 = valueList[12]
        inheritance2 = valueList[13]
        inheritance3 = valueList[14]

        if inheritance1 in inheritanceDict.keys():
            inheritance1=inheritanceDict[inheritance1]
        else: 
            inheritance1=5

        if inheritance2 in inheritanceDict.keys():
            inheritance2=inheritanceDict[inheritance2]
        else: 
            inheritance2=5

        if inheritance3 in inheritanceDict.keys():
            inheritance3=inheritanceDict[inheritance3]
        else: 
            inheritance3=5

        inheritance.add(inheritance1)
        inheritance.add(inheritance2)
        inheritance.add(inheritance3)

    if approvedGeneSymbol not in phenoDict.keys():
        phenoDict[approvedGeneSymbol]=[approvedGeneSymbol, inheritance]
    else:   
        phenoDict[approvedGeneSymbol][1]=phenoDict[approvedGeneSymbol][1].union(inheritance)


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

for key in sorted(phenoDict.keys()):
    stringList = [str(x) for x in sorted(list(phenoDict[key][1]))]
    print(phenoDict[key][0]+"\t"+''.join(stringList))

exit()