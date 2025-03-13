#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# This takes the homo sapiens gene file from gene and from HGNC, and converts a list of genes to the appropriate format. 
# It will print the original gene name, the converted gene name. If a gene is not in the list, it will print the original gene name, and flag with an "x".
#

# Imports
import sys
import re
import argparse

parser = argparse.ArgumentParser()

#-i INTERVAL_LIST -g GENEFILE -o OUTPUT
parser.add_argument("-i", "--Gene", help="Gene Input")
parser.add_argument("-g", "--HGNC", help="HGNC Input")
parser.add_argument("-f", "--file", help="Gene List Input")
parser.add_argument("-o", "--output", help="Output File")

args = parser.parse_args()

#Make list of input gene names
Gene_File = open(args.Gene, mode = 'r')
GeneLines = Gene_File.readlines()
Gene_File.close()

HGNC_File = open(args.HGNC, mode = 'r')
HGNCLines = HGNC_File.readlines()
HGNC_File.close()

input_File = open(args.file, mode = 'r')
inputLines = input_File.readlines()
input_File.close()

output_File = open(args.output, mode = 'w')

GeneDict={}
HGNCDict={}

#Read in both gene lists

#Creates a NLM dictionary with key approved gene name

counter = 0
for line in GeneLines:
    if counter == 0:
        counter = counter+1
        continue
    line = line.strip('\n')
    valueList = line.split('\t')

    geneSymbol = valueList[10].strip()
    synonym_list = valueList[4].strip()

    if geneSymbol in GeneDict.keys():
        if synonym_list=="-":
            GeneDict[geneSymbol].add(geneSymbol)
        else:
            GeneDict[geneSymbol].union(set(synonym_list.split("|")))
            GeneDict[geneSymbol].add(geneSymbol)
    else:
        if synonym_list=="-":
            GeneDict[geneSymbol]=set([geneSymbol])
        else:
            GeneDict[geneSymbol]=set(synonym_list.split("|"))
            GeneDict[geneSymbol].add(geneSymbol)

#Creates an HGNC dictionary with key approved gene name
counter = 0
for line in HGNCLines:
    if counter == 0:
        counter = counter+1
        continue
    line = line.strip('\n')
    valueList = line.split('\t')

    geneSymbol = valueList[1].strip()
    previousName = valueList[3].strip()
    aliasSymbol = valueList[4].strip()

    if geneSymbol in HGNCDict.keys():
        HGNCDict[geneSymbol].add(geneSymbol)
        if previousName!="":
            HGNCDict[geneSymbol].add(previousName)
        if aliasSymbol!="":
            HGNCDict[geneSymbol].add(aliasSymbol)
    else:
        HGNCDict[geneSymbol]=set([geneSymbol])
        if previousName!="":
            HGNCDict[geneSymbol].add(previousName)
        if aliasSymbol!="":
            HGNCDict[geneSymbol].add(aliasSymbol)

#Sorted master gene list, dictionary of approved gene symbol and synonyms
master_keys = sorted(list(set(GeneDict.keys()).union(set(HGNCDict.keys()))))
masterDict = {}

for key in master_keys:
    if key in GeneDict.keys() and key in HGNCDict.keys():
        masterDict[key]=GeneDict[key].union(HGNCDict[key])
    if key in GeneDict.keys() and key not in HGNCDict.keys():
        masterDict[key]=GeneDict[key]
    if key in HGNCDict.keys() and key not in GeneDict.keys():
        masterDict[key]=HGNCDict[key]

#Create reverse dictionary, synonym key will return approved gene symbol
reverseDict = {}

for key in master_keys:
    #Entries for each master key
    dictList = sorted(list(masterDict[key]))

#Ambiguous if multiple aliases map to different validated symbols. 
    for item in dictList:
        if item in reverseDict.keys():
            reverseDict[item]="ambiguous"
        else:
            reverseDict[item]=key

#Compare gene list with reverseDict; if ambiguous assume gene name if in master_keys
for line in inputLines:
    line = line.strip('\n')
    valueList = line.split('\t')
    geneName = valueList[0].strip()
    if geneName in reverseDict.keys():
        if reverseDict[geneName]=="ambiguous" and geneName in master_keys:
            output_File.write(geneName + "\t" + geneName + "\t" + "\n")
        else:
            output_File.write(geneName + "\t" + reverseDict[geneName] + "\t" + "\n")
    else: 
        output_File.write(geneName + "\t" + "X" + "\n")

output_File.close()
input_File.close()
exit()