#!/usr/bin/env python
# -*- coding: utf-8 -*-


#
# This is a simple script to parse the genemap2.txt file that
# can be downloaded from https://omim.org/
#
# The file can downloaded from https://omim.org/downloads
# (registration required).
#


# Imports
import sys
import re

# Read from stdin

print("chromosome"+"\t"+"genomicPositionStart" + "\t" + "genomicPositionEnd" + "\t" + "mimNumber"+"\t"+"geneSymbols"+"\t"+"geneName"+"\t"+"approvedGeneSymbol"+"\t"+"entrezGeneID"+"\t"+"ensemblGeneID"+"\t"+"phenotype"+"\t"+"phenotypeMim"+"\t"+"phenotypeMapping"+"\t"+"inheritance")

for line in sys.stdin:

    # Skip comments
    if line.startswith('#'):
        continue

    # Strip trailing new line
    line = line.strip('\n')

    # Get the values
    valueList = line.split('\t')

    # Get the fields
    chromosome = valueList[0]
    genomicPositionStart = valueList[1]
    genomicPositionEnd = valueList[2]
    cytoLocation = valueList[3]
    computedCytoLocation = valueList[4]
    mimNumber = valueList[5]
    geneSymbols = valueList[6]
    geneName = valueList[7]
    approvedGeneSymbol = valueList[8]
    entrezGeneID = valueList[9]
    ensemblGeneID = valueList[10]
    comments = valueList[11]
    phenotypeString = valueList[12]
    mouse = valueList[13]

    # Skip empty phenotypes
    if not phenotypeString:
        print(chromosome+"\t"+ genomicPositionStart + "\t" + genomicPositionEnd + "\t" + mimNumber+"\t"+geneSymbols+"\t"+geneName+"\t"+approvedGeneSymbol+"\t"+entrezGeneID+"\t"+ensemblGeneID+"\n", end="")
        continue

    # Parse the phenotypes
    for phenotype in phenotypeString.split(';'):

        print(chromosome+"\t"+ genomicPositionStart + "\t" + genomicPositionEnd + "\t" + mimNumber+"\t"+geneSymbols+"\t"+geneName+"\t"+approvedGeneSymbol+"\t"+entrezGeneID+"\t"+ensemblGeneID, end="")

        # Clean the phenotype
        phenotype = phenotype.strip()

        # Long phenotype
        matcher = re.match(r'^(.*),\s(\d{6})\s\((\d)\)(|, (.*))$', phenotype)
        if matcher:
            # Get the fields
            phenotype = matcher.group(1)
            phenotypeMimNumber = matcher.group(2)
            phenotypeMappingKey = matcher.group(3)
            inheritances = matcher.group(5)

            print("\t"+phenotype+"\t"+phenotypeMimNumber+"\t"+phenotypeMappingKey, end="")

            # Get the inheritances, may or may not be there
            if inheritances:
                for inheritance in inheritances.split(','):
                    inheritance = inheritance.strip()
                    print("\t"+inheritance, end="")
        # Short phenotype
        else:

            # Short phenotype
            matcher = re.match(r'^(.*)\((\d)\)(|, (.*))$', phenotype)
            if matcher:

                # Get the fields
                phenotype = matcher.group(1)
                phenotypeMappingKey = matcher.group(2)
                inheritances = matcher.group(3)

                print("\t"+phenotype+"\t"+"NA"+"\t"+phenotypeMappingKey, end="")

                # Get the inheritances, may or may not be there
                if inheritances:
                    for inheritance in inheritances.split(','):
                        inheritance = inheritance.strip()
                        print("\t"+inheritance, end="")
        print("\n", end="")