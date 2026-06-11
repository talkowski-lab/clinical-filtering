#!/usr/bin/env python
# -*- coding: utf-8 -*-


#
# This is a script that takes a list of genes, adds inheritances (combined inheritance list), adds max(evidence) from genCC file (genCC.tsv) or OMIM file (geneMap), 
# and adds 6 if inheritance not in inheritance sheet. 
#
# Inputs: 
# Master list of genes file
# Combined OMIM_genCC inheritance
# genCC raw file
# OMIM genemap file
# 
# Outputs: 
# Master list file
# (OMIM genes)
# (GenCC genes)
#
# 11 July 2025
#
# The file can downloaded from https://omim.org/downloads
# (registration required).
#


# Imports
import sys
from pathlib import Path
import re
import argparse
import pandas as pd
from io import StringIO

parser = argparse.ArgumentParser()

parser.add_argument("-m", "--masterList", help="Path to folder of gene lists")
parser.add_argument("-g", "--genCCfile", help="genCC File")
parser.add_argument("-u", "--OMIMfile", help="OMIM geneMap File")
parser.add_argument("-i", "--inheritanceFile", help="Inheritance File")
parser.add_argument("-o", "--outputFile", help="Output File")

parser.add_argument("--omim_output", help="Flag: Output OMIM Gene List", action='store_true')
parser.add_argument("--genCC_output", help="Flag: Output genCC Gene List", action='store_true')

parser.set_defaults(omim_output=False)
parser.set_defaults(genCC_output=False)

args = parser.parse_args()

#GenCC Classification Order
classification_order = [
    "Definitive",
    "Strong",
    "Moderate",
    "Limited",
    "Supportive",
    "Disputed Evidence",
    "No Known Disease Relationship",
    "Refuted Evidence"
]

order_rank = {c: i for i, c in enumerate(classification_order)}

##Function to get max classification

def get_max_classification(group):
    # Remove empty entries for comparison, but keep them if all are empty
    non_blank = [x for x in group if x]
    if non_blank:
        return min(non_blank, key=lambda x: order_rank.get(x, float('inf')))
    else:
        return ""  # All blank

##Generate the master gene list
folder_path = Path(args.masterList)

# Set to hold gene names (to avoid duplicates)
gene_names = set()

# Loop through all files in the folder
for file in folder_path.glob("*"):
    with file.open("r") as f:
        for line in f:
            gene = line.strip()
            if gene:  # skip empty lines
                gene_names.add(gene)

# Sort the unique gene names in the master list
sorted_genes_master = sorted(gene_names)

##Get the maximum classification for each gene
#Load genCC file
df = pd.read_csv(args.genCCfile, sep="\t", dtype=str)

# Fill missing classification_title with empty string
df["classification_title"] = df["classification_title"].fillna("")

#genCC dictionary, key gene_symbol, entry max classification title
genCC_dict = (
    df.groupby("gene_symbol")["classification_title"]
      .apply(get_max_classification)
      .to_dict()
)

##Print all genCC genes if needed
if args.genCC_output==True:
    with open("genCC_Genes_all.txt", "w") as out:
        for gene in genCC_dict.keys():
            out.write(gene + "\n")

##Open the inheritance file and make a dictionary, gene name/inheritance
df = pd.read_csv(args.inheritanceFile, sep="\t", dtype=str)

# Make the dictionary
inheritance_dict = dict(zip(df["approvedGeneSymbol"], df["inheritance_code"]))

##Open the OMIM file and make a unique gene list. Output if desired. 
with open(args.OMIMfile) as f:
    for idx, line in enumerate(f):
        if line.startswith("# Chromosome"):
            header_line_idx = idx
            break

with open(args.OMIMfile) as f:
    lines = f.readlines()

# Remove the leading '#' from header
lines[header_line_idx] = lines[header_line_idx].lstrip("#").strip() + "\n"

# Convert to dataframe (tab-separated)
df = pd.read_csv(StringIO("".join(lines[header_line_idx:])), sep="\t", dtype=str)

# Get non-missing Approved Gene Symbols, deduplicate, and sort
genes = df["Approved Gene Symbol"].dropna()
unique_OMIM_genes = sorted(set(genes))

#Print OMIM genes if needed
if args.omim_output==True: 
    with open("OMIM_genes_all.txt", "w") as out:
        for gene in unique_OMIM_genes:
            out.write(gene + "\n")

##For each gene in the master list, output final file, approvedGeneSymbol / inheriance_code / genCC_classification

outputFile = open(args.outputFile, mode = 'w')
outputFile.write("approvedGeneSymbol" + "\t" + "inheritance_code" + "\t" + "genCC_classification" + "\n")

for gene in sorted_genes_master:
    inheritance_code = "6"
    if gene in inheritance_dict.keys():
        inheritance_code = inheritance_dict[gene]

    genCC_classification = ""
    if gene in genCC_dict.keys():
        genCC_classification = genCC_dict[gene]

    if gene not in genCC_dict.keys() and gene in unique_OMIM_genes: 
        genCC_classification = "OMIM"

    outputFile.write(gene + "\t" + inheritance_code + "\t" + genCC_classification + "\n")

outputFile.close()