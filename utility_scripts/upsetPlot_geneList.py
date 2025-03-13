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
import argparse
import pandas as pd

parser = argparse.ArgumentParser()

parser.add_argument("-i", "--list", help="List of Gene Names; <name of gene list>  <path to genelist>")
parser.add_argument("-o", "--output", help="R input for upset plot")

args = parser.parse_args()

geneLists = open(args.list, mode = 'r')
geneLines = geneLists.readlines()
geneLists.close()

#Dictionary with keys list name, and entries dataframes (series) of genes
listDict = {}

for line in geneLines:
	lineList = line.strip().split("\t")
	name = lineList[0]
	file = lineList[1]

	listDict[name]=pd.read_csv(file, sep="\t", usecols=[0], header=None, dtype=str)

#Make master series of genes to use as column names
masterdf = pd.Series(dtype=str)
for key in listDict.keys():
	masterdf = pd.concat([masterdf, listDict[key]], axis=0)

#remove duplicates
masterdf = masterdf.sort_values(by=[0]).drop_duplicates().reset_index(drop=True)

#clone masterdf to use as rownames; colnames are listDict.keys()
bool_df = masterdf
colnames = list(listDict.keys())

#Make dataframe of bools, with first column as row names
for key in listDict.keys():
	bool_df = pd.concat([bool_df, masterdf[0].isin(listDict[key][0])], axis=1)

#Renumber the columns and remove gene names
bool_df.columns = range(len(bool_df.columns))
bool_df=bool_df.drop(bool_df.columns[0], axis=1)

#Assign gene names as the row index, after first renumbering column index
#bool_df = bool_df.set_index([0], drop=True)

#Assign column names
bool_df = bool_df.set_axis(colnames, axis=1)

#Change to 0/1 for R input
for col in bool_df.columns:
    bool_df[col] = bool_df[col].astype(int)

#write to CSV
bool_df.to_csv(args.output, index=False, sep="\t")