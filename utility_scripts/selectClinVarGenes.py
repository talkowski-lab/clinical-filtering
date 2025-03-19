#!/usr/bin/env python
# -*- coding: utf-8 -*-

#18 March 2025 -- M. Duyzend
#
# Counts the number of ClinVar variants in the categories below, and outputs the set of variants as bed files. 
#
# Inputs: 
#-g GeneList of Interest (single column)
#-i Input ClinVar.tsv
#Usage: preprocess using, bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/ALLELEID\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\t%INFO/CLNSIGCONF\t%INFO/GENEINFO\n' -H clinvar_20250202.vcf.gz | python selectClinVarGenes.py -g GeneList.txt -o output -c Counts.txt
#-o output_root --> will make several outputs (.tsv): 
#	-All Variants (output_all.txt)
#	-Not Benign/Likely Benign **(takes "maximum" classification; for example would include uncertain|benign/likely benign)**
#	-Any P/LP variant count **(takes any variants that has a P/LP classification; even if other classifications)**
#	-All P/LP 1*+ (output_PLP_1.txt) [including conflicting]
#	-All P/LP 2*+ (output_PLP_2.txt)
#-c counts_output.tsv (outputs .tsv with the following columns)
#	1. Gene Name
#	2. All Variant Count
#	3. Not Benign/Likely Benign Variant Count (takes "maximum" classification; for example would include uncertain|benign/likely benign)**
#	4. Any P/LP variant count **(takes any variants that has a P/LP classification; even if other classifications)**
#	5. All P/LP 1*+ Variant Count
#	6. All P/LP 2*+ Variant Count
#   

#NB: If multiple genes, will output counts for variants associated with either gene

# Imports
import sys
import re
import argparse
import gzip
import os

def is_gzipped(filepath):
    """Checks if a file is gzipped by attempting to read its header.

    Args:
        filepath: The path to the file.

    Returns:
        True if the file appears to be gzipped, False otherwise.
        Returns False if the file does not exist.
    """
    if not os.path.exists(filepath):
        return False

    try:
        with gzip.open(filepath, 'rb') as f:
            f.peek(1)  # Read at least one byte without consuming it
        return True
    except gzip.BadGzipFile:
        return False
    except Exception:  # Handle other potential errors like file not found or permission issues
        return False

def check_any_member_present(list1, list2):
    """
    Checks if any element of list1 is present in list2.

    Args:
        list1: The first list.
        list2: The second list.

    Returns:
        True if at least one element of list1 is in list2, False otherwise.
    """
    for element in list1:
        if element in list2:
            return True
    return False

def pathogenicity(criteria):
	"""Takes in a list object of pathogenicity criteria and provides a determination of maximum pathogenicity. 

	Args:
		criteria: list object of criteria

	Returns: 
		String object of the pathogenicity. Outputs: PLP, VUS, BLB, other
	"""

	path_list = ["Pathogenic", "Pathogenic/Likely_pathogenic", "Pathogenic/Likely_pathogenic/Likely_risk_allele", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|other", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|risk_factor", "Pathogenic/Likely_pathogenic|association", "Pathogenic/Likely_pathogenic|other", "Pathogenic/Likely_pathogenic|risk_factor", "Pathogenic/Likely_risk_allele", "Pathogenic/Likely_risk_allele|risk_factor", "Pathogenic/Pathogenic,_low_penetrance|other", "Pathogenic/Pathogenic,_low_penetrance|other|risk_factor", "Pathogenic|Affects", "Pathogenic|association", "Pathogenic|association|protective", "Pathogenic|confers_sensitivity", "Pathogenic|drug_response", "Pathogenic|other", "Pathogenic|protective", "Pathogenic|risk_factor", "Likely_pathogenic", "Likely_pathogenic,_low_penetrance", "Likely_pathogenic/Likely_risk_allele", "Likely_pathogenic|Affects", "Likely_pathogenic|association", "Likely_pathogenic|drug_response", "Likely_pathogenic|other", "Likely_pathogenic|risk_factor"]
	benign_list = ["Benign", "Benign/Likely_benign", "Benign/Likely_benign|association", "Benign/Likely_benign|drug_response", "Benign/Likely_benign|drug_response|other", "Benign/Likely_benign|other", "Benign/Likely_benign|other|risk_factor", "Benign/Likely_benign|risk_factor", "Benign|Affects", "Benign|Affects|association|other", "Benign|association", "Benign|confers_sensitivity", "Benign|drug_response", "Benign|other", "Benign|protective", "Benign|risk_factor", "Likely_benign", "Likely_benign|Affects|association", "Likely_benign|association", "Likely_benign|drug_response", "Likely_benign|drug_response|other", "Likely_benign|other", "Likely_benign|protective", "Likely_benign|risk_factor"]
	vus_list = ["Conflicting_classifications_of_pathogenicity", "Conflicting_classifications_of_pathogenicity|Affects", "Conflicting_classifications_of_pathogenicity|association", "Conflicting_classifications_of_pathogenicity|association|risk_factor", "Conflicting_classifications_of_pathogenicity|drug_response", "Conflicting_classifications_of_pathogenicity|drug_response|other", "Conflicting_classifications_of_pathogenicity|other", "Conflicting_classifications_of_pathogenicity|other|risk_factor", "Conflicting_classifications_of_pathogenicity|protective", "Conflicting_classifications_of_pathogenicity|risk_factor", "Uncertain_risk_allele", "Uncertain_risk_allele|risk_factor", "Uncertain_significance", "Uncertain_significance/Uncertain_risk_allele", "Uncertain_significance|association", "Uncertain_significance|drug_response", "Uncertain_significance|other", "Uncertain_significance|risk_factor"]

	if check_any_member_present(criteria, path_list): 
		return "PLP"
	elif check_any_member_present(criteria, vus_list): 
		return "VUS"
	elif check_any_member_present(criteria, benign_list):
		return "BLB"
	else: 
		return "other"

def clinsig(significance):
	"""Takes in a string object of review criteria and provides a determination. 

	Args:
		significance: string object of significance

	Returns: 
		String object of the significance. Outputs: two_star_plus, one_star, no_star
	"""

	two_star_list = ["practice_guideline", "reviewed_by_expert_panel", "criteria_provided,_multiple_submitters,_no_conflicts"]
	one_star_list = ["criteria_provided,_single_submitter", "criteria_provided,_conflicting_classifications"]
	
	if significance in two_star_list:
		return "two_star_plus"
	elif significance in one_star_list:
		return "one_star"
	else: 
		return "no_star"

parser = argparse.ArgumentParser()

#-i INTERVAL_LIST -g GENEFILE -o OUTPUT
parser.add_argument("-g", "--Gene", help="Gene List Input (Optional)")
parser.add_argument("-i", "--input", help="ClinVar Input TSV (optional, otherwise takes from STDIN)")
parser.add_argument("-o", "--output", help="Output Root")
parser.add_argument('--genelist', action='store_true', help='Gene List Present')
parser.add_argument('--outputbed', action='store_true', help='Output Bed Files')

args = parser.parse_args()

#True if gene list
flag = args.genelist
bedflag = args.outputbed

#Take input from stdin (uncompressed) or compressed or uncompressed tsv if input is given
if args.input:
	if is_gzipped(args.input):
		inputFile = gzip.open(args.input, mode = 'rt')
	else: 
		inputFile = open(args.input, mode = 'r')
else: 
	inputFile = sys.stdin

#Make dictionary of genes
#Value is counters of the columns of counts
#
#	1. Gene Name
#	2. All Variant Count
#	3. Not Benign/Likely Benign Variant Count (takes "maximum" classification; for example would include uncertain|benign/likely benign)**
#	4. Any P/LP variant count **(takes any variants that has a P/LP classification; even if other classifications)**
#	5. All P/LP 1*+ Variant Count
#	6. All P/LP 2*+ Variant Count
#   

geneDict = {}

if flag == True: 
	with open(args.Gene, mode = 'r') as file:
		for line in file:
			geneDict[line.strip()]=[0,0,0,0,0]

if bedflag == True: 
	ending = ".bed"
else: 
	ending = ".txt"

#Set of all genes in ClinVar (to check missed genes)
geneSet = set()

#Open variant level output files
All_File = open(args.output + "_all.txt", mode = 'w')
NB_File = open(args.output + "_non_benign" + ending, mode = 'w')
Any_PLP_File = open(args.output + "_PLP_any" + ending, mode = 'w')
All_PLP_1_File = open(args.output + "_PLP_1" + ending, mode = 'w')
All_PLP_2_File = open(args.output + "_PLP_2+" + ending, mode = 'w')

#counters for variants types

for line in inputFile: 
	line = line.strip()
	if line.startswith("#"):
		headerList=line.strip('\n').split('\t')

		#Remove numbers from header
		for i in range(0, len(headerList)):
			headerList[i]=headerList[i].strip().split("]")[1]
			header = "\t".join(headerList)

		if bedflag == True: 
			header = "chrom" + "\t" + "chromStart" + "\t" + "chromEnd"

		All_File.write(header + "\n")
		NB_File.write(header + "\n")
		Any_PLP_File.write(header + "\n")
		All_PLP_1_File.write(header + "\n")
		All_PLP_2_File.write(header + "\n")

		#add other files
	else:
		entryList=line.strip('\n').split('\t')
		clnrev = entryList[headerList.index("CLNREVSTAT")].strip()
		clnsig = entryList[headerList.index("CLNSIG")].strip().split("|")
		genes = entryList[headerList.index("GENEINFO")].strip().split("|")

		path = pathogenicity(clnsig)
		stars = clinsig(clnrev)

	#Count numbers of variants

		#Loop over the genes that the variant falls in
		gene_flag = "N"
		for i in range(len(genes)):
			genes[i] = genes[i].split(":")[0]
			gene = genes[i]
			geneSet.add(gene)

			if flag==False and gene not in geneDict.keys(): 
				geneDict[gene]=[0,0,0,0,0]

		#If a variant is present in more than one gene, count that variant in each gene
			if gene in geneDict.keys():
				gene_flag = "Y"
				#count all variants
				geneDict[gene][0]=geneDict[gene][0]+1
				if path!="BLB":
					#count all non-benign variants
					geneDict[gene][1]=geneDict[gene][1]+1
				if path=="PLP":
					#Count any P/LP variants
					geneDict[gene][2]=geneDict[gene][2]+1
					#Count all 1*+ P/:P variants
					if stars=="two_star_plus" or stars=="one_star":
						geneDict[gene][3]=geneDict[gene][3]+1
					#Write all 2*+ P/LP variants
					if stars=="two_star_plus":
						geneDict[gene][4]=geneDict[gene][4]+1

		#If a variant is in more than one gene, only output that variant once

		if bedflag == True: 
			line_list = line.split("\t")
			chrom = "chr" + str(line_list[0].strip())
			chromStart = str(int(line_list[1].strip())-1)
			chromEnd = line_list[1].strip()
			line = chrom + "\t" + chromStart + "\t" + chromEnd

		if gene_flag=="Y":
			#write all variants
			All_File.write(line + "\n")
			if path!="BLB":
				#write all non-benign variants
				NB_File.write(line + "\n")
			if path=="PLP":
				#Write any P/LP variants
				Any_PLP_File.write(line + "\n")
				#Write all P/LP 1*+ variants
				if stars=="two_star_plus" or stars=="one_star":
					All_PLP_1_File.write(line + "\n")
				#Write all P/LP 2*+ variants
				if stars=="two_star_plus":
					All_PLP_2_File.write(line + "\n")

All_File.close()
NB_File.close()
Any_PLP_File.close()
All_PLP_1_File.close()
All_PLP_2_File.close()

#Open counts output file
counts_file  = open(args.output + "_counts.txt", mode = 'w')
error_file = open(args.output + "_no_count_error.txt", mode = 'w')

counts_file.write("Gene" + "\t" + "All" + "\t" + "Not_BLB" + "\t" "PLP_Any" + "\t" + "PLP_1*+" + "\t" + "PLP_2*+" + "\n")

#All genes in clinVar
geneList = list(geneSet)

#write counts file for gene of interest
#error if not in total set of clinvar genes

for gene in geneDict.keys():
	counts_file.write(gene + "\t" + str(geneDict[gene][0]) + "\t" + str(geneDict[gene][1]) + "\t" + str(geneDict[gene][2]) + "\t" + str(geneDict[gene][3]) + "\t" + str(geneDict[gene][4]) + "\n")
	if gene not in geneList:
		error_file.write(gene + "\n")

counts_file.close()
error_file.close()