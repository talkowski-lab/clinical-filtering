#!/usr/bin/env python
# -*- coding: utf-8 -*-

#28 March 2025 -- M. Duyzend
#
# Counts the total number of ClinVar variants in the categories below, and outputs the set of variants as VCF files. 
# Now uses hail for VCF operations
#
# Inputs: 
#-g GeneList of Interest (single column text file)
#-i Input clinvar.vcf
#
#-o output_root --> will make several outputs (.vcf.gz): 
#	-All Variants (output_all.vcf.gz)
#	-Not Benign/Likely Benign (output_notBenign.vcf.gz) **(takes "maximum" classification; for example would include uncertain|benign/likely benign)**
#	-Any P/LP variant count (output_anyPLP.vcf.gz) **(takes any variants that has a P/LP classification; even if other classifications)**
#	-All P/LP 1*+ (output_PLP_1.vcf.gz) [including conflicting]
#	-All P/LP 2*+ (output_PLP_2.vcf.gz)
#
#   -counts_output.txt (outputs .tsv with the following columns)
#	1. Gene Name
#	2. All Variant Count
#	3. Not Benign/Likely Benign Variant Count (takes "maximum" classification; for example would include uncertain|benign/likely benign)**
#	4. Any P/LP variant count **(takes any variants that has a P/LP classification; even if other classifications)**
#	5. All P/LP 1*+ Variant Count
#	6. All P/LP 2*+ Variant Count
#   

#NB: If multiple genes, will output counts for variants associated with either gene
#NB: The counts the variant type, for example if different variant at the same position, will count twice. 

#Note, in hail: 
#CLNREVSTAT: array<str> --> review status
#CLINSIGCONF: array<str> --> present only if conflicting to say what conflicts are
#CLINSIG: array<str> --> clinical significants
#GENEINFO: str

# Imports
import sys
import re
import argparse
import hail as hl
import subprocess

parser = argparse.ArgumentParser()

#Inputs: -i INTERVAL_LIST -g GENEFILE -o OUTPUT
parser.add_argument("-g", "--genelist", help="Gene List Input (Optional)", nargs='?', default=None)
parser.add_argument("-i", "--input", help="ClinVar Input VCF")
parser.add_argument("-o", "--output", help="Output Root")
parser.add_argument("--vcf", help="VCF output for each category (Optional)", action='store_true')

args = parser.parse_args()

vcfflag = args.vcf

#CLNSIG; Pathogenic, benign, VUS
#Include parsing by comma as hail parses by this. 
path_list = ["Pathogenic", "Pathogenic/Likely_pathogenic", "Pathogenic/Likely_pathogenic/Likely_risk_allele", "Pathogenic/Likely_pathogenic/Pathogenic", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|other", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|risk_factor", "Pathogenic/Likely_pathogenic|association", "Pathogenic/Likely_pathogenic|other", "Pathogenic/Likely_pathogenic|risk_factor", "Pathogenic/Likely_risk_allele", "Pathogenic/Likely_risk_allele|risk_factor", "Pathogenic/Pathogenic", "Pathogenic/Pathogenic,_low_penetrance|other", "Pathogenic/Pathogenic,_low_penetrance|other|risk_factor", "Pathogenic|Affects", "Pathogenic|association", "Pathogenic|association|protective", "Pathogenic|confers_sensitivity", "Pathogenic|drug_response", "Pathogenic|other", "Pathogenic|protective", "Pathogenic|risk_factor", "Likely_pathogenic", "Likely_pathogenic", "Likely_pathogenic,_low_penetrance", "Likely_pathogenic/Likely_risk_allele", "Likely_pathogenic|Affects", "Likely_pathogenic|association", "Likely_pathogenic|drug_response", "Likely_pathogenic|other", "Likely_pathogenic|risk_factor"]
benign_list = ["Benign", "Benign/Likely_benign", "Benign/Likely_benign|association", "Benign/Likely_benign|drug_response", "Benign/Likely_benign|drug_response|other", "Benign/Likely_benign|other", "Benign/Likely_benign|other|risk_factor", "Benign/Likely_benign|risk_factor", "Benign|Affects", "Benign|Affects|association|other", "Benign|association", "Benign|confers_sensitivity", "Benign|drug_response", "Benign|other", "Benign|protective", "Benign|risk_factor", "Likely_benign", "Likely_benign|Affects|association", "Likely_benign|association", "Likely_benign|drug_response", "Likely_benign|drug_response|other", "Likely_benign|other", "Likely_benign|protective", "Likely_benign|risk_factor"]
vus_list = ["Conflicting_classifications_of_pathogenicity", "Conflicting_classifications_of_pathogenicity|Affects", "Conflicting_classifications_of_pathogenicity|association", "Conflicting_classifications_of_pathogenicity|association|risk_factor", "Conflicting_classifications_of_pathogenicity|drug_response", "Conflicting_classifications_of_pathogenicity|drug_response|other", "Conflicting_classifications_of_pathogenicity|other", "Conflicting_classifications_of_pathogenicity|other|risk_factor", "Conflicting_classifications_of_pathogenicity|protective", "Conflicting_classifications_of_pathogenicity|risk_factor", "Uncertain_risk_allele", "Uncertain_risk_allele|risk_factor", "Uncertain_significance", "Uncertain_significance/Uncertain_risk_allele", "Uncertain_significance|association", "Uncertain_significance|drug_response", "Uncertain_significance|other", "Uncertain_significance|risk_factor"]

#CLNREVSTAT; One and two star list
#Include parsing by comma as hail parses by this. 
two_star_list = ["practice_guideline", "reviewed_by_expert_panel", "criteria_provided,_multiple_submitters,_no_conflicts", "_multiple_submitters"]
one_star_list = ["criteria_provided,_single_submitter", "_single_submitter", "criteria_provided,_conflicting_classifications", "_conflicting_classifications"]
clnrevstat_two_star_plus = two_star_list
clnrevstat_one_star_plus = one_star_list+two_star_list

output_root = args.output 

#Import VCF into a hail matrix table, recoding if necessary. 
if vcfflag == True:
	recode = {f"{i}":f"chr{i}" for i in (list(range(1, 23)) + ['X', 'Y'])}
	ds = hl.import_vcf(args.input, reference_genome='GRCh38', contig_recoding=recode, force_bgz=True, skip_invalid_loci=True, array_elements_required=False)

#Make gene dict for counts
geneDict = {}

#If a gene list is provided, read through the gene list
if args.genelist is not None: 
	with open(args.genelist, mode = 'r') as file:
		for line in file:
			geneDict[line.strip()]=[0,0,0,0,0]

	if vcfflag == True:
		geneList = geneDict.keys()
		#filter vcf so only contains genes of interest
		ds = ds.filter_rows(hl.any([ds.info.GENEINFO.contains(item) for item in geneList]))


##Use hail to create filtered VCFs as output

if vcfflag == True: 
	##All variants##
	hl.export_vcf(ds, output + "_all.vcf.bgz", tabix=True)

	##Not Benign/Likely Benign Variants##
	N_BLB_ds = ds.filter_rows(hl.any([ds.info.CLNSIG.contains(category) for category in benign_list]), keep=False)
	hl.export_vcf(N_BLB_ds, output + "_not_BLB.vcf.bgz", tabix=True)

	##Any P/LP variant## 
	#Select any CLNSIGCONF or CLNSIG pathogenic in path list

	#Select any CLNSIGCONF with a pathogenic entry (true) otherwise false
	#a few entries of CLNSIGCONF with commas, but this will still capture)
	anyPLP_Filt = hl.any([ds.info.CLNSIGCONF[0].contains(category) for category in path_list])
	anyPLP_Filt = hl.if_else(hl.is_missing(anyPLP_Filt), False, anyPLP_Filt)

	#Select any CLNSIG with a pathogenic entry otherwise false
	anyPLP_Filt2 = hl.any([ds.info.CLNSIG.contains(category) for category in path_list])

	#Keep CLNSIG pathogenic or CLNSIGCONF pathogenic
	anyPLP_Filt = anyPLP_Filt | anyPLP_Filt2

	#Filter
	any_PLP_ds = ds.filter_rows(anyPLP_Filt, keep=True)
	hl.export_vcf(any_PLP_ds, output + "_any_PLP.vcf.bgz", tabix=True)


	##All P/LP 1*+ variant##
	PLP_1_filt = hl.any([ds.info.CLNREVSTAT.contains(category) for category in clnrevstat_one_star_plus])
	PLP_1_filt2 = hl.any([ds.info.CLNSIG.contains(item) for item in path_list])

	PLP_1_filt = PLP_1_filt & PLP_1_filt2

	PLP_1_ds = ds.filter_rows(PLP_1_filt)
	hl.export_vcf(PLP_1_ds, output + "_1_plus_PLP.vcf.bgz", tabix=True)

	#All P/LP 2*+ variant
	PLP_2_filt = hl.any([ds.info.CLNREVSTAT.contains(category) for category in clnrevstat_two_star_plus])
	PLP_2_filt2 = hl.any([ds.info.CLNSIG.contains(item) for item in path_list])

	PLP_2_filt = PLP_2_filt & PLP_2_filt2

	PLP_2_ds = ds.filter_rows(PLP_2_filt)
	hl.export_vcf(PLP_2_ds, output + "_2_plus_PLP.vcf.bgz", tabix=True)


##Counts
#Approach; loop over each gene in gene list and get hail counts based on those search criteria

#Open counts output file
counts_file  = open(args.output + "_counts.txt", mode = 'w')
error_file = open(args.output + "_no_count_error.txt", mode = 'w')

counts_file.write("Gene" + "\t" + "All" + "\t" + "Not_BLB" + "\t" "PLP_Any" + "\t" + "PLP_1*+" + "\t" + "PLP_2*+" + "\n")

#Get the GENE INFO, CLNREVSTAT, CLNSIG, 
command = ["bcftools", "query", "-f", "%INFO/GENEINFO\t%INFO/CLNREVSTAT\t%INFO/CLNSIG\t%INFO/CLNSIGCONF", input]
process = subprocess.Popen(command, stdout=subprocess.PIPE, text=True)

#Set of all genes in ClinVar (to check missed genes)
geneSet = set()

while True: 
	line = process.stdout.readline()
	if line.startswith("#")
		continue
	if not line:
		break
	line_array = line.strip().split("\t")

	gene_array = line_array[0].strip().split("|")
	cln_rev_stat = line_array[1].strip()
	cln_sig = line_array[2].strip()
	cln_sig_conf = line_array[3].strip()

	#Loop over the genes that the variant falls in
	for i in range(len(gene_array)):
		gene_array[i] = gene_array[i].split(":")[0]
		gene = gene_array[i]
		geneSet.add(gene)

		if args.genelist is None and gene not in geneDict.keys(): 
			geneDict[gene]=[0,0,0,0,0]

		#If a variant is present in more than one gene, count that variant in each gene
		if gene in geneDict.keys():
			#count all variants
			geneDict[gene][0]=geneDict[gene][0]+1
			if cln_sig not in benign_list:
				#count all non-benign variants
				geneDict[gene][1]=geneDict[gene][1]+1
			if cln_sig in path_list or [item in cln_sig_conf for item in path_list]:
				geneDict[gene][2]=geneDict[gene][2]+1
			if cln_sig in path_list and cln_rev_stat in clnrevstat_one_star_plus:
				geneDict[gene][3]=geneDict[gene][3]+1
			if cln_sig in path_list and cln_rev_stat in clnrevstat_two_star_plus:
				geneDict[gene][4]=geneDict[gene][4]+1

# Check for errors
if process.returncode != 0:
  	print(f"Error running bcftools: {process.stderr}")
  	exit()


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