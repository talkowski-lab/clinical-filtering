#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
# This is a simple script to: 
# Filter the ClinVar TSV small variant output from NIFS
# Based on user-based criteria
# Also strips extra spaces and text

#Adds tiers to outputs

#Additional Filters
#No SEGDUP, STR, SIMPLEREP
#Flag CA0 and CA2 only (except XLD)
#Genotype quality >20
#ECNT<=5
#NCOUNT<=10 (divided by depth, percentage less than 5%)

#Tier 1
#ClinVar only P/LP, AD/XLD only, CA0 and CA2
#GeneLists (AD genes from these lists)
#No inframe changes
#No low/modifier
# [x] GenomicsEngland Fetal Anomalies Panel v5
# [x] HPO fetal anomaly gene list
# [x] HPO Stillbirth 
# [x] GUARDIAN v1
# [x] -1,023 penetrant genes (Ceyhan-Birsoy, GiM, 2017)
# [x] DDD strong/definitive
# [x] fetal demise
# OR clinvar 2* path (any variant)

#Tier 2:
#Same as Tier 1 but include all clusters

#Tier 3: 
#Same as above but not restricted to gene lists (including clinVar)

#(Same outputs but with no filters)

## CHANGE LOG:
'''
2/10/2025:
- merged parseCompHetTSV_NIFS.py and parseDomTSV_NIFS.py into this script (only difference is OMIM_inheritance_code filtering)
- added inheritance_map to parse dominant/recessive outputs using same script
2/25/2025:
- added 'other' category for inheritance_type and inheritance_map
2/27/2025:
- removed keeping not in ClinVar ('or clin_sig==""')
'''
###

import sys
import gzip
import re
import argparse
import os

def is_gzipped(file_path):
    try:
        with gzip.open(file_path, 'rb') as f:
            f.read(1)  # Read a small amount of data to test
        return True
    except gzip.BadGzipFile:
        return False

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

parser.add_argument("-i", "--input", help="Input_TSV")
parser.add_argument("-o", "--output", help="Output_File")
parser.add_argument("-l", "--list", help="Gene_List")
parser.add_argument("-t", "--type", help="Type (dominant, recessive, other)")

args = parser.parse_args()

#tries to open binary file, if not, then assumes uncompressed
if is_gzipped(args.input):
    input = gzip.open(args.input, mode = 'rt')
else: 
    input = open(args.input, mode = 'r')

#uncompressed output
output = open(args.output, mode="w")

# NEW 2/10/2025: added inheritance_map to parse dominant/recessive outputs using same script
# NEW 2/25/2025: added 'other' category for inheritance_type and inheritance_map
inheritance_type = args.type
inheritance_map = {'dominant': [1,3], 'recessive': [2,4], 'other': [1,2,3,4,5]}

#clinSig
#Include ClinVar P/LP
#clinSigPossible
clinSigPath = ["Likely_pathogenic", "Likely_pathogenic,_low_penetrance", "Likely_pathogenic/Likely_risk_allele", "Likely_pathogenic/Pathogenic,_low_penetrance", "Likely_pathogenic|Affects", "Likely_pathogenic|association", "Likely_pathogenic|drug_response", "Likely_pathogenic|other", "Likely_pathogenic|protective", "Likely_pathogenic|risk_factor", "Pathogenic", "Pathogenic/Likely_pathogenic", "Pathogenic/Likely_pathogenic/Likely_risk_allele", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance/Established_risk_allele", "Pathogenic/Likely_pathogenic/Pathogenic,_low_penetrance|other", "Pathogenic/Likely_pathogenic|drug_response", "Pathogenic/Likely_pathogenic|other", "Pathogenic/Likely_pathogenic|risk_factor", "Pathogenic/Likely_risk_allele", "Pathogenic/Likely_risk_allele|risk_factor", "Pathogenic/Pathogenic,_low_penetrance|other", "Pathogenic/Pathogenic,_low_penetrance|other|risk_factor", "Pathogenic|Affects", "Pathogenic|association", "Pathogenic|association|protective", "Pathogenic|confers_sensitivity", "Pathogenic|drug_response", "Pathogenic|other", "Pathogenic|protective", "Pathogenic|risk_factor"]
clinSigUncertain = ["Conflicting_classifications_of_pathogenicity", "Conflicting_classifications_of_pathogenicity|Affects", "Conflicting_classifications_of_pathogenicity|association", "Conflicting_classifications_of_pathogenicity|association|risk_factor", "Conflicting_classifications_of_pathogenicity|drug_response", "Conflicting_classifications_of_pathogenicity|drug_response|other", "Conflicting_classifications_of_pathogenicity|other", "Conflicting_classifications_of_pathogenicity|other|risk_factor", "Conflicting_classifications_of_pathogenicity|protective", "Conflicting_classifications_of_pathogenicity|risk_factor", "Uncertain_significance", "Uncertain_significance/Uncertain_risk_allele", "Uncertain_significance|association", "Uncertain_significance|drug_response", "Uncertain_significance|other", "Uncertain_significance|risk_factor"]
clinSigIncludeP = clinSigPath
clinSigIncludeV = clinSigPath + clinSigUncertain

#clinRev 
#Include 2*+ ClinVar
#clinRevPossible = ["criteria_provided,_conflicting_classifications", "criteria_provided,_multiple_submitters,_no_conflicts", "criteria_provided,_single_submitter", "no_assertion_criteria_provided", "no_classification_for_the_single_variant", "no_classification_provided", "no_classifications_from_unflagged_records", "practice_guideline", "reviewed_by_expert_panel"]
clinRevOnePlus = ["criteria_provided,_multiple_submitters,_no_conflicts", "criteria_provided,_single_submitter", "practice_guideline", "reviewed_by_expert_panel"]
clinRevUncertain = ["criteria_provided,_conflicting_classifications"]
clinRevIncludeP = clinRevOnePlus
clinRevIncludeV = clinRevOnePlus + clinRevUncertain

#geneList
#geneLists to include

geneListInclude = ["DDD_Strong_Definitive_v6_24_C", "Fetal_Anomalies_NHS_v5_C", "GenCC_DefStrong_6_1_24_C", "Guardian_v1_List_C", "0003826_0005268_0034241_Combined_4_24_C", "HP_0034057_1_3_2025_C", "Newborn_Screen_Rehm_954_C"]

#Genotype Include
proGTInclude = ["0/1", "1/1", "0|1", "1|1"]
matGTInclude = ["0/0", "0/1"]
matGTIncludeX = ["0/1", "0|1", "1|1"]

#Truth possibility
truthList = ["TRUE", "True", True]

#Flags
clinVarFlagP = "N"
clinVarFlagV = "N"
geneListFlag = "N"
filterFlag = "N"
filterFlagNoGQ = "N"
mainFilterFlag = "N"
GTFlag = "N"
GTFlagMat = "N"
impactFlag = "N"
CSQFlag = "N"
IDFlag = "N"
privateFlag = "N"

Tier1Flag = "N"
Tier2FlagQ ="N"
Tier2FlagP = "N"
Tier2FlagM = "N"
Tier2Flag = "N"
Tier3Flag = "N"
Tier4Flag = "N"

#Extra gene lists
eGeneListFlag = "N"
counter = 0
geneDict = {}
gene_List = open(args.list, 'r')
for line in gene_List.readlines():
    line = line.strip("\n")
    if counter==0:
        counter=1
        continue
    else: 
        columnList = line.split("\t")
        geneDict[columnList[0]]=columnList[1]

counter=0

#Read through variants and filter/assign Tiers
for line in input.readlines(): 
    #Assumes header and makes header list
    if counter==0:
        counter=1
        headerList=line.strip('\n').split('\t')
        ##
        ##output.write("locus" + "\t" + "clinVarFlag" + "\t" + "geneListFlag" + "\t" + "filterFlag" + "\t" + "mainFilterFlag" + "\t" + "GTFlag" + "\t" + "impactFlag" + "\t" + "CSQFlag" + "\t" + "eGeneListFlag" + "\t" + "privateFlag" + "\n")
        ##
        #output.write(line.strip('\n') + "\t" + "Tier" + "\t" + "Tier Group" + "\n")
        output.write("locus" + "\t" + "alleles" + "\t" + "Tier" + "\t" + "CLNSIG" + "\t" + "CLNREVSTAT" + "\t" + "Consequence" + "\t" + "IMPACT" + "\t"+ "SYMBOL" + "\t" + "HGVSc" + "\t" + "HGVSp" + "\t" + "\t" + "\t" + line.strip('\n') + "\t" + "Tier" + "\t" + "Tier Group" + "\n")
        continue
    else: 
        entryList=line.strip('\n').split('\t')
        #strip quotes etc. from every entry
        for i in range(0, len(entryList)):
            entryList[i]=entryList[i].strip('\n').replace('\"','').replace('[','').replace(']','')

        #Overall Filters
        segdup = entryList[headerList.index("info.SEGDUP")]
        shortTandem = entryList[headerList.index("info.STR")]
        simpleRep = entryList[headerList.index("info.SIMPLEREP")]
        GQ = entryList[headerList.index("proband_entry.GQ")]
        ECNT = entryList[headerList.index("info.ECNT")]
        NCount = entryList[headerList.index("info.NCount")]
        probandDP = entryList[headerList.index("proband_entry.DP")]

        try:
            float(GQ)
        except:
            GQ =100

        try: 
            float(ECNT)
        except: 
            ECNT=0

        try: 
            float(NCount)
        except:
            NCount=1

        try:
            float(probandDP)
        except:
            probandDP=1000
        
        if segdup in truthList or float(GQ)<20 or float(ECNT)>5 or float(NCount)/float(probandDP)>0.05:
            filterFlag="N"
        else: 
            filterFlag="Y"

        if segdup in truthList or float(ECNT)>5 or float(NCount)/float(probandDP)>0.05:
            filterFlagNoGQ="N"
        else: 
            filterFlagNoGQ="Y"

        #Preset Filters
        mainFilter = entryList[headerList.index("filters")]

        #if mainFilter != "":
        if mainFilter != "" and mainFilter != "ABBINOM":
            mainFilterFlag="N"
        else: 
            mainFilterFlag="Y"

        #Impact
        impact = entryList[headerList.index("vep.transcript_consequences.IMPACT")]
        if impact == "HIGH" or impact == "MODERATE":
            impactFlag = "Y"
        else:
            impactFlag = "N"

        #Consequence
        consequence = entryList[headerList.index("vep.transcript_consequences.Consequence")]
        CSQList = consequence.strip().split("_")

        if "inframe" in CSQList:
            CSQFlag = "N"
        else: 
            CSQFlag = "Y"

        #Genotypes

        probandGT = entryList[headerList.index("proband_entry.GT")]
        motherGT = entryList[headerList.index("mother_entry.GT")]
        chrom = entryList[headerList.index("locus")].strip().split(":")[0]

        if (probandGT in proGTInclude and motherGT in matGTInclude) or (probandGT in proGTInclude and motherGT in matGTIncludeX and chrom=="chrX"): 
            GTFlag = "Y"
        else: 
            GTFlag = "N"
        ##Enter desired values to filter on here##

        if (probandGT in proGTInclude or motherGT in matGTInclude) or (probandGT in proGTInclude or motherGT in matGTIncludeX and chrom=="chrX"): 
            GTFlagMat = "Y"
        else: 
            GTFlagMat = "N"

        #CLINVAR

        clin_sig=entryList[headerList.index("info.CLNSIG")]
        clin_revstat=entryList[headerList.index("info.CLNREVSTAT")]

        # NEW 2/27/2025: Removed keeping not in ClinVar ('or clin_sig==""')
        #Clinvar P/LP with 1* evidence or not present
        if (clin_sig in clinSigIncludeP and clin_revstat in clinRevIncludeP):
            clinVarFlagP="Y"
        else:
            clinVarFlagP="N"

        #Clinvar VUS
        if (clin_sig in clinSigIncludeV and clin_revstat in clinRevIncludeV):
            clinVarFlagV="Y"
        else:
            clinVarFlagV="N"

        #Gene_List
        
        geneList=entryList[headerList.index("vep.transcript_consequences.gene_list")].strip().split("&")
        geneListFlag = "N"
        for i in range(0, len(geneList)):
            if geneList[i] in geneListInclude: 
                geneListFlag = "Y"
                continue

        #Extra gene list
        symbol = entryList[headerList.index("vep.transcript_consequences.SYMBOL")]

        if symbol in geneDict.keys():
            # NEW 2/10/2025: added inheritance_map to parse dominant/recessive outputs using same script
            if inheritance_map[inheritance_type][0] in separate_digits(int(geneDict[symbol])) or inheritance_map[inheritance_type][1] in separate_digits(int(geneDict[symbol])):
                eGeneListFlag = "Y"
            else: 
                eGeneListFlag = "N"
        else:
            eGeneListFlag = "N"

        #Private flag
        #Set to 3 given allele count
        AC = entryList[headerList.index("info.cohort_AC")]
        if int(AC)<=3:
            privateFlag = "Y"
        else: 
            privateFlag = "N"

        output_line = "\t".join(entryList)

        #Filter criteria
        #GeneList Flag Seperate
        if clinVarFlagP=="Y" and filterFlag=="Y" and mainFilterFlag=="Y" and GTFlag=="Y" and impactFlag=="Y" and eGeneListFlag=="Y" and CSQFlag=="Y" and privateFlag == "Y": 
                Tier1Flag = "Y"
        else: 
                Tier1Flag = "N"

        if clinVarFlagV=="Y" and filterFlagNoGQ == "Y" and mainFilterFlag=="Y" and GTFlagMat=="Y" and impactFlag=="Y" and eGeneListFlag=="Y" and CSQFlag=="Y" and privateFlag == "Y": 
                if filterFlag == "Y":
                    Tier2FlagQ = "Y" #Yes if high quality
                else: 
                    Tier2FlagQ = "N" 

                if clinVarFlagP == "Y":
                    Tier2FlagP == "Y" #Yes if pathogenic
                else: 
                    Tier2FlagP == "N"

                if GTFlag == "Y":
                    Tier2FlagM == "Y" #Yes if only in proband
                else: 
                    Tier2FlagM == "N"            

                Tier2Flag = "Y" #In Tier 2
        else: 
                Tier2Flag = "N"

#Tier 3 no post-hoc filters except AC. 

        if mainFilterFlag=="Y" and eGeneListFlag=="Y" and privateFlag == "Y": 
                Tier3Flag = "Y"
        else: 
                Tier3Flag = "N"

#Tier 4, no filter, everything else

        #Logic for tiers
        tierSet = set()
        modSet = set()
        tier = "3"

        if Tier1Flag == "Y":
            tierSet.add(1)

        if Tier2Flag == "Y":
            tierSet.add(2)

        if Tier2FlagQ == "N":
            modSet.add("LQ")

        if Tier2FlagM == "N" and GTFlag!="Y":
            modSet.add("Mat")

        if Tier2FlagP == "N" and clin_sig!="":
            modSet.add("VUS")

        if Tier3Flag == "Y":
            tierSet.add(3)

        tierList = list(tierSet)

        geneListFlag=="Y" 

#Consider if want to change to non-definitive inheritance too. 

        if 1 in tierList: 
            tier = "1*"
            if geneListFlag=="Y":
                tier = "1"
        elif 2 in tierList: 
            mod = ";".join(sorted(list(modSet)))
            tier ="2*" + ";" + mod
            if geneListFlag=="Y":
                tier = "2" + ";" + mod
        elif 3 in tierList:
            tier = "3*"
            if geneListFlag=="Y":
                tier = "3"
        else:
            tier = "4*"
            if geneListFlag=="Y":
                tier = "4"

        #if tier!="4" and tier!="4*":
        #    tier_group = tier
        #    output.write(output_line + "\t" + tier + "\t" + tier_group + "\n")

        tier_group = tier
        output.write(str(entryList[headerList.index("locus")]) + "\t" + str(entryList[headerList.index("alleles")]) + "\t" + tier + "\t" + str(entryList[headerList.index("info.CLNSIG")]) + "\t" + str(entryList[headerList.index("info.CLNREVSTAT")]) + "\t" + str(entryList[headerList.index("vep.transcript_consequences.Consequence")]) + "\t" + str(entryList[headerList.index("vep.transcript_consequences.IMPACT")]) + "\t"+ str(entryList[headerList.index("vep.transcript_consequences.SYMBOL")]) + "\t" + str(entryList[headerList.index("vep.transcript_consequences.HGVSc")]) + "\t" + str(entryList[headerList.index("vep.transcript_consequences.HGVSp")]) + "\t" + "\t" + "\t" + output_line + "\t" + tier + "\t" + tier_group + "\n")

'''
        #original 
        if clinVarFlagP=="Y" and geneListFlag=="Y" and filterFlag=="Y" and mainFilterFlag=="Y" and GTFlag=="Y" and impactFlag=="Y" and CSQFlag=="Y" and eGeneListFlag =="Y" and privateFlag == "Y":  
        #Test
        #if clinVarFlagP=="Y" and geneListFlag=="Y" and filterFlag=="Y" and mainFilterFlag=="Y" and GTFlag=="Y" and impactFlag=="Y" and CSQFlag=="Y" and eGeneListFlag =="Y" and privateFlag == "Y": 
            output.write(output_line + "\n")
            ##
            #output.write(entryList[headerList.index("locus")] + "\t" + clinVarFlag + "\t" + geneListFlag + "\t" + filterFlag + "\t" + mainFilterFlag + "\t" + GTFlag + "\t" + impactFlag + "\t" + CSQFlag + "\t" + eGeneListFlag + "\t" + privateFlag + "\n")
            ##
        #    continue
        else:
            ##for testing
            #output.write(entryList[headerList.index("locus")] + "\t" + clinVarFlag + "\t" + geneListFlag + "\t" + filterFlag + "\t" + mainFilterFlag + "\t" + GTFlag + "\t" + impactFlag + "\t" + CSQFlag + "\t" + eGeneListFlag + "\t" + privateFlag + "\n")
            ##
            continue
'''

input.close()
output.close()
exit()