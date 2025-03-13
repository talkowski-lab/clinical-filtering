#!/usr/bin/env RScript
# -*- coding: utf-8 -*-

#Code for making the plot

#Load upset library
library(UpSetR)

#Set working directory where you have "R_output.tsv"
setwd("working directory")
genes <- read.csv("R_output.tsv", header = TRUE, sep = "\t")
upset(genes, sets = colnames(genes), sets.bar.color = "#56B4E9", order.by="freq", empty.intersection = "on")