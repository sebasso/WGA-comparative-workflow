#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os

"""line[2]=snp
	line[3]=position
	line[4]=strand
	line[5]=filename
	line[6]=indeitifier"""

def parseline(line):
	newline = ""
	sep="\t"
	line=line.split(sep)
	placeholder=line[3].split(" ")
	line[5]=line[5][0:len(line[5])-1]

	newline=line[4]+sep+placeholder[0]+sep+line[2]+sep+line[5]+sep+placeholder[1]
	return newline


def parsesnps():

	inputSnps=""
	with open(sys.argv[1],"r") as f:
		inputSnps = f.readlines()


	outputFormattedSnps = []
	for line in inputSnps:
		if line in '\n':
			continue
		else:
			outputFormattedSnps.append(parseline(line))


	currpath=os.getcwd()
	outputfile=currpath+"/"+"kSNP_SNPs_POS_formatted.tsv"

	#sort on SNP position
	outputFormattedSnps = sorted(outputFormattedSnps, key=lambda x: int(x.split("\t")[1]))

	with open(outputfile,"w") as f:
		f.write("%s\n" %"kSNP")
		for lineoutput in outputFormattedSnps:
			f.write("%s\n" %lineoutput)



if __name__ == '__main__':
	print sys.argv
	#usage: sys.argv[1] = SNPs_all file from kSNP output
	parsesnps()
