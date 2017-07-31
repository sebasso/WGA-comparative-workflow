#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os

"""
	line[2]=snp
	line[3]=position
	line[4]=strand
	line[5]=filename
	line[6]=identifier
	known position:
	0       AAAAAC.GTACTA   C       3544 F  EEE_FL93-939.fasta      gb|EF151502.1|gi|119633046|Eastern
	0       CATT.GCTA       C       616396 F        seq30_DC2_Shuffled_1.fna        Shuffled-1-DC2

	unknown position:
	0       CATT.GCTA       G       x       seq24_CD2_Shuffled_1.fna

	to:
	EEE_FL93-939.fasta      9       C       gb|EF151502.1|gi|119633046|Eastern      F
"""

def parseline(line):
	newline = ""
	sep="\t"
	line=line.split(sep)
	if len(line) == 6:
		placeholder=line[3].split(" ")
		line[5]=line[5][0:len(line[5])-1]
		newline=line[4]+sep+placeholder[0]+sep+line[2]#+sep+line[5]+sep+placeholder[1]
	elif len(line) == 5:
		newline=line[4]+sep+line[3]+sep+line[2]#+sep+line[4]

	return newline


def parsesnps():
	inputSnps=""
	with open(sys.argv[1],"r") as f:
		inputSnps = f.readlines()

	outputFormattedSnps = []

	for i in xrange(0,len(inputSnps)):
		if inputSnps[i] in '\n':
			continue
		else:
			outputFormattedSnps.append(parseline(inputSnps[i]))

	snp_groups = inputSnps[len(inputSnps)-1].split("\t")[0]

	print "num_snp_groups: "+str(int(snp_groups)+1)#starts at 0

	currpath=os.getcwd()
	outputfile=currpath+"/"+"kSNP_SNPs_POS_formatted.tsv"
	#sort on SNP position
	outputFormattedSnps = sorted(outputFormattedSnps, key=lambda x: int(x.split("\t")[1]))
	output = "kSNP\n"+str(snp_groups)+"\n"
	for lineoutput in outputFormattedSnps:
		output += lineoutput+"\n"

	with open(outputfile,"w") as f:
		f.write(output)

	exit(0)

if __name__ == '__main__':
	print sys.argv
	#usage: sys.argv[1] = SNPs_all file from kSNP output
	parsesnps()



