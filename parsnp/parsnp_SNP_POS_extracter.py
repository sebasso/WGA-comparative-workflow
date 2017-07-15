#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os

"""
while line startswith ##
	build mapping between sequenceindex and sequence file && capture heading
"""

def parsesnps():

	inputSnps=""
	with open(sys.argv[1],"r") as f:
		inputSnps = f.readlines()

	iters = int(inputSnps[1].split(" ")[1])
	seq_mapping = {}

	#sequencid mapping with filename and identifier
	for i in xrange(2,iters*4,4):
		seq_mapping[inputSnps[i].split(" ")[1].rstrip()] = inputSnps[i+1].split(" ")[1].rstrip()+"\t"+inputSnps[i+2].split(" ")[1].rstrip()

	startpos = 0
	#endpos = 0
	snpid = None
	snps = {}
	poscount = 0
	#snps["1"] = list([POS, SNP],[])

	for j in xrange(i+5,len(inputSnps)):
		if inputSnps[j].startswith(">"):
			poscount = 0
			header = inputSnps[j].split(" ")[0]
			header = header.split("-")
			snpid = header[0].split(":")[0][1:]
			startpos = header[0].split(":")[1]

			"""sys.stderr.write("\nsnpID:\n")
			sys.stderr.write(snpid)
			sys.stderr.write("\n")
			sys.stderr.write("\startpos:\n")
			sys.stderr.write(startpos)
			sys.stderr.write("\n")"""

			if not snpid in snps:
				snps[snpid] = []

		elif inputSnps[j].startswith("="):#done with xmfa file
			break
		else:
			for char in inputSnps[j]:#linewise
				poscount += 1
				if char.isupper():
					position = poscount+int(startpos)
					mock = []
					mock.append(str(position))
					mock.append(char)
					snps[snpid].append(mock)
				else:
					continue


	outputFormattedSnps = []
	sep="\t"
	for key,val in snps.items():
		for entry in val:
			placeholder=seq_mapping[key].split("\t")
			if ".ref" in placeholder[0]:
				placeholder[0] = placeholder[0][:-4]
			output = placeholder[0]+sep+entry[0]+sep+entry[1]+sep+placeholder[1]
			outputFormattedSnps.append(output)

	outputfile=sys.argv[2]+"/parsnp_SNPs_POS_formatted.tsv"

	#sort on SNP position
	outputFormattedSnps = sorted(outputFormattedSnps, key=lambda x: int(x.split("\t")[1]))

	with open(outputfile,"w") as f:
		f.write("%s\n" %"parsnp")
		for lineoutput in outputFormattedSnps:
			f.write("%s\n" %lineoutput)



if __name__ == '__main__':
	print sys.argv
	#usage: sys.argv[1] = xmfa file from parsnp output
	parsesnps()
