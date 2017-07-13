#!/usr/bin/env python
# -*- coding: utf-8 -*-

# usage: python snp_comparator.py comma separated list of SNP formatted files
import sys
import os
from collections import Counter
import json

"""
    further development:
    1. create own for snp_formatted files so galaxy only provides those as options
    - modify config/datatypes_conf.xml

    x. make extensible for > 3 tool comparisons

    output example:
{
  "total_snps": [
    705,
    376
  ],
  "same_loci_filename": [

  ],
  "tool_names": [
    "dataset_589.dat",
    "dataset_598.dat"
  ],
  "SNPs_per_genome": [
    {
      "dataset_4.dat": 208,
      "dataset_3.dat": 208,
      "dataset_2.dat": 204,
      "dataset_1.dat": 85
    },
    {
      "dataset_78.dat": 94,
      "dataset_77.dat": 94,
      "dataset_76.dat": 94,
      "dataset_131.dat": 94
    }
  ],
  "same_loci_snps": [

  ]
}

"""


def readfiles(files):
    all_files = []
    tool_names = []
    for f in files:
        with open(f,"r") as f:
            tmp = f.readlines()
            filename = tmp[0]
            sys.stderr.write("\n"+filename+"\n")
            all_files.append(tmp[1:])
            tool_names.append(filename)

    return all_files, tool_names

def find_SNPs_in_same_position(files):
    comparisons = len(files)-1

    SNPs_same_position = []
    SNPs_same_position_and_file = []

    if len(files)>2:
        #multiple loop
        print "not implemented yet"
        exit(0)
    else:
        for lines in files[0]:
            for lines2 in files[1]:
                l1 = lines.split("\t")
                l2 = lines2.split("\t")
                if  int(l1[1]) == int(l2[1]):
                    SNPs_same_position.append(l1[0])
                    SNPs_same_position.append(l2[0])
                    if l1[0] == l2[0]:
                        SNPs_same_position_and_file.append(l1[0])
                        SNPs_same_position_and_file.append(l2[0])
                elif int(l1[1]) > int(l2[1]):
                    break

    return (SNPs_same_position, SNPs_same_position_and_file)


def count_snps_per_genome(snp_lists_files):
    SNPs_per_genome = []
    for f in snp_lists_files:
        genome_SNP_count = Counter(map(lambda x: x.split("\t")[0], f)) # dict response, should consider this on the count of position as well
        #genome_loci_occurence = Counter(map(lambda x: x.split("\t")[1], f))
        SNPs_per_genome.append(genome_SNP_count)

    return SNPs_per_genome


def compare_snps(SNP_files):
    files = SNP_files
    #files = SNP_files
    numfiles = len(files)

    if numfiles <= 1:
        sys.stderr.write("snp_comparator requires a minimum of 2 formatted SNP lists with their respective position and file name")
        exit(1)

    snp_lists_files, tool_names = readfiles(files)

    stats = {}

    stats["tool_names"] = tool_names
    stats["total_snps"] = []

    #total number of snps per tool = len(all_files[0]), len(all_files[1]) ... len(all_files[x])
    for i in xrange(0,len(snp_lists_files)):
        stats["total_snps"].append(len(snp_lists_files[i]))

    stats["SNPs_per_genome"] = count_snps_per_genome(snp_lists_files)

    SNPs_loci, SNPs_loci_filename = find_SNPs_in_same_position(snp_lists_files)
    stats["same_loci_snps"] = SNPs_loci
    stats["same_loci_filename"] = SNPs_loci_filename

    r = json.dumps(stats)

    for entry, value in stats.items():
        print 'key: ',entry,'\t value: ',value

    with open("snps_stats.json","w") as f:
        f.write(r)

    print "json stats: \n", r
    #SNPSs in total     #total number of snps per tool = len(all_files[0]), len(all_files[1]) ... len(all_files[x])
    #collections.counter can easily find how many snps occur at one position, however not so relevant in comparison

    """
    NB:
    1. parsnp doesnt compare all fasta files, will be differen output there
    2. parsnp -> core genome always same amount of snps?
    """
if __name__ == '__main__':
    compare_snps(sys.argv[1:])
