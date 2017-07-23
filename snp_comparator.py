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
            all_files.append(tmp[1:])
            tool_names.append(filename.rstrip())

    sys.stdout.write("Comparing tools: "+str(tool_names)+"\n")

    return all_files, tool_names


def find_SNPs_in_same_position(files, tool_names):
    comparisons = len(files)-1

    SNPs_same_position = {}

    if len(files)>2:
        #multiple loop
        print "not implemented yet"
        exit(0)
    else:
        print "\nstat check\n"
        for lines in files[0]:
            for lines2 in files[1]:
                same_file = 0
                l1 = lines.split("\t")
                l2 = lines2.split("\t")
                key = int(l1[1])
                if  key == int(l2[1]):

                    if key in SNPs_same_position: # key= position in genome
                        if l1[0] == l2[0]:
                            SNPs_same_position[key][tool_names[0]][l1[0]] = []
                            SNPs_same_position[key][tool_names[0]][l1[0]].append(l1[2])
                            SNPs_same_position[key][tool_names[0]][l1[0]].append("common\t"+tool_names[1])
                            SNPs_same_position[key][tool_names[1]][l2[0]] = []
                            SNPs_same_position[key][tool_names[1]][l2[0]].append(l2[2])
                            SNPs_same_position[key][tool_names[1]][l2[0]].append("common\t"+tool_names[0])
                        else:
                            SNPs_same_position[key][tool_names[0]][l1[0]] = l1[2]
                            SNPs_same_position[key][tool_names[1]][l2[0]] = l2[2]

                    else:
                        SNPs_same_position[key] = {}
                        SNPs_same_position[key][tool_names[0]] = {}
                        SNPs_same_position[key][tool_names[1]] = {}
                        if l1[0] == l2[0]:#TODO, handle this in js?
                            SNPs_same_position[key][tool_names[0]][l1[0]] = []
                            SNPs_same_position[key][tool_names[0]][l1[0]].append(l1[2])
                            SNPs_same_position[key][tool_names[0]][l1[0]].append("common\t"+tool_names[1])
                            SNPs_same_position[key][tool_names[1]][l2[0]] = []
                            SNPs_same_position[key][tool_names[1]][l2[0]].append(l2[2])
                            SNPs_same_position[key][tool_names[1]][l2[0]].append("common\t"+tool_names[0])

                        else:
                            SNPs_same_position[key][tool_names[0]][l1[0]] = l1[2]
                            SNPs_same_position[key][tool_names[1]][l2[0]] = l2[2]

                    # same filename && position

                elif int(l1[1]) > int(l2[1]):
                    break

    print "\nstat check done ---->\n"
    return SNPs_same_position

def count_snps_per_genome(snp_lists_files):
    SNPs_per_genome = []
    for f in snp_lists_files:
        genome_SNP_count = Counter(map(lambda x: x.split("\t")[0], f)) # dict response, should consider this on the count of position as well
        #genome_loci_occurence = Counter(map(lambda x: x.split("\t")[1], f))
        SNPs_per_genome.append(genome_SNP_count)

    return SNPs_per_genome


def compare_snps(outputdir, SNP_files):
    files = SNP_files

    snp_lists_files, tool_names = readfiles(files)

    stats = {}
    stats["tool_names"] = tool_names
    stats["total_snps"] = []

    ### SETTING stats:
    #total number of snps per tool = len(all_files[0]), len(all_files[1]) ... len(all_files[x])
    for i in xrange(0,len(snp_lists_files)):
        stats["total_snps"].append(len(snp_lists_files[i]))

    stats["SNPs_per_genome"] = count_snps_per_genome(snp_lists_files)

    SNPs_loci = find_SNPs_in_same_position(snp_lists_files, tool_names)
    stats["same_loci_snps"] = SNPs_loci

    r = json.dumps(stats, indent=4, encoding="utf-8", sort_keys=True)

    with open(outputdir+"/snps_stats.json","w") as f: #warning will write this relative to exection path - sys.executables
        f.write(r)

    ### Printing stats
    #for entry, value in stats.items():
    #    print 'key: ',entry,'\t value: ',value
    print "\nSnps comparator OUTPUT: \n", r, "\n"
    #SNPSs in total     #total number of snps per tool = len(all_files[0]), len(all_files[1]) ... len(all_files[x])
    #collections.counter can easily find how many snps occur at one position, however not so relevant in comparison

    """
    NB:
    1. parsnp doesnt compare all fasta files, will be different output there
    2. parsnp -> core genome always same amount of snps?
    """

if __name__ == '__main__':
    files = sys.argv[1:]
    if len(files) <= 2:
        sys.stderr.write("snp_comparator requires a minimum of 2 formatted SNP lists with their respective position and file name")
        exit(1)
    compare_snps(files[0],files[1:])

