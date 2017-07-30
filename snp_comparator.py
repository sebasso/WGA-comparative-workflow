#!/usr/bin/env python
# -*- coding: utf-8 -*-

# usage: python snp_comparator.py comma separated list of SNP formatted files
import sys
import os
from collections import Counter
import json
import re

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


def readfiles(files):#files[0]=ksnp, files[1]=parsnp
    all_files = []
    tool_names = []
    snp_group_count = []

    for f in files:
        with open(f,"r") as f:
            tmp = f.readlines()
            tool_names.append(tmp[0].rstrip())
            snp_group_count.append(tmp[1].rstrip())
            all_files.append(tmp[2:])

    sys.stdout.write("Comparing tools: "+str(tool_names)+"\n")
    print "snp group count(snp_comparator): ",snp_group_count,"\n"
    return all_files, tool_names


def find_SNPs_in_same_position(files, tool_names):
    SNPs_same_position = {}

    for i, snp_file in enumerate(files):
        for j in xrange(0,len(snp_file)):
            l1 = snp_file[j].rstrip().split("\t")
            if l1[1] != "x": #unknown position snp from ksnp, still relevant when constructing tree from alignment there it is here
                key = int(l1[1])
                if key in SNPs_same_position:
                    if tool_names[i] in SNPs_same_position[key]:
                        if l1[0] in SNPs_same_position[key][tool_names[i]]:
                            SNPs_same_position[key][tool_names[i]][l1[0]] = l1[2]
                        else:
                            SNPs_same_position[key][tool_names[i]][l1[0]] = {}
                            SNPs_same_position[key][tool_names[i]][l1[0]] = l1[2]
                    else:
                        SNPs_same_position[key][tool_names[i]] = {}
                        SNPs_same_position[key][tool_names[i]][l1[0]] = {}
                        SNPs_same_position[key][tool_names[i]][l1[0]] = l1[2]
                else:
                    SNPs_same_position[key] = {}
                    SNPs_same_position[key][tool_names[i]] = {}
                    SNPs_same_position[key][tool_names[i]][l1[0]] = {}
                    SNPs_same_position[key][tool_names[i]][l1[0]] = l1[2]
            else:
                SNPs_same_position[key][tool_names[i]][l1[0]] = l1[2]

    return SNPs_same_position



def count_snps_per_genome(snp_lists_files):
    SNPs_per_genome = []
    grouped_snp_pos = []
    for f in snp_lists_files:
        genome_SNP_count = Counter(map(lambda x: x.split("\t")[0], f))
        snp_group_count = Counter(map(lambda x: x.split("\t")[1], f)).keys()
        #genome_loci_occurence = Counter(map(lambda x: x.split("\t")[1], f))
        SNPs_per_genome.append(genome_SNP_count)
        grouped_snp_pos.append(len(snp_group_count))

    return SNPs_per_genome, grouped_snp_pos


def compare_snps(outputdir, SNP_files):
    files = SNP_files

    snp_lists_files, tool_names = readfiles(files)

    stats = {}
    stats["tool_names"] = tool_names
    stats["total_snps"] = []
    #SNP_files[0] = ".../ksnp/ksnp_output/kSNP_SNPs_POS_formatted.tsv"
    #TODO: change this for universal use, same in parsnp
    ksnp_output = SNP_files[0][:-27] #kSNP_SNPs_POS_formatted.tsv
    # grep for num_snp_groups:

    ksnp_output += "stdout"
    print "ksnp path: ",ksnp_output,"\n"
    with open(ksnp_output,"r") as f:
        match = re.search('num_snp_groups: (.*)', f.read())

    if match:
        print match.group(0)
        ksnp_snp_groups = int(match.group(0).split(" ")[1])
        
    ### SETTING stats:
    #total number of snps per tool = len(all_files[0]), len(all_files[1]) ... len(all_files[x])
    for i in xrange(0,len(snp_lists_files)):
        stats["total_snps"].append(len(snp_lists_files[i])) #total_snps is the sum of snps even in the same position
        #TODO: incorporate SNP_group which will give unique snps after position

    stats["SNPs_per_genome"], grouped_snp_pos = count_snps_per_genome(snp_lists_files)
    print "\n"
    print "grouped_snp_pos: ",grouped_snp_pos,"\n"
    SNPs_loci = find_SNPs_in_same_position(snp_lists_files, tool_names)
    stats["same_loci_snps"] = SNPs_loci

    r = json.dumps(stats, indent=4, encoding="utf-8", sort_keys=True)

    with open(outputdir+"/snps_stats.json","w") as f: #warning will write this relative to exection path - sys.executables
        f.write(r)

    ### Printing stats
    #for entry, value in stats.items():
    #    print 'key: ',entry,'\t value: ',value
    #print "\nSnps comparator OUTPUT: \n", r, "\n"
    #SNPSs in total     #total number of snps per tool = len(all_files[0]), len(all_files[1]) ... len(all_files[x])
    #collections.counter can easily find how many snps occur at one position, however not so relevant in comparison

    """
    TODO: count how many snp groups in ksnp?
        new metric: SNPs_per_genome / num_snp_groups
    """

if __name__ == '__main__':
    print "Snp comparator @args:"
    print str(sys.argv)+"\n"
    files = sys.argv[1:]
    if len(files) <= 2:
        sys.stderr.write("snp_comparator requires a minimum of 2 formatted SNP lists with their respective position and file name")
        exit(1)
    compare_snps(files[0],files[1:])
