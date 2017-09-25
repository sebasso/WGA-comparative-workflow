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
"""

def readfiles(files):#files[0]=ksnp, files[1]=parsnp
    print files
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


def remove_most_frequent_base_per_position_ksnp_if_all_genomes_have_snp(num_genomes, SNPs_loci): #ksnp only
    # ONLY remove most frequent base IFF all genomes are listed
    for pos,value in SNPs_loci.iteritems():
        for toolname,snpdict in value.iteritems(): #value= dict where key= toolname and value = list of snips
            if toolname == "kSNP" and len(snpdict) == num_genomes:
                 # snpdict.iteritems() key=filname value= snps
                snp_mapping = snpdict.items()
                snp_counts = Counter(map(lambda x: x[1], snp_mapping)) # count snps
                most_frequent_snp = snp_counts.most_common(1)[0][0] # fetch most frequent snp

                for filename, snp in snp_mapping: #remove most frequent in genomes
                    if snp == most_frequent_snp:
                        del snpdict[filename]

    return SNPs_loci


def find_SNPs_in_same_position(num_genomes, files, tool_names):
    SNPs_same_position = {}

    for i, snp_file in enumerate(files):
        for j in xrange(0,len(snp_file)):
            l1 = snp_file[j].rstrip().split("\t")
            if l1[1] != "x": # unknown position snp from ksnp, still relevant when constructing tree from alignment there it is here
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

    # dict{ position {  tool_name{  filename: SNP  }}}
    if "kSNP" in tool_names:
        SNPs_same_position = remove_most_frequent_base_per_position_ksnp_if_all_genomes_have_snp(num_genomes, SNPs_same_position)

    return SNPs_same_position


def snps_per_genome_total(SNPs_loci, tool_names):
    SNPs_per_genome = Counter()
    for pos,tooldict in SNPs_loci.iteritems():
        for toolname, snpdict in tooldict.iteritems():#snpdict key=filename, value = snp
            for filename in snpdict.iterkeys():
                SNPs_per_genome[filename] += 1

    return SNPs_per_genome


def snps_per_genome_per_tool(SNPs_loci, tool_names):
    Snps_per_tool = {}
    for tool in tool_names:
        SNPs_per_genome = Counter()
        for pos,tooldict in SNPs_loci.iteritems():
            for toolname, snpdict in tooldict.iteritems():#snpdict key=filename, value = snp
                if toolname == tool: #TODO: Change this to be valid for all
                    for filename in snpdict.iterkeys():
                        SNPs_per_genome[filename] += 1

        Snps_per_tool[tool] = SNPs_per_genome
    return SNPs_per_genome


def snp_postional_count_stats(SNPs_same_position, tool_names):

    tool_snp_count = {}
    num_snp_positions = len(SNPs_same_position)
    tool_snp_count["total_snp_positons"] = num_snp_positions

    for tools in tool_names:
        if not tools in tool_snp_count:
            tool_snp_count[tools] = 0

    common_positions = "common_snp_positions"
    tool_snp_count[common_positions] = 0 #finds common snp positions
    for key,value in SNPs_same_position.iteritems():
        keys = value.keys()
        if len(keys) == 2:
            tool_snp_count[common_positions] += 1
        for tool in keys:
            tool_snp_count[tool] += 1

    tool_snp_percent_of_common_snps = []
    for tool in tool_names:
        tool_snp_percent = float(tool_snp_count[common_positions])/tool_snp_count[tool]*100
        tool_snp_percent_of_common_snps.append(tool+" "+str(tool_snp_percent))

    tool_snp_count["common_snp_percentage_of_common_snps"] = tool_snp_percent_of_common_snps

    return tool_snp_count


def snp_positional_plus_allgenomes_count_stats(SNPs_same_position):
    common_snps_results = {}

    common_snps = 0
    for pos,tooldict in SNPs_same_position.iteritems():
        tmp = []
        toolsnpdict = tooldict.items()
        if len(toolsnpdict) == 2:
            #tool1 = set(toolsnpdict[0][1].items())
            tool2 = toolsnpdict[1][1].items()
            for entry in toolsnpdict[0][1].iteritems():
                if entry in tool2:
                    common_snps += 1
                    tmp.append(dict(entry))
        if len(tmp) > 0:
            common_snps_results[pos] = tmp

    return common_snps, common_snps_results


def find_total_snps(SNPs_loci, tool_names):

    tool_snp_count = {}
    for tools in tool_names:
        if not tools in tool_snp_count:
            tool_snp_count[tools] = 0

    for pos,tool_snp_dict in SNPs_loci.iteritems():
        for toolname,snpdict in  tool_snp_dict.iteritems():
            tool_snp_count[toolname] += len(snpdict)

    return tool_snp_count


def compare_snps(num_genomes, outputdir, SNP_files):

    snp_lists_files, tool_names = readfiles(SNP_files)

    stats = {}
    stats["tool_names"] = tool_names

    SNPs_loci = find_SNPs_in_same_position(num_genomes, snp_lists_files, tool_names)
    stats["SNP_positions"] = SNPs_loci

    stats["total_SNPs_per_genome"] = snps_per_genome_total(SNPs_loci, tool_names)

    stats["SNPs_per_genome_per_tool"] = snps_per_genome_per_tool(SNPs_loci, tool_names)#TODO what is this?:

    stats["snp_position_counts"] = snp_postional_count_stats(SNPs_loci, tool_names)

    total_snps = find_total_snps(SNPs_loci,tool_names)
    total_snps["common_snps_all"], common_snps_overview = snp_positional_plus_allgenomes_count_stats(SNPs_loci)

    stats["total_snps"] = total_snps
    
    stats["_common_snps_overview_"] = common_snps_overview

    r = json.dumps(stats, indent=4, encoding="utf-8", sort_keys=True)

    with open(outputdir+"/snps_stats.json","w") as f: #warning will write this relative to exection path - sys.executables
        f.write(r)

    """
    TODO: count how many snp groups in ksnp?
        new metric: SNPs_per_genome / num_snp_groups
    """

if __name__ == '__main__':
    print "Snp comparator @args:"
    print str(sys.argv)+"\n"

    if len(sys.argv) < 4:
        sys.stderr.write("snp_comparator requires a minimum of 2 formatted SNP lists with their respective position and file name")
        exit(1)
    num_genomes = sys.argv[1]
    outputdir = sys.argv[2]
    files = sys.argv[3:]
    compare_snps(int(num_genomes), outputdir, files)
