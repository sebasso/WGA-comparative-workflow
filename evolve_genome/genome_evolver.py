#!/usr/bin/env python
# -*- coding: utf-8 -*-

#usage: python genome_evolver.py path_to_genome percent identitiy number_of_genomes
"""
This script evolves a genome until a % identitiy is reached.
% identitiy is defined as number of bases equal in the same position between two genomes.


reason for tree structure is to simulate evolution properly i.e. for each generation each new genome is based upon its parents genome. 
"""
from math import log
from math import pow
import sys
import Tree as t
#from guppy import hpy


def create_


#consider global values 
#NB: Must ha a treshold of length on refgenome towards percent identity
def evolve_genome(refgenome, percent, num_genomes):
	#h = hpy() 
	genome = ""

	with open(refgenome,"r") as f:
		meta = f.readline()
		genome = f.read()

	startdiff = 0
	maxdiff = int(float(len(genome)) * (int(percent)/100.0)) # in characters
	generations = int(log(float(num_genomes))/log(2.0))

	if pow(2,generations) < num_genomes: # 2^generations must be bigger than num_genomes
		generations += 1 

	wait_generations = 0
	if maxdiff < generations:
		wait_generations = maxdiff - generations # pass this on so the mutation doesnt start until it makes the most divergent.
		if wait_generations < 0:
			wait_generations = 0
		#dont mutate until maxdiff == generations_left
		#each node has a generation counter

	print "generations: {} max_snp_diff: {} wait_generations: {}".format(generations, maxdiff, wait_generations)

	tree = t.Tree(startdiff, maxdiff, generations, wait_generations, genome) #*args

	leaves_list = tree.get_leaves()
	print "done making structure"
	print "num leaves: ", len(leaves_list)
	
	identities = []
	for leave in leaves_list:
		identities.append(leave.calc_identitiy(genome)[0])
		#print leave.snp_positions
		print "num snps: ",len(leave.snp_positions)
		for key,val in leave.snp_positions.items():
			print key,val
		exit(0)


	# select correct amount of genomes:
	#for i in xrange(0, num_genomes):
		#print leaves_list[i]	
	#	leaves_list[i] 


#test
#	for identity in identities:
#		print identity, identity < 100.0, identity > 98.0

	exit(0)
	#print h.heap() # profiler

#TODO: should support different probabilities for mutation rates
if __name__ == '__main__':

	if len(sys.argv) > 2:
		evolve_genome(sys.argv[1], sys.argv[2], sys.argv[3])
	else:
		exit(0)



"""
98 % identitiy
2 % unequal
parsnp requires > 97.5 %

"""
