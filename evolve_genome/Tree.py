#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import numpy as np

class Tree(object):
	"""docstring for Tree"""
	def __init__(self, *args):
		self.start_snp_diff = args[0]
		self.max_snp_diff = args[1]
		self.max_generations = args[2]
		self.wait_generations = args[3]
		self.genome = args[4]

		self.generation = 0
		self.snp_positions = {}
		self.leaves = []
		self.root = Node(self.start_snp_diff, self.max_snp_diff, self, self.generation, self.max_generations, self.wait_generations, self.genome, self.snp_positions, self.leaves)

	def traverse_tree(self):
		print "todo"

	def get_leaves(self):
		return self.leaves


class Node(object):
	"""docstring for Node"""
	def __init__(self, *args):
		#print "generation: ", args[3]+1
		self.snp_diff = args[0]
		self.max_snp_diff = args[1]
		self.parent = args[2]#node
		self.generation = args[3] + 1
		self.max_generations = args[4]
		self.wait_generations = args[5]
		#self.genome = args[6] #string
		#self.snp_positions = args[7]  #self.positions should only be kept in leaves, same with genomes.

		if self.generation <= self.max_generations: # as long as more generations are needed to achienge number of genomes as desired	
			#print "selfgeneration: ", self.generation, " max generation: ", self.max_generations, " waitgen: ",self.wait_generations
			#print self.wait_generations == 0
			#print self.snp_diff < self.max_snp_diff
			#print self.snp_diff, self.max_snp_diff
			#this enables mutation to start WHEN there is only mutation per generation if percent identitiy is very high and sequence very short.
			if self.wait_generations == 0 and self.snp_diff <= self.max_snp_diff: # if maxdiff is achieved, simply do nothing
				self.mutations, self.snp_diff = self.mutate_genome(args[6], self.snp_diff, self.max_snp_diff, args[7]) 
				args[7].update(self.mutations)
			else:
				if self.snp_diff == 0: # make sure it hasnt started, because doesnt want it to start counting if snp diff is reached
					self.wait_generations -= 1 # counts down until it shall start producing snps

			# always two childern because thats how it happends in nature
			self.left = Node(self.snp_diff, self.max_snp_diff, self, self.generation, self.max_generations, self.wait_generations, args[6], args[7], args[8])	
			self.right = Node(self.snp_diff, self.max_snp_diff, self, self.generation, self.max_generations, self.wait_generations, args[6], args[7], args[8])

		else:
			#doesnt need more generations
			args[8].append(self)
			self.genome = args[6]
			self.snp_positions = args[7]
			self.left = None
			self.rigth = None
			#should add all leaves to a list for easy retrieval

	"""
	1. decide how many snps
	2. then decide locus
	"""
	def mutate_genome(self, genome, snp_diff, maxdiff, snp_positions): 
		#positions are those positions already in genome, they are changed and will not be changed again.  
		nucleotides = ['A','T','C','G']
		mutations = {} # position,SNP
		nucleotide_probabilities = [0.25, 0.25, 0.25, 0.25]
		snp_interval = [1,2,3]
		num_snps = np.random.choice(snp_interval, p=[0.7, 0.21, 0.09])

		for i in range(0,num_snps):
			
			position = random.randrange(0, len(genome)) # pick position -> assumption: all positions in genome are equally prone to a snp (should be changed to differentiate core genome)
			while position in snp_positions:
				position = random.randrange(0, len(genome))

			snp = np.random.choice(nucleotides, p=nucleotide_probabilities)
			
			while snp == genome[position]:
				snp = np.random.choice(nucleotides, p=nucleotide_probabilities)
			
			mutations[position] = snp
			genome_modable = list(genome)
			genome_modable[position] = snp
			genome = "".join(genome_modable)
		
		return (mutations, (snp_diff + num_snps))

	
	def calc_identitiy(self, refgenome):
		c = 0
		for i in xrange(0,len(self.genome)):
			if self.genome[i] == refgenome[i]:
				c += 1

		return (c/float(len(self.genome))*100), c

	
	def create_fasttree_snps(self):
		fasta_format = ""
		gap="-"
		
		sorted(self.snp_positions.iterkeys())

		for key in sorted(self.snp_positions.iterkeys()):
			key, self.snp_positions[key] # sorted



		return fasta_format






