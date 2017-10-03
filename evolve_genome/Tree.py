#!/usr/bin/env python
# -*- coding: utf-8 -*-

import random
import numpy as np

class Tree(object):
	"""docstring for Tree"""
	def __init__(self, *args):
		self.start_snp_diff = args[0]
		self.max_snp_diff = args[1]
		self.max_generations = args[2] + 1 # because node under starts at 1 shoul be zero but but easier to have + 1 here, if not to few genomes are produced
		self.wait_generations = args[3]
		self.genome = args[4]
		self.maxdiff_generation = args[5]
		self.base_probabilities = args[6]
		self.generation = 0
		self.snp_positions = {}
		self.leaves = []
		self.root = Node(self.start_snp_diff, self.max_snp_diff, self, self.generation, self.max_generations, \
			self.wait_generations, self.genome, self.snp_positions, self.leaves, self.maxdiff_generation, self.base_probabilities)

	def get_leaves(self):
		return self.leaves


	"""
	This works for evolvement without recombination
	"""
	def create_fasttree_file(self, filename, outputdir, num_genomes):

		# output UNIQUE snp positions in for self.leaves.snp_positions
		leaves = self.get_leaves()[0:num_genomes]

		mainsnps = {} # main snp
		total_number_of_snps = 0
		for i in xrange(0,num_genomes):
			snp_positions = leaves[i].snp_positions
			print "Num snps for genome: ",i," : ", len(snp_positions)
			total_number_of_snps += len(snp_positions)
			mainsnps.update(snp_positions) # group all positions
		tmp = mainsnps.copy()
		mainsnps = sorted(mainsnps.iteritems()) #list
		print "Total number of snps: ", len(mainsnps)

		open(outputdir + "/fasttreeoutput.fasta", 'w').close()
		max_snp_length = len(mainsnps)

		with open(outputdir + "/fasttreeoutput.fasta", "a") as f:
			for i in xrange(0, num_genomes):
				curr_genome = [None] * max_snp_length
				for k in xrange(0, max_snp_length):
					pos = mainsnps[k][0]
					snp = mainsnps[k][1]
					if pos in leaves[i].snp_positions:  # hashlookup
						curr_genome[k] = snp

				for j in xrange(0, max_snp_length):
					if curr_genome[j] is None:
						curr_genome[j] = "-" #what fasttree requires as gap ref: http://www.microbesonline.org/fasttree/

				f.write(">" + filename + str(i) + "\n" + ''.join(curr_genome) + "\n")

		"""
		This generates the each a whole genome with the snps. It is too memory consuming to be used as for large
		sequences it causes fasttree to use above > 60 GB of ram due the scaling of large sequences.
		Instead a strategy based on snps solely is chosen, this solution only outputs the all snp each
		genome contains and adds gaps between snps. It doesnt care if the two snps are 100 bases apart
		it treats it as one as all these 100 bases are similar it wont effect fasttree and greatly reduce
		memory consumtion and computation time. 
		
		2. It also makes it more realistic when comapring to ksnp and parsnp as it defines
			the output to fasttree in the same manner
		
		common_snps = {} # all snp positions with their counts, remove?
		for pos, snp in mainsnps: # for snp position there is check which leaves it is in, if its not in the leave append base from reference_genome
			for count, leave in enumerate(leaves): #must have enumeration here
				if pos in leave.snp_positions:
					if not pos in common_snps:
						common_snps[pos] = 1
					else:
						common_snps[pos] += 1

		print "Shared snps among genomes: ",len(common_snps.values())

		open(outputdir+"/fasttreeoutput2.fasta", 'w').close()

		genome_output = list(self.genome)
		with open(outputdir+"/fasttreeoutput2.fasta", "a") as f:
			for count, leave in enumerate(leaves): #must have enumeration here
				for pos, snp in leave.snp_positions.iteritems():
					genome_output[pos] = snp
				f.write( ">"+filename+str(count)+"\n"+''.join(genome_output)+"\n")
		"""

		return tmp, total_number_of_snps

	"""
	max snpp diff er lengden per genome
	- ha med snippene og - ellers.	
	"""


class Node(object):
	"""Represents internal and leaf nodes"""
	def __init__(self, *args):
		self.snp_diff = args[0]
		self.max_snp_diff = args[1]
		self.parent = args[2]#node
		self.generation = args[3] + 1
		self.max_generations = args[4]
		self.wait_generations = args[5]
		#self.genome = args[6] #string
		self.snp_positions = args[7].copy() # every mutation must be local

		if self.generation <= self.max_generations: # as long as more generations are needed to achienge number of genomes as desired
			#this enables mutation to start WHEN there is only mutation per generation if percent identitiy is very high and sequence very short.

			if self.wait_generations == 0 and self.snp_diff <= self.max_snp_diff: # if maxdiff is achieved, simply do nothing
				self.mutations, new_snps = self.mutate_genome(args[6], self.max_snp_diff, self.snp_positions, args[9], args[10])
				self.snp_diff += new_snps
				self.snp_positions.update(self.mutations)
			else:
				if self.snp_diff == 0: # make sure it hasnt started, because doesnt want it to start counting if snp diff is reached
					self.wait_generations -= 1 # counts down until it shall start producing snps

			if self.generation == self.max_generations:
				#doesnt need more generations
				args[8].append(self) # appends leave nodes
				#print "\nLeaf:\nsnip_diff: {} total snps: {}\nall snps: {}\nparent generation: {}".format(self.snp_diff, len(self.snp_positions), self.snp_positions, self.parent.generation)
				self.left = None
				self.rigth = None
			else:
				# always two childern because thats how it happends in nature
			#	if self.generation == 6:
					#print "num snps: {} snps before left: {}".format(len(self.snp_positions), self.snp_positions)
				self.left = Node(self.snp_diff, self.max_snp_diff, self, self.generation, self.max_generations, self.wait_generations, args[6], self.snp_positions, args[8], args[9], args[10])
			#	if self.generation == 6:
					#print "num snps: {} snps before right: {}".format(len(self.snp_positions), self.snp_positions)
				self.right = Node(self.snp_diff, self.max_snp_diff, self, self.generation, self.max_generations, self.wait_generations, args[6], self.snp_positions, args[8], args[9], args[10])


	"""
	1. decide how many snps
	2. then decide locus
	"""
	def mutate_genome(self, genome , maxdiff, snp_positions, maxdiff_generation, nucleotide_probabilities):
		#positions are those positions already in genome, they are changed and will not be changed again.
		nucleotides = ['A','T','C','G']
		mutations = {} # position,SNP
		nucleotide_probabilities = nucleotide_probabilities

		snp_interval = [x for x in xrange(1,maxdiff_generation+1)] #maxdiff_generation is needed to reach max dissimilarity
		num_snps = random.choice(snp_interval) #uniform distribution, if not use np.random.choide(seq, list_of_probaility distribution)
		#snp_interval = [1,2,3]
		#num_snps = np.random.choice(snp_interval, p=[0.7, 0.21, 0.09])
		for i in range(0,num_snps):

			position = random.randrange(0, len(genome)) # pick position -> assumption: all positions in genome are equally prone to a snp (should be changed to differentiate core genome)
			while position in snp_positions or genome[position] == "\n":
				position = random.randrange(0, len(genome))

			snp = np.random.choice(nucleotides, p=nucleotide_probabilities)
			#snp = random.choice(nucleotides)
			while snp == genome[position]: #make sure its not the same as already are there
				snp = np.random.choice(nucleotides, p=nucleotide_probabilities)
				#snp = random.choice(nucleotides)

			mutations[position] = snp

		return (mutations, num_snps)


	def calc_identitiy(self, refgenome):
		num_snps = len(self.snp_positions)
		reflen = len(refgenome)
		common_bases = reflen - num_snps
		return ((float((common_bases)) / float(reflen)) * 100.0), num_snps, common_bases


	def get_mutated_genome(self, refgenome):
		refgenome_list = list(refgenome)

		for pos, snp in self.snp_positions.items():
			refgenome_list[pos] = snp

		tmp = "".join(refgenome_list)

		a, b = self.check_identitiy(refgenome,tmp) #TODO: remove when not in testing phase
		ans, bens, common_bases =  self.calc_identitiy(refgenome) #TODO: remove when not in testing phase
		assert a == ans
		assert b == common_bases

		return tmp


	def check_identitiy(self, refgenome, genome):
		c = 0
 		for i in xrange(0,len(genome)):
 			if genome[i] == refgenome[i]:
 				c += 1

		return ( (c / float(len(genome))) * 100), c


	def create_fasttree_snps(self, refgenome):
		fasta_format = ""
		gap="-"
		sorted(self.snp_positions.iterkeys())

		for key in sorted(self.snp_positions.iterkeys()):
			key, self.snp_positions[key] # sorted

		return fasta_format


	def get_snp_stats(self):
		return self.snp_positions
