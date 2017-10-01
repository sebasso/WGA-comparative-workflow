#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
usage: python genome_evolver.py path_to_genome percent identitiy(%) number_of_genomes nucleotide_probabilities(ATCG)
example: python genome_evolver.py /Users/sebastiansoberg/Downloads/Example1/EEE_NJ-60.fasta 2 200 0.25,0.25,0.25,0.25
Percent identitiy is only a outer bound, not an absolute metric. This means if percent=20 the distribution range of probabilities is: (100 - percent) < x < 100

When assigining different nucleotide probabilities np.random.choide is used, this is much slower than the built in random.choice, so use it with caution.
	* numpy.random.choice is used when assigining selective probabilities to nucleotides else pythons random.choise is used.

profiling:
python -m cProfile -s tottime genome_evolver.py /Users/sebastiansoberg/Downloads/Example1/EEE_NJ-60.fasta 20 100
	sorted by tottime -> for the total time spent in the given function (and excluding time made in calls to sub-functions)
"""
from math import log
from math import pow
import sys
import os
import argparse
import Tree as t
import subprocess

"""
3. output snps per leave as done in ksnp and parsnp, check how they use position
	- use gaps or not?

future:
	nb: all of this requires running ksnp_chooser for ksnp and modifying and printing out the core-genome from parsnp.1
	1. ksnp testing with snp_density, have extra constraint on snp density > (kmersize/2)-1
	2. parsnp -> find core-genome and mutate only that section of the genome.
"""

def evolve_genome(args):
	refgenome = args.refgenome
	filename = os.path.basename(args.refgenome)
	percent = args.percent
	num_genomes = args.num_genomes
	base_probabilities = [float(x) for x in args.nucleotide_probabilities.split(",")]
	assert len(base_probabilities) == 4
	assert sum(base_probabilities) == 1
	outputdir = args.outputdir

	genome = ""
	with open(refgenome,"r") as f:
		meta = f.readline()
		genome = f.read()
	num_newlines = -1

	for char in genome:
		if char == '\n':
			num_newlines += 1

	print "newlines: ",num_newlines
	# sanitize input:
	startdiff = 0
	maxdiff = int(float(len(genome)-num_newlines) * (percent/100.0)) # in characters
	generations = int(log(float(num_genomes))/log(2.0))
	assert generations > 0

	if pow(2,generations) < num_genomes: # 2^generations must be bigger than num_genomes, if not add one generation
		generations += 1

	wait_generations = 0
	if maxdiff < generations:
		wait_generations = maxdiff - generations # pass this on so the mutation doesnt start until it makes the most divergent.
		if wait_generations < 0:
			wait_generations = 0
		#dont mutate until maxdiff == generations_left
		#each node has a generation counter

	"""
	if maxdiff > generations*3,4,5,6,7 etc...
	need to adjust have som with more than "minimum change"
	maxdiff = max snp
	max snps per generation = int(maxdiff/generations)
	however there is a built in check on max snps.

	must have snp strategy for each node if max snps per generation > 1 || or all the time?
	1. low
	2. medium
	3. high(aggressive)
	"""
	print "refgenmelength: {} generations: {} max_snp_diff: {} wait_generations: {}".format(len(genome)-num_newlines, generations, maxdiff, wait_generations)
	#check that max % change is in bounds:
	#identitiy = (100 - (float(maxdiff)/float(((len(genome)-num_newlines) * 100))))

	identitiy = 100 - (float(maxdiff)/float((len(genome)-num_newlines)) *100)
	#print "id2: ", id2
	#-(maxdiff/refgenmelength)*100
	#100-(82074/1641482)*100
	print "identitiy: ", identitiy
	print "percent: ", percent
	assert identitiy < 100.0
	assert identitiy >=	 (100.0 - percent)

	maxdiff_generation = int(maxdiff/generations)
	print "maxdiff_generation: ", maxdiff_generation
	assert maxdiff_generation >= 1
	tree = t.Tree(startdiff, maxdiff, generations, wait_generations, genome, maxdiff_generation, base_probabilities)
	leaves_list = tree.get_leaves()[0:num_genomes]

	print "done making structure"
	print "num leaves: ", len(leaves_list)

	if not os.path.exists(outputdir):
		os.makedirs(outputdir)
	outputdirgenomes = outputdir+os.sep+"genomes"

	if not os.path.exists(outputdirgenomes):
		os.makedirs(outputdirgenomes)
	statsfolders = outputdir+os.sep+"snp_stats"

	if not os.path.exists(statsfolders):
		os.makedirs(statsfolders)

	mainsnps_grouped, total_number_of_snps  = tree.create_fasttree_file(filename, outputdir, num_genomes)

	identities = []
	max_snp_diff = 0
	low = sys.maxint
	high = 0

	snp_positions = []

	outputfilename = outputdirgenomes + os.sep + filename

	snp_formatted_file = "reference\n"
	snp_formatted_file += str(len(mainsnps_grouped))+"\n"

	snp_formatted_file_alt = "reference\n"
	snp_formatted_file_alt += str(len(mainsnps_grouped))+"\n"

	sep = "\t"
	statslist = []
	side_bases = 6
	statslist_alt = []

	for i in xrange(0, num_genomes):
		statsMummerInput = ""
		leave_genome = leaves_list[i].get_mutated_genome(genome)

		with open(outputfilename+str(i),"w") as f:
			f.write(meta)
			f.write(leave_genome)

		stats = leaves_list[i].get_snp_stats()

		with open(statsfolders+os.sep+filename+str(i),"w") as f:
			f.write(str(stats))

		"""
		Map snp_stat position back on to genomes with mummer
		then extract it && write to file
		"""
		genome_max_lengt = len(leave_genome)
		for position, snp in stats.iteritems():
			print "snp position: ", position, " leftposition: ", position-1-side_bases, " rightposition: ", position+1+side_bases
			print "actual snp: ", snp, " == ",leave_genome[position], " \tleftmostbase: ",leave_genome[position-1-side_bases], " rightmostbase: ",leave_genome[position+1+side_bases]
			print "Bases-> " , leave_genome[position-1-side_bases: position+1+side_bases], " num bases: ", len(leave_genome[position-1-side_bases: position+1+side_bases])
			print ""

			#TODO: PUT assertion on edge positions here

			"""
				if position - side_bases -1 < 0 
				
				if position + side_bases + 1 > genome_max_length  
			"""

			left = leave_genome[ ( (position) - side_bases) : position ]
			right = leave_genome[ position + 1 : ( (position) + side_bases + 1) ]

			if "\n" in left:
				left = ''.join(leave_genome[((position) - side_bases - 1) : position].splitlines())
			if "\n" in right:
				right = ''.join(leave_genome[ position + 1 : ((position) + side_bases + 2 ) ].splitlines())
			print len(left), len(right)
			print left, "\t", right
			#assert len(left) == side_bases
			#assert len(right) == side_bases

			surronding_bases = left + '.' + right

			if position < side_bases+1 or position+side_bases == len(leave_genome):
				print "snp pos: ",position
			else:
				assert len(leave_genome[((position)-side_bases):position]) == len(leave_genome[position:((position)+side_bases)])

			if len(surronding_bases) < 12:
				print "surrounding ONE\n"
				print surronding_bases
				print "len: ",len(leave_genome)
				print "snp pos: ", position, " diff ", len(leave_genome)-position
			if "\n" in surronding_bases:
				print "surrounding two\n"
				print surronding_bases
				print "len: ", len(leave_genome)
				print "snp pos: ", position, " diff ", len(leave_genome) - position
				exit(1)

			statsMummerInput += ">" + surronding_bases + "_" + snp + "\n" + surronding_bases.replace(".",snp) + "\n"
		print "statsMummerInput: "
		print statsMummerInput
		print "\n\n"
		print outputfilename + str(i)

		#with open(outputfilename + str(i),"r") as f:
		#	print f.read()
		#print "finished reading: ", outputfilename + str(i)
		outputMummerInput = outputdir+"/statsMummerinput.fasta"

		with open(outputMummerInput,"w") as f:
			f.write(statsMummerInput)

		#print "mummer args: ", args.mummer
		mummer_args = " -maxmatch -l 13 -b -c"
		newmum = args.mummer + mummer_args + " "+ outputMummerInput + " " + outputfilename+str(i)
		print newmum
		for position, snp in stats.iteritems():
			print position, snp

		p = subprocess.Popen([newmum], stdout=subprocess.PIPE,
							 stderr=subprocess.PIPE, shell=True)
		stdout, stderr = p.communicate()
		sys.stderr.write("mummer genome evolver stderr:\n" + stderr)
		sys.stdout.write("mummer genome evolver:\n" + stdout + "\n\n\n")
		# now stdout should be filled with mummer.out
		# it is

		c = 0
		for line in stdout.splitlines(): #TODO: can REVERSE strands arise here?
			if c == 2:
				print "reversestrand: ", line
				exit(0)
			if line[0] != ">":
				snp = line[14] # snp
				formatted_line = line.split()
				position = formatted_line[2]
				statslist_alt.append([filename + str(i), str(int(position)), snp+"ALT", "\n"])
			else:
				c += 1

		#TODO: doublecheck this position, either plus or minus 1 or nothing
		#ownversion
		for position, snp in stats.iteritems():
			statslist.append([filename+str(i), str(position + 1), snp, "\n"])# pos+1 as it starts at 0 here but not in ksnp/parsnp


	for snp_entry in sorted(statslist_alt, key=lambda entry: int(entry[1])):
		snp_formatted_file_alt += sep.join(snp_entry)

	for snp_entry in sorted(statslist, key=lambda entry: int(entry[1])):
		snp_formatted_file += sep.join(snp_entry)

	snp_formatted_file_both = ""
	print statslist[0]
	print statslist_alt[0]
	print statslist[0]+statslist_alt[0]
	print "nonextended:\t ",len(statslist_alt)
	statslist_alt.extend(statslist)
	print "extended:\t ",len(statslist_alt)
	print "\n\n\n\n"
	print statslist[0]
	print statslist_alt[0]
	print "\n\n\n\n"
	print "\n\n\n\n"

	for snp_entry in sorted(statslist_alt, key=lambda entry: int(entry[1])):
		snp_formatted_file_both += sep.join(snp_entry)

	with open(statsfolders+os.sep+"reference_formatted_snps.tsv","w") as f:#TODO: create switch here for when simulation and not
		f.write(snp_formatted_file)


	for leave in leaves_list:
		identitiy, num_snps, common_bases = leave.calc_identitiy(genome)

		if identitiy < low:
			low = identitiy
		if identitiy > high:
			high = identitiy
		identities.append(identitiy)

		snp_positions.append(" ".join(str(v) for v in sorted(leave.snp_positions.keys())))

		if leave.snp_diff > max_snp_diff:
			max_snp_diff = leave.snp_diff


	print "Stats:"
	print "Total number of snps: ", total_number_of_snps, "Number of snps unique snps: ", len(mainsnps_grouped)
	print "Lowest identitiy: ",low, " Highest identitiy: ",high
	print "Max snps per leave: ", max_snp_diff, " of: ", maxdiff
	print "Unique probability distributions based on % identitity(bases) between reference and leaves: ",len(set(identities)), " of distributions: ", len(identities)
	print "There are: ", abs(len(set(identities)) - len(set(snp_positions))) , " distributions with equal % identitiy(num bases) but unequal positions on these bases"
	print "Unique distributions of leaves based on snp unique positions: ", len(set(snp_positions)), " of ", len(identities), " distributions. I.e: ", ((float(len(set(snp_positions))) / len(identities)) * 100.0), " % unique genome leaves"
	exit(0)



if __name__ == '__main__':
	#python genome_evolver.py /Users/sebastiansoberg/Downloads/Example1/EEE_NJ-60.fasta 5 200
	"""
	basic:
	python genome_evolver.py --refgenome /Users/sebastiansoberg/Downloads/Example1/EEE_NJ-60.fasta --percent 2 --num_genomes 200
	with probs:
	python genome_evolver.py --refgenome /Users/sebastiansoberg/Downloads/Example1/EEE_NJ-60.fasta --percent 2 --num_genomes 200 --nucleotide_probabilities 0.25,0.25,0.4,0.1
	"""
	parser = argparse.ArgumentParser()
	parser.add_argument("--refgenome", required=True)
	parser.add_argument("--percent", type=float, required=True)
	parser.add_argument("--num_genomes", type=int, required=True)
	parser.add_argument("--nucleotide_probabilities", required=False, default="0.25,0.25,0.25,0.25")
	parser.add_argument("--outputdir", required=False, default=os.getcwd()+"/output")
	parser.add_argument("--mummer", required=False)
	args = parser.parse_args()

	evolve_genome(args)
