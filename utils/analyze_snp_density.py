import sys

def find_snp_density(k, snp_file):
	
	snps = ""
	with open(snp_file,'r') as f:
		snps=eval(f.read())
	sorted_snps = sorted(snps.items())

	k = k
	closerthan_k_dividedby2 = []
	closerthan_k = []

	for i in xrange(0,len(sorted_snps)-1):
		temp = sorted_snps[i+1][0] - sorted_snps[i][0]
		if(temp < k/2):
			closerthan_k_dividedby2.append(tuple(sorted_snps[i+1]+sorted_snps[i]))
		if(temp < k):
			closerthan_k.append(tuple(sorted_snps[i+1]+sorted_snps[i]))

	print "Density < k/2: ", len(closerthan_k_dividedby2), " Density < k: ",len(closerthan_k)


if __name__ == '__main__':
	print len(sys.argv)
	if len(sys.argv) < 3:
		print "python analyze_snp_density.py 4 somesnpfiledict.fasta1"
		exit(1)
	else:
		find_snp_density(int(sys.argv[1]), sys.argv[2])

