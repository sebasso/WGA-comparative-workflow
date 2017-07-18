import sys
import subprocess
import os
from os import listdir
from os.path import isfile, join

#NOTE; works only on FASTA

def execute_kChooser_on_all(in_file,path,os_path):
	files = []
	with open(in_file,"r") as f:
		files = f.readlines()

	filepaths = []
	for file in files:
		filepaths.append(file.split("\t")[0])

	allthis = 0
	perl_exec = path+'/Kchooser-modded.pl'
	sys.stderr.write("Perl path:\n")
	sys.stderr.write(perl_exec)
	sys.stderr.write("\n")
	sys.stderr.write(os.getcwd())
	sys.stderr.write("\n")
	sys.stderr.write("path:\n"+path)
	sys.stderr.write("\nos_path:\n"+os_path)
	sys.stderr.write("\n")
	for file in filepaths:
		proc = subprocess.Popen(['perl', perl_exec, file, path, os_path], stdout=subprocess.PIPE)
		allthis += int(proc.communicate()[0])

	answer = allthis/len(filepaths)
	if (answer%2)==0:
		print (answer-1)
	else:
		print answer

	#cleanup
	cwd = os.getcwd()
	if os.path.isfile(cwd+"/output0"):
		os.remove(cwd+"/output0")

	if os.path.isfile(cwd+"/shortSeq.fasta"):
		os.remove(cwd+"/shortSeq.fasta")

	if os.path.isfile(cwd+"/medSeq.fasta"):
		os.remove(cwd+"/medSeq.fasta")

	if os.path.isfile(cwd+"/jellyout.txt"):
		os.remove(cwd+"/jellyout.txt")



if __name__ == '__main__':
	if sys.argv[1]:
		execute_kChooser_on_all(sys.argv[1], sys.argv[2], sys.argv[3])
	else:
		exit(0)
