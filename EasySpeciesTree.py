#!/usr/bin/env python
# -*- coding: UTF-8 -*-
# @Author: "Wei Dong"
# @Contanct: "1369852697@qq.com"
# @Time: "2018-06-04"
# @Version = "Version 1.0"
# @Discripton: "This script was designed to easily construct the species tree based on the single-copy gene obtained from the OrthoFinder results."

#import sys
#import os.path as op
import os
from collections import OrderedDict
import argparse
import datetime
import time
import subprocess

############### MODIFY THE FOLLWINGS PATHS FOR ALL DEPENDENT PROGRAMS ###############
MAFFT = '/usr/local/bin/mafft'
RAxML = '/Users/Davey/Desktop/RAxML-master/raxmlHPC-PTHREADS'
ASTRAL = '/Users/Davey/Desktop/ASTRAL/astral.5.6.2.jar'
TRIMAL = '/usr/local/bin/trimal'
#####################################################################################

#SpeciesID = ["Aco","Musa","VIT","Potri","AT"] #offer all abbreviated species id information
SpeciesID = []
SingleOrtho = []
SingleGeneID = []
mydict1 = OrderedDict()
mydict2 = OrderedDict()
mydict3 = OrderedDict()
OrthoDict = OrderedDict()

# get every species abbreviated id information
def get_SpeciesID(speciesID):
	with open(speciesID) as fh:
		for line in fh:
			SpeciesID.append(line.strip())

# execute the shell command
def run_command(cmd):
    # print("INFO: Running command: {0}".format(cmd), flush=True)
    print(cmd)
    return_code = subprocess.call(cmd, shell=True)
    if return_code != 0:
        print("ERROR: [{2}] Return code {0} when running the following command: {1}".format(return_code, cmd, datetime.datetime.now()))

def get_SingleGeneID(SingleOrthoGroup, OrthoGroups):
	'''
	Get the single-copy gene's id from the SingleOrthoGroup
	'''
	with open(SingleOrthoGroup) as fh:
		for line in fh:
			SingleOrtho.append(line.strip())

	with open(OrthoGroups) as fh:
		for line in fh:
			if line.startswith("\s"):
				continue
			else:
				lines = line.strip().split()
				OrthoDict[lines[0]] = "-".join(lines[1:])

	for item in SingleOrtho:
		if item in OrthoDict:
			SingleGeneID.append(item + ":" + OrthoDict[item])
	return SingleGeneID

def get_SingleGeneSeq(SingleGeneIDs, allseqs):
	'''
	Get the single-copy gene's sequence from the all species protein sequences
	''' 
	with open(allseqs) as fh:
		for line in fh:
			if line.startswith(">"):
				seqid = line.strip(">").split()[0]
				mydict1[seqid] = []
			else:
				mydict1[seqid].append(line.strip())

	os.mkdir("SingleGene")
	os.chdir("SingleGene")
	for item in SingleGeneIDs:
		items = item.strip().split(":")
		seqfile = items[0] + ".fas"
		geneids = items[1].split("-")
		with open(seqfile,"w") as fh:
			for item in geneids:
				if item in mydict1:
					for ids in SpeciesID:
						if item.startswith(ids):
							fh.write(">" + ids + "\n" + "\n".join(mydict1[item]) + "\n")
	os.chdir("../")

def MSA_SingleCopyGene(SingleGeneIDs):
	'''
	Multiple sequences alignment for each single-copy gene in every species
	'''
	os.mkdir("SingleGene_MSA")
	os.chdir("SingleGene_MSA")
	for item in SingleGeneIDs:
		items = item.strip().split(":")
		seqfile = items[0] + ".fas"
		seq_aln_file = items[0] + "_aln.fas"
		seq_aln_trimed_file = items[0] + "_aln_trimed.fas"
		#os.system(MAFFT + ' ../SingleGene/' + seqfile + ' >' + seq_aln_file)
		cmd1 = MAFFT + ' ../SingleGene/' + seqfile + ' >' + seq_aln_file
		run_command(cmd1)
		cmd2 = TRIMAL + ' -in ' + seq_aln_file + ' -out ' + seq_aln_trimed_file + ' -automated1'
		run_command(cmd2)

def merge_SingleCopyGene(SingleGeneIDs):
	'''
	Merge every aligned single-copy gene into one super-gene matirx
	'''
	for item in SingleGeneIDs:
		items = item.strip().split(":")
		seq_aln_file = items[0] + "_aln.fas"
		seq_aln_trimed_file = items[0] + "_aln_trimed.fas"
		with open(seq_aln_trimed_file) as fh:
		#with open(seq_aln_file) as fh:
			for line in fh:
				if line.startswith(">"):
					seqid = line.strip(">").split()[0] + "." + items[0]
					mydict2[seqid] = []
				else:
					mydict2[seqid].append(line.strip())

			for ids in SpeciesID:
				mydict3[ids] = []
				for key,value in mydict2.items():
					if key.startswith(ids):
						mydict3[ids].append("".join(value))

	with open("merged_allSingleGenes.fas","w") as fh:
		for key,value in mydict3.items():
			fh.write(">" + key + "\n" + "".join(value) + "\n")

def MLtree_Concatenation():
	'''
	Construct the ML species tree with the concatenate methold
	'''
	os.mkdir("Concatenation")
	os.chdir("Concatenation")
	cmd = RAxML + ' -T ' + thread + ' -f a -N ' + nb + ' -m ' + model + ' -x 123456 -p 123456 -s ../SingleGene_MSA/merged_allSingleGenes.fas -n concatenation_out.nwk'
	run_command(cmd)
	os.chdir("../")

def MLtree_Coalescence(SingleGeneIDs):
	'''
	Construct the ML species tree with the coalescent methold
	'''
	os.mkdir("Coalescence")
	os.chdir("Coalescence")
	for item in SingleGeneIDs:
		items = item.strip().split(":")
		seq_aln_file = items[0] + "_aln.fas"
		seq_aln_trimed_file = items[0] + "_aln_trimed.fas"
		out_tree_file = items[0] + "_out.nwk"
		cmd1 = RAxML + ' -T ' + thread + ' -f a -N ' + nb + ' -m ' + model + ' -x 123456 -p 123456 -s ../SingleGene_MSA/' + seq_aln_trimed_file + ' -n ' + out_tree_file
		run_command(cmd1)
		cmd2 = 'cat RAxML_bipartitions.' + out_tree_file + ' >>allSingleGenes_tree.nwk'
		#os.system('cat ' + out_tree_file + ' >>allSingleGenes_tree.nwk')
		#cmd = 'cat ' + out_tree_file + ' >>allSingleGenes_tree.nwk'
		run_command(cmd2)
		cmd3 = 'echo ./RAxML_bootstrap.' + out_tree_file + ' >>allSingleGenes_bootstrap.txt'
		#os.system('echo ./' + out_tree_file + ' >>allSingleGenes_bootstrap.txt')
		#cmd = 'echo ./' + out_tree_file + ' >>allSingleGenes_bootstrap.txt'
		run_command(cmd3)

	cmd4 = 'java -jar ' + ASTRAL +' -i allSingleGenes_tree.nwk -b allSingleGenes_bootstrap.txt -r ' + nb + ' -o Astral.coalescent_out.result'
	run_command(cmd4)
	os.system('tail -n 1 Astral.coalescent_out.result >Astral.coalescence_tree.nwk')
	os.chdir("../")

def main(args):
	reader1 = args.input1
	reader2 = args.input2
	reader3 = args.input3
	reader4 = args.input4
	global thread,nb,model
	thread = args.thread
	nb = args.bootstrap
	model = args.model
	get_SpeciesID(reader1)
	print("Step1: Get each single-copy gene id and its sequences.\n")
	single_gene_ids = get_SingleGeneID(reader2,reader3)
	get_SingleGeneSeq(single_gene_ids,reader4)

	print("Step2: Conduct multiple sequences alignment for each single-copy gene.\n")
	MSA_SingleCopyGene(single_gene_ids)

	print("Step3: Merge all aligned single-copy gene into one file.\n")
	merge_SingleCopyGene(single_gene_ids)
	os.chdir("../")

	print("Step4: Construct the ML species tree with the concatenate method.\n")
	MLtree_Concatenation()

	print("Step5: Construct the ML species tree with the coalescent method.\n")
	MLtree_Coalescence(single_gene_ids)

	print("Congratulations! All tasks were finished!\n")

if __name__ == "__main__":
	begin = time.clock()
	parser = argparse.ArgumentParser(
		prog="EasySpeciesTree",
		formatter_class=argparse.RawTextHelpFormatter,
		description='''
-------------------------------------------------------------------------------------------------------
%(prog)s <SpeciesID prefix> <SingleCopyOrtho> <Orthogroups> <protein file> [thread] [bootstrap] [model]
Author: Wei Dong <1369852697@qq.com>, FAFU
Version: v1.0
Easily construct the ML species tree with all single-copy gene's protein sequences
-------------------------------------------------------------------------------------------------------
''')
	parser.add_argument(
		'-in1', '--input1', 
		required=True, 
#	 	dest="Prefix",
		help="offer the prefix of all abbreviated species id ")
	parser.add_argument(
		'-in2', '--input2', 
		required=True, 
#	 	dest="SingleCopyOrtho",
		help="offer the Single-copy Orthogroups file, SingleCopyOrthogroups.txt")
	parser.add_argument(
		'-in3', '--input3', 
		required=True, 
#	 	dest="Orthogroups",
		help="offer the all Orthogroups file, Orthogroups.csv")
	parser.add_argument(
		'-in4', '--input4',
		required=True, 
#	 	dest="Proteins",
		help="offer all species protein sequences")
	parser.add_argument(
		'-t', '--thread',
		default="10", 
		help="set the number of thread, default=10")
	parser.add_argument(
		'-nb', '--bootstrap',
		default="100", 
		help="set the number of bootstrap, default=100")
	parser.add_argument(
		'-m', '--model',
		default="PROTGAMMAJTT", 
		help="set the model of amino acid substitution, default=PROTGAMMAJTT")
	args = parser.parse_args()
	main(args)
	end = time.clock()
	print("All tasks used time: %ss" % (end - begin))

