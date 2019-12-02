# Phillip Richmond
# 2019-07-23
# Updated 2019-10-28 for simulating N and complex STRs

# The purpose of this script is to generate expanded fasta files from STR loci
# Steps:
# 1. 


# Imports
import argparse
import sys
import os
import pybedtools
import random


def GetArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-T","--Table",help="Input Table of STRs, taken from STR_Analysis Github", required=True)
	parser.add_argument("-F","--Fasta",help="Input Fasta file to extract regions from", required=True)
	parser.add_argument("-D","--Directory",help="Output directory to store simulated fasta files", required=True)
	parser.add_argument("-W","--Window",help="Window size for simulating around the candidate STR", required=True)
	parser.add_argument("-S","--Sizes",help="Sizes for simulation, in addition to the lower pathogenic range: e.g. 0,10,100,1000",default="0,10,100,1000")
	parser.add_argument("-C","--Complex",help="These are complex STRs")
	args = parser.parse_args()
	return args

# This function will get the upstream sequence, not including the first repeat nucleotide
# It will create the extra sequence, given a motif and size
# and it will fetch the downstream sequence, defined as the reference genomic repeat itself + window downstream
# NOTE: This will simulate as size + refLen
# And then it will paste the three together
# Fetching sequence with this
#seq = pybedtools.BedTool.seq((chromosome,start,stop),fasta)
def BuildSimulatedSequence(size,window,chrom,start,end,motif,fasta):
	upstreamSeq = pybedtools.BedTool.seq((chrom,int(start)-window,int(start)),fasta)
	if 'N' in motif:
		NewBase = random.choice(['A','C','T','G'])
		newmotif = motif.replace('N',NewBase)
		extraSeq = newmotif * (size)
	else:
		extraSeq = motif * (size)
	
	downstreamSeq = pybedtools.BedTool.seq((chrom,int(start),int(end) + window),fasta)
	returnSeq = upstreamSeq + extraSeq + downstreamSeq
	return returnSeq


# This function will parse an input table and make fasta files for each of the simulated sizes for
# each of the repeats in the table.
# It will also make a fasta for the reference without any extra size
# Table file looks like this:
# #Gene PMID;Source chrom   start(0)	end(0)  Repeat_Motif	Pathogenic_Lower_Bound
def Simulate(tableinfile,outputdirectory,sizes,window,fasta):
	infile = open(tableinfile,'r')
	header = infile.readline() 
	# parse the infile, fetch the sequence, print the fasta with the information
	for line in infile:
		cols = line.strip('\n').split('\t')
		gene = cols[0]
		chrom = cols[2]
		repeat_start = int(cols[3])
		repeat_end = int(cols[4])
		motif = cols[5]
		pathogenic_lower_bound = int(cols[6])
		
		full_start = repeat_start - window
		full_end = repeat_end + window

		# Define outfile identifier
		# Output files will have informative names, split by _
		# gene_motif_repeatlocus_fullwindowlocus_SimSize
		# Where repeatlocus and fullwindowlocus are chrom-start-end
		repeatlocus = '%s-%d-%d'%(chrom,repeat_start,repeat_end)
		fullwindowlocus = '%s-%d-%d'%(chrom,full_start,full_end)

		# iterate over window sizes within the array SizeRanges
		for size in sizes:
			SizeToSim = pathogenic_lower_bound + int(size)
			# Here I wrap Egor's tool to build a JSON then have it produce a fasta file



			seq = BuildSimulatedSequence(SizeToSim,window,chrom,repeat_start,repeat_end,motif,fasta)

			outfile_identifier = '%s_%s_%s_%s_%d.fasta'%(gene,motif,repeatlocus,fullwindowlocus,SizeToSim)
			outfile = open("%s/%s"%(outputdirectory,outfile_identifier),'w')
			outfile.write(">%s\n"%outfile_identifier)
			outfile.write("%s\n"%seq)

		# additionally simulate a reference only version
		
		SizeToSim = 0
		seq = BuildSimulatedSequence(SizeToSim,window,chrom,repeat_start,repeat_end,motif,fasta)
		outfile_identifier = '%s_%s_%s_%s_%d.fasta'%(gene,motif,repeatlocus,fullwindowlocus,SizeToSim)
		outfile = open("%s/%s"%(outputdirectory,outfile_identifier),'w')
		outfile.write(">%s\n"%outfile_identifier)
		outfile.write("%s\n"%seq)


def Main():
	ARGS = GetArgs()
	SizeRanges = ARGS.Sizes.split(',')
	print("You are going to simulate with these sizes:")
	for each in SizeRanges:
		print(each)
	WindowSize = int(ARGS.Window)
	
	# test fasta file is valid
	#seq = pybedtools.BedTool.seq(('2',191745598-WindowSize,191745646+WindowSize),ARGS.Fasta)
	#print(seq)

	# Test BuildSimulatedSequence(size,window,chrom,start,end,motif,fasta)
	#testSeq = BuildSimulatedSequence(10,10,'2','191745598','191745646','GCA',ARGS.Fasta)
	#print(testSeq)

	Simulate(ARGS.Table,ARGS.Directory,SizeRanges,WindowSize,ARGS.Fasta)



if __name__=="__main__":
	Main()






