import sys
import os
import argparse


def GetArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-C","--Complex",help="These are complex STRs")
	parser.add_argument("-T","--Table",help="The STR Table for which you are going to collate outputs",required=True)
	parser.add_argument("-D","--Directory",help="The directory which has the output data from the EHdn analysis",required=True)
	parser.add_argument("-O","--Output",help="Output table for results",required=True)
	args = parser.parse_args()
	return args


# This function parses a STRetch tsv output
# And with a given locus it will search for that locus 
# And extract the rank and genotyped size
def ParseSTRetch(stretchInfilename,locus):
	return


# This function parses an EHDN tsv output from the outlier script
# And with a given locus it will look for that locus and give the rank
# and Z-score
# NOTE: Assumes Z-score sorted file with header line
def ParseEHDNLocus(ehdnInfilename,locus_chrom,locus_start,locus_end):
	try:
		infile = open(ehdnInfilename,'r')
	except(FileNotFoundError):
		print("This file is missing, returning NA,NA")
		rank=0
		zscore=0
		return rank,zscore
	Hits=infile.readlines()
	for counter,line in enumerate(Hits,0):
		if line[0]=='#':
			continue
		#print(counter)
		#print(line)
		cols = line.strip('\n').split('\t')
		chrom = cols[0]
		start = int(cols[1])
		end = int(cols[2])
		zscore = float(cols[4])
		motif = cols[3]
		# see if the locus is within the window which EHdn called
		if chrom == locus_chrom:
			if (locus_start >= start) and (locus_end <= end):
				# because I skip the headerline, counter == rank
				rank = counter
				return rank,zscore
	# If you complete the file and don't find the hit, then your z-score and rank are 0
	rank = 0
	zscore = 0
	return rank,zscore

# Runs reverse compliment
def reverse_compliment(seq):
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	letters = list(seq) 
	letters = [basecomplement[base] for base in letters] 
	comp = ''.join(letters)
	revcomp = comp[::-1]
	return revcomp

def rotate(input,d): 
    # slice string in two parts for left and right 
	Lfirst = input[0 : d] 
	Lsecond = input[d :] 

# This function take in a motif, and can produce all possible motifs which mirror it
# e.g. CGG --> CCG --> GCG --> CGC
def MotifToPossibleMotifs(motif):
	Motifs = [motif]
	reverseMotif = reverse_compliment(motif)
	Motifs.append(reverseMotif)
	for i in range(len(motif)+1):
		newMotif = motif[len(motif)-i : len(motif)] + motif[0:len(motif) - i]
		Motifs.append(newMotif)
	for i in range(len(reverseMotif)+1):
		newMotif = reverseMotif[len(reverseMotif)-i : len(reverseMotif)] + reverseMotif[0:len(motif) - i]
		Motifs.append(newMotif)
	finalMotifs = list(set(Motifs))
	sortedMotifs = sorted(finalMotifs)
	#print(sortedMotifs)
	return sortedMotifs


# This function parses an EHDN Motif output file
# And produces the rank and z-score of the matching
# motif to the locus query.
def ParseEHDNMotif(ehdnInfilename,locus_motif):
	try:
		infile = open(ehdnInfilename,'r')
	except(FileNotFoundError):
		print("This file is missing, returning 0,0")
		rank=0
		zscore=0
		return rank,zscore
	Hits=infile.readlines()
	#print(locus_motif)
	for counter,line in enumerate(Hits):
		cols = line.strip('\n').split('\t')
		if cols[0]=='motif':
			continue
		zscore = float(cols[1])
		motif = cols[0]

		# Expand the motif to all possible similar motifs
		AllMotifs = MotifToPossibleMotifs(motif)
		# see if the motif is in the possible motifs
		if locus_motif in AllMotifs:
			rank = counter
			return rank,zscore
	# If you complete the file and don't find the hit, then your z-score and rank are 0
	rank = 0
	zscore = 0
	return rank,zscore

def ParseTableCollateOutput(strTableFilename,ehdnOutputDir,outFilename):
	outfile = open(outFilename,'w')
	# write header
	outfile.write("#Gene\tchrom\tstart\tend\tPathogenic_Lower_Bound\tMotif\tSize\tMotif_Z-score\tMotif_Rank\tLocus_Z-score\tLocus_Rank\n")
	strfile = open(strTableFilename,'r')
	
	# ignore header
	header=strfile.readline()
	# Loop through the STR table which we simulated earlier
	# For each STR, create a bed file, and compare it against the locus file
	# as well as the motif ranked file
	for line in strfile:
		cols=line.strip('\n').split('\t')
		gene = cols[0]
		chrom = cols[2]
		start = int(cols[3])
		end = int(cols[4])
		motif = cols[5]
		pathogenicLowerBound = int(cols[6])
		for size in [0, 10, 20, 50, 100, 1000]:
			Size = pathogenicLowerBound + size
		
			# Parse the locus, retrieve rank and z-score
			LocusFilename='%s/%s_%d_Het_EHDN.str_profile_OutlierAnalysis_Locus.txt.sorted.tsv'%(ehdnOutputDir,gene,Size)
			Locusrank,Locuszscore = ParseEHDNLocus(LocusFilename,chrom,start,end)
			
			# Parse the motif, retrieve rank and z-score
			MotifFilename='%s/%s_%d_Het_EHDN.str_profile_OutlierAnalysis_Motif.txt.sorted.tsv'%(ehdnOutputDir,gene,Size)
			Motifrank,Motifzscore = ParseEHDNMotif(MotifFilename,motif)

			outfile.write("%s\t%s\t%d\t%d\t%d\t%s\t%d\t%.3f\t%d\t%.3f\t%d\n"%(gene,chrom,start,end,pathogenicLowerBound,motif,Size,Motifzscore,Motifrank,Locuszscore,Locusrank))


def Main():
	ARGS = GetArgs()
	ParseTableCollateOutput(ARGS.Table,ARGS.Directory,ARGS.Output)
	
	# test motif expansion
	#print(MotifToPossibleMotifs('GGCCTG'))
	#print(MotifToPossibleMotifs('CAG'))
	#print(MotifToPossibleMotifs('CGG'))

if __name__=="__main__":
    Main()







