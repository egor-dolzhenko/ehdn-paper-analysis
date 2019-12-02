import sys
import os
import argparse

def GetArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-C","--Complex",help="These are complex STRs")
	parser.add_argument("-T","--Table",help="The STR Table for which you are going to collate outputs",required=True)
	parser.add_argument("-EM","--EHDN_Motifs",help="The EHDN motif files (z-score sorted)",nargs="*")
	parser.add_argument("-EL","--EHDN_Locuses",help="The EHDN locus files (z-score sorted)",nargs="*")
	parser.add_argument("-S","--STRetches",help="The output STRetch_STRs.tsv files",nargs="*")
	parser.add_argument("-O","--Output",help="Output table for results",required=True)
	parser.add_argument("-L","--Lower",help="Only use lower-bound",action='store_true')
	args = parser.parse_args()
	
	print("You have selected this as your input STR table:\t%s"%args.Table)

	print("You have selected this as your output file:\t%s"%args.Output)

	print("You have selected these as your EHDN_Locus files:")
	if args.EHDN_Locuses is not None:
		for each in args.EHDN_Locuses:
			print(each)
	else:
		EHDN_Locuses = []

	print("You have selected these as your EHDN_Motif files:")
	if args.EHDN_Motifs is not None:
		for each in args.EHDN_Motifs:
			print(each)
	else:
		args.EHDN_Motifs = []

	print("You have selected these as your STRetch files:")
	if args.STRetches is not None:
		for each in args.STRetches:
			print(each)
	else:
		args.STRetches = []

	return args

# This function takes two genomic regions and checks to see if they overlap
# covers scenarios:
# 1)
# A1-----------A2
#           B1--------B2
# 2)
# B1--------B2
#       A1---------A2
# 3)
# A1-----------A2
#      B1--B2
# 4)
# B1--------------B2
#       A1---A2
def Overlap(Astart,Aend,Bstart,Bend):
	if (Aend < Bend) and (Aend < Bstart):
		return False
	if (Astart > Bend) and (Aend > Bstart):
		return False
	array = [Astart,Aend,Bstart,Bend]
	sortedarray=sorted(array)
	Inner=int(sortedarray[2])-int(sortedarray[1])
	if (Inner ):
		return True
	else:
		return False

# This function parses a STRetch tsv output
# And with a given locus it will search for that locus 
# And extract the rank and genotyped size
def ParseSTRetch(stretchInfilename,locus_chrom,locus_start,locus_end,locus_motif):
	try:
		infile = open(stretchInfilename,'r')
	except(FileNotFoundError):
		print("This file is missing:%s\n Returning 0,0"%stretchInfilename)
		rank=0
		pval='.'
		size=0
		return size,rank,pval
	Hits=infile.readlines()
	
	for counter,line in enumerate(Hits,0):
		cols = line.strip('\n').split('\t')
		if cols[0]=='chrom':
			continue
		chrom = cols[0]
		start = int(cols[1])
		end = int(cols[2])
		motif = cols[4]
		pval = cols[8]
		size = float(cols[10])

		# Expand the motif to all possible similar motifs
		AllMotifs = MotifToPossibleMotifs(motif)

		# See if there is overlap between the locus we are querying and the STRetch hit
		if (chrom==locus_chrom) or (chrom==('chr'+locus_chrom)):
			#print(chrom,locus_chrom)
			if ( Overlap(start,end,locus_start,locus_end) ):
		# see if the motif is in the possible motifs
				if locus_motif in AllMotifs:
					rank=counter
					#print(size,rank,pval)
					return size,rank,pval
	rank = 0
	size = 0
	pval = '.'
	return size,rank,pval

# This function parses an EHDN tsv output from the outlier script
# And with a given locus it will look for that locus and give the rank
# and Z-score
# NOTE: Assumes Z-score sorted file with header line
def ParseEHDNLocus(ehdnInfilename,locus_chrom,locus_start,locus_end):
	try:
		infile = open(ehdnInfilename,'r')
	except(FileNotFoundError):
		print("This file is missing:%s\n Returning 0,0"%ehdnInfilename)
		rank=0
		zscore=0
		return rank,zscore
	Hits=infile.readlines()
	strcounter = 0
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
		# keep a counter for STR separate from longer motifs
		if len(motif)<=6:
			strcounter+=1
		# see if the locus is within the window which EHdn called
		if chrom == locus_chrom:
			if (locus_start >= start) and (locus_end <= end):
				# because I skip the headerline, counter == rank
				rank = counter
				strrank=strcounter
				return rank,strrank,zscore
	# If you complete the file and don't find the hit, then your z-score and rank are 0
	rank = -1
	zscore = -1
	strrank = -1
	return rank,strrank,zscore

# Runs reverse compliment
def reverse_compliment(seq):
	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
	letters = list(seq) 
	letters = [basecomplement[base] for base in letters] 
	comp = ''.join(letters)
	revcomp = comp[::-1]
	return revcomp

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
		print("This file is missing:%s\n Returning 0,0"%ehdnInfilename)
		rank=-1
		zscore=-1
		return rank,zscore
	Hits=infile.readlines()
	#print(locus_motif)
	strcounter=0
	for counter,line in enumerate(Hits):
		cols = line.strip('\n').split('\t')
		if cols[0]=='motif':
			continue
		zscore = float(cols[1])
		motif = cols[0]
		# keep a counter for STR separate from longer motifs
		if len(motif)<=6:
			strcounter+=1
	
		# Expand the motif to all possible similar motifs
		AllMotifs = MotifToPossibleMotifs(motif)
		# see if the motif is in the possible motifs
		if locus_motif in AllMotifs:
			rank = counter
			strrank=strcounter
			return rank,strrank,zscore
	# If you complete the file and don't find the hit, then your z-score and rank are 0
	rank = -1
	zscore = -1
	strrank = -1
	return rank,strrank,zscore


# The primary function which parses the STR_table
def ParseTableCollateOutput(strTableFilename,EHdn_Locus_FileList,EHdn_Motif_FileList,STRetch_FileList,outFilename,LowerBoundOnly):
	outfile = open(outFilename,'w')
	# write header
	outfile.write("#Gene\tchrom\tstart\tend\tPathogenic_Lower_Bound\tMotif\n")
	strfile = open(strTableFilename,'r')
	
	# ignore header
	header=strfile.readline()
	# Loop through the STR table which we simulated earlier
	# For each STR, create a bed file, and compare it against the locus file
	# as well as the motif ranked file
	outfile.write("GeneID,PathogenicLowerBound,Motif,SimulatedSize,STR_Tool,Rank,Zscore||Pval\n")
	for line in strfile:
		cols=line.strip('\n').split('\t')
		gene = cols[0]
		chrom = cols[2]
		start = int(cols[3])
		end = int(cols[4])
		motif = cols[5]
		# skip motifs with N
		if 'N' in motif:
			continue
		pathogenicLowerBound = int(cols[6])
		
		# Process EHDN Locus files for matches
		if len(EHdn_Locus_FileList) > 0:
			for ehdnlocusfile in EHdn_Locus_FileList:
				# split filename to get rid of the filepath, and then split the name into columns
				filenamecols = ehdnlocusfile.split('/')[-1].split('_')
				file_gene = filenamecols[0]
				# ignore the complex STRs
				if '-' in filenamecols[1]:
					continue
				file_strsize = int(filenamecols[1])
				# Check to make sure we have the right gene
				if file_gene != gene:
					continue
				if LowerBoundOnly:
					if file_strsize != pathogenicLowerBound:
						continue
				EHDNLocusRank,EHDNLocusSTRRank,EHDNLocuszscore = ParseEHDNLocus(ehdnlocusfile,chrom,start,end)
				outfile.write("%s,%d,%s,%d,EHDN_Locus,%d,%d,%.3f\n"%(gene,pathogenicLowerBound,motif,file_strsize,EHDNLocusRank,EHDNLocusSTRRank,EHDNLocuszscore))

		if len(EHdn_Motif_FileList) > 0:
			for ehdnmotiffile in EHdn_Motif_FileList:
				# split filename to get rid of the filepath, and then split the name into columns
				filenamecols = ehdnmotiffile.split('/')[-1].split('_')
				file_gene = filenamecols[0]
				# ignore the complex STRs
				if '-' in filenamecols[1]:
					continue
				file_strsize = int(filenamecols[1])
				# Check to make sure we have the right gene
				if file_gene != gene:
					continue
				if LowerBoundOnly:
					if file_strsize != pathogenicLowerBound:
						continue
				EHDNMotifRank,EHDNMotifSTRRank,EHDNMotifzscore = ParseEHDNMotif(ehdnmotiffile,motif)
				outfile.write("%s,%d,%s,%d,EHDN_Motif,%d,%d,%.3f\n"%(gene,pathogenicLowerBound,motif,file_strsize,EHDNMotifRank,EHDNMotifSTRRank,EHDNMotifzscore))

		if len(STRetch_FileList) > 0:
			for stretchfile in STRetch_FileList:
				# split filename to get rid of the filepath, and then split the name into columns
				filenamecols = stretchfile.split('/')[-1].split('_')
				file_gene = filenamecols[0]
				# ignore the complex STRs
				if '-' in filenamecols[1]:
					continue
				file_strsize = int(filenamecols[1])
				# Check to make sure we have the right gene
				if file_gene != gene:
					continue
				if LowerBoundOnly:
					if file_strsize != pathogenicLowerBound:
						continue
				# Parse STRetch, retrieve size, rank and p-value
				STRetchSize,STRetchRank,STRetchPval = ParseSTRetch(stretchfile,chrom,start,end,motif)
				outfile.write("%s,%d,%s,%d,STRetch,%d,%d,%s\n"%(gene,pathogenicLowerBound,motif,file_strsize,STRetchRank,STRetchRank,STRetchPval))
			

# The main function
def Main():
	ARGS = GetArgs()
	ParseTableCollateOutput(ARGS.Table,ARGS.EHDN_Locuses,ARGS.EHDN_Motifs,ARGS.STRetches,ARGS.Output,ARGS.Lower)
	
	# test motif expansion
	#print(MotifToPossibleMotifs('GGCCTG'))
	#print(MotifToPossibleMotifs('CAG'))
	#print(MotifToPossibleMotifs('CGG'))

if __name__=="__main__":
    Main()




