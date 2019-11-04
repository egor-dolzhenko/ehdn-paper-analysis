import sys
import argparse


def GetArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-C","--Complex",help="These are complex STRs")
	parser.add_argument("-T","--Table",help="The STR Table for which you are going to collate outputs")
	parser.add_argument("-D","--Directory",help="The directory which has the output data from the EHdn analysis")
	parser.add_argument("-O","--Output",help="Output table for results")
	args = parser.parse_args()
	return args


# This function parses a STRetch tsv output
# And with a given locus it will search for that locus 
# And extract the rank and genotyped size
def ParseSTRetch(stretchInfile,locus):
	return


# This function parses an EHDN tsv output from the outlier script
# And with a given locus it will look for that locus and 
def ParseEHDNLocus(ehdnInfile,locus):
	
	return


def ParseEHDNMotif(ehdnInfile,motif):
	return


def ParseTableCollateOutput(strTableFilename,ehdnOutputDir,outFilename):
	outfile = open(outFilename,'w')
	strfile = open(strTableFilename,






def Main():
	ARGS = GetArgs()
	

if __name__=="__main__":
    Main()







