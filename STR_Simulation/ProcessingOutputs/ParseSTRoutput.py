import sys
import argparse


def GetArgs():
	parser = argparse.ArgumentParser()
	parser.add_argument("-C","--Complex",help="These are complex STRs")
	args = parser.parse_args()
	return args


# This function parses a STRetch tsv output
# And with a given locus it will search for that locus 
# And extract the rank and genotyped size
def ParseSTRetch(stretchInfile,locus):

	return


# This function parses an EHDN tsv output from the outlier script
# And with a given locus it will look for that locus and 
def ParseEHDN(ehdnInfile,locus):
	return










def Main():
	ARGS = GetArgs()


if __name__=="__main__":
    Main()







