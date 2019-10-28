import sys,os,argparse


def GetArgs():
	# Redone with Argparse
	parser = argparse.ArgumentParser()
	parser.add_argument("-T","--Type",help="GangSTR || STRetch", required=True)
	parser.add_argument("-I","--Infile",help="Input File",required=True)
	parser.add_argument("-O","--Outfile",help="Output File",required=True)
	parser.add_argument("--Min",help="Min Motif Length",required=True,type=int)
	parser.add_argument("--Max",help="Max Motif Length",required=True,type=int)
	args = parser.parse_args()
	return args



#STRetch 
# chr1	10000	10468	TAACCC	77.2
#GangSTR 
# chr1	14070	14081	4	CCTC

def ParseFile(infilename,outfilename,maxlength,minlength,Type):
	infile = open(infilename,'r')
	outfile = open(outfilename,'w')
	for line in infile:
		cols = line.strip('\n').split('\t')
		if Type == 'GangSTR':
			motifCol = cols[4]
		elif Type == 'STRetch':
			motifCol= cols[3]
		if (len(motifCol) < minlength) or (len(motifCol) > maxlength):
			continue
		outfile.write(line)


def Main():
	# get args
	args = GetArgs()
	ParseFile(args.Infile,args.Outfile,args.Max,args.Min,args.Type)


if __name__=="__main__":
	Main()

