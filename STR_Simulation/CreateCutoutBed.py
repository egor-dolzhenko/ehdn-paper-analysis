# Phillip Richmond
# 2019-07-23


# The purpose of this script is to generate expanded fasta files from STR loci
# Steps:
# 1. 


# Imports
import argparse
import sys


def GetArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-O","--OutputBed",help="Output Filename", required=True)
    parser.add_argument("-S","--SampleID",help="Sample ID, this is a fasta identifier in the format: \
            DIP2B_GGC_12-50898784-50898805_12-50896784-50900805 \
            <gene>_<repeatMotif>_<repeatPosition>_<repeatRegionToCutOut>",required=True)
    parser.add_argument("-C","--ChromosomeSizes",help="Genome chromosome sizes in the format of chr\tendPos",required=True)
    args = parser.parse_args()
    return args



# For the variables here:
# p == print
# g == genomic (e.g. info from chromFile)
# r == region, extracted from the sample identifier
def CreateCutoutBed(chromFilename,outputFilename,sampleID):
    chromFile = open(chromFilename,'r')
    outfile = open(outputFilename,'w')
    regionInfo = sampleID.split('_')[3]
    # region chrom, region start, region end from the sample identifier
    rChrom,rStart,rEnd = regionInfo.split('-')
   
    # genome chrom,end information
    for line in chromFile:
        cols = line.strip('\n').split('\t')
        # The chrom file stores the chrom and genomic end of each chromosome
        gChrom = cols[0]
        gEnd = cols[1]
        # If the cutout chrom is the same as the genomic chrom, then the start is 1, and the end is the region start
        # and then there is another line which is the 
        if gChrom == rChrom:
            pChrom = gChrom
            # write the first half, e.g. the chromosome up to the beginning of the repeat region
            pStart = 0
            pEnd = rStart
            outfile.write('%s\t%s\t%s\n'%(pChrom,pStart,pEnd))
            # write the second half, which is the repeat region end to the end of the chromosome
            pStart = rEnd
            pEnd = gEnd
            outfile.write('%s\t%s\t%s\n'%(pChrom,pStart,pEnd))
        else:
            pStart = 0
            pEnd = gEnd
            pChrom = gChrom
            outfile.write('%s\t%s\t%s\n'%(pChrom,pStart,pEnd))


def Main():
    ARGS = GetArgs()
    CreateCutoutBed(ARGS.ChromosomeSizes,ARGS.OutputBed,ARGS.SampleID)



if __name__=="__main__":
    Main()






