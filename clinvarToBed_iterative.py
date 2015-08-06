#coding:utf-8
#!/usr/bin/python


# Converts the ClinVar DB XML file to bed format
#Requires:
#-elementtree: http://effbot.org/zone/element-index.htm#installation
# Usage: python clinvarToBed_iterative.py -i ClinVarFullRelease_2015-04.xml
#Args:
#-i: input file, clinVar XML file
#Output:
#'clinvar.bed': clinVar DB in bed format

import getopt
import sys
from elementtree.ElementTree import iterparse
 
def get_params(argv):
    try:
        opts, args = getopt.getopt(argv, "i:", ["infile"])
    except getopt.GetoptError:
        sys.exit("Converts the ClinVar DB XML file to bed format\nUsage: python clinvarToBed_iterative.py -i ClinVarFullRelease_2015-04.xml\nArgs:\n-i: input file, clinVar XML file\nOutput:\n'clinvar.bed': clinVar DB in bed format")
    infilename = "none"
    for opt,arg in opts:
        if opt =='-i':
            infilename = arg
    return infilename


def writeSequenceLocations(infileName):
    outFile = open('clinvar.bed', 'w')
    seqCount = 0
    print "Parsing file and writing sequence locations..."
    context = iterparse(infileName, events=("start", "end"))
    context = iter(context)
    event, root = context.next()
    for event, elem in context:
        if elem.tag == "MeasureRelationship":
            elem.clear()
            root.clear()
            continue
        for access in elem.findall("ClinVarAccession"):
            acc = access.get("Acc")
        for seqLoc in elem.findall("SequenceLocation"):
            if seqLoc.get("Assembly") == "GRCh37" and event == "end":
                chrom = "chr"+seqLoc.get("Chr")
                start = seqLoc.get("start")
                stop = seqLoc.get("stop")
                if seqLoc.get("Strand") == "+" or seqLoc.get("Strand") == "-":
                    strand = seqLoc.get("Strand")
                else:
                    strand = "."
                seqTuple = (chrom, start, stop, acc, 0, strand)
                if str(seqTuple[1]) != "None" and str(seqTuple[2]) != "None":
                    line = str(seqTuple[0])+"\t"+str(seqTuple[1])+"\t"+str(seqTuple[2])+"\t"+str(seqTuple[3])+"\t"+str(seqTuple[4])+"\t"+str(seqTuple[5])+"\n"
                    outFile.write(line)
                    seqCount= seqCount+1
                elem.clear()
                root.clear()
    print "Parsing done!"
    outFile.flush()
    outFile.close()
    print "Output file 'clinvar.bed': "+str(seqCount)+" locations found."
    print "Finished!"
    
if __name__ == '__main__':
    infileName = get_params(sys.argv[1:])
    if infileName=="none":
        sys.exit("Converts the ClinVar DB XML file to bed format\nUsage: python clinvarToBed_iterative.py -i ClinVarFullRelease_2015-04.xml\nArgs:\n-i: input file, clinVar XML file\nOutput:\n'clinvar.bed': clinVar DB in bed format")
    writeSequenceLocations(infileName)
