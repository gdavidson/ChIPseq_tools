'''
Created on 15 mai 2015

@author: gdavidson
'''

# Changes first field of a GTF file to match UCSC standards (chr1, chr2 ...).
# Arguments:
# -i <file.gtf>, the GTF file to edit.
# -a <assembly>, name of the assembly (hg19 or mm9).
# Optional:
# -p <mm9_patches.tsv>, the file containing patches and loci for your assembly (alt_scaffold_placement.txt in genebank genome DB).

import sys
import getopt

progHelp = "Changes first field of a GTF file to match UCSC standards (chr1, chr2 ...).\nArguments:\n-i <file.gtf>, the GTF file to edit.\n-a <assembly>, name of the assembly (hg19 or mm9).\nOptional:\n-p <mm9_patches.tsv>, the file containing patches and loci for your assembly (alt_scaffold_placement.txt in genebank genome DB)."

def get_params(argv):
    try:
        opts, args = getopt.getopt(argv, "i:p:a:", ["infile"])
    except getopt.GetoptError:
        sys.exit("Invalid argument:\n"+progHelp)
    infilename = "none"
    patchfile = "none"
    assembly = "none"
    for opt,arg in opts:
        if opt =='-i':
            infilename = arg
        if opt =="-p":
            patchfile = arg
        if opt == "-a":
            assembly = str(arg).lower()
    return infilename, patchfile, assembly

def getGTFList(infileName):
    inputFile = open(infileName, "r")
    gtfList = []
    for line in inputFile:
        gtfList.append(line)
    inputFile.flush()
    inputFile.close()
    print "Reading Input file "+infileName+" : "+str(len(gtfList))+" lines."
    return gtfList            

def getPatchMap(patchfile):
    pFile = open(patchfile, "r")
    patchMap = {}
    for line in pFile:
        if line.startswith("#"):
            continue
        fullID = str(line).split("\t")[0]
        contigID = fullID.split("|")[1]
        truncID = contigID.split(".")[0]
        chrom = "chr"+str(line).split("\t")[3]
        patchMap[truncID] = chrom
    print "Reading patches and alternate loci file "+patchfile+" : "+str(len(patchMap))+" entries."
    return patchMap
             
def getNewGTFlocations(gtfList, patchMap = "none"):
    newGTFList = []
    for line in gtfList:
        newLine = ""
        featNum = 0
        newLoc = ""
        for feat in line.split("\t"):
            if featNum == 0:
                if feat.startswith("NT") or feat.startswith("NW"):
                    if feat in patchMap:
                        newLoc = patchMap[feat]
                    else:
                        newLoc = "NA"
                else:
                    newLoc = "chr"+feat
                newLine=newLoc+"\t"
                featNum = featNum+1
            else:
                newLine = newLine+feat+"\t"
        newGTFList.append(newLine.strip())
    return newGTFList

def getPatchMapHuman(patchfile):
    pFile = open(patchfile, "r")
    patchMap = {}
    for line in pFile:
        if line.startswith("#"):
            continue
        fullID = str(line).split("\t")[2].upper().strip()
        chrom = "chr"+str(line).split("\t")[5]
        patchMap[fullID] = chrom
    print "Reading patches and alternate loci file "+patchfile+" : "+str(len(patchMap))+" entries."
    return patchMap

def getNewGTFlocationsHuman(gtfList, patchMap):
    newGTFList = []
    for line in gtfList:
        newLine = ""
        featNum = 0
        newLoc = ""
        for feat in line.split("\t"):
            if featNum == 0:
                feat = str(feat).upper().strip()
                if feat in patchMap:
                    newLoc = patchMap[feat]
                elif not feat in patchMap and feat.startswith("HSCHR"):
                    chrom = feat.split("_")[0]
                    chrnum = chrom.split("CHR")[1]
                    newLoc = "chrom"+str(chrnum)
                elif not feat in patchMap and feat.startswith("GL"):
                    newLoc = "NA"
                else:
                    newLoc = "chrom"+str(feat)
                newLine=newLoc+"\t"
                featNum = featNum+1
            else:
                newLine = newLine+feat+"\t"
        newGTFList.append(newLine.strip())
    return newGTFList

def writeList(inList, outName):
    outFile = open(outName, "w")
    for elem in inList:
        outFile.write(str(elem)+"\n")
    print "Writing file '"+str(outName)+"': "+str(len(inList))+" lines."
    outFile.flush()
    outFile.close()

if __name__ == '__main__':
    infilename, patchfile, assembly = get_params(sys.argv[1:])
    if infilename == "none" or assembly == "none":
        sys.exit(progHelp)
    gtfList = getGTFList(infilename)
    patchMap = "none"
    if not patchfile == "none":
        if assembly =="mm9":
            patchMap = getPatchMap(patchfile)
        elif assembly == "hg19":
            patchMap = getPatchMapHuman(patchfile)
    if assembly == "mm9":
        newGTF = getNewGTFlocations(gtfList, patchMap)
    elif assembly == "hg19":
        newGTF = getNewGTFlocationsHuman(gtfList, patchMap)
    writeList(newGTF, infilename.split(".")[0]+"_locations_edited.gtf")
    print "Finished."
        
    