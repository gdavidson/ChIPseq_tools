'''
Created on 20 mai 2015

@author: gdavidson
'''
#Writes the genomic locations of a MEME-found motif in bed format.
# Arguments:
#-m <meme_motif.txt>, copy paste the locations from the MEME HTML output into a .txt file (example line:'1684.    hg19_ct_macs2summits_2357_sox10_peak_3016    -    49    1.24e-8    ACAACAACAC    AAAAGGCCCCTTTGT    TACGGCCCTG')
#-p <macs_peaks.bed>, MACS predicted peaks in bed format.
#-w <motif width(integer)>, width of the motif (you can fiddle with the size)
#-n <experiment name>, name of the experiment, must be the '-n' argument you used in MACS.
# Optional:
#-f true, if your '-m' file is an output from FiMo

import sys
import getopt

progHelp = "Writes the genomic locations of a MEME-found motif in bed format.\nArguments:\n-m <meme_motif.txt>, copy paste the locations from the MEME HTML output into a .txt file (example line:'1684.    hg19_ct_macs2summits_2357_sox10_peak_3016    -    49    1.24e-8    ACAACAACAC    AAAAGGCCCCTTTGT    TACGGCCCTG')\n-p <macs_peaks.bed>, MACS predicted peaks in bed format.\n-w <motif width(integer)>, width of the motif (you can fiddle with the size)\n-n <experiment name>, name of the experiment, must be the '-n' argument you used in MACS.\nOptional:\n-f true, if your '-m' file is an output from FiMo"

def get_params(argv):
    try:
        opts, args = getopt.getopt(argv, "m:p:w:n:f:", ["memeFile", "peakFile", "motifWidth", "expName","fimo"])
    except getopt.GetoptError:
        sys.exit("Invalid argument:\n"+progHelp)
    memeFile = "none"
    peakFile = "none"
    motifWidth = "none"
    expName = "none"
    fimo = False
    for opt,arg in opts:
        if opt =='-m':
            memeFile = arg
        if opt =="-p":
            peakFile = arg
        if opt == "-w":
            motifWidth = arg
        if opt == "-n":
            expName = arg
        if opt == "-f":
            fimo = True  
    return memeFile, peakFile, motifWidth, expName, fimo

def parsePeakNumber(peakID, expName):
    peakNum = peakID.split(expName+"_")[1]
    return peakNum

# returns hashmap: key = peak number, value = list of motif starts     
def getMotifPeakLocationMap(memeFile, expName, fimo = False):
    mFile = open(memeFile, "r")
    memeMap = {}
    printBool = True
    for line in mFile:
        if line.startswith("#"):
            continue
        peakID = str(line).split("\t")[1]
        if fimo == True:
            motifStart = str(line).split("\t")[2]
            if printBool == True:
                print "File from FiMo, reading motif start in column 3."
                printBool = False
        else:
            motifStart = str(line).split("\t")[3]
            if printBool == True:
                print "File from MEME, reading motif start in column 4."
                printBool = False
        peakNum = str(parsePeakNumber(peakID, expName)).strip()
        if not peakNum in memeMap:
            memeMap[peakNum] = [motifStart]
        else:
            starts = memeMap[peakNum]
            starts.append(motifStart)
            memeMap[peakNum] = starts
            
    print "Reading meme output '"+str(memeFile)+"': motif found in "+str(len(memeMap))+" sequences."
    return memeMap

def getPeaksWithMotif(peakFile, memeMap, expName):
    pFile = open(peakFile, "r")
    peakList = []
    for line in pFile:
        if str(line).startswith("#") or str(line).startswith("track"):
            continue
        peakID = str(line).split("\t")[3]
        peakNum = str(parsePeakNumber(peakID, expName)).strip()
        if peakNum in memeMap:
            peakList.append(line)
    print "Reading peakfile '"+str(peakFile)+"': found "+str(len(peakList))+" peaks with motif."
    return peakList

def getMotifList(peakList, memeMap, expName, motifWidth):
    motifList = []
    for line in peakList:
        peakID = str(line).split("\t")[3].strip()
        peakNum = str(parsePeakNumber(peakID, expName)).strip()
        starts = memeMap[peakNum]
        for motifStart in starts:
            peakStart= str(line).split("\t")[1]
            newStart = int(peakStart)+int(motifStart)
            newEnd = int(newStart)+int(motifWidth)
            newLine =  str(line).split("\t")[0]+"\t"+str(newStart)+"\t"+str(newEnd)+"\t"+peakID
            motifList.append(newLine)
    return motifList

def writeBed(bedList, outfileName):
    outfile = open(outfileName, "w")
    for line in bedList:
        outfile.write(line+"\n")
    print "Writing: '"+outfileName+"', "+str(len(bedList))+" lines."
    outfile.flush()
    outfile.close()
           
if __name__ == '__main__':
    memeFile, peakFile, motifWidth, expName, fimo = get_params(sys.argv[1:])
    if memeFile == "none" or peakFile == "none" or motifWidth == "none" or expName == "none":
        sys.exit(progHelp)
    memeMap = getMotifPeakLocationMap(memeFile, expName, fimo)
    peakList = getPeaksWithMotif(peakFile, memeMap, expName)
    motifList = getMotifList(peakList, memeMap, expName, motifWidth)
    writeBed(motifList, expName+"_motif.bed")
    print "Finished."
     