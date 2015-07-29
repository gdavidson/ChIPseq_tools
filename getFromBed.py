'''
Created on 7 mai 2015

@author: gdavidson
'''

# Allows operations on BED Files.
# Arguments:
#-i <file1.bed>, input file in bed format
# Options (choose one):
#-o <offset(integer)>, changes the coordinates of each entry of the input file by the offset amount (start-offset, end+offset)
#-b <file2.bed>, finds all entries in file2 which ID/name matches any ID/name in file1 (4th column). Can be a simple list of IDs (one ID per line) then set '-f' to 1.
# Optional:
#-f <field number(integer)>, only with '-b', field (column) which contains IDs for file2.bed (default:4)


import sys
import getopt

progHelp = "Allows operations on BED Files.\nArguments:\n-i <file1.bed>, input file in bed format\nOptions (choose one):\n-o <offset(integer)>, changes the coordinates of each entry of the input file by the offset amount (start-offset, end+offset)\n-b <file2.bed>, finds all entries in file2 which ID/name matches any ID/name in file1(4th column). Can be a simple list of IDs (one ID per line) then set '-f' to 1.\nOptional:\n-f <field number(integer)>, only with '-b', field (column) which contains IDs for file2.bed (default:4)"

def get_params(argv):
    try:
        opts, args = getopt.getopt(argv, "i:o:b:f:", ["infilename", "offset", "secondBedFile","idField"])
    except getopt.GetoptError:
        sys.exit("Invalid argument:\n"+progHelp)
    infilename = "none"
    queryType = "none"
    offset = "none"
    secondBedFile = "none"
    idField = 4
    for opt,arg in opts:
        if opt =='-i':
            infilename = arg
        if opt =="-o":
            offset = int(arg)
            queryType = "Offset"
        if opt == "-b":
            secondBedFile = arg
            queryType = "CommonLinesID"
        if opt == "-f":
            idField = int(arg)
    return infilename, offset, queryType, secondBedFile, idField


def getBedList(infileName):
    inputFile = open(infileName, "r")
    bedList = []
    for line in inputFile:
        bedList.append(line)
    inputFile.flush()
    inputFile.close()
    print "Reading Input file "+infileName+" : "+str(len(bedList))+" lines."
    return bedList

def changeCoord(bedList, offset):
    newList = []
    for entry in bedList:
        fnum=0
        newEntry=""
        for field in entry.split("\t"):
            if fnum == 1:
                newStart = int(field)-offset
                newEntry = newEntry+str(newStart)+"\t"
            elif fnum == 2:
                newEnd = int(field)+offset
                newEntry = newEntry+str(newEnd)+"\t"
            else:
                newEntry = newEntry+str(field)+"\t"
            fnum = fnum+1
        newList.append(newEntry.strip())
    return newList
                
def writeBed(bedList, outfileName):
    outfile = open(outfileName, "w")
    for line in bedList:
        outfile.write(line+"\n")
    print "Writing '"+outfileName+"': "+str(len(bedList))+" lines."
    outfile.flush()
    outfile.close()
    
def getIDMap(bedList):
    idMap = {}
    for line in bedList:
        bID = str(line).split("\t")[3].strip().upper()
        idMap[bID] = line
    return idMap

def getLinesWithIDs(bedList, idMap, idField = 4):
    foundLines = []
    for line in bedList:
        sbID = str(line).split("\t")[(idField-1)].strip().upper()
        if sbID in idMap:
            foundLines.append(str(idMap[sbID]).strip())
    return foundLines
    
if __name__ == '__main__':
    infileName, offset, queryType, secondBedFile, idField = get_params(sys.argv[1:])
    if infileName == "none" or queryType == "none":
        sys.exit(progHelp)
    if queryType == "Offset":
        bedList = getBedList(infileName)
        newBedList = changeCoord(bedList, offset)
        outfileName = str(infileName).split(".")[0]+"_offset"+str(offset)+".bed"
        writeBed(newBedList, outfileName)
    if queryType == "CommonLinesID":
        bedList = getBedList(infileName)
        secondBedList = getBedList(secondBedFile)
        idMap = getIDMap(bedList)
        commonLineList = getLinesWithIDs(secondBedList, idMap, idField)
        writeBed(commonLineList, "commonLines.bed")
    print "Finished !"