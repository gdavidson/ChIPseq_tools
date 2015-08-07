'''
Created on 22 mai 2015

@author: gdavidson
'''

# Writes a table summarizing single-nucleotide variants (SNVs) found in a motif.
# Arguments:
# -i <intersect_motif_peaks.bed>, output file from intersectBed with '-a': clinvar Bed file, '-b': the motif locations (from getMotifLocations.py) and '-wo'.
# Optional:
# -v <clinvarMain.bed>, clinvar bed file (with all columns) containing more information on the variants. The output will contain more details on each SNV.
# -h <homer_annotation.txt>, homer annotation  of the peak file containing the motif. The output will contain annotation details of the peak containing the variant.
# -d <disease_name.txt>, clinvar disease file (ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/disease_names)

import sys
import getopt


progHelp = "Writes a table summarizing single-nucleotide variants (SNVs) found in a motif.\n Arguments:\n-i <intersect_motif_peaks.bed>, output file from intersectBed with '-a': clinvar Bed file, '-b': the motif locations (from getMotifLocations.py) and '-wo'.\n Optional:\n-v <clinvarMain.bed>, clinvar bed file (with all columns) containing more information on the variants. The output will contain more details on each SNV.\n-h <homer_annotation.txt>, homer annotation  of the peak file containing the motif. The output will contain annotation details of the peak containing the variant.\n-d <disease_name.txt>, clinvar disease file (ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/disease_names)"

def get_params(argv):
    try:
        opts, args  = getopt.getopt(argv, "i:v:h:d:", ["infile"])
    except getopt.GetoptError:
        sys.exit("Invalid argument:\n"+progHelp)
    infilename = "none"
    variantFile = "none"
    homerFile = "none"
    diseaseFile = "none"
    for opt,arg in opts:
        if opt =='-i':
            infilename = arg
        if opt =='-v':
            variantFile = arg
        if opt == '-h':
            homerFile = arg
        if opt == '-d':
            diseaseFile = arg

    return infilename, variantFile, homerFile, diseaseFile

def parseBed(bedFile):
    bFile = open(bedFile,"r")
    motifWidth = 0
    variantTupleList = []
    for line in bFile:
        vStart = int(str(line).split("\t")[1])
        vEnd = int(str(line).split("\t")[2])
        if (vEnd - vStart)>2:
            continue
        mStart = int(str(line).split("\t")[5])
        mEnd = int(str(line).split("\t")[6])
        
        if motifWidth == 0:
            motifWidth = mEnd-mStart
        peakID = str(line).split("\t")[7]
        varID = str(line).split("\t")[3]
        varLocation = (vStart - mStart)+1
        varTuple = (varLocation, varID, peakID)
        variantTupleList.append(varTuple)
    bFile.close()
    bFile.close()
    print "Reading '"+infilename+"', found "+str(len(variantTupleList))+" SNVs, motif length: "+str(motifWidth)+"pb."
    return variantTupleList, motifWidth
    
def mapTSVcolumns(tsvFile, col1, col2):
    tFile = open(tsvFile, "r")
    colMap = {}
    for line in tFile:
        key = str(line).split("\t")[(col1-1)].strip()
        val = str(line).split("\t")[(col2-1)].strip()
        colMap[key] = val
    print "Reading '"+tsvFile+"', mapping column "+str(col1)+" to column "+str(col2)+"."
    tFile.flush()
    tFile.close()
    return colMap

def getSNVRepresentation(seqLength, snvPosition, nucleotide):
    snvRep = []
    for i in range(0,seqLength):
        if i == (snvPosition-1):
            snvRep.append("["+nucleotide+"]")
        else:
            snvRep.append("N")
    snvRep = "".join(snvRep)
    return snvRep
        
def parseMedGenID(condition):
    medID = "NA"
    for fullID in condition.split(","):
        if fullID.startswith("MedGen"):
            medID = fullID.split(":")[1].strip()
    return medID
                 
def getTableList(variantTupleList, motifWidth, idToNameMap = "none", idToClinSinMap = "none",  homerMap = "none", idToConditionMap = "none", diseaseMap = "none"):
    tableList = []
    for vTuple in variantTupleList:
        varLoc = vTuple[0]
        varID = str(vTuple[1]).strip()
        peakID = vTuple[2]
        varName = ""
        nucleotide = "X"
        if not idToNameMap == "none":
            varName = idToNameMap[varID]
            if varName.startswith("A") or varName.startswith("T") or varName.startswith("G") or varName.startswith("C"):
                nucleotide = varName.split(">")[1]
            clinSin = idToClinSinMap[varID]
        distance, geneName, geneDesc = "NA", "NA", "NA"
        if not homerMap == "none":
            geneTuple = homerMap[peakID]
            distance = geneTuple[0]
            geneName = geneTuple[1]
            geneDesc = geneTuple[2]
        if not idToConditionMap == "none" and not diseaseMap == "none":
            condition = idToConditionMap[varID]
            medGenID = parseMedGenID(condition)
            if medGenID in diseaseMap:
                diseaseName = diseaseMap[medGenID]
            else:
                diseaseName = medGenID
        snvMotif = getSNVRepresentation(motifWidth, varLoc, nucleotide)
        tableRow = snvMotif+"\t"+varID+"\t"+peakID+"\t"+str(varLoc)+"\t"+varName+"\t"+clinSin+"\t"+diseaseName+"\t"+geneName+"\t"+geneDesc+"\t"+distance
        tableList.append(tableRow)
    return tableList

def writeTable(tableList):
    header = "#Motif\tSNV AccVar\tFound in Peak\tPosition in motif\tSNV Name\tClinical Significance\tCondition\tNearest Gene Name\tGene Description\tDistance to TSS"
    outFile = open("variantTable.xls","w")
    outFile.write(header+"\n")
    outFile.write(""+"\n")
    for row in tableList:
        outFile.write(row+"\n")
    outFile.flush()
    outFile.close()
    print "Writing 'variantTable.xls: "+str((len(tableList)+2))+" lines."
        
def getHomerMap(homerFile):
    hFile = open(homerFile, "r")
    hMap = {}
    for line in hFile:
        pID = str(line).split("\t")[0]
        distance = str(line).split("\t")[9]
        geneName = str(line).split("\t")[15]
        geneDesc = str(line).split("\t")[17]
        geneTuple = (distance, geneName, geneDesc)
        hMap[pID] = geneTuple
    print "Reading '"+str(homerFile)+"', "+str(len(hMap))+" annotations."
    return hMap
        
def getDiseaseMap(diseaseFile):
    dFile = open(diseaseFile, "r")
    diseaseMap = {}
    for line in dFile:
        dName = line.split("\t")[0]
        medID = line.split("\t")[2]
        diseaseMap[medID] = dName
    dFile.flush()
    dFile.close()
    return diseaseMap
                            
if __name__ == '__main__':
    infilename, variantFile, homerFile, diseaseFile = get_params(sys.argv[1:])
    if infilename == "none":
        sys.exit(progHelp)
    variantTupleList, motifWidth = parseBed(infilename)
    idToNameMap, idToClinSinMap, idToConditionMap, homerMap, diseaseMap = "none", "none", "none", "none", "none"
    if not variantFile == "none":
        idToNameMap = mapTSVcolumns(variantFile, 20, 4)
        idToClinSinMap = mapTSVcolumns(variantFile, 20, 17)
        idToConditionMap = mapTSVcolumns(variantFile, 20, 22)
    if not homerFile == "none":
        homerMap = getHomerMap(homerFile)
    if not diseaseFile == "none":
        diseaseMap = getDiseaseMap(diseaseFile)
    varList = getTableList(variantTupleList, motifWidth, idToNameMap, idToClinSinMap, homerMap, idToConditionMap, diseaseMap)
    writeTable(varList)
    print "Finished"