'''
Created on 11 mai 2015

@author: gdavidson
'''

# Compares a gene list to a homer annotatePeaks.py output and finds common genes.
# Arguments:
#-i <genes.txt>, gene list with two tab separated fields per line (ENSEMBL Gene ID, gene common name).
#-h <homerOut.tsv>, annotated file from homer with ensembl transcriptID as reference.
#-e <biomartOut.tsv>, biomart output for your ensembl release and organism (two columns: geneID, transcriptID).
# Optional:
#-b true, if your '-i' file is an homer annotation file. 

import sys
import getopt
import pylab as P

progHelp = "Compares a gene list to a homer annotatePeaks.py output and finds common genes.\nArguments:\n-i <genes.txt>, gene list with two tab separated fields per line (ENSEMBL Gene ID, gene common name).\n-h <homerOut.tsv>, annotated file from homer with ensembl transcriptID as reference.\n-e <biomartOut.tsv>, biomart output for your ensembl release and organism (two columns: geneID, transcriptID).\n Optional:\n-b true, if your '-i' file is an homer annotation file."

def get_params(argv):
    try:
        opts, args = getopt.getopt(argv, "i:h:e:b:", ["infilename", "homerFile", "ensemblFile", "bothHomer"])
    except getopt.GetoptError:
        sys.exit("Invalid argument:\n"+progHelp)
    infilename = "none"
    homerFile = "none"
    ensemblFile = "none"
    bothHomer = False
    for opt,arg in opts:
        if opt =='-i':
            infilename = arg
        if opt =="-h":
            homerFile = arg
        if opt =="-e":
            ensemblFile = arg
        if opt == '-b':
            bothHomer = True
    return infilename, homerFile, ensemblFile, bothHomer

def getHomerDictionary(homerFile, ensemblTransToGeneMap):
    homerDic = {}
    distanceMap = {}
    hFile = open(homerFile, "r")
    for line in hFile:
        transID = str(line).split("\t")[10].upper().strip()
        if transID in ensemblTransToGeneMap:
            geneID = ensemblTransToGeneMap[transID]
            peakID = str(line).split("\t")[0]
            homerDic[geneID] = peakID
            distanceToTSS = str(line).split("\t")[9]
            distanceMap[geneID] = distanceToTSS
    hFile.flush()
    hFile.close()
    print "Reading homer output '"+str(homerFile)+"': "+str(len(homerDic))+" genes."
    return homerDic, distanceMap

def getCommonGenesList(refList, geneDic):
    commonGenes = []
    refIDList = []
    for gene in refList:
        if gene in geneDic:
            commonGenes.append(gene)
            refID = geneDic[gene]
            refIDList.append(refID)
    print "Comparing gene lists: "+str(len(commonGenes))+" genes in common."
    return commonGenes, refIDList
    
def getGeneList(infilename):
    iFile = open(infilename, "r")
    geneList = []
    geneIDToNameMap = {}
    for line in iFile:
        geneID = str(line).split("\t")[0].upper().strip()
        geneList.append(geneID)
        geneName = str(line).split("\t")[1].upper().strip()
        geneIDToNameMap[geneID] = geneName
    iFile.flush()
    iFile.close()
    print "Reading geneID file '"+str(infilename)+"': "+str(len(geneList))+" genes."
    return geneList, geneIDToNameMap

def writeList(inList, outName, nameMap = "none"):
    outFile = open(outName, "w")
    for elem in inList:
        name = "NA"
        elemUp = str(elem).strip().upper()
        if elemUp in nameMap:
            name = nameMap[elemUp]
        outFile.write(str(elem).strip()+"\t"+str(name)+"\n")
    print "Writing file '"+str(outName)+"': "+str(len(inList))+" lines."
    outFile.flush()
    outFile.close()
    
def getEnsemblTransToGeneMap(ensemblFile):
    eFile = open(ensemblFile,"r")
    geneMap = {}
    for line in eFile:
        geneID = str(line).split("\t")[0].upper().strip()
        transID = str(line).split("\t")[1].upper().strip()
        geneMap[transID] = geneID
    eFile.flush()
    eFile.close()
    print "Reading biomart output '"+str(ensemblFile)+"': "+str(len(geneMap))+" transcripts."
    return geneMap

def getDistanceList(commonGenesList, distanceMap):
    distanceList = []
    for gene in commonGenesList:
        distance = distanceMap[gene]
        distanceList.append(int(distance))
    return distanceList

def histDistances(commonGenesList, distanceMap):
    distanceList = getDistanceList(commonGenesList, distanceMap)
    dbin = [-12000,-10000,-8000,-6000,-4000,-2000,0, 2000, 4000,6000,8000,10000,12000]
    P.hist(distanceList, bins=dbin, color = 'blue')
    P.xlabel('Distance to TSS (pb)')
    P.ylabel('Count')
    P.xticks(dbin)
    P.show()
 
def compareGeneMaps(geneMap, refGeneMap):
    commonGenes = []
    for geneID in geneMap.keys():
        if geneID in refGeneMap:
            geneAnn = geneMap[geneID]
            commonGenes.append(geneAnn)
    print "Comparing gene lists: "+str(len(commonGenes))+" genes in common."
    return commonGenes
    
def getGeneMap(infilename, ensemblTransToGeneMap):
    geneMap = {}
    hFile = open(infilename, "r")
    for line in hFile:
        transID = str(line).split("\t")[10].upper().strip()
        if transID in ensemblTransToGeneMap:
            geneID = ensemblTransToGeneMap[transID]
            geneMap[geneID] = line
    hFile.flush()
    hFile.close()
    print "Reading homer output '"+str(infilename)+"': "+str(len(geneMap))+" genes."
    return geneMap
    
if __name__ == '__main__':
    infilename, homerFile, ensemblFile, bothHomer = get_params(sys.argv[1:])
    if infilename == "none" or homerFile == "none":
        sys.exit(progHelp)
    if bothHomer == False:   
        ensemblMap = getEnsemblTransToGeneMap(ensemblFile)
        homerDic,distanceMap = getHomerDictionary(homerFile, ensemblMap)
        refList, geneIDToNameMap = getGeneList(infilename)
        commonGenesList, peakList = getCommonGenesList(refList, homerDic)
        writeList(commonGenesList, "common_genes.txt", geneIDToNameMap)
        writeList(peakList,"common_genes_reference_sequences.txt")
        histDistances(commonGenesList, distanceMap)
    else:
        ensemblMap = getEnsemblTransToGeneMap(ensemblFile)
        geneMap = getGeneMap(infilename, ensemblMap)
        refGeneMap = getGeneMap(homerFile, ensemblMap)
        commonGenes = compareGeneMaps(geneMap, refGeneMap)
        writeList(commonGenes, "common_genes.xls")       
    print "Finished."  
    