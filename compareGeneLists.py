'''
Created on 11 mai 2015

@author: gdavidson
'''

# Compares a gene list to a homer annotatePeaks.py output and finds common genes.
# Requires:
# -matplotlib: http://matplotlib.org/
# -matplotlib-venn: https://pypi.python.org/pypi/matplotlib-venn
# Arguments:
#-i <genes.txt>, gene list with two tab separated fields per line (ENSEMBL Gene ID, gene common name).
#-h <homerOut.tsv>, annotated file from homer with ensembl transcriptID as reference.
#-e <biomartOut.tsv>, biomart output for your ensembl release and organism (two columns: geneID, transcriptID). If -g, first column has to be geneID then other fields are optional. 
# Optional:
#-b true, if your '-i' file is also a homer annotation file.
#-g true, if your '-h' file is also a list with gene IDs in the first column.
#-q <N (INTEGER)>, only with '-g', divides the '-i' gene list into N sublists and compares them to the '-h' gene list.
#
#Examples:
# Compares a gene list to a homer annotation file:
#python compareGeneLists.py -i siSOX10_upreg_genes.txt -h chip_sox10_peaks_annotations.xls -e mart_export_hg19_release80.tsv
# Compares two annotation files:
#python compareGeneLists.py -i chip_sox10_peaks_annotations.xls -h chip_mitf_peaks_annotations.xls -e mart_export_hg19_release80.tsv -b true
# Compares two gene lists:
#python compareGeneLists.py -i siSOX10_upreg_genes.txt -h siMITF_upreg_genes.txt -e mart_export_hg19_release80.tsv -g true
 

import sys
import getopt
import pylab as P
from matplotlib_venn import venn2
import math

progHelp = "Compares a gene list to a homer annotatePeaks.py output and finds common genes.\nArguments:\n-i <genes.txt>, gene list with two tab separated fields per line (ENSEMBL Gene ID, gene common name).\n-h <homerOut.tsv>, annotated file from homer with ensembl transcriptID as reference.\n-e <biomartOut.tsv>, biomart output for your ensembl release and organism (two columns: geneID, transcriptID). If -g, first column has to be geneID then other fields are optional.\n Optional:\n-b true, if your '-i' file is also a homer annotation file.\n-g true, if your '-h' file is also a list with gene IDs in the first column.\n-q <N (INTEGER)>, only with '-g', divides the '-i' gene list into N sublists and compares them to the '-h' gene list.\n\nExamples:\nCompares a gene list to a homer annotation file:\npython compareGeneLists.py -i siSOX10_upreg_genes.txt -h chip_sox10_peaks_annotations.xls -e mart_export_hg19_release80.tsv\nCompares two annotation files:\npython compareGeneLists.py -i chip_sox10_peaks_annotations.xls -h chip_mitf_peaks_annotations.xls -e mart_export_hg19_release80.tsv -b true\nCompares two gene lists:\npython compareGeneLists.py -i siSOX10_upreg_genes.txt -h siMITF_upreg_genes.txt -e mart_export_hg19_release80.tsv -g true"

def get_params(argv):
    try:
        opts, args = getopt.getopt(argv, "i:h:e:b:g:q:", ["infilename", "homerFile", "ensemblFile", "bothHomer"])
    except getopt.GetoptError:
        sys.exit("Invalid argument:\n"+progHelp)
    infilename = "none"
    homerFile = "none"
    ensemblFile = "none"
    bothHomer = False
    bothLists = False
    quantiles = "none"
    for opt,arg in opts:
        if opt =='-i':
            infilename = arg
        if opt =="-h":
            homerFile = arg
        if opt =="-e":
            ensemblFile = arg
        if opt == '-b':
            bothHomer = True
        if opt == '-g':
            bothLists = True
        if opt == '-q':
            quantiles = int(arg)     
    return infilename, homerFile, ensemblFile, bothHomer, bothLists, quantiles

def getHomerDictionary(homerFile, ensemblTransToGeneMap):
    homerDic = {}
    distanceMap = {}
    hFile = open(homerFile, "r")
    for line in hFile:
        if str(line).isspace() or line.startswith("#"):
            continue
        transID = str(line).split("\t")[10].upper().strip()
        if transID in ensemblTransToGeneMap:
            geneID = ensemblTransToGeneMap[transID]
            peakID = str(line).split("\t")[0]
            if geneID in homerDic:
                peakIDs = homerDic[geneID]
                peakIDs.append(peakID)
                homerDic[geneID] = peakIDs
            else:
                homerDic[geneID] = [peakID]
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
            refIDs = geneDic[gene]
            for refID in refIDs:
                refIDList.append(refID)
    print "Comparing gene lists: "+str(len(commonGenes))+" genes in common."
    return commonGenes, refIDList

def getGeneList(infilename):
    iFile = open(infilename, "r")
    geneList = []
    geneIDToNameMap = {}
    for line in iFile:
        if str(line).isspace() or line.startswith("#"):
            continue   
        geneID = str(line).split("\t")[0].upper().strip()
        geneList.append(geneID)
        geneName = str(line).split("\t")[1].upper().strip()
        geneIDToNameMap[geneID] = geneName
    iFile.flush()
    iFile.close()
    print "Reading gene list '"+str(infilename)+"': "+str(len(geneList))+" genes."
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
        if str(line).isspace():
            continue
        transID = str(line).split("\t")[10].upper().strip()
        if transID in ensemblTransToGeneMap:
            geneID = ensemblTransToGeneMap[transID]
            geneMap[geneID] = line
    hFile.flush()
    hFile.close()
    print "Reading homer output '"+str(infilename)+"': "+str(len(geneMap))+" genes."
    return geneMap

def getEnsemblMap(ensemblFile):
    eFile = open(ensemblFile,"r")
    geneMap = {}
    for line in eFile:
        if str(line).isspace():
            continue
        geneID = str(line).split("\t")[0].upper().strip()
        geneMap[geneID] = line
    eFile.flush()
    eFile.close()
    print "Reading biomart output '"+str(ensemblFile)+"': "+str(len(geneMap))+" genes."
    return geneMap

def getIDList(fileName):
    iFile = open(fileName, "r")
    idList = []
    for line in iFile:
        if str(line).isspace():
            continue
        geneID = str(line).split("\t")[0].upper().strip()
        if not geneID in idList:
            idList.append(geneID)
        else:
            continue
    iFile.flush()
    iFile.close()
    print "Reading gene list '"+str(fileName)+"': "+str(len(idList))+" IDs."
    return idList

def compareLists(list1, list2):
    commonGenes = []
    uniqueGenesL1 = []
    uniqueGenesL2 = []
    for geneID in list1:
        if geneID in list2:
            commonGenes.append(geneID)
        else:
            uniqueGenesL1.append(geneID)
            
    for geneID in list2:
        if geneID in list1:
            continue
        else:
            uniqueGenesL2.append(geneID)
    print "Comparing lists ..."
    print str(len(commonGenes))+" genes in common. ("+str(int((float(len(commonGenes))/len(list1))*100))+"% of list1, "+str(int((float(len(commonGenes))/len(list2))*100))+"% of list2)."
    print str(len(uniqueGenesL1))+" genes unique to the first list"
    print str(len(uniqueGenesL2))+" genes unique to the second list"
    return commonGenes, uniqueGenesL1, uniqueGenesL2        

def writeLines(inList, outName, lineMap):
    oFile = open(outName, "w")
    for geneID in inList:
        if geneID in lineMap:
            line = str(lineMap[geneID]).strip()
        else:
            line = str(geneID).strip()
            print str(geneID)+" not in biomart file."
        oFile.write(line+"\n")
    oFile.flush()
    oFile.close()
    print "Writing '"+outName+"': "+str(len(inList))+" lines."
 
def drawVennDiagram(list1, list2):
    print "Displaying Venn Diagram ..."
    venn2([set(list1), set(list2)], ('list1', 'list2'))
    P.show()
    
def divideIntoSublists(inputList, subListNumber):
    incr = int(math.ceil(float(len(inputList))/subListNumber))
    sublists = []
    for x in xrange(0, len(inputList), incr):
        if x+incr < len(inputList):
            sublists.append(inputList[x:x+incr])
        else:
            sublists.append(inputList[x:len(inputList)])
    return sublists
                                  
if __name__ == '__main__':
    infilename, homerFile, ensemblFile, bothHomer, bothLists, quantiles = get_params(sys.argv[1:])
    if infilename == "none" or homerFile == "none" or ensemblFile == "none":
        sys.exit(progHelp)
    if bothHomer == False and bothLists == False:   
        ensemblMap = getEnsemblTransToGeneMap(ensemblFile)
        homerDic,distanceMap = getHomerDictionary(homerFile, ensemblMap)
        refList, geneIDToNameMap = getGeneList(infilename)
        commonGenesList, peakList = getCommonGenesList(refList, homerDic)
        writeList(commonGenesList, "common_genes.txt", geneIDToNameMap)
        writeList(peakList,"common_genes_reference_sequences.txt")
        histDistances(commonGenesList, distanceMap)
    elif bothHomer == True and bothLists == False:
        ensemblMap = getEnsemblTransToGeneMap(ensemblFile)
        geneMap = getGeneMap(infilename, ensemblMap)
        refGeneMap = getGeneMap(homerFile, ensemblMap)
        commonGenes = compareGeneMaps(geneMap, refGeneMap)
        writeList(commonGenes, "common_genes.xls")
    elif bothHomer == False and bothLists == True:
        ensemblMap = getEnsemblMap(ensemblFile)
        list1 = getIDList(infilename)
        list2 = getIDList(homerFile)
        commonGenes, uniqueL1, uniqueL2 = compareLists(list1, list2)
        writeLines(commonGenes, "common_genes.txt",ensemblMap)
        writeLines(uniqueL1, infilename.split(".")[0]+"_unique.txt", ensemblMap)
        writeLines(uniqueL2, homerFile.split(".")[0]+"_unique.txt", ensemblMap)
        drawVennDiagram(list1, list2)
        if not quantiles == "none":
            print "Dividing '"+infilename+"' into "+str(quantiles)+" quantiles."
            sublists = divideIntoSublists(list1, quantiles)
            index = 1
            for sublist in sublists:
                print "Quantile Number "+str(index)+":"
                index = index+1
                compareLists(sublist, list2)
    print "Finished."      