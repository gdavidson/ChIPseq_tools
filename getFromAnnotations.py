'''
Created on 5 juin 2015

@author: gdavidson
'''

# Allows operations on homer 'annotatePeaks.pl' output file.
# Requires:
# -pylab module from package matplotlib (http://matplotlib.org/)
# Arguments:
# -i <homer_annotation.txt>, homer output file.
# Options:
# -h <histogram title>, makes a histogram of distances to nearest TSS.
# -p <piechart title>, makes a chart summarizing annotations (number of peaks annotated as 'promoter-TSS', 'exon, 'intron' ect...).
# -r <keyword>, retrieves lines annotated as the keyword (<promoter-tss>, <exon>, <intron>, <TTS>, <intergenic>, <non-coding>, <3'UTR>, <5'UTR>)
# -b <biomart_export.txt>, completes annotation with a biomart output file containing the following fields (in that order): transcript ID, gene ID, gene name, description. Adds these three fields to the homer file.
# -l <peaks.txt>, retrieves lines corresponding to a list of peaks.
#
#Examples:
# Makes a chart and a histogram summarizing annotations:
#python getFromAnnotations.py -i chip_mitf_peaks_annotations.xls -h 'MITF distances to nearest TSS' -p 'MITF Annotations'
# Adds ENSEMBL gene IDs, gene names and descriptions:
#python getFromAnnotations.py -i chip_mitf_peaks_annotations.xls -b mart_export_hg19_release69.tsv
# Retrieves lines corresponding to TSS:
#python getFromAnnotations.py -i chip_mitf_peaks_annotations.xls -r 'promoter-tss'
# Retrieves lines corresponding to a list of peaks:
#python getFromAnnotations.py -i chip_mitf_peaks_annotations.xls -l mitf_sox_commonpeaks.txt
  
import sys
import getopt
import pylab as P
import matplotlib.pyplot as plt

progHelp = "Allows operations on homer 'annotatePeaks.pl' output file.\nRequires:\n -pylab module from package matplotlib (http://matplotlib.org/)\nArguments:\n -i <homer_annotation.txt>, homer output file.\nOptions:\n -h <histogram title>, makes a histogram of distances to nearest TSS.\n -p <piechart title>, makes a chart summarizing annotations (number of peaks annotated as 'promoter-TSS', 'exon, 'intron' ect...).\n -r <keyword>, retrieves lines annotated as the keyword (<promoter-tss>, <exon>, <intron>, <TTS>, <intergenic>, <non-coding>, <3'UTR>, <5'UTR>).\n -b <biomart_export.txt>, completes annotation with a biomart output file containing the following fields (in that order): transcript ID, gene ID, gene name, description. Adds these three fields to the homer file.\n-l <peaks.txt>, retrieves lines corresponding to a list of peaks.\n\nExamples:\nMakes a chart and a histogram summarizing annotations:\npython getFromAnnotations.py -i chip_mitf_peaks_annotations.xls -h 'MITF distances to nearest TSS' -p 'MITF Annotations'\nAdds ENSEMBL gene IDs, gene names and descriptions:\npython getFromAnnotations.py -i chip_mitf_peaks_annotations.xls -b mart_export_hg19_release69.tsv\nRetrieves lines corresponding to TSS:\npython getFromAnnotations.py -i chip_mitf_peaks_annotations.xls -r 'promoter-tss'\nRetrieves lines corresponding to a list of peaks:\npython getFromAnnotations.py -i chip_mitf_peaks_annotations.xls -l mitf_sox_commonpeaks.txt"

def get_params(argv):
    try:
        opts, args  = getopt.getopt(argv, "i:h:p:r:b:l:", ["infile"])
    except getopt.GetoptError:
        sys.exit("Invalid argument:\n"+progHelp)
    infilename = "none"
    queryType = "none"
    plotDist = False
    plotLoc = False
    chartTitle = "pieChart"
    histTitle = "histogram"
    typeToRetrieve = "none"
    biomartFile = "none"
    completeFile = False
    peakListFile = "none"
    for opt,arg in opts:
        if opt =='-i':
            infilename = arg
        if opt == '-h':
            if queryType == "none":
                queryType = "distances"
            plotDist = True
            histTitle = arg
        if opt == '-p':
            if queryType == "none":
                queryType = "locations"
            plotLoc = True
            chartTitle = arg
        if opt == '-r':
            typeToRetrieve = str(arg).lower().strip()
            queryType = "retrieveLines"
        if opt == "-b":
            biomartFile = arg
            completeFile = True
            if queryType == "none":
                queryType = "biomart" 
        if opt == "-l":
            peakListFile = arg
            if queryType == "none":
                queryType = "peakList"
                           
    return infilename, queryType, plotDist, plotLoc, chartTitle, histTitle, typeToRetrieve, biomartFile, completeFile, peakListFile

def getAnnotationList(annotationFile):
    aFile = open(annotationFile, "r")
    annList = []
    for line in aFile:
        annList.append(line)
    aFile.flush()
    aFile.close()
    print "Reading homer annotation file '"+str(annotationFile)+"': "+str(len(annList))+" lines."
    return annList

def getDistanceList(annotationList):
    distances = []
    countMap = {}
    for annotation in annotationList:
        if str(annotation).startswith("Peak"):
            continue
        distance = int(annotation.split("\t")[9])
        distances.append(distance)
        if distance in countMap:
            countMap[distance]=countMap[distance]+1
        else:
            countMap[distance] = 1
    return distances, countMap

def histDistances(distanceList, hisTitle):
    print "Displaying histogram 'distances from nearest TSS' ..."
    dbin = [-50000,-40000, -30000,-20000,-10000, 0, 10000, 20000, 30000, 40000, 50000]
    P.hist(distanceList, bins=dbin, color = 'blue')
    P.xlabel('Distance to TSS (pb)')
    P.ylabel('Count')
    P.xticks(dbin)
    P.title(hisTitle, bbox={'facecolor':'0.8', 'pad':5})
    P.show()
    
def getPieChartMap(annotationList):
    countMap = {"3'UTR": 0, "5'UTR": 0, "Exon":0, "Intergenic": 0, "Intron": 0, "promoter-TSS": 0, "TTS": 0, "non-coding": 0}
    for ann in annotationList:
        location = str(str(ann).split("\t")[7])
        if location.startswith("3'"):
            countMap["3'UTR"] = countMap["3'UTR"]+1
        elif location.startswith("5'"):
            countMap["5'UTR"] = countMap["5'UTR"]+1
        elif location.startswith("exon"):
            countMap["Exon"] = countMap["Exon"]+1
        elif location.startswith("Intergenic"):
            countMap["Intergenic"] = countMap["Intergenic"]+1
        elif location.startswith("intron"):
            countMap["Intron"] = countMap["Intron"]+1
        elif location.startswith("promoter"):
            countMap["promoter-TSS"] = countMap["promoter-TSS"]+1
        elif location.startswith("TTS"):
            countMap["TTS"] = countMap["TTS"]+1
        elif location.startswith("non-coding"):
            countMap["non-coding"] = countMap["non-coding"]+1     
    return countMap

def getSumFromMap(aMap):
    total = 0
    for key in aMap.keys():
        val = int(aMap[key])
        total = int(total+val)
    return total
    
def pieChart(pieChartMap, chartTitle):
    print "Displaying Annotation Chart ..."
    P.figure(1, figsize=(15,15))
    labels = []
    fracs = []
    total = getSumFromMap(pieChartMap)
    colors = ['yellowgreen', 'gold', 'lightskyblue', 'cyan', 'lightcoral', 'lightgreen', 'gray', 'pink']
    for key in pieChartMap.keys():
        count = pieChartMap[key]
        countPercentage = (float(count)/total)*100
        label = str(key)
        labels.append(label+"("+str(count)+" - "+str("{:10.2f}".format(countPercentage)).strip()+"%)")
        fracs.append(countPercentage)
    
    P.pie(fracs,colors=colors, autopct='%1.1f%%', shadow=False, startangle=90)
    P.title(chartTitle, bbox={'facecolor':'0.8', 'pad':5})
    P.legend(labels, bbox_to_anchor=(1, 1))
    P.show()

def getLinesFromType(annotationList, typeToRetrieve):
    rLines = []
    rLines.append(annotationList[0])
    for ann in annotationList:
        location = str(str(ann).split("\t")[7]).lower().strip()
        if location.startswith(typeToRetrieve):
            rLines.append(ann)
    print "Looking for lines with annotation '"+typeToRetrieve+"': "+str(len(rLines))+" lines found."
    return rLines

def getLinesFromPeaks(annotationList, peakList):
    rLines = []
    for ann in annotationList:
        peakID = str(ann).split("\t")[0].upper().strip()
        if peakID in peakList:
            rLines.append(ann)
    return rLines

def getEnsemblMap(biomartFile):
    bFile = open(biomartFile, "r")
    ensemblMap = {}
    for line in bFile:
        transID = str(line).split("\t")[0].strip().upper()
        geneID = str(line).split("\t")[1]
        geneName = str(line).split("\t")[2]
        geneDesc = str(line).split("\t")[3].strip()
        ensemblMap[transID] = geneID+"\t"+geneName+"\t"+geneDesc
    print "Reading biomart file '"+str(biomartFile)+"': "+str(len(ensemblMap))+" transcripts."
    bFile.flush()
    bFile.close()
    return ensemblMap
    
def writeList(inList, outName):
    outFile = open(outName, "w")
    for elem in inList:
        outFile.write(str(elem))
    print "Writing file '"+str(outName)+"': "+str(len(inList))+" lines."
    outFile.flush()
    outFile.close()
    
def fillLine(line, featureNumber):
    fieldsToAdd = int(featureNumber-len(line.split("\t")))
    filledLine = line
    while fieldsToAdd>0:
        filledLine = filledLine+"\tNA"
        fieldsToAdd = fieldsToAdd -1
    return filledLine

def completeAnnotations(annotationList, ensemblMap):
    detailedLines = []
    for line in annotationList:
        if line.startswith("Peak"):
            firstline = str(line).strip()+"\tGene ID\tGene Name\tGene Description\n"
            detailedLines.append(firstline)
            featureNumber = len(line.split("\t"))
            continue
        transID = str(line).split("\t")[10].upper().strip()
        addedAnnotations = "NA"
        if transID in ensemblMap:
            addedAnnotations = ensemblMap[transID]
        else:
            print "ID not found in biomart file: "+str(transID)+" !"
        newline = str(line).strip()
        newline = fillLine(newline, featureNumber)
        newline = newline+"\t"+addedAnnotations+"\n"
        detailedLines.append(newline)
    return detailedLines

def getIDList(fileName):
    iFile = open(fileName, "r")
    idList = []
    for line in iFile:
        if str(line).isspace():
            continue
        peakID = str(line).split("\t")[0].upper().strip()
        idList.append(peakID)
    iFile.flush()
    iFile.close()
    print "Reading peak list '"+str(fileName)+"': "+str(len(idList))+" IDs."
    return idList

def plotDistances(countMap):
    distances = []
    distCount = []
    keys = range(-250,250,1)
    for key in keys:
        if key in countMap:
            distances.append(key)
            distCount.append(countMap[key])
    print "Displaying TSS plot at base pair resolution..."
    plt.plot(distances,distCount, 'ro-')
    plt.xticks(range(-250,250,20))
    plt.grid(True)
    plt.xlabel("Distance to TSS (pb)")
    plt.ylabel("Peak Count")
    plt.show()
                 
if __name__ == '__main__':
    infilename, queryType, plotDist, plotLoc, chartTitle, histTitle, typeToRetrieve, biomartFile, completeFile, peakListFile = get_params(sys.argv[1:])
    if infilename == "none" or queryType == "none":
        sys.exit(progHelp)
    annotationList = getAnnotationList(infilename)
    if plotDist == True:
        distanceList,countMap = getDistanceList(annotationList)
        histDistances(distanceList, histTitle)
        plotDistances(countMap)
    if plotLoc == True:
        pieChartMap = getPieChartMap(annotationList)
        pieChart(pieChartMap, chartTitle)
    if queryType == "retrieveLines":
        retLines = getLinesFromType(annotationList, typeToRetrieve)
        writeList(retLines, "retrievedAnnotations_"+typeToRetrieve+".xls")
    if completeFile == True:
        ensemblMap = getEnsemblMap(biomartFile)
        detailedLines = completeAnnotations(annotationList, ensemblMap)
        writeList(detailedLines, str(infilename).split(".")[0]+"_detailed.xls")
    if not peakListFile == "none":
        peakList = getIDList(peakListFile)
        retLines = getLinesFromPeaks(annotationList, peakList)
        writeList(retLines, "foundLines.txt")
    print "Finished."