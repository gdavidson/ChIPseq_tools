#coding:utf-8
#!/usr/bin/python
'''
Created on 4 mai 2015

@author: gdavidson
'''
# Allows different operations on fasta files.
# Requires:
# -biopython: http://biopython.org
# Arguments:
# -i <file.fasta>, sequence(s) file in fasta format
# Options (choose one):
# -r <regexp>,  looks for the regular expression within input sequences.
# -u true, generates a new multi fasta file with new uniques IDs for each sequence.
# -f <ids.txt>, file containing a list of sequences IDs (1 per line) to retrieve from '-i' file.
# Optional:
# -p <separator>, only if '-f', parses the IDs of both files (ex: sequenceID(file.fasta): MITF_peak_201, ID(ids.txt): CT_HG19_MITF_peak_201, then use either 'MITF' or 'peak' as separator).
#
#Examples:
#Retrieves sequences from a list of IDs:
#python getFromFasta.py -i chip_mitf_peaks.fasta -f peak_ids.txt 
#Generates new IDs for each sequences:
#python getFromFasta.py -i unnamed_sequences.fasta -u true
#Looks for TATATA or TATATATA in the input sequences:
#python getFromFasta.py -i chip_mitf_peaks.fasta -r 'TA{3,4}'

import sys
import getopt
import Bio.SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

progHelp = "Allows different operations on fasta files.\nArguments:\n -i <file>, sequence(s) file in fasta format\nOptions (choose one):\n -r <regexp>,  looks for the input regular expression within input sequences.\n -u true, generates a new multi fasta file with new uniques IDs for each sequence.\n -f <ids.txt>, file containing a list of sequences IDs (1 per line) to retrieve from '-i' file.\nOptional:\n -p <separator>, only if '-f', parses the IDs of both files (ex: sequenceID(file.fasta): MITF_peak_201, ID(ids.txt): CT_HG19_MITF_peak_201, then use either 'MITF' or 'peak' as separator).\n\nExamples:\nRetrieves sequences from a list of IDs:\npython getFromFasta.py -i chip_mitf_peaks.fasta -f peak_ids.txt\nGenerates new IDs for each sequences:\npython getFromFasta.py -i unnamed_sequences.fasta -u true\nLooks for TATATA or TATATATA in the input sequences:\npython getFromFasta.py -i chip_mitf_peaks.fasta -r 'TA{3,4}'"

def get_params(argv):
    try:
        opts, args = getopt.getopt(argv, "i:r:u:f:p:", ["infile", "regexp","uniqueIDs","idFile","parseIDs"])
    except getopt.GetoptError:
        sys.exit("Invalid argument:\n"+progHelp)
    infilename = "none"
    queryType = "none"
    pattern = "none"
    idFile = "none"
    parseIDs = False
    expName = "none"
    for opt,arg in opts:
        if opt =='-i':
            infilename = arg
        if opt =="-r":
            pattern = str(arg)
            queryType = "regexp"
        if opt == "-u":
            queryType = "uniqueIDs"
        if opt == "-f":
            idFile = arg
            queryType = "getSequences"
        if opt == "-p":
            parseIDs = True
            expName = arg      
    return infilename, pattern, queryType, idFile, parseIDs, expName

def getRecordDic(infileName):
    recordDic = Bio.SeqIO.index(infileName, "fasta")
    print "Reading Input file '"+infileName+"': "+str(len(recordDic))+" sequence(s)."
    return recordDic

def searchRegexp(pattern, recordDic):
    prog = re.compile(pattern)
    idList = []
    for key in recordDic.keys():
        record = recordDic[key]
        seq = str(record.seq)
        rID = record.id
        result = prog.findall(seq)
        if len(result) > 0:
            print "Found "+str(len(result))+" matches in '"+str(rID)+"':"
            for match in result:
                print match
            idList.append(rID)
    print "Found "+str(len(idList))+" sequences with match(es)."
    return idList        
                    
def getUniqueRecords(recordList):
    newRecordList = []
    idCount = 0
    for record in recordList:
        newID = "sequence_"+str(idCount)
        idCount = idCount+1
        newRecord = SeqRecord(seq=record.seq, id=newID, name=record.name, description=record.description)
        newRecordList.append(newRecord)
    return newRecordList

def getRecordList(infileName):
    records = list(Bio.SeqIO.parse(infileName, "fasta"))
    print "Reading Input file '"+str(infileName)+"': "+str(len(records))+" sequence(s)."
    return records

def parsePeakNumber(peakID, expName):
    peakNum = peakID.split(expName+"_")[1]
    return peakNum

def getIDList(idFile):
    iFile = open(idFile)
    idList = []
    for line in iFile:
        seqID = str(line).split("\t")[0]
        idList.append(seqID)
    iFile.flush()
    iFile.close()
    print "Reading sequences ID '"+str(idFile)+"': "+str(len(idList))+" sequences to retrieve."
    return idList
    
def writeRecords(recordList, outName = "output.fasta"):
    print "Writing Output File:'"+outName+"'."
    Bio.SeqIO.write(recordList, outName, "fasta")

def getRecordListFromIDs(recordList, idList, parseIDs = False, expName = "peak"):
    records = []
    displayID = True
    for record in recordList:       
        if displayID == True:
            print "Sequence file ID is like: "+str(record.id)
            print "List ID is like: "+str(idList[0])
            if parseIDs == True:
                print "Parsing IDs using '"+expName+"_' as separator."
            displayID = False
        for seqID in idList:
            if parseIDs == True:
                id1 = parsePeakNumber(record.id, expName)
                id2 = parsePeakNumber(seqID, expName)
            else:
                id1 = record.id
                id2 = seqID
            if id1 == id2:
                records.append(record)
    print "Found "+str(len(records))+" sequences matching IDs."
    return records

def writeList(inList, outName):
    outFile = open(outName, "w")
    for elem in inList:
        outFile.write(str(elem)+"\n")
    print "Writing file '"+str(outName)+"': "+str(len(inList))+" lines."
    outFile.flush()
    outFile.close()
                       
if __name__ == '__main__':
    infileName, pattern, queryType, idFile, parseIDs, expName = get_params(sys.argv[1:])
    if infileName == "none" or queryType == "none":
        sys.exit(progHelp)  
    if queryType == "regexp":
        recordDic = getRecordDic(infileName)
        idList = searchRegexp(pattern, recordDic)
        writeList(idList, "sequences_with_pattern.txt")
    if queryType == "uniqueIDs":
        oldRecordList = getRecordList(infileName)
        recordList = getUniqueRecords(oldRecordList)
        writeRecords(recordList)
    if queryType == "getSequences":
        recordList = getRecordList(infileName)
        idList = getIDList(idFile)
        records = getRecordListFromIDs(recordList, idList, parseIDs, expName)
        writeRecords(records, "retrieved_sequences.fasta")      
    print "Finished."     