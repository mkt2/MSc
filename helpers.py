#coding:utf8
import numpy as np
import matplotlib.pyplot as plt
import dbg
import collections
from math import log
import os.path
import Graph, Bloom

def returnTime(timeInSeconds):
    s = timeInSeconds%60
    m = timeInSeconds/60
    return "Time in seconds: "+str(timeInSeconds)+"\n" \
    + str(m)+" minutes and "+str(s)+" seconds\n"
    assert m*60+s==timeInSeconds

def printResultsToFile(BF,G,outFolder="defaultOutFolder",timeInSeconds=-1,nextLine=-1):
    print "printResultsToFile(BF, G, outFolder="+str(outFolder)+", timeInSeconds="+str(timeInSeconds)+")"
    assert isinstance(outFolder, str)
    assert isinstance(timeInSeconds, int)
    f = open('Output/'+outFolder+"/G.txt", 'w')
    b = open('Output/'+outFolder+'/other_info.txt', 'w')
    G.printToFile(f)
    f.close()
    b.write(str(BF))
    b.write(BF.bitarray_str())
    if nextLine!=-1:
        b.write("\nWe were about to add the segment in line nr "+str(nextLine)+" from the input file\n")
    if timeInSeconds!=-1:
        b.write("\n"+str(returnTime(int(timeInSeconds))))
    b.close()


def printAllInfoFromFiles(fn,k,BF):
    print "helpers.printAllInfoFromFiles(fn,k,BF)"
    kd = collections.defaultdict(int)       #Stores all kmers in the files
    kd_BF = collections.defaultdict(int)    #Stores all kmers in the files
                                            #seen for the second time according to BF
    numberOfKmers = 0
    print "Running through the files:"
    for f in fn:
        h = open(f, "rU")
        for lineNr,line in enumerate(h):
            if (lineNr%4 == 1):
                segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
                for s in segments:
                    for km in dbg.kmers(s,k):
                        numberOfKmers += 1
                        rep_km = min(km,dbg.twin(km))
                        #kd[rep_km] += 1
                        if not rep_km in BF:
                            BF.add(rep_km)
                        else:
                            kd_BF[rep_km] += 1
    print "Done running through the files"
    G_twice = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
    numberOfAtLeastTwice = 0
    for km,num in kd.iteritems():
        if num > 1:
            G_twice.addSegmentToGraph(km)
            numberOfAtLeastTwice += 1
    G_BF = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
    for km in kd_BF:
        G_BF.addSegmentToGraph(km)
    
    print "\n"
    print "Total number of kmers in the files (with duplicates): %9i" % numberOfKmers
    print "Number of unique kmers in the files:                  %9i" % len(kd)
    print "Number of kmers in the BF:                            %9i" % len(BF)
    print "Number of kmers occuring at least twice in the files: %9i" % numberOfAtLeastTwice
    print "Same number according to BF:                          %9i" % len(kd_BF)
    print "Number of kmers in G_twice:                           %9i" % len(G_twice.kmers)
    print "Number of contigs in G_twice:                         %9i" % len(G_twice)
    print "Number of kmers in G_BF:                              %9i" % len(G_BF.kmers)
    print "Number of contigs in G_BF:                            %9i" % len(G_BF)
    print BF

def createKmerDictFromSegmentList(SL,k):
    kd = collections.defaultdict(int)
    for segment in SL:
        for km in dbg.kmers(segment,k):
            kd[km] += 1
        for km in dbg.kmers(dbg.twin(segment),k):
            kd[km] += 1
    return kd

def createNaiveFromKmerDict(kd,k):
    G,cs = dbg.all_contigs(kd,k)
    G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
    dbg.createGraphObject(G,cs,k,G_naive,pfn=False,ps=False)
    return G_naive

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

if __name__ == "__main__":
    k = 31
    fn = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
    fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
    BF = Bloom.Bloom(0.01,8000000,pfn=True)
    printAllInfoFromFiles(fn,k,BF)