#coding:utf8
import numpy as np
import matplotlib.pyplot as plt
import dbg
import collections
from math import log
import os.path
import Graph, Bloom
import time
import sys

def returnTime(timeInSeconds):
    timeInSeconds = int(timeInSeconds)
    if timeInSeconds<60:
        return "Time in seconds: "+str(timeInSeconds)+"\n"
    else:
        s = timeInSeconds%60
        m = timeInSeconds/60
        assert m*60+s==timeInSeconds
        return str(m)+" minutes and "+str(s)+" seconds\n"
        
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

def readKmersFromFileToDict(kmerFile):
    kmerDict = collections.defaultdict(int)
    f = open(kmerFile, 'rU')
    for km in f:
        km = km.rstrip('\n')
        kmerDict[km] = 1
    f.close()
    return kmerDict

def readGenomeInfoFromFile(infoFile):
    def nextValue(f):
        return f.readline().strip().split(",")[-1]
    f = open(infoFile, 'rU')
    genomeName = nextValue(f)
    k = nextValue(f)
    numKmersInGenome = nextValue(f)
    numReadsPerFile = nextValue(f)
    numKmersPerRead = nextValue(f)
    numKmersInReads = nextValue(f)
    numKmersInReads_twice = nextValue(f)
    numKmersInReads_twice_BF = nextValue(f)
    BF_ratio_1_vs_0 = nextValue(f)
    return genomeName,k,numKmersInGenome,numReadsPerFile,numKmersPerRead, \
    numKmersInReads,numKmersInReads_twice,numKmersInReads_twice_BF,BF_ratio_1_vs_0

#Prints the run times to a file in a format which can be copied 
#directly into a latex tabular environment
def printRuntimesToFile(genomeName, maxCovs, runTimes):
    print "printRunTimesToFile()"
    #modify runTimes so it becomes a string with the format:
        #x min y sec    <--- if time>1 min
        #y sec          <----if time<=1 min
    def modRunTimes(runTimes):
        runTimeStrings = [""]*len(runTimes)
        for i, rt in enumerate(runTimes):
            if rt > 60:
                rt_s = int(rt)%60
                rt_m = int(rt)/60
                runTimeStrings[i] = str(rt_m) + " min, " + str(rt_s) + " sec"
            else:
                runTimeStrings[i] = str(rt) + " sec"
        return runTimeStrings

    runTimes = modRunTimes(runTimes)
    timeFile = "Output/runTimes_t.txt"
    if genomeName=="t":
        tf = open(timeFile, 'w')
        for i in range(0,len(maxCovs)):
            tf.write(str(maxCovs[i])+" & "+runTimes[i]+" & \\\\\n")
        tf.close
    elif genomeName=="sa":
        tf_old = open(timeFile, 'r')
        tf = open("Output/runTimes_tAndSA.txt", 'w')
        i=0
        for line in tf_old:
            line=line[0:-3]
            lineCov = int(line.strip().split(" & ")[0])
            if lineCov==maxCovs[i]:
                tf.write(str(line)+runTimes[i]+" \\\\\n")
            i+=1
        for j in range(i-1,len(maxCovs)):
            tf.write(str(maxCovs[j])+" & - & " + runTimes[i] + " \\\\\n")
    else:
        raise Exception("Illegal value for genomeName="+str(genomeName)+". It must be either t or sa")

#Input:
#   kmersInGenome:      A dict storing all kmers actually occurring in the genome (stores rep_km)
#   numKmersInGenome:   The number of k-mers in the genome
#   Object:             Can either be a graph or a bloom filter
#   isGraph:            Tells us whether object is a graph or a bloom filter (can't have both)
#       isGraph==1: Object is a Graph
#       isGraph==0: Object is a BF
#Returns:
#   The fraction of kmers from the genome that occur in Object    
def fractionFromGenomeInObject(kmersInGenome,numKmersInGenome,Object,isGraph=1):
    assert(isGraph==0 or isGraph==1), "We must either be working with a graph or a bloom filter"
    count = 0
    for km in kmersInGenome:
        rep_km = min(km,dbg.twin(km))
        if isGraph==1:
            #Check whether rep_km occurs in the Graphs k-merdict
            if rep_km in Object.kmers:
                count += 1
        if isGraph==0:
            #Check whether rep_km occurs in the bloom filter
            if rep_km in Object:
                count += 1
    fraction = float(count) / numKmersInGenome
    assert (fraction>=0) and (fraction<=1), "A fraction must be between 0 and 1"
    return fraction

#Skilar öllum k-merum sem koma fyrir í kd1 en ekki
#í kd2. I.e. skilar mis"mengi" dictanna tveggja
def difference(kd1,kd2):
    return { x : kd1[x] for x in set(kd1) - set(kd2) }

def printKmerdictToFile(kd,outFile):
    f = open(outFile, 'w')
    for km in kd:
        f.write(km+"\n")
    f.close()

#Tekur lista af maxCovs og prentar mismengi af kd_-1 og kd_maxCov
#fyrir hvert maxCov
#maxCov = [-1,mc1,mc2,...]
#maxCovFiles = [outDirectory+"/kmersFromG_maxCov_"+str(maxCov)+".txt"...]
def diff_GNoMaxCov_and_GMaxCov(maxCovs,outDirectory):
    kd1 = readKmersFromFileToDict(outDirectory+"/kmersFromG_maxCov_-1.txt")
    print len(kd1)
    for maxCov in maxCovs[1:]:
        kd2 = readKmersFromFileToDict(outDirectory+"/kmersFromG_maxCov_"+str(maxCov)+".txt")
        print len(kd2)
        printKmerdictToFile(difference(kd1,kd2),outDirectory+"/diff_GNoMaxCov_and_GMaxCov_"+str(maxCov)+".txt")

#Before: G is a graph and kd is a subset of the k-mers from G
#After:  Return True if kd form a bubble <----- hvað á ég við? Er nóg að þeir séu t.d. hluti af bubble
#        False otherwise
def isBubble(G,kd):
    return

#Before: G is a graph and kd is a subset of the k-mers from G
#After:  Return True if kd form a tip
#        False otherwise
def isTip():
    return

#Before: G is a graph and kd is a subset of the k-mers from G
#After:  Return True if kd are isolated. I.e. not connected to the rest of G
#        False otherwise
def isIsolated():
    return

#Gera einnig fall sem finnur öll bubble, tip og isolated.
#Spurning hvort ég láti duga að skila hver mörg prósent allra 
#k-mera í G eru bubble etc
    
if __name__ == "__main__":
    #kmersInGenome = readKmersFromFileToDict("Output/t/kmers_genome.txt")
    #print readGenomeInfoFromFile("Output/t/genome_info.csv")
    kd1 = collections.defaultdict(int)
    kd2 = collections.defaultdict(int)
    kd1["A"] = "a"
    kd1["B"] = "b"
    kd1["C"] = "c"
    kd2["D"] = "d"
    kd2["A"] = "A"
    #kd1 = A B C
    #kd2 = A D
    #kd = B C
    kd = difference(kd1,kd2)
    assert(("B" in kd) and ("C" in kd) and not ("A" in kd) and not ("D" in kd))
    print kd
    #printKmerdictToFile(kd,"thisIsATestFile.txt")
    diff_GNoMaxCov_and_GMaxCov([-1, 5, 10, 15, 20, 30],"Output/t")

