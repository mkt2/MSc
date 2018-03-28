#coding:utf8
import collections, sys
import dbg, Bloom, Graph, helpers
import os.path
from shutil import copyfile
import fileinput
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log, exp, pow
import time
import pathlib2
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

alphabet = ["A","C","G","T","N"]

def printBFStatus(lineNr,f,BF):
    print "We are reading the segment from line "+str(lineNr)+" from the file "+str(f)
    ratio = BF.computeRatio()[0]
    B = BF.hasAcceptableRatio(ratio)
    if B:
        print "The BF is not full. Ratio="+str(ratio)
    else:
        print "The BF is full. Ratio="+str(ratio)

def newDictValues(COV,G,BF,maxAddCov,genome_dict,genome_len):
    G_len = G.numKmerPairs()
    G_frac = helpers.fractionFromGenomeInObject(genome_dict,genome_len,Object=G,isGraph=1)
    if maxAddCov==float('inf'):
        BF_len = len(BF)
        BF_frac = helpers.fractionFromGenomeInObject(genome_dict,genome_len,Object=BF,isGraph=0)
        return [G_len, G_frac, BF_len, BF_frac]
    else:
        return [G_len, G_frac]

def updateDict(COV,covDict,maxAddCov,values,i):
    #values = [G_len,G_frac,BF_len,BF_frac]
    #covDict[maxAddCov][0][i] = COV
    covDict[maxAddCov][0].append(COV)
    assert(len(covDict[maxAddCov])>1)
    if maxAddCov==float('inf'):
        assert(len(values)==4)
        for j in [1,2,3,4]:
            #covDict[maxAddCov][j][i] = values[j-1]
            covDict[maxAddCov][j].append(values[j-1])
    else:
        assert(len(values)==2)
        for j in [1,2]:
            #covDict[maxAddCov][j][i] = values[j-1]
            covDict[maxAddCov][j].append(values[j-1])

def naive_createDict(fn,k):
    #Creates a dict, kd, storing every k-mer from the reads
    kd = collections.defaultdict(int)
    for s in helpers.segments(fn,k):
        for km in dbg.kmers(s,k):
            kd[km] += 1
            kd[dbg.twin(km)] += 1
    return kd

def naive_createDict_nonSingleton(fn,k):
    #Creates a dict, kd, storing every non singleton k-mer from the reads
    #by filtering out singletons from kd
    kd = naive_createDict(fn,k)
    d1 = [x for x in kd if kd[x] <= 1]
    for x in d1:
        del kd[x]
    return kd

def naive_createDict_nonSingleton_BF(fn,k,BF):
    #Filters out singletons from the reads using a BF, before adding them to kd
    #Creates a dict, kd, storing every non singleton (some will be False Positives from the BF)
    kd = collections.defaultdict(int)
    for s in helpers.segments(fn,k):
        for km in dbg.kmers(s,k):
            twin_km = dbg.twin(km)
            rep_km = min(km,twin_km)
            if not (rep_km in BF):
                BF.add(rep_km)
            else:
                kd[km] += 1
                kd[twin_km] += 1
    return kd

def naive_GFA(fn,k,fileName):
    #Creates a dict, kd, storing every k-mer from the reads
    #Adds the k-mers from kd to the DBG, one k-mer at a time      
    kd = naive_createDict(fn,k)
    G,cs = dbg.all_contigs(kd,k)
    dbg.print_GFA_to_file(G,cs,k,fileName)

def naive_GFA_nonSingleton(fn,k,fileName):
    #print "naive_GFA_nonSingleton(fn,k,fileName)"
    #Creates a dict, kd, storing every k-mer from the reads
    #Filters out singletons from kd
    #Adds the k-mers, from kd, to the DBG, one k-mer at a time
    kd = naive_createDict_nonSingleton(fn,k)
    G,cs = dbg.all_contigs(kd,k)
    dbg.print_GFA_to_file(G,cs,k,fileName)

def naive_GFA_nonSingleton_BF(fn,k,BF,fileName,pfn=True,printProgress=False):
    #Creates kd with no singletons. Filters using a BF
    #Adds the k-mers, from kd, to the DBG, one k-mer at a time
    kd = naive_createDict_nonSingleton_BF(fn,k,BF)
    G,cs = dbg.all_contigs(kd,k)
    dbg.print_GFA_to_file(G,cs,k,fileName)

def BFAdder(fn,k,BF,G,maxAddCov,covDict,read_len,genome_len,sampleAtReadsNumber,pfn=False,printProgress=False):
    if pfn:
        print "BFAdder",locals()
    if not isinstance(fn, list):
        print "fn:",fn
        raise Exception('fn has to be a list')
    MSC = float('inf')    #Max Split Cov
    alpha = 0.01
    covFactor = (float(read_len) / genome_len) * pow(1-alpha,k) #Til að flýta fyrir reikningnum á COV
    COV = 0		#Hversu oft við höfum séð hvern k-mer að meðaltali
    measureIndex = 0
    for segmentCount, s in enumerate(helpers.segments(fn,k)):
        print segmentCount
        if segmentCount==500:
            break
        assert(G.assertLegal==False)
        if segmentCount in sampleAtReadsNumber:
            if printProgress:
                print "Taking sample nr. "+str(measureIndex)+" out of "+str(len(sampleAtReadsNumber)-1)+" at segment nr. "+str(segmentCount)
            values = newDictValues(COV,G,BF,maxAddCov,genome_dict,genome_len)
            updateDict(COV,covDict,maxAddCov,values,measureIndex)
            measureIndex+=1

        start = 0
        for i, km in enumerate(dbg.kmers(s,k)):
            rep_km = min(km,dbg.twin(km))
            if not (rep_km in BF):
                if COV<=maxAddCov:
                    BF.add(rep_km)
                if i>start:
                    goodSequence = s[start:i+k-1]
                    G.addSegmentToGraph(goodSequence,CS=(COV>MSC))
                start = i+1

        #If we reach the end of s we add the current goodSequence
        if len(s)-start>=k:
            G.addSegmentToGraph(s[start:],CS=(COV>MSC))
        COV = segmentCount * covFactor

    #Final measurement when finished running through the last file:
    if printProgress:
        print "measureIndex="+str(measureIndex)+" after finishing the loop"
        print "Taking sample nr. "+str(measureIndex)+" out of "+str(len(sampleAtReadsNumber)-1)+" at segment nr. "+str(segmentCount)
    values = newDictValues(COV,G,BF,maxAddCov,genome_dict,genome_len)
    updateDict(COV,covDict,maxAddCov,values,measureIndex)
    assert(G.carefulSplit==False), "Á meðan ég er ekki byrjaður að kalla á BFAdder með maxSplitCov ætti G.carefulSplit aldrei að verða True"

def selectGenome(genomeName):
	if genomeName=="t":
		maxAddCovs = [5, 10, 15, 20, 30,float('inf')]
		fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
		numberOfKmers = 8000000
	elif genomeName=="sa":
		maxAddCovs = [5, 10, 15, 20, 30,float('inf')]
		fn = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
		numberOfKmers = 100000000
	else:
		raise Exception("The genomeName must be either 't' or 'sa'!")
	return maxAddCovs, fn, numberOfKmers

#Need to run preprocessInput.py before running this file
if __name__ == "__main__":
    genomeName = sys.argv[1]
    maxAddCovs, fn, numberOfKmers = selectGenome(genomeName)
    #maxAddCovs = maxAddCovs[-1:]
    runTimes = [-1]*len(maxAddCovs)
    pfn = False
    printProgress = True

    #Lesum úr preprocess skránum:
    p = 0.01
    outDir = "Output/"+genomeName
    assert(pathlib2.Path(outDir+"/genome_info.csv").is_file()), "We have to create this file using preprocessInput.py before running BF_counter.py"
    assert(pathlib2.Path(outDir+"/kmers_genome.txt").is_file()), "We have to create this file using preprocessInput.py before running BF_counter.py"
    assert(pathlib2.Path(outDir+"/kmers_reads.txt").is_file()), "We have to create this file using preprocessInput.py before running BF_counter.py"
    preprocFiles = [outDir+"/genome_info.csv",outDir+"/kmers_genome.txt",outDir+"/kmers_reads.txt"]
    genomeInfo = helpers.readGenomeInfoFromFile(preprocFiles[0])
    k = int(genomeInfo[1])
    genomeLen = int(genomeInfo[2])   #Number of k-mers in the genome
    numReadsPerFile = int(genomeInfo[3])
    numKmersPerRead = int(genomeInfo[4])
    genome_dict = helpers.readKmersFromFileToDict(preprocFiles[1])  #A dictionary storing all k-mers from the genome

    #Select where we want to sample:
    NoM = 10    #numberOfMeasurements
    totNumberOfReads = numReadsPerFile*2
    sampleEvery = int(totNumberOfReads / (NoM-1))
    sampleAtReadsNumber = range(0,totNumberOfReads-sampleEvery,sampleEvery)
    if (totNumberOfReads%(NoM-1))!=0:
        print 'adding an extra value'
        sampleAtReadsNumber.append(totNumberOfReads+1)
    print "totNumberOfReads", totNumberOfReads
    print sampleAtReadsNumber, len(sampleAtReadsNumber)

    #covDict geymir upplýsingar fyrir myndirnar fyrir sérhvert maxAddCov
    #   Geymir COV, lenG, fracG, lenBF og fracBF fyrir G
    #   Geymir COV, lenG og fracG fyrir sérhvert Gx
    #   ATH: Þetta eru allt listar með NoM gildum hver
    t = len(sampleAtReadsNumber)    #annaðhvort NoM eða NoM+1
    print "t="+str(t)+", NoM="+str(NoM)
    #t = NoM+1   #XX var að bæta við +1 til að prófa fyrir SA
    covDict = collections.defaultdict(list)
    for maxAddCov in maxAddCovs:
        if maxAddCov==float('inf'):
            #covDict[maxAddCov] = [[0]*t,[0]*t,[0]*t,[0]*t,[0]*t]
            covDict[maxAddCov] = [[],[],[],[],[]]
        else:
            #covDict[maxAddCov] = [[0]*t,[0]*t,[0]*t]
            covDict[maxAddCov] = [[],[],[]]
    
    start = time.time()
    #Keyrum BF_counter fyrir sérhvert maxAddCov og tökum tímana
    start = time.time()
    for i, maxAddCov in enumerate(maxAddCovs):
        start_i = time.time()
        print "Starting on maxAddCov="+str(maxAddCov)
        print "Initializing an empty bloom filter BF and Graph G"+str(maxAddCov)
        BF = Bloom.Bloom(p,numberOfKmers,pfn=False)
        theGraph = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=False)
        assert(theGraph.isEmpty())

        print "Running BFAdder:"
        BFAdder(fn,k,BF,theGraph,maxAddCov,covDict,numKmersPerRead,genomeLen,sampleAtReadsNumber,pfn,printProgress)
        timeInSeconds = int(time.time()-start_i)
        runTimes[i] = timeInSeconds
        print "Done creating G"+str(maxAddCov)
        break   #Tímabundið bara að keyra BFAdder fyrir eitt maxCov og ekki save-a niðurstöður eða neitt
"""
        print "Saving G"+str(maxAddCov)+" to G"+str(maxAddCov)+".txt"
        theGraph.printToFile(outDir+"/G"+str(maxAddCov)+".txt")

        print "Finished creating G"+str(maxAddCov)+". Time: "+helpers.returnTime(timeInSeconds)
        print "-------------------------------------------------------\n"
	
    helpers.printRuntimesToFile(genomeName, maxAddCovs, runTimes)
    print "About to print covDict to fileName="+outDir+"/covDict.txt"
    helpers.printCovDictToFile(covDict,fileName=outDir+"/covDict.txt")
    print "Finished running BFAdder.py on genome="+str(genomeName)
    print "Total time: "+helpers.returnTime(int(time.time()-start))"""

