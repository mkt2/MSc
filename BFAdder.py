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

def updateDict(COV,covDict,MAC,MSC,values):
    #values = [G_len,G_frac,BF_len,BF_frac]
    #covDict[maxAddCov][0][i] = COV
    key = (MAC,MSC)
    covDict[key][0].append(COV)
    assert(len(covDict[key])>1)
    if MAC==float('inf'):
        assert(len(values)==4)
        for j in [1,2,3,4]:
            covDict[key][j].append(values[j-1])
    else:
        assert(len(values)==2)
        for j in [1,2]:
            covDict[key][j].append(values[j-1])

def printCurrentKeyToFile(covDict,key,fileName):
    f = open(fileName, 'a')
    #for key, values in covDict.iteritems():
    values = covDict[key]
    MAC = key[0]
    MSC = key[1]
    #values = [COV, G_len,G_frac,BF_len,BF_frac] fyrir G (maxAddCov = float('inf'))
    #values = [COV, G_len,G_frac] fyrir sérhvert Gx
    if MAC == float('inf'):
        assert(len(values)==5)
        f.write(str(key)+";"+str(values[0])+";"+str(values[1])+";"+str(values[2])+";"+str(values[3])+";"+str(values[4])+"\n")
    else:
        assert(len(values)==3)
        f.write(str(key)+";"+str(values[0])+";"+str(values[1])+";"+str(values[2])+"\n")
    f.close()


def printKeyToFile(key,COV,values,fileName):
    f = open(fileName, 'a')
    #for key, values in covDict.iteritems():
    MAC = key[0]
    MSC = key[1]
    #values = [G_len,G_frac,BF_len,BF_frac] fyrir G (maxAddCov = float('inf'))
    #values = [G_len,G_frac] fyrir sérhvert Gx
    if MAC == float('inf'):
        assert(len(values)==4)
        printString = str(key)+";"+str(COV)+";"+str(values[0])+";"+str(values[1])+";"+str(values[2])+";"+str(values[3])+"\n"
    else:
        assert(len(values)==2), "len(values)="+str(len(values))
        printString = str(key)+";"+str(COV)+";"+str(values[0])+";"+str(values[1])+"\n"
    f.write(printString)
    print printString
    f.close()

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

def BFAdder(fn,k,BF,G,maxAddCov,maxSplitCov,read_len,genome_len,sampleAtReadsNumber,pfn=False,printProgress=False):
    #assert(maxSplitCov == float('inf'))
    if pfn:
        print "BFAdder",locals()
    if not isinstance(fn, list):
        print "fn:",fn
        raise Exception('fn has to be a list')
    alpha = 0.01
    covFactor = (float(read_len) / genome_len) * pow(1-alpha,k) #Til að flýta fyrir reikningnum á COV
    COV = 0		#Hversu oft við höfum séð hvern k-mer að meðaltali
    measureIndex = 0
    #for segmentCount, s in enumerate(helpers.segments(fn,k)):
    for s, segmentCount in helpers.segments(fn,k):
        #print segmentCount
        #if segmentCount==10000:
        #    break
        #assert(G.assertLegal==False)
        if segmentCount in sampleAtReadsNumber:
            if printProgress:
                print "Taking sample nr. "+str(measureIndex)+" out of "+str(len(sampleAtReadsNumber)-1)+" at segment nr. "+str(segmentCount)
            values = newDictValues(COV,G,BF,maxAddCov,genome_dict,genome_len)
            yield COV, values
            measureIndex+=1

        start = 0
        for i, km in enumerate(dbg.kmers(s,k)):
            rep_km = min(km,dbg.twin(km))
            if not (rep_km in BF):
                if COV<=maxAddCov:
                    BF.add(rep_km)
                if i>start:
                    goodSequence = s[start:i+k-1]
                    G.addSegmentToGraph(goodSequence,CS=(COV>maxSplitCov))
                    #assert(G.isLegalDBG())
                start = i+1

        #If we reach the end of s we add the current goodSequence
        if len(s)-start>=k:
            G.addSegmentToGraph(s[start:],CS=(COV>maxSplitCov))
            #assert(G.isLegalDBG())
        COV = segmentCount * covFactor

    print 'segmentCount='+str(segmentCount)+" after the loop"
    print 'read_len='+str(read_len)
    print 'read_len*2='+str(read_len*2)
    #Final measurement when finished running through the last file:
    if printProgress:
        print "measureIndex="+str(measureIndex)+" after finishing the loop"
        print "Taking sample nr. "+str(measureIndex)+" out of "+str(len(sampleAtReadsNumber)-1)+" at segment nr. "+str(segmentCount)
    values = newDictValues(COV,G,BF,maxAddCov,genome_dict,genome_len)
    yield COV, values

def selectGenome(genomeName):
    I = float('inf')
    if genomeName=="t":
        #Gömul gildi:
        '''
        filters = [ \
        (float('inf'),2),(10,11),(20,10), \
        (30,10),(float('inf'),5),(float('inf'),10),(float('inf'),20), #(float('inf'),2)\ 
        (5,float('inf')), \
        (10,12),(10,15),(10,float('inf')), #fæ error í merge fyrir (10,11)\
		(15,16),(15,20),(15,float('inf')), \
		(20,float('inf')), \
		(30,float('inf')), \
        (15,10),(20,15),(30,20), \
		(float('inf'),float('inf'))]
        '''
        #Ný gildi:
        filters = [ \
        (I,2),(I,3), \
        (10,2),(10,5),(10,11),(10,I), \
        (15,2),(15,5),(15,10),(15,I), \
        (20,5),(20,I), \
        (I,I), \
        ]
        fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
        numberOfKmers = 8000000
    elif genomeName=="sa":
        #Gömul gildi:
        '''
        filters = [ \
        (15,16),(15,20),(15,float('inf')), \
		(20,21),(20,25),(20,float('inf')), \
		(30,float('inf')), \
        (float('inf'),float('inf')), \
        (15,10),(20,15),(30,20), \
		]
        '''
        #Ný gildi:
        filters = [ \
        (5,2),(5,5),(5,10),(5,I), \
        (I,5), \
        (10,5),(10,11),(10,I), \
        (15,5),(15,I), \
        (I,I) \
        ]
        #vil eiga en þarf ekki að keyra aftur svo lengi sem ég hendi ekki skránum
        #(15,10), (15,I)
        #(I,I)
        #Pæling að bæta við:
        #(20,5)
        fn = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
        numberOfKmers = 100000000
    else:
        raise Exception("The genomeName must be either 't' or 'sa'!")
    return filters, fn, numberOfKmers

#Need to run preprocessInput.py before running this file
if __name__ == "__main__":
    genomeName = sys.argv[1]
    filters, fn, numberOfKmers = selectGenome(genomeName)
    numGraphs = len(filters)
    runTimes = [-1]*numGraphs
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
        sampleAtReadsNumber.append(totNumberOfReads+1)
    print "totNumberOfReads:", totNumberOfReads
    print sampleAtReadsNumber, len(sampleAtReadsNumber)

    #covDict geymir upplýsingar fyrir myndirnar fyrir sérhvert maxAddCov
    #   Geymir COV, lenG, fracG, lenBF og fracBF fyrir G
    #   Geymir COV, lenG og fracG fyrir sérhvert Gx
    #   ATH: Þetta eru allt listar með NoM gildum hver
    covDict = collections.defaultdict(list)
    for key in filters:
        (MAC,MSC) = key
        if MAC==float('inf'):
            covDict[key] = [[],[],[],[],[]]
        else:
            covDict[key] = [[],[],[]]
    
    start = time.time()
    #Keyrum BF_counter fyrir sérhvert maxAddCov og tökum tímana
    start = time.time()
    fileName=outDir+"/covDict.txt"
    try:
        os.remove(fileName)
    except OSError:
        pass
    for i, key in enumerate(filters):
        (MAC,MSC) = key
        gn = helpers.createGraphName(MAC,MSC)
        print "Starting on gn="+gn+":"
        start_i = time.time()
        print "Initializing an empty bloom filter and graph"
        BF = Bloom.Bloom(p,numberOfKmers,pfn=False)
        theGraph = Graph.Graph(k,pfn=False,ps=False,al=False,pil=True,printInit=False)
        #theGraph = Graph.Graph(k,pfn=False,ps=False,al=True,pil=True,printInit=False)  #<----to test
        assert(theGraph.isEmpty())
        print "Running BFAdder:"
        for COV, dictValues in BFAdder(fn,k,BF,theGraph,MAC,MSC,numKmersPerRead,genomeLen,sampleAtReadsNumber,pfn,printProgress):
            updateDict(COV,covDict,MAC,MSC,dictValues)
            #printKeyToFile(key,COV,dictValues,fileName)
        printCurrentKeyToFile(covDict,key,fileName)
        timeInSeconds = int(time.time()-start_i)
        runTimes[i] = timeInSeconds
        #print "Done creating "+gn
        gn_f = gn+".txt"
        #print "Saving "+gn+" to "+gn_f
        theGraph.printToFile(outDir+"/"+gn_f)

        print "Finished creating "+gn+" and saving it to "+gn_f+". Time: "+helpers.returnTime(timeInSeconds)
        print "-------------------------------------------------------\n"
	
    #print "About to print covDict to fileName="+outDir+"/covDict.txt"
    #helpers.printCovDictToFile(covDict,fileName=outDir+"/covDict.txt")

    print "About to print runTimes to file"
    for i, t in enumerate(runTimes):
        print i,filters[i],t
    helpers.printRuntimesToFile(genomeName, filters, runTimes)
    print "Finished running BFAdder.py on genome="+str(genomeName)
    print "Total time: "+helpers.returnTime(int(time.time()-start))

