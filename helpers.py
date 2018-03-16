#coding:utf8
import dbg
import collections
import os.path
import Graph
import time
import sys

alphabet = ["A","C","G","T","N"]

#------------------------Helper functions------------------------
def splitString(s,i,k,ps=False):
    #Before:    s is a string representing a segment
    #               0    i    L
    #           s: [    |    ]
    #After:
    #   s has been split up into
    #     0     i               L
    #   s[0:i-1]  ->  s[i-k:L-1]
    if ps:
        print "splitString(s="+str(s)+", i="+str(i)+", k="+str(k)+")"
    L = len(s)
    #print "\ns: " + str(s) + ". i: " + str(i) + ". k: " + str(k) + ". L: " + str(L)
    assert(i>=k)    #i=k is the lowest i allowed
    #assert(i<=L-k)
    assert(i<=L-1)  #L-1 is the highest i allowed
    s0 = s[0:i]
    s1 = s[i-k+1:]
    if ps:
        print s0,s1
    return s0,s1

def reverseList(L):
    for i, (ID,B) in enumerate(L):
        L[i] = (ID,not B)

#Before:    c is a string representing a contig
#After:     Returns the number of k-mers in c
def cLen(c,k):
    assert(len(c)>=k), "c must have a minimum length k"
    num_kmers = len(list(dbg.kmers(c,k)))
    num_kmers_fast = len(c)-k+1
    assert(num_kmers==num_kmers_fast), "num_kmers_fast should return the number of k-mers."+str(num_kmers)+", "+str(num_kmers_fast)
    return num_kmers_fast

#Skilar öllum k-merum sem koma fyrir í kd1 en ekki í kd2 sem mengi
#I.e. skilar mismengi kd1 og kd2
def difference(kd1,kd2):
    return { x : kd1[x] for x in set(kd1) - set(kd2) }

def getIDsFromSetOfKmers(G,kmerSet):
    ids = set()
    for km in kmerSet:
        ids.add(G.kmers[km][0])
    return ids


#Returns a set containing all nodes (as tuples) from DCorA which neither occur in visited nor skipNodes
def DCorA_toAdd(DCorA,visited,skipNodes):
    toAdd = set()
    for (x_ID,x_B) in DCorA:
        if not ((x_ID in visited) or (x_ID in skipNodes)):
            toAdd.add((x_ID,x_B))
    return toAdd

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

def createNaiveFromKmerDict(kd,k):
    G,cs = dbg.all_contigs(kd,k)
    G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
    dbg.createGraphObject(G,cs,k,G_naive,pfn=False,ps=False)
    return G_naive

def createNaiveFromReads(fn,k,GraphObject):
    d = dbg.build(fn,k,1)
    G,cs = dbg.all_contigs(d,k)
    dbg.createGraphObject(G,cs,k,GraphObject)

def printKmerdictToFile(kd,outFile):
    f = open(outFile, 'w')
    for km in kd:
        f.write(km+"\n")
    f.close()

def createKmerDictFromSegmentList(SL,k):
    kd = collections.defaultdict(int)
    for segment in SL:
        for km in dbg.kmers(segment,k):
            kd[km] += 1
        for km in dbg.kmers(dbg.twin(segment),k):
            kd[km] += 1
    return kd



#---------------------------------------------------------------------------
#------------------------Functions to work with time------------------------
#---------------------------------------------------------------------------
def returnTime(timeInSeconds):
    timeInSeconds = int(timeInSeconds)
    if timeInSeconds<60:
        return "Time in seconds: "+str(timeInSeconds)+"\n"
    else:
        s = timeInSeconds%60
        m = timeInSeconds/60
        assert m*60+s==timeInSeconds
        return str(m)+" minutes and "+str(s)+" seconds\n"
        
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
        tf.close()
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



#----------------------------------------------------------------------------
#--------Functions to create k-merdicts from genome and read files-----------
#----------------------------------------------------------------------------
#Creates a k-merdict from a .fa or .fasta file
def createKmerDictFromGenomeFile(k,genomeFile):
    kmersInGenome = collections.defaultdict(int)
    h = open(genomeFile, "rU")
    genomeFileExtension = os.path.splitext(genomeFile)[1]
    if genomeFileExtension==".fa":
        genome = h.readline()
        genome = h.readline().rstrip('\n')
        h.close()
        for km in dbg.kmers(genome,k):
            rep_km = min(km,dbg.twin(km))
            kmersInGenome[rep_km] += 1
    elif genomeFileExtension==".fasta":
        genome = h.readline()
        for line in h:
            line = line.strip()
            line = line if all([c in alphabet for c in line]) else ""
            tl = dbg.twin(line)
            for km in dbg.kmers(line,k):
                rep_km = min(km,dbg.twin(km))
                kmersInGenome[rep_km] += 1
    numKmersInGenome = len(kmersInGenome)
    return kmersInGenome,numKmersInGenome

#Creates a k-merdict from two .fastq files
def createKmerDictFromReadFiles(k,readFiles):
    kmersInReads = collections.defaultdict(int)
    numReadsPerFile = 0
    for f in readFiles:
        h = open(f, "rU")
        for lineNr,line in enumerate(h):
            if (lineNr%4 == 1):
                numReadsPerFile += 1
                segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
                for s in segments:
                    for km in dbg.kmers(s,k):
                        rep_km = min(km,dbg.twin(km))
                        kmersInReads[rep_km] += 1
        h.close()
    numReadsPerFile = numReadsPerFile/2
    numKmersInReads = len(kmersInReads)
    return kmersInReads,numKmersInReads,numReadsPerFile



#----------------------------------------------------------------------------
#--------------------Functions to work with genome_info.csv------------------
#----------------------------------------------------------------------------
def writeGenomeInfoToFile(infoFile,variables):
    [genomeName,k,numKmersInGenome,numReadsPerFile,numKmersPerRead,numKmersInReads,numKmersInReads_twice,numKmersInReads_twice_BF,BF_ratio_1_vs_0,preprocessTime] = variables
    f = open(infoFile, 'w')
    f.write("genomeName,"+str(genomeName)+"\n")
    f.write("k,"+str(k)+"\n")
    f.write("numKmersInGenome,"+str(numKmersInGenome)+"\n")
    f.write("numReadsPerFile,"+str(numReadsPerFile)+"\n")
    f.write("numKmersPerRead,"+str(numKmersPerRead)+"\n")
    f.write("numKmersInReads,"+str(numKmersInReads)+"\n")
    f.write("numKmersInReads_twice,"+str(numKmersInReads_twice)+"\n")
    f.write("numKmersInReads_twice_BF,"+str(numKmersInReads_twice_BF)+"\n")
    f.write("BF_ratio_1_vs_0,"+str(BF_ratio_1_vs_0)+"\n")
    f.write("preprocessTime,"+str(preprocessTime)+"\n")
    f.close()

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


#----------------------------------------------------------------------------
#------------------------Other functions------------------------
#----------------------------------------------------------------------------
def readKmersFromFileToGraph(kmerFile,G):
    #Tekur skrá með lista af k-merum sem inntak ásamt tómu grafi G
    #Les alla k-mera úr skránni og bætir hverjum og einum við grafið
    assert(G.isEmpty())
    f = open(kmerFile, 'rU')
    for km in f:
        km = km.rstrip('\n')
        G.addSegmentToGraph(km)
    f.close()

def writeTotalsAndPercToFile(fileName,totals,percs):
	#totals = [tot_iso,tot_tip,tot_bub3,tot_bub4,tot_non]
	f = open(fileName, 'w')
	assert(len(totals)==len(percs))
	for i in xrange(len(totals)):
		f.write(str(totals[i])+";"+str(percs[i])+"\n")
	f.close()

def readTotalsAndPercFromFile(fileName):
    f = open(fileName, 'rU')
    totals = []
    percs = []
    for i, line in enumerate(f):
        line = line.rstrip('\n').split(";")
        totals.append(int(line[0]))
        percs.append(float(line[1]))
    f.close()
    return totals, percs

def fractionFromGenomeInObject(kmersInGenome,numKmersInGenome,Object,isGraph=1):
    #Input:
    #   kmersInGenome:      A dict storing all kmers actually occurring in the genome (stores rep_km)
    #   numKmersInGenome:   The number of k-mers in the genome
    #   Object:             Can either be a graph or a bloom filter
    #   isGraph:            Tells us whether object is a graph or a bloom filter (can't have both)
    #       isGraph==1: Object is a Graph
    #       isGraph==0: Object is a BF
    #Returns:
    #   The fraction of kmers from the genome that occur in Object    
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

#Tekur lista af maxCovs (-1 fremst) og prentar mismengi af kd_-1 og kd_maxCov fyrir hvert maxCov
def diff_GNoMaxCov_and_GMaxCov(maxCovs,outDirectory):
    assert(maxCovs[0]==-1), "-1 must be first in the list"
    kd1 = readKmersFromFileToDict(outDirectory+"/kmersFromG_maxCov_-1.txt")
    for maxCov in maxCovs[1:]:
        kd2 = readKmersFromFileToDict(outDirectory+"/kmersFromG_maxCov_"+str(maxCov)+".txt")
        printKmerdictToFile(difference(kd1,kd2),outDirectory+"/diff_GNoMaxCov_and_GMaxCov_"+str(maxCov)+".txt")

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

