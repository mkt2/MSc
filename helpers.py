#coding:utf8
import dbg
import collections
import os.path
import Graph
import time
import sys
from math import ceil
#sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

alphabet = ["A","C","G","T","N"]

#------------------------Helper functions------------------------
def fix_i(i,L,k):
    #L is the length of a segment
    #i is the first index of the k-mer
    #k is the k-mer length
    #Takes the index in the twin and returns the index in the sequence
    return L-k-i

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

def splitCov(s0,s1,k,COV):
    #pre:   s0 and s1 are outputs of splitString(s,i,k)
    L0 = float(cLen(s0,k))
    L1 = cLen(s1,k)
    COV0 = int(ceil(float(COV)*(L0/(L0+L1))))
    COV1 = COV - COV0
    assert(COV0>=2), "COV0="+str(COV0)+". COV="+str(COV)
    assert(COV1>=2), "COV1="+str(COV1)
    return COV0, COV1

def reverseList(L):
    for i, (ID,B) in enumerate(L):
        L[i] = (ID,not B)

def assert_isKmer(km,k):
    assert(isinstance(km,str)), "The k-mer must be a str"
    assert(len(km)==k), "The k-mer must be of length k"

"""
def splitOnSeen(s,k,kmersIn_s,out=[]):
    #kmersIn_s=collections.defaultdict(int)
    #splits a sequence s up into sub sequences where we skip k-mers previously seen in s
    #if no k-mers occur twice in s then we return [s]
    if len(s)<=k:
        return out
    for i, km in enumerate(dbg.kmers(s,k)):
        rep_km = min(km,dbg.twin(km))
        if not (rep_km in kmersIn_s):
            kmersIn_s[rep_km] += 2
        else:   #rep_km in kmersIn_s
            kmersIn_s[rep_km] += 1
            #We have seen s[i:i+k] before
            #split s into s[0:i+k-1]    s[i:i+k]    s[i+1:]
            #             add to out     throw     recursion
            if i==0:
                return splitOnSeen(s[1:],k,kmersIn_s,out)
            else:
                out.append(s[0:i+k-1])
                if i==len(s)-k:
                    return out
                else:
                    return splitOnSeen(s[i+1:],k,kmersIn_s,out)
    return out+[s]
"""

def splitOnConnToSelf(s,s_twin,L,k,i,j,start):
    #    i                          j         L
    #s: |                          |  trash  |

    #            i   i+k            j         L
    #s: | added |       unseen     |  trash  |
    #           | a |
    #                
    #            start  i           j         L
    #s: | added | good |   unseen  |  trash  |
    #                  | a |
    #                       i+k
    #print "splitOnConnToSelf(s="+str(s)+", s_twin="+str(s_twin)+", k="+str(k)+", i="+str(i)+", L="+str(L)+", start="+str(start)+")"
    if (j-start<k) or (j-i<k):
        #print "case -1"
        return
    if j-start==k:
        #print "case 0"
        #yield s[start:j]
        yield start, j
        return
    a = s[i:i+k]
    a_twin = s_twin[L-k-i:L-i]
    assert(a_twin==dbg.twin(a))
    #assert(len(a_twin)==k), a_twin
    for y in dbg.bw(a):
        #a->a
        if (y==a):
            #print "case 1"
            if i>0 and i>start:
                yield start, i+k-1
                #yield s[start:i+k-1]
            #yield a
            yield i,i+k
            for temp1,temp2 in splitOnConnToSelf(s,s_twin,L,k,i+1,j,i+1):
                yield temp1,temp2
            return
        #twin(a)->a
        if (y==a_twin):
            #print "case 2"
            if i>0 and i>start:
                yield start, i+k-1
                #yield s[start:i+k-1]
            start = i
        #b->a
        if i>0 and i>start:
            index = s.find(s, i+1)
            if index!=-1:
                #print "case 3"
                yield start, i+k-1
                #yield s[start:i+k-1]
                start = i
        #twin(b)->a
        if i>0 and i>start:
            index = s_twin.find(y,(L-k-i)+1)
            if index!=-1:
                #print "case 4"
                yield start, i+k-1
                #yield s[start:i+k-1]
                start = i
    #a->twin(a)
    for y in dbg.fw(a):
        if y==a_twin:
            #print "case 5"
            #yield s[start:i+k]
            yield start, i+k
            for temp1,temp2 in splitOnConnToSelf(s,s_twin,L,k,i+1,j,i+1):
                yield temp1,temp2
            return
        #a->b
        index = s.find(y, i+2)
        if index!=-1:
            #print "case 6, y="+str(y)
            yield start, i+k
            #yield s[start:i+k]
            for temp1,temp2 in splitOnConnToSelf(s,s_twin,L,k,i+1,j,i+1):
                yield temp1,temp2
            return
    #a->twin(b) = b->twin(a)
    for y in dbg.bw(a_twin):
        index = s.find(y, i+1)
        if index!=-1:
            #print "case 7"
            #yield s[start:i+k]
            yield start, i+k
            for temp1,temp2 in splitOnConnToSelf(s,s_twin,L,k,i+1,j,i+1):
                yield temp1,temp2
            return
    #If a is the last kmer we return s[start:]
    if i+k==j:
        #print "case 8"
        #yield s[start:j]
        yield start, j
        return
    #Recursively call to get the next a
    #print "case 9"
    for temp1,temp2 in splitOnConnToSelf(s,s_twin,L,k,i+1,j,start):
        yield temp1,temp2

def canConn(a,b,A,B,k,assertStuff=False):
    #Returns True if k-mers a and b cann connect together
    #according to booleans A and B
    if assertStuff:
        assert_isKmer(a,k)
        assert_isKmer(b,k)
        assert(isinstance(A,bool)), "A must be a boolean"
        assert(isinstance(B,bool)), "B must be a boolean"
        assert(isinstance(k,int)), "k must be an integer"
    if (A==True) and (B==True):
        for y in dbg.fw(a):
            if y==b:
                return True
        return False
        #return a[1:]==b[:-1]
    if (A==True) and (B==False):
        for y in dbg.fw(a):
            if y==dbg.twin(b):
                return True
        return False
        #return a[1:]==dbg.twin(b)[:-1]
    if (A==False) and (B==True):
        for y in dbg.bw(b):
            if y==dbg.twin(a):
                return True
        return False
        #return dbg.twin(a)[1:]==b[:-1]
    if (A==False) and (B==False):
        for y in dbg.fw(b):
            if y==a:
                return True
        return False
        #return canConn(b,a,True,True,k,False)

def canConnFromTo(a,b,k,assertStuff=True):
    #returns True if a->b or a->twin(b) possible
    if assertStuff:
        assert_isKmer(a,k)
        assert_isKmer(b,k)
        assert(isinstance(k,int)), "k must be an integer"
    return canConn(a,b,True,True,k,False) or canConn(a,b,True,False,k,False)

def canConnToFrom(a,b,k,assertStuff=True):
    #returns True if b->a or twin(b)->a possible
    if assertStuff:
        assert_isKmer(a,k)
        assert_isKmer(b,k)
        assert(isinstance(k,int)), "k must be an integer"
    return canConn(b,a,True,True,k,False) or canConn(b,a,False,True,k,False)

def canConnToSelf(a,k,assertStuff=True):
    if assertStuff:
        assert_isKmer(a,k)
        assert(isinstance(k,int)), "k must be an integer"
    return canConn(a,a,True,True,k,False) or canConn(a,a,True,False,k,False) or canConn(a,a,False,True,k,False)

def canConnToSelfFront(a,k,assertStuff=True):
    #Returns True if the front of a can connect to a og twin(a)
    if assertStuff:
        assert_isKmer(a,k)
        assert(isinstance(k,int)), "k must be an integer"
    return canConn(a,a,True,True,k,False) or canConn(a,a,True,False,k,False)

def canConnToSelfBack(a,k,assertStuff=True):
    if assertStuff:
        assert_isKmer(a,k)
        assert(isinstance(k,int)), "k must be an integer"
    return canConn(a,a,True,True,k,False) or canConn(a,a,False,True,k,False)

#Before:    c is a string representing a contig
#After:     Returns the number of k-mers in c
def cLen(c,k):
    assert(len(c)>=k), "c must have a minimum length k"
    num_kmers = len(list(dbg.kmers(c,k)))
    num_kmers_fast = len(c)-k+1
    assert(num_kmers==num_kmers_fast), "num_kmers_fast should return the number of k-mers."+str(num_kmers)+", "+str(num_kmers_fast)
    return num_kmers_fast

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

def createGraphName(maxAddCov,maxSplitCov,withSpaces=False,latexFormat=False):
    latex_inf = "$\\infty$"
    if latexFormat:
        s = "G"
    else:
        s = "G_"
    if withSpaces:
        s+=" "
    if maxAddCov == float('inf'):
        if latexFormat:
            s += latex_inf
        else:
            s += 'inf'
    else:
        s += str(maxAddCov)
    if latexFormat:
        pass
    else:
        s += "_"
    if maxSplitCov==float('inf'):
        if latexFormat:
            s += latex_inf
        else:
            s += 'inf'
    else:
        if withSpaces:
            s += str(maxSplitCov)
        else:
            s += str(maxSplitCov)
    #s = G10s15
    return s
    #G_10_15
    #G_inf_inf

#---------------------------------------------------------------------------
#------------------------Functions to work with time------------------------
#---------------------------------------------------------------------------
def returnTime(timeInSeconds):
    timeInSeconds = int(timeInSeconds)
    val = secTo_m_s(timeInSeconds)
    if timeInSeconds<60:
        return str(timeInSeconds)+" seconds"
    else:
        return str(val[0])+" minutes and "+str(val[1])+" seconds"

def secTo_m_s(timeInSeconds):
    timeInSeconds = int(timeInSeconds)
    if timeInSeconds<60:
        return timeInSeconds
    else:
        s = timeInSeconds%60
        m = timeInSeconds/60
        assert m*60+s==timeInSeconds
        return m, s
        
#Prints the run times to a file in a format which can be copied 
#directly into a latex tabular environment
#def printRuntimesToFile(genomeName, maxAddCovs, maxSplitCovs, runTimes):
def printRuntimesToFile(genomeName, filters, runTimes):
    print "printRunTimesToFile()"
    def modRunTimes(runTimes):
        runTimeStrings = [""]*len(runTimes)
        for i, rt in enumerate(runTimes):
            val = secTo_m_s(rt)
            if rt < 60:
                runTimeStrings[i] = str(val) + " sec"
            else:
                runTimeStrings[i] = str(val[0]) + " min, " + str(val[1]) + " sec"
        return runTimeStrings
    
    #maxAddCovs = [5, 10, 15, 20, 30,float('inf')]
    runTimes = modRunTimes(runTimes)
    timeFile = "Output/runTimes_t.txt"
    latex_inf = "$\\\\infty$"
    latex_newline = "\\\\\n"
    if genomeName=="t":
        tf = open(timeFile, 'w')
        for i,(MAC,MSC) in enumerate(filters):
            #for maxSplitCov in maxSplitCovs[i]:
            s = createGraphName(MAC,MSC)
            s += " & "+runTimes[i]+" & "+latex_newline
            #s = G10s15 & 59 sec & newline
            tf.write(s)
        tf.close()
    elif genomeName=="sa":
        #raise Exception("printRunTimesToFile er ekki tilbúið fyrir sa")
        tf_old = open(timeFile, 'r')
        tf = open("Output/runTimes_tAndSA.txt", 'w')
        #graphNames = []
        #rest = []
        nameDict = collections.defaultdict(tuple)
        #G20->[t_time, sa_time]
        for i, line in enumerate(tf_old):
            line=line[0:-3]
            #print line
            temp = line.strip().split(" & ")
            temp[1] = temp[1][0:-2]
            #print temp
            if temp[0]=="G$\\infty$":
                nameDict["G"+latex_inf] = [temp[1],-1]
            else:
                nameDict[temp[0]] = [temp[1],-1]
        runTimeCounter = -1
        for i,(MAC,MSC) in enumerate(filters):
            #for maxSplitCov in maxSplitCovs[i]:
            runTimeCounter += 1
            gn = createGraphName(MAC,MSC)
            if gn in nameDict:
                nameDict[gn][1] = runTimes[runTimeCounter]
            else:
                nameDict[gn] = ["",runTimes[runTimeCounter]]
                #tf.write(gn+" & "+rest[j]+" "+runTimes[index]+latex_newline)
        
        for key,value in nameDict.iteritems():
            print key,value
        Ginf_key = "G_inf_inf"
        temp = nameDict[Ginf_key]
        Ginf_val = "G$\\\\infty$ & "+temp[0]+" & "+temp[1]+latex_newline
        del nameDict[Ginf_key]
        nameList = [k+" & "+v[0]+" & "+str(v[1])+latex_newline for k,v in nameDict.iteritems()]
        print nameList
        print nameList[0].split("_")

        nameList.sort(key=lambda x: (int(x.split("_")[1]), float(x.split("_")[2][0:3])))
        for i, v in enumerate(nameList):
            nameList[i] = nameList[i][0]+nameList[i][2:]
            if "s " in nameList[i]:
                index = nameList[i].find("s ")
                nameList[i] = nameList[i][0:index-1]+nameList[i][index]+nameList[i][index+2:]
        nameList.append(Ginf_val)
        for v in nameList:
            #print v
            tf.write(v)
    else:
        raise Exception("Illegal value for genomeName="+str(genomeName)+". It must be either t or sa")



#----------------------------------------------------------------------------
#-------------Functions to work with raw genome and read files---------------
#----------------------------------------------------------------------------
def segments(fn,k):
    #yields segments from the reads, one segment at a time
    segmentCount = -1
    for f in fn:
        h = open(f, "rU")
        for lineNr,line in enumerate(h):
            if (lineNr%4 == 1):
                segmentCount+=1
                segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
                for s in segments:
                    yield s, segmentCount

def createKmerDictFromGenomeFile(k,genomeFile):
    #Creates a k-merdict from a .fa or .fasta file
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
            for km in dbg.kmers(line,k):
                rep_km = min(km,dbg.twin(km))
                kmersInGenome[rep_km] += 1
    return kmersInGenome

def createKmerDictFromReadFiles(k,readFiles):
    #Creates a k-merdict from two .fastq files
    kmersInReads = collections.defaultdict(int)
    for s in segments(readFiles,k):
        for km in dbg.kmers(s,k):
            rep_km = min(km,dbg.twin(km))
            kmersInReads[rep_km] += 1
    return kmersInReads


#---------------------------------------------------------------------------
#---------------------Functions to work with kmerDicts----------------------
#---------------------------------------------------------------------------
#Skilar öllum k-merum sem koma fyrir í kd1 en ekki í kd2 sem mengi
#I.e. skilar mismengi kd1 og kd2
def difference(kd1,kd2):
    return { x : kd1[x] for x in set(kd1) - set(kd2) }

def createNaiveFromKmerDict(kd,k):
    G,cs = dbg.all_contigs(kd,k)
    G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
    dbg.createGraphObject(G,cs,k,G_naive,pfn=False,ps=False)
    return G_naive

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


def createDictOfNonSingletons(fn,k):
    kmerDict_correct = collections.defaultdict(int)
    for segment in segments(fn,k):
        for km in dbg.kmers(segment,k):
            kmerDict_correct[km] += 1
        for km in dbg.kmers(dbg.twin(segment),k):
            kmerDict_correct[km] += 1

    #Throw all singletons:
    kmerDict_correct = {k:v for k,v in kmerDict_correct.items() if v != 0}
    return kmerDict_correct

def dictEqualsOther(dict1,dict2):
    if not len(dict1)==len(dict2):
        return False
    for km in dict1:
        if not km in dict2:
            return False
    for km in dict2:
        if not km in dict1:
            raise Exception("this shouldn't happen")
    return True

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

def printCovDictToFile(covDict,fileName):
    f = open(fileName, 'w')
    #for maxAddCov, values in covDict.iteritems():   #<---virkar ekki því covDict.keys() er núna (maxAddCov,maxSplitCov)
    for key, values in covDict.iteritems():
        maxAddCov = key[0]
        maxSplitCov = key[1]
        #values = [COV, G_len,G_frac,BF_len,BF_frac] fyrir G (maxAddCov = float('inf'))
        #values = [COV, G_len,G_frac] fyrir sérhvert Gx
        if maxAddCov == float('inf'):
            assert(len(values)==5)
            #[COV, G_len, G_frac, BF_len, BF_frac] = values
            f.write(str(key)+";"+str(values[0])+";"+str(values[1])+";"+str(values[2])+";"+str(values[3])+";"+str(values[4])+"\n")
        else:
            assert(len(values)==3)
            f.write(str(key)+";"+str(values[0])+";"+str(values[1])+";"+str(values[2])+"\n")
    f.close()

def readCovDictFromFile(fileName):
    f = open(fileName, 'rU')
    covDict = collections.defaultdict(list)
    inf = float('inf')  #Used in the eval statement to convert inf to float('inf)
    for line in f:
        line = line.rstrip('\n').split(";")
        key = eval(line[0])
        maxAddCov = float(key[0])
        maxSplitCov = float(key[1])
        #print maxAddCov, maxSplitCov
        if maxAddCov!=float('inf'):
            maxAddCov = int(maxAddCov)
        covDict[key] = [[]]*(len(line)-1)
        for j in range(0,len(line)-1):
            covDict[key][j] = eval(line[j+1])
    f.close()
    return covDict

'''
#Tekur lista af maxAddCovs (-1 fremst) og prentar mismengi af kd_-1 og kd_maxAddCov fyrir hvert maxAddCov
def diff_GNoMaxCov_and_GMaxCov(maxAddCovs,outDirectory):
    assert(maxAddCovs[0]==-1), "-1 must be first in the list"
    kd1 = readKmersFromFileToDict(outDirectory+"/kmersFromG_maxAddCov_-1.txt")
    for maxAddCov in maxAddCovs[1:]:
        kd2 = readKmersFromFileToDict(outDirectory+"/kmersFromG_maxAddCov_"+str(maxAddCov)+".txt")
        printKmerdictToFile(difference(kd1,kd2),outDirectory+"/diff_GNoMaxCov_and_GMaxCov_"+str(maxAddCov)+".txt")
'''

def readSkipPrint(L):
    assert(len(sys.argv)<=L), "L is either the length of sys.argv or the index of it's last value"
    if len(sys.argv)==L:
        if sys.argv[L-1]=="False":
            skipPrint = False
        elif sys.argv[L-1]=="True":
            skipPrint = True
        else:
            assert(False), "Illegal value for skipPrint"
    else:
        skipPrint = False
    return skipPrint

def readKmersFromFileToDict(kmerFile):
    kmerDict = collections.defaultdict(int)
    f = open(kmerFile, 'rU')
    for km in f:
        km = km.rstrip('\n')
        kmerDict[km] = 1
    f.close()
    return kmerDict

def createNaiveFromReads(fn,k,GraphObject):
    d = dbg.build(fn,k,1)
    G,cs = dbg.all_contigs(d,k)
    dbg.createGraphObject(G,cs,k,GraphObject)



if __name__ == "__main__":
    genomeName = "sa"
    filters = [ \
    (15,20),(15,30),(15,float('inf')), \
    (20,30),(20,float('inf')), \
    (30,float('inf')), \
    (float('inf'),float('inf'))]
    numGraphs = len(filters)
    runTimes = [-1]*numGraphs
    for i in range(0,len(runTimes)):
        runTimes[i] = i

    print runTimes
    printRuntimesToFile(genomeName, filters, runTimes)