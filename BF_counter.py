#coding:utf8
import collections, sys
import dbg, Bloom, Graph, helpers

#Inputs:
#	k: kmerLength
#	fn: a list of fastq files containing reads
#returns all substrings from segments whose kmers we've seen 2 or more times
#(with a false positive rate of p)
#whatToRun=0: Keyra BF_counter_naive
#whatToRun=1: Keyra BF_counter
def geraAllt(fn,k,BF,G,pfn=False,printProgress=False,startAtLine=0,sizeOfGenome=-1,whatToRun=1):
	if pfn:
		print "generateSegments_from_fastq(fn="+str(fn)+", k="+str(k)+", BF, pfn="+str(pfn)+", printProgress="+str(printProgress)+", startAtLine="+str(startAtLine)+")"
	if not isinstance(fn, list):
		print "fn:",fn
		raise Exception('fn has to be a list')
	#Open the file and get ready for gathering data for the Figures
	if sizeOfGenome==-1:
		sizeOfGenome = 70000
	sampleDensity, readsPerLine = getInfoFromInputFile(fn,sizeOfGenome)
	lineNumber = [0]
	cov = [0]
	kmersInGraph = [0]
	kmersInBF = [0]
	G_ratio = [0]
	BF_ratio = [0]
	count = 0
	kmersInGenome = helpers.createKmerDictFrom_fa("Input/t/t.fa",k)
	count = 0
	for f in fn:
		h = open(f, "rU")
		for lineNr,line in enumerate(h,start=startAtLine):
			if (lineNr%4 == 1):
				segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
				for s in segments:
					if not isinstance(s, basestring):
						print "s:",s
						raise Exception('Each segment has to be a string')

					#Prentum út stöðuna annað slagið:
					if printProgress and (lineNr%5000==1):
						print "We are reading the segment from line "+str(lineNr)+" from the file "+str(f)
						if lineNr%40000==1:
							ratio = BF.computeRatio()[0]
							B = BF.hasAcceptableRatio(ratio)
							if B:
								print "The BF is not full. Ratio="+str(ratio)
							else:
								print "The BF is full. Ratio="+str(ratio)
					#Bætum öllum k-merum úr s við BF
					#Bætum öllum k-merum úr s sem við höfum séð tvisvar eða 
					#oftar við G (bætum segmenti sem inniheldur k-mera)
					L = len(s)
					start = 0
					for i, kmer in enumerate(dbg.kmers(s,k)):
						#print i,kmer,start,end
						if not (kmer in BF):
							BF.add(kmer)
							BF.add(dbg.twin(kmer))
							if i-start>0:
								assert len(s[start:i+k-1])>=k
								#yield s[start:i+k-1],lineNr
								G.addSegmentToGraph(s)
							start = i+1
					#If we reach the end of segment we add the current sequence
					if L-start>=k:
						assert len(s[start:])>=k
						#yield s[start:],lineNr
						G.addSegmentToGraph(s)

					#á "sampleDensity" segmenta fresti þá mælum við fjölda k-mera í
					#G og BF ásamt coverage og geymum niðurstöðurnar
					if count%sampleDensity==0:
						print "count%sampleDensity==0. sampleDensity="+str(sampleDensity)+", count="+str(count)
						#cov = numberOfReads*readsPerLine/sizeOfGenome
						lineNumber.append(lineNr)
						cov.append(int(count/2) * readsPerLine / sizeOfGenome)
						kmersInGraph.append(len(G))
						kmersInBF.append(len(BF))
						g,b = helpers.ratioInGenome(kmersInGenome,G,BF)
						G_ratio.append(g)
						BF_ratio.append(b)
				count += 1
		#Geymum einnig niðurstöðurnar í lokin
		lineNumber.append(lineNr)
		cov.append(int(count/2) * readsPerLine / sizeOfGenome)
		kmersInGraph.append(len(G))
		kmersInBF.append(len(BF))
		g,b = helpers.ratioInGenome(kmersInGenome,G,BF)
		G_ratio.append(g)
		BF_ratio.append(b)
	#Prentum niðurstöðurnar í skrá sem við getum notað síðar til að búa til mynd
	printFigureFile(lineNumber,cov,kmersInGraph,kmersInBF)
	printFigureFile2(lineNumber,cov,G_ratio,BF_ratio)
			


#fn:			A list of 1 or more .fastq files
#					fn = [file1.fastq, file2.fastq, ...]
#BF:			Bloom filter
#k:				kmer length
#G:				Graph
#startAtLine:	The first line in the input file we read from 
#				Useful when we have already read the earlier lines from the
#				file and want to continue from there instead of starting over
def BF_counter(fn,BF,k,G,pfn=False,printProgress=False,startAtLine=0,sizeOfGenome=-1):
	geraAllt(fn,k,BF,G,pfn,printProgress,startAtLine,sizeOfGenome,whatToRun=1)

def BF_counter_naive(fn,BF,k,G_naive,pfn=True,printProgress=False):
	if pfn:
		print "BF_counter_naive(fn, BF, k="+str(k)+", G_naive, pfn="+str(pfn)+", printProgress="+str(printProgress)+")"
	#Generate the kmerDict of all seen kmers
	kmerDict = createKmerDict(fn,BF,k,pfn,printProgress)
	createNaiveFrom_kmerDict(kmerDict,k,G_naive,pfn)

def createKmerDict(fn,BF,k=31,pfn=False,printProgress=False):
	if pfn:
		print createKmerDict.__name__+": "+str(locals())
	kmerDict = collections.defaultdict(int)
	count = 0
	for s,lineNr in generateSegmentsFrom_fq(fn,k,BF,pfn,printProgress,startAtLine=0):
		assert isinstance(s,str)
		assert isinstance(lineNr,int)
		for kmer in dbg.kmers(s,k):
			kmerDict[kmer] = 1
		for kmer in dbg.kmers(dbg.twin(s),k):
			kmerDict[kmer] = 1
	return kmerDict

def createKmerDictFromSegmentList(segmentList,k):
    kmerDict = collections.defaultdict(int)
    for s in segmentList:
        for kmer in dbg.kmers(s,k):
            kmerDict[kmer] = 1
        for kmer in dbg.kmers(dbg.twin(s),k):
            kmerDict[kmer] = 1
    return kmerDict

#Generate the DBG graph from the kmers using dbg.py
def createNaiveFrom_kmerDict(kmerDict,k,G_naive,pfn=False):
	if pfn:
		print createNaiveFrom_kmerDict.__name__
	G,cs = dbg.all_contigs(kmerDict,k)
	dbg.createGraphObject(G,cs,k,G_naive,pfn,ps=False)

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def getInfoFromInputFile(fn,sizeOfGenome):
	f = fn[0]
	numberOfLines  = file_len(f)
	numberOfReads = int(numberOfLines/4)
	h = open(f, "rU")
	line = h.readline()
	line = h.readline()
	readsPerLine = len(line)
	h.close()
	numberOfMeasurements = 11
	sampleDensity = numberOfReads/(numberOfMeasurements-1)
	return sampleDensity, readsPerLine

def printFigureFile(lineNumber,cov,kmersInGraph,kmersInBF):
	print "printFigureFile(lineNumber,cov,kmersInGraph,kmersInBF)"
	h = open("Output/defaultOutFolder/figureData.csv","w")
	#print lineNumber
	#print cov
	#print kmersInGraph
	#print kmersInBF
	for i in xrange(0,len(cov)):
		h.write(str(lineNumber[i])+","+str(cov[i])+","+str(kmersInGraph[i])+","+str(kmersInBF[i])+"\n")
	h.close()

def printFigureFile2(lineNumber,cov,G_ratio,BF_ratio):
	print "printFigureFile2(lineNumber,cov,G_ratio,BF_ratio)"
	h = open("Output/defaultOutFolder/figureData2.csv","w")
	#print lineNumber
	#print cov
	#print G_ratio
	#print BF_ratio
	for i in xrange(0,len(cov)):
		h.write(str(lineNumber[i])+","+str(cov[i])+","+str(G_ratio[i])+","+str(BF_ratio[i])+"\n")
	h.close()

def BF_counter_and_naive(fn,BF,k,G,G_naive,printProgress=False,pfn=False):
	if pfn:
		print "BF_counter_and_naive(fn, BF, k="+str(k)+", G, G_naive, printProgress="+str(printProgress)+", pfn="+str(pfn)+")"
	kmerDict = collections.defaultdict(int)	#for naive
	count = 0
	for segment,lineNr in generateSegmentsFrom_fq(fn,k,BF,pfn):
		if printProgress==True:
			print "Adding segment", count#, segment
			count+=1
		G.addSegmentToGraph(segment)

		#Add every kmer to a dictionary (naive):
		for kmer in dbg.kmers(segment,k):
			kmerDict[kmer] = 1
		for kmer in dbg.kmers(dbg.twin(segment),k):
			kmerDict[kmer] = 1

	#Generate the DBG graph from the kmers using dbg.py (naive)
	createNaiveFrom_kmerDict(kmerDict,k,G_naive,pfn)
	return G, G_naive


if __name__ == "__main__":
	segments = ["AAAAA","TCTATCCTATA","ATATATANNN","ATTTATANNNNNNCCCCTTATATA"]
	#segments = ["AAANCCCNAAAN"]
	y = onlyACGT(segments,3)
	print segments
	print y