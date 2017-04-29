#coding:utf8
import collections, sys
import dbg, Bloom, Graph, helpers
import os.path


#Inputs:
#	k: kmerLength
#	fn: a list of fastq files containing reads
#returns all substrings from segments whose kmers we've seen 2 or more times
#(with a false positive rate of p)
#whatToRun=0: Keyra BF_counter_naive
#whatToRun=1: Keyra BF_counter
def geraAllt(fn,k,BF,G,pfn=False,printProgress=False,startAtLine=0,sizeOfGenome=-1,maxCov=-1,whatToRun=1):
	if pfn:
		print "geraAllt",locals().keys(),"\n"
		#print "geraAllt(fn="+str(fn)+", k="+str(k)+", BF, G, pfn="+str(pfn)+", printProgress="+str(printProgress)+", startAtLine="+str(startAtLine)+")"
	if not isinstance(fn, list):
		print "fn:",fn
		raise Exception('fn has to be a list')
	assert (whatToRun==0 or whatToRun==1)
	#Open the file and get ready for gathering data for the Figures
	sampleDensity, readsPerLine = getInfoFromInputFile(fn,sizeOfGenome)
	lineNumber = []
	cov = []
	kmersInGraph = []
	kmersInBF = []
	G_ratio = []
	BF_ratio = []
	kmersInGenome = helpers.createKmerDictFrom_fa("Input/t/t.fa",k)
	count = 0
	if not maxCov==-1:
		assert os.path.isfile("Output/defaultOutFolder/Without_BF_stopper/figureData.csv")
		assert os.path.isfile("Output/defaultOutFolder/Without_BF_stopper/figureData2.csv") 
		print "The input files exist"
	print "maxCov="+str(maxCov)
	for f in fn:
		h = open(f, "rU")
		for lineNr,line in enumerate(h,start=startAtLine):
			if (lineNr%4 == 1):
				segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
				if len(segments)>1:
					print "lineNr="+str(lineNr)+", segments:",segments, "len(segments)="+str(len(segments))
				for s in segments:
					assert isinstance(s, basestring)

					#Prentum út stöðuna annað slagið:
					if printProgress and (lineNr%10000==1):
						print "We are reading the segment from line "+str(lineNr)+" from the file "+str(f)
						if lineNr%100000==1:
							ratio = BF.computeRatio()[0]
							B = BF.hasAcceptableRatio(ratio)
							if B:
								print "The BF is not full. Ratio="+str(ratio)
							else:
								print "The BF is full. Ratio="+str(ratio)

					#á "sampleDensity" segmenta fresti þá mælum við fjölda k-mera í
					#G og BF ásamt coverage og geymum niðurstöðurnar
					if count%sampleDensity==0:
						print "count%sampleDensity==0. sampleDensity="+str(sampleDensity)+", count="+str(count)
						lineNumber.append(lineNr)
						cov.append(count * readsPerLine / sizeOfGenome)
						kmersInGraph.append(len(G.kmers))
						kmersInBF.append(len(BF))
						g,b = helpers.ratioInGenome(kmersInGenome,G,BF)
						G_ratio.append(g)
						BF_ratio.append(b)

					#Bætum öllum k-merum úr s við BF
					#Bætum öllum k-merum úr s sem við höfum séð tvisvar eða 
					#oftar við G (bætum segmenti sem inniheldur k-mera)
					L = len(s)
					start = 0
					for i, kmer in enumerate(dbg.kmers(s,k)):
						#print lineNr,i,kmer,start
						if not (kmer in BF):
							if (maxCov==-1) or (cov[-1] <= maxCov):
								#print maxCov,cov
								BF.add(kmer)
								BF.add(dbg.twin(kmer))
							#else:
							#	if not maxCov==-1:
							#		print "We don't add to BF due to the cov stopper"
							if i-start>0:
								goodSequence = s[start:i+k-1]
								assert len(goodSequence)>=k
								for temp_km in dbg.kmers(goodSequence,k):
									assert temp_km in BF, "temp_km="+str(temp_km)
								G.addSegmentToGraph(goodSequence)
							start = i+1
					#If we reach the end of segment we add the current sequence
					if L-start>=k:
						assert len(s[start:])>=k
						G.addSegmentToGraph(s[start:])

				count += 1
				#print "lineNr="+str(lineNr)+", count="+str(count)+", len(BF)="+str(len(BF))
				#print len(s),len(list(dbg.kmers(s,k)))
				#break
		#Geymum einnig niðurstöðurnar í lokin
		lineNumber.append(lineNr)
		cov.append(int(count) * readsPerLine / sizeOfGenome)
		kmersInGraph.append(len(G.kmers))
		kmersInBF.append(len(BF))
		g,b = helpers.ratioInGenome(kmersInGenome,G,BF)
		G_ratio.append(g)
		BF_ratio.append(b)
	
	#Prentum niðurstöðurnar í skrá sem við getum notað síðar til að búa til mynd
	if maxCov==-1:
		outFolder="Output/defaultOutFolder/Without_BF_stopper"
		printFigureFiles(lineNumber,cov,kmersInGraph,kmersInBF,G_ratio,BF_ratio,outFolder)
		helpers.printFigureFromFile("Output/defaultOutFolder/Without_BF_stopper/figureData.csv",sizeOfGenome,"#k-mers as a function of cov (log scale below)",outFolder)
		helpers.printFigureFromFile2("Output/defaultOutFolder/Without_BF_stopper/figureData2.csv","#Ratio of k-mers in G and BF",outFolder)
	else:
		#Assume that the figureData.csv and figureData2.csv have
		#already been created without using a BF-stopper
		#i.e. there should be 4 columns in the file, not 6
		#before we call changeFile
		outFolder="Output/defaultOutFolder/BF_stopper"
		helpers.changeFile("Output/defaultOutFolder/Without_BF_stopper/figureData.csv",kmersInGraph,kmersInBF,outFolder)
		helpers.changeFile("Output/defaultOutFolder/Without_BF_stopper/figureData2.csv",G_ratio,BF_ratio,outFolder)
		helpers.printFigureFromFile("Output/defaultOutFolder/BF_stopper/figureData.csv",sizeOfGenome,"#k-mers as a function of cov (log scale below) with BF-stopper at cov="+str(maxCov),outFolder)
		helpers.printFigureFromFile2("Output/defaultOutFolder/BF_stopper/figureData2.csv","#Ratio of k-mers in G and BF with BF-stopper at cov="+str(maxCov),outFolder)

#fn:			A list of 1 or more .fastq files
#					fn = [file1.fastq, file2.fastq, ...]
#BF:			Bloom filter
#k:				kmer length
#G:				Graph
#startAtLine:	The first line in the input file we read from 
#				Useful when we have already read the earlier lines from the
#				file and want to continue from there instead of starting over
def BF_counter(fn,BF,k,G,pfn=False,printProgress=False,startAtLine=0,sizeOfGenome=-1,maxCov=-1):
	geraAllt(fn,k,BF,G,pfn,printProgress,startAtLine,sizeOfGenome,maxCov,whatToRun=1)

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
	numberOfMeasurements = 10
	sampleDensity = numberOfReads/(numberOfMeasurements-1)
	return sampleDensity, readsPerLine

def printFigureFiles(lineNumber,cov,kmersInGraph,kmersInBF,G_ratio,BF_ratio,outFolder=""):
	print "printFigureFile(lineNumber,cov,kmersInGraph,kmersInBF)"
	if outFolder=="":
		h1 = open("Output/defaultOutFolder/Without_BF_stopper/figureData.csv","w")
	else:
		h1 = open(outFolder+"/figureData.csv","w")
	#print lineNumber
	#print cov
	#print kmersInGraph
	#print kmersInBF
	for i in xrange(0,len(cov)):
		h1.write(str(lineNumber[i])+","+str(cov[i])+","+str(kmersInGraph[i])+","+str(kmersInBF[i])+"\n")
	h1.close()

	if outFolder=="":
		h2 = open("Output/defaultOutFolder/Without_BF_stopper/figureData2.csv","w")
	else:
		h2 = open(outFolder+"/figureData2.csv","w")
	for i in xrange(0,len(cov)):
		h2.write(str(lineNumber[i])+","+str(cov[i])+","+str(round(G_ratio[i],4))+","+str(round(BF_ratio[i],4))+"\n")
	h2.close()
"""
def printFigureFile2(lineNumber,cov,G_ratio,BF_ratio):
	print "printFigureFile2(lineNumber,cov,G_ratio,BF_ratio)"
	h = open("Output/defaultOutFolder/figureData2.csv","w")
	for i in xrange(0,len(cov)):
		h.write(str(lineNumber[i])+","+str(cov[i])+","+str(G_ratio[i])+","+str(BF_ratio[i])+"\n")
	h.close()
"""

def printAllInfoFromFiles(fn,k):
	kd = collections.defaultdict(int)		#Stores all kmers in the files
	kd_BF = collections.defaultdict(int)	#Stores all kmers in the files
											#seen for the second time according to BF
	numberOfKmers = 0
	BF = Bloom.Bloom(0.01,6000000,pfn=True)
	print "Running through the files:"
	for f in fn:
		h = open(f, "rU")
		for lineNr,line in enumerate(h):
			if (lineNr%4 == 1):
				segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
				if len(segments)>1:
					print "lineNr="+str(lineNr)+", segments:",segments, "len(segments)="+str(len(segments))
				for s in segments:
					for km in dbg.kmers(s,k):
						kd[km] += 1
						numberOfKmers += 1
						if not km in BF:
							BF.add(km)
						else:
							kd_BF[km] += 1
					for km in dbg.kmers(dbg.twin(s),k):
						kd[km] += 1
						numberOfKmers +=1
						if not km in BF:
							BF.add(km)
						else:
							kd_BF[km] += 1

	G_twice = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
	numberOfAtLeastTwice = 0
	for km,num in kd.iteritems():
		if num > 1:
			G_twice.addSegmentToGraph(km)
			numberOfAtLeastTwice += 1
	G_BF = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
	for km,num in kd_BF.iteritems():
		G_BF.addSegmentToGraph(km)
	
	print "Total number of kmers in the files:                  ", numberOfKmers
	print "Number of unique kmers in the files:                 ", len(kd)
	print "Number of kmers in the BF:                           ", len(BF)
	print "Number of kmers occuring at least twice in the files: ", numberOfAtLeastTwice
	print "Same number according to BF:                          ", len(kd_BF)
	print "Number of kmers in G_twice:                           ", len(G_twice.kmers)
	print "Number of contigs in G_twice:                          ", len(G_twice)
	print "Number of kmers in G_BF:                              ", len(G_BF.kmers)
	print "Number of contigs in G_BF:                            ", len(G_BF)
	print BF
	BF.print_bitarray()

def BF_counter_naive(fn,BF,k,G_naive,pfn=True,printProgress=False):
	if pfn:
		print "BF_counter_naive(fn, BF, k="+str(k)+", G_naive, pfn="+str(pfn)+", printProgress="+str(printProgress)+")"
	kd = collections.defaultdict(int)
	numberOfKmers = 0
	for f in fn:
		h = open(f, "rU")
		for lineNr,line in enumerate(h):
			if (lineNr%4 == 1):
				segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
				if len(segments)>1:
					print "lineNr="+str(lineNr)+", segments:",segments, "len(segments)="+str(len(segments))
				for s in segments:
					for km in dbg.kmers(s,k):
						numberOfKmers += 1
						if not km in BF:
							BF.add(km)
						else:
							kd[km] += 1

					for km in dbg.kmers(dbg.twin(s),k):
						numberOfKmers += 1
						if not km in BF:
							BF.add(km)
						else:
							kd[km] += 1

	print "Total number of kmers in the files:                  ", numberOfKmers
	print "Number of kmers occuring at least twice in the files: ", len(kd)

	G,cs = dbg.all_contigs(kd,k)
	dbg.createGraphObject(G,cs,k,G_naive,pfn,ps=False)
	print "Number of kmers in G_naive:",len(G_naive.kmers)

if __name__ == "__main__":
	fn = ['Input/t/r1.fastq', 'Input/t/r2.fastq']
	k = 31
	printAllInfoFromFiles(fn,k)
	#G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
	#BF = Bloom.Bloom(0.01,6000000,pfn=True)
	#BF_counter_naive(fn,BF,k,G_naive,pfn=True,printProgress=False)