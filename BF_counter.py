#coding:utf8
import collections, sys
import dbg, Bloom, Graph, helpers
import os.path

alphabet = ["A","C","G","T","N"]

#dasdf
#Inputs:
#	k: kmerLength
#	fn: a list of fastq files containing reads
#returns all substrings from segments whose kmers we've seen 2 or more times
#(with a false positive rate of p)
#whatToRun=0: Keyra BF_counter_naive
#whatToRun=1: Keyra BF_counter
def BF_counter(fn,k,BF,G,pfn=False,printProgress=False,startAtLine=0,sizeOfGenome=-1,maxCov=-1,genomeFile="",outDirectory=""):
	if pfn:
		#print "BF_counter",locals().keys(),"\n"
		print "BF_counter",locals()
		#print "geraAllt(fn="+str(fn)+", k="+str(k)+", BF, G, pfn="+str(pfn)+", printProgress="+str(printProgress)+", startAtLine="+str(startAtLine)+")"
	if not isinstance(fn, list):
		print "fn:",fn
		raise Exception('fn has to be a list')
	#assert (whatToRun==0 or whatToRun==1)
	#Open the file and get ready for gathering data for the Figures
	if genomeFile=="":
		print "Using the default genomeFile"
		genomeFile = "Input/t/t.fa"
	kmersInGenome,num_kmers_in_genome,sampleDensity, readsPerLine = getInfoFromInputFile(fn,k,genomeFile,numberOfMeasurements=10)
	print "sizeOfGenome:",sizeOfGenome
	print "num_kmers_in_genome:",num_kmers_in_genome
	lineNumber = []
	cov = []
	kmersInGraph = []
	kmersInBF = []
	G_ratio = []
	BF_ratio = []
	count = 0
	noMaxCovFile = "figureData_maxCov_-1.csv"
	if outDirectory=="":
		print "Using default outDirectory inside BF_counter"
		outFolder = "Output/defaultOutFolder"
	else:
		if not os.path.exists(outDirectory):
			raise Exception("The outDirectory doesn't exist. outDirectory="+str(outDirectory))
	if not os.path.exists(outDirectory+"/Without_BF_stopper"):
		if maxCov==-1:
			os.makedirs(outDirectory+"/Without_BF_stopper")
		else:
			raise Exception("We should have already created this folder")
	noMaxCovFilePath = outDirectory+"/Without_BF_stopper/"+noMaxCovFile
	if maxCov==-1:
		outFolder = outDirectory+"/Without_BF_stopper"
	else:
		outFolder = outDirectory+"/BF_stopper"
		if not os.path.exists(outFolder):
			os.makedirs(outFolder)
		assert os.path.isfile(noMaxCovFilePath),"The file with data from maxCov=-1 doesn't exist"

	genomeName = os.path.dirname(genomeFile).split("/")[-1]
	print "genomeName:",genomeName
	for f in fn:
		h = open(f, "rU")
		for lineNr,line in enumerate(h,start=startAtLine):
			if (lineNr%4 == 1):
				segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
				#if len(segments)>1:
				#	print "lineNr="+str(lineNr)+", segments:",segments, "len(segments)="+str(len(segments))
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
						if printProgress:
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
		printFigureFile(lineNumber,cov,kmersInGraph,kmersInBF,G_ratio,BF_ratio,outFolder,outFile=noMaxCovFile)
		#inputFile=outFolder+"/figureData.csv"
		titles=[str(genomeName)+". No maximum coverage","","","",""]
		helpers.createFigureFromFile(noMaxCovFilePath,sizeOfGenome,titles,outFolder,maxCov,genomeName)
		#helpers.createFigureFromFile2("Output/defaultOutFolder/Without_BF_stopper/figureData2.csv","Ratio of k-mers in G and BF",outFolder)
	else:
		#Assume that the figureData.csv and figureData2.csv have
		#already been created without using a BF-stopper
		#i.e. there should be 4 columns in the file, not 6
		#before we call changeFile
		#Til minnis til að þurfa ekki að skrolla upp:
		#	noMaxCovFilePath = outDirectory+"/Without_BF_stopper/"+noMaxCovFile
		#	outFolder = outDirectory+"/BF_stopper"
		newCols = [kmersInGraph,kmersInBF,G_ratio,BF_ratio]
		newFileName = "figureData_maxCov_"+str(maxCov)+".csv"
		helpers.changeFile(noMaxCovFilePath,newCols,outFolder,newFileName)
		titles=[str(genomeName)+". maxCov="+str(maxCov),"","","",""]
		newFilePath = outFolder+"/"+newFileName
		print "newFilePath="+str(newFilePath)
		helpers.createFigureFromFile(newFilePath,sizeOfGenome,titles,outFolder,maxCov,genomeName)
		#helpers.createFigureFromFile2("Output/defaultOutFolder/BF_stopper/figureData2.csv","Ratio of k-mers in G and BF with BF-stopper at cov="+str(maxCov),outFolder)

#fn:			A list of 1 or more .fastq files
#					fn = [file1.fastq, file2.fastq, ...]
#BF:			Bloom filter
#k:				kmer length
#G:				Graph
#startAtLine:	The first line in the input file we read from 
#				Useful when we have already read the earlier lines from the
#				file and want to continue from there instead of starting over
#def BF_counter(fn,BF,k,G,pfn=False,printProgress=False,startAtLine=0,sizeOfGenome=-1,maxCov=-1):
	         #geraAllt(fn,k,BF,G,pfn,      printProgress,      startAtLine,  sizeOfGenome,   maxCov,   whatToRun=1)
	       #BF_counter(fn,k,BF,G,pfn=False,printProgress=False,startAtLine=0,sizeOfGenome=-1,maxCov=-1,whatToRun=1)
#BF_counter.BF_counter(fn,k,BF,G,pfn,      printProgress,      startAtLine,  sizeOfGenome,   maxCov)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def getInfoFromInputFile(fn,k,genomeFile,numberOfMeasurements):
	f = fn[0]
	numberOfLines  = file_len(f)
	numberOfReads = int(numberOfLines/4)
	h = open(f, "rU")
	line = h.readline()
	line = h.readline()
	readsPerLine = len(line)
	h.close()
	if numberOfMeasurements==-1:
		numberOfMeasurements = 10
	sampleDensity = numberOfReads/(numberOfMeasurements-1)


	#Now create kmersInGenome
	genomeFileExtension = os.path.splitext(genomeFile)[1]
	kmersInGenome = collections.defaultdict(int)
	g = open(genomeFile, "rU")
	if genomeFileExtension==".fa":
		line = g.readline()
		line = g.readline().strip()
		for km in dbg.kmers(line,k):
			kmersInGenome[km] += 1
		for km in dbg.kmers(dbg.twin(line),k):
			kmersInGenome[km] += 1
	elif genomeFileExtension==".fasta":
		line = g.readline()
		for line in g:
			line = line.strip()
			line = line if all([c in alphabet for c in line]) else ""
			tl = dbg.twin(line)
			for km in dbg.kmers(line,k):
				kmersInGenome[km] += 1
			for km in dbg.kmers(tl,k):
				kmersInGenome[km] += 1
	else:
		raise Exception("genomeFile does not have a legal extension. genomeFile="+str(genomeFile))
	num_kmers_in_genome = len(kmersInGenome)
	g.close()
	return kmersInGenome,num_kmers_in_genome,sampleDensity, readsPerLine

def printFigureFile(lineNumber,cov,kmersInGraph,kmersInBF,G_ratio,BF_ratio,outFolder="",outFile=""):
	print "printFigureFile"
	if outFolder=="":
		print "Using default outFolder inside printFigureFile"
		outFolder = "Output/defaultOutFolder/Without_BF_stopper"
	if outFile=="":
		print "Using default outFile inside printFigureFile"
		outFile = "figureData.csv"
	h1 = open(outFolder+"/"+outFile,"w")

	for i in xrange(0,len(cov)):
		h1.write(str(lineNumber[i])+","+str(cov[i])+","+str(kmersInGraph[i])+","+str(kmersInBF[i])+","+str(round(G_ratio[i],4))+","+str(round(BF_ratio[i],4))+"\n")
	h1.close()

	"""
	if outFolder=="":
		h2 = open("Output/defaultOutFolder/Without_BF_stopper/figureData2.csv","w")
	else:
		h2 = open(outFolder+"/figureData2.csv","w")
	for i in xrange(0,len(cov)):
		h2.write(str(lineNumber[i])+","+str(cov[i])+","+str(round(G_ratio[i],4))+","+str(round(BF_ratio[i],4))+"\n")
	h2.close()
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