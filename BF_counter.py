#coding:utf8
import collections, sys
import dbg, Bloom, Graph, helpers
import os.path

alphabet = ["A","C","G","T","N"]

class infoKeeper:
	def __init__(self,fn,k,maxCov,sizeOfGenome,outDirectory,genomeFile):
		#outFolder, noMaxCovFile,genomeName,noMaxCovFilePath,
		print "Initializing the infoKeeper"
		self.lineNumber = []
		self.cov = []
		self.kmersInGraph = []
		self.kmersInBF = []
		self.G_ratio = []
		self.BF_ratio = []
		self.readsPerLine = -1
		self.sampleDensity = -1
		self.kmersInGenome = collections.defaultdict(int)
		self.num_kmers_in_genome = -1		#The computed size of the genome
		self.sizeOfGenome = sizeOfGenome	#The sizeOfGenome given as input
											#Why isn't this the same number?
		self.noMaxCovFolder = outDirectory+"/Without_BF_stopper"
		self.noMaxCovFile = "figureData_maxCov_-1.csv"
		self.noMaxCovFilePath = self.noMaxCovFolder+"/"+self.noMaxCovFile
		#Gather info about the genome and create the files where to save the output
		if not os.path.exists(outDirectory):
			raise Exception("The outDirectory doesn't exist. outDirectory="+str(outDirectory))
		if maxCov==-1:
			#self.outFolder = outDirectory+"/Without_BF_stopper"
			if not os.path.exists(self.noMaxCovFolder):
				os.makedirs(self.NoMaxCovFolder)
		else:
			if not os.path.isfile(self.noMaxCovFilePath):
				raise Exception("This file is supposed to store the noMaxCovFile already created when running BF_counter with maxCov=-1")
			self.outFolder = outDirectory+"/BF_stopper"
			if not os.path.exists(outFolder):
				os.makedirs(self.outFolder)
		self.genomeName = os.path.dirname(genomeFile).split("/")[-1]

		#Create values for:
		#	self.readsPerLine
		#	self.sampleDensity
		numberOfMeasurements=10
		f = fn[0]
		numberOfLines  = self.file_len(f)
		numberOfReads = int(numberOfLines/4)
		h = open(f, "rU")
		line = h.readline()
		line = h.readline()
		self.readsPerLine = len(line)
		h.close()
		if numberOfMeasurements==-1:
			numberOfMeasurements = 10
		self.sampleDensity = numberOfReads/(numberOfMeasurements-1)

		#Create values for:
		#	self.kmersInGenome
		#	self.num_kmers_in_genome
		genomeFileExtension = os.path.splitext(genomeFile)[1]
		g = open(genomeFile, "rU")
		if genomeFileExtension==".fa":
			line = g.readline()
			line = g.readline().strip()
			for km in dbg.kmers(line,k):
				self.kmersInGenome[km] += 1
			for km in dbg.kmers(dbg.twin(line),k):
				self.kmersInGenome[km] += 1
		elif genomeFileExtension==".fasta":
			line = g.readline()
			for line in g:
				line = line.strip()
				line = line if all([c in alphabet for c in line]) else ""
				tl = dbg.twin(line)
				for km in dbg.kmers(line,k):
					self.kmersInGenome[km] += 1
				for km in dbg.kmers(tl,k):
					self.kmersInGenome[km] += 1
		else:
			raise Exception("genomeFile does not have a legal extension. genomeFile="+str(genomeFile))
		self.num_kmers_in_genome = len(self.kmersInGenome)
		g.close()

		print "genomeName:",self.genomeName
		print "sizeOfGenome:",self.sizeOfGenome
		print "num_kmers_in_genome:",self.num_kmers_in_genome
	
	def file_len(self,fname):
		with open(fname) as f:
			for i, l in enumerate(f):
				pass
		return i + 1

	def append(self, lineNr,count,G,BF):
		self.lineNumber.append(lineNr)
		self.cov.append(count * self.readsPerLine / self.sizeOfGenome)
		self.kmersInGraph.append(len(G.kmers))
		self.kmersInBF.append(len(BF))
		g,bf = self.ratioInGenome(G,BF)
		self.G_ratio.append(g)
		self.BF_ratio.append(bf)

	#Input:
	#   kmerdict:  A dict storing all kmers actually occurring in the genome. Also stores twins
	#   G:         Our Graph
	#   BF:        Our Bloom filter
	#Returns:
	#   The ratio of kmers in the genome that occur in G
	#   The ratio of kmers in the genome that occur in BF
	def ratioInGenome(self,G,BF):
		#kmerdict = self.kmersInGenome
		G_count = 0
		BF_count = 0
		L = len(self.kmersInGenome)
		for km in self.kmersInGenome:
			if km in G.kmers:
				G_count += 1
			if km in BF:
				BF_count += 1
		G_ratio = float(G_count) / L
		BF_ratio = float(BF_count) / L
		return G_ratio,BF_ratio

	def printResults(self,maxCov):
		if maxCov==-1:
			h1 = open(self.noMaxCovFilePath,"w")
			for i in xrange(0,len(self.cov)):
				h1.write(str(self.lineNumber[i])+","+str(self.cov[i])+","+str(self.kmersInGraph[i])+","+str(self.kmersInBF[i])+","+str(round(self.G_ratio[i],4))+","+str(round(self.BF_ratio[i],4))+"\n")
			h1.close()
			titles=[str(self.genomeName)+". No maximum coverage","","","",""]
			helpers.createFigureFromFile(self.noMaxCovFilePath,self.sizeOfGenome,titles,self.noMaxCovFolder,maxCov,self.genomeName)
		else:
			newFileName = "figureData_maxCov_"+str(maxCov)+".csv"
			newFilePath = self.outFolder+"/"+newFileName
			copyfile(self.noMaxCovFilePath, newFilePath)
			for i,line in enumerate(fileinput.input(newFile, inplace=True)):
				assert line.count(",")==5
				print "%s,%s,%s,%s,%s" % (line.strip(),str(self.kmersInGraph[i]),str(self.kmersInBF[i]),str(self.G_ratio[i]),str(self.BF_ratio[i])+"\n"),
			titles=[str(self.genomeName)+". maxCov="+str(maxCov),"","","",""]
			helpers.createFigureFromFile(newFilePath,self.sizeOfGenome,titles,self.outFolder,maxCov,self.genomeName)

#dasdf
#Inputs:
#	k: kmerLength
#	fn: a list of fastq files containing reads
#returns all substrings from segments whose kmers we've seen 2 or more times
#(with a false positive rate of p)
#whatToRun=0: Keyra BF_counter_naive
#whatToRun=1: Keyra BF_counter
#def BF_counter(fn,k,BF,G,pfn=False,printProgress=False,startAtLine=0,sizeOfGenome=-1,maxCov=-1,genomeFile="",outDirectory=""):
def BF_counter(fn,k,BF,G,IK,maxCov,pfn=False,printProgress=False,startAtLine=0):
	if pfn:
		#print "BF_counter",locals().keys(),"\n"
		print "BF_counter",locals()
		#print "geraAllt(fn="+str(fn)+", k="+str(k)+", BF, G, pfn="+str(pfn)+", printProgress="+str(printProgress)+", startAtLine="+str(startAtLine)+")"
	if not isinstance(fn, list):
		print "fn:",fn
		raise Exception('fn has to be a list')

	#genomeFile, noMaxCovFile, outDirectory, maxCov, outFolder,noMaxCovFilePath
	#The actual work is done here:
	#IK = infoKeeper(fn,k,genomeFile,sizeOfGenome)	#A class object to keep track of values needed for our figures
	sampleDensity = IK.sampleDensity
	#maxCov = IK.maxCov
	count = 0
	for f in fn:
		h = open(f, "rU")
		for lineNr,line in enumerate(h,start=startAtLine):
			if (lineNr%4 == 1):
				segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
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
						IK.append(lineNr,count,G,BF)

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
		#Geymum einnig niðurstöðurnar í lokin
		IK.append(lineNr,count,G,BF)
	
	#Prentum niðurstöðurnar í skrá sem við getum notað síðar til að búa til mynd
	IK.printResults(maxCov)


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