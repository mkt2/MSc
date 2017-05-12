#coding:utf8
import collections, sys
import dbg, Bloom, Graph, helpers
import os.path
from shutil import copyfile
import fileinput
import numpy as np
import matplotlib.pyplot as plt
from math import log
import time

alphabet = ["A","C","G","T","N"]

class infoKeeper:
	def __init__(self,fn,k,outDirectory,genomeFile):
		print "Initializing the infoKeeper"
		self.cov = []
		self.num_kmers_InGraph_noMaxCov = []
		self.num_kmers_InBF_noMaxCov = []
		self.G_ratio_noMaxCov = []
		self.BF_ratio_noMaxCov = []
		self.num_kmersPerLine = -1
		self.sampleDensity = -1
		self.kmersInGenome = collections.defaultdict(int)
		self.num_kmers_in_genome_counted = -1		#The computed size of the genome
											#Why isn't this the same number?
		self.outDirectory = outDirectory
		self.noMaxCovFile = "figureData_maxCov_-1.csv"
		self.noMaxCovFilePath = self.outDirectory+"/"+self.noMaxCovFile
		self.genomeName = os.path.dirname(genomeFile).split("/")[-1]
		self.maxCov=-1
		
		#Create the files where to save the output
		if not os.path.exists(self.outDirectory):
			os.makedirs(self.outDirectory)

		#Create values for:
		#	self.num_kmersPerLine
		#	self.sampleDensity
		numberOfMeasurements=10
		f = fn[0]
		numberOfLines  = helpers.file_len(f)
		numberOfReads = int(numberOfLines/4)
		h = open(f, "rU")
		line = h.readline()
		line = h.readline()
		self.num_kmersPerLine = len(list(dbg.kmers(line,k)))
		h.close()
		if numberOfMeasurements==-1:
			numberOfMeasurements = 10
		self.sampleDensity = numberOfReads/(numberOfMeasurements-1)

		#Create values for:
		#	self.kmersInGenome
		#	self.num_kmers_in_genome_counted
		genomeFileExtension = os.path.splitext(genomeFile)[1]
		g = open(genomeFile, "rU")
		if genomeFileExtension==".fa":
			line = g.readline()
			line = g.readline().strip()
			for km in dbg.kmers(line,k):
				self.kmersInGenome[min(km,dbg.twin(km))] += 1
		elif genomeFileExtension==".fasta":
			line = g.readline()
			for line in g:
				line = line.strip()
				line = line if all([c in alphabet for c in line]) else ""
				tl = dbg.twin(line)
				for km in dbg.kmers(line,k):
					self.kmersInGenome[min(km,dbg.twin(km))] += 1
		else:
			raise Exception("genomeFile does not have a legal extension. genomeFile="+str(genomeFile))
		self.num_kmers_in_genome_counted = len(self.kmersInGenome)
		g.close()

		print "genomeName:",self.genomeName
		print "length_of_t:                   70000"
		print "Number of bp in S_aureus:    2903081"
		print "num_kmers_in_genome_counted:",self.num_kmers_in_genome_counted
	
	def setMaxCov(self,newMaxCov):
		self.maxCov = maxCov
		#self.maxCovFile = "figureData_maxCov_"+str(maxCov)+".csv"
		#self.maxCovFilePath = self.outDirectory+"/"+self.maxCovFile
		#if not os.path.isfile(self.noMaxCovFilePath):
		#	raise Exception("This file is supposed to store the noMaxCovFile already created when running BF_counter with maxCov=-1")
		self.num_kmers_InGraph_maxCov = []
		self.G_ratio_maxCov = []

	def append(self,G,BF,cov):
		#print "append(len(G)="+str(len(G))+", len(BF)="+str(len(BF))+", cov="+str(cov)+")"
		g,bf = self.ratioInGenome(G,BF)
		if self.maxCov==-1:
			self.cov.append(cov)
			self.num_kmers_InGraph_noMaxCov.append(G.numKmerPairs())
			self.num_kmers_InBF_noMaxCov.append(len(BF))
			self.G_ratio_noMaxCov.append(g)
			self.BF_ratio_noMaxCov.append(bf)
		else:
			self.num_kmers_InGraph_maxCov.append(G.numKmerPairs())
			self.G_ratio_maxCov.append(g)
		

	#Input:
	#   kmerdict:  A dict storing all kmers actually occurring in the genome. Also stores twins
	#   G:         Our Graph
	#   BF:        Our Bloom filter
	#Returns:
	#   The ratio of kmers in the genome that occur in G
	#   The ratio of kmers in the genome that occur in BF
	def ratioInGenome(self,G,BF):
		G_count = 0
		BF_count = 0
		for km in self.kmersInGenome:
			if km in G.kmers:
				G_count += 1
			if km in BF:
				BF_count += 1
		G_ratio = float(G_count) / self.num_kmers_in_genome_counted
		BF_ratio = float(BF_count) / self.num_kmers_in_genome_counted
		assert G_ratio>=0 and G_ratio<=1, "A ratio must be between 0 and 1"
		assert BF_rati>=0 and BF_ratio<=1, "A ratio must be between 0 and 1"
		return G_ratio,BF_ratio

	def printResults(self):
		self.createFigure()
		"""
		if self.maxCov==-1:
			h1 = open(self.noMaxCovFilePath,"w")
			for i in xrange(0,len(self.cov)):
				h1.write(str(self.lineNumber[i])+","+str(self.cov[i])+","+str(self.kmersInGraph_noMaxCov[i])+","+str(self.kmersInBF[i])+","+str(round(self.G_ratio_noMaxCov[i],4))+","+str(round(self.BF_ratio[i],4))+"\n")
			h1.close()
			self.createFigure()
		else:
			copyfile(self.noMaxCovFilePath, self.maxCovFilePath)
			for i,line in enumerate(fileinput.input(self.noMaxCovFilePath, inplace=True)):
				assert line.count(",")==5
				print "%s,%s,%s,%s,%s" % (line.strip(),str(self.kmersInGraph_maxCov[i]),str(self.kmersInBF[i]),str(self.G_ratio_maxCov[i]),str(self.BF_ratio[i])+"\n"),
			self.createFigure(maxCov)
		"""

	def createFigure(self):
		#print "IK.createFigure()"
		#Initialize the plot:
		fig, (ax1,ax2) = plt.subplots(1,2)
		#fig1.subplots_adjust(hspace=.5)
		fig.subplots_adjust(wspace=.5)

		#Plot on the subplot to the left
		ax1.axhline(y=self.num_kmers_in_genome_counted,label="Genome")
		ax1.plot(self.cov,self.num_kmers_InBF_noMaxCov,"k-",label="BF. noMaxCov")
		ax1.plot(self.cov,self.num_kmers_InGraph_noMaxCov,"k--",label='G. noMaxCov')
		x1,x2,y1,y2 = ax1.axis()
		ax1.axis((x1,max(self.cov),0,max(self.num_kmers_InGraph_noMaxCov)*1.3))

		#Plot on the subplot to the right:
		ratio_1_log_scale = -log(1-0.9999)
		ax2.axhline(y=ratio_1_log_scale,label="99.99% ratio")
		BF_ratio_noMaxCov_log = [ratio_1_log_scale if x==1 else -log(1-x) for x in self.BF_ratio_noMaxCov]
		G_ratio_noMaxCov_log = [ratio_1_log_scale if x==1 else -log(1-x) for x in self.G_ratio_noMaxCov]
		ax2.plot(self.cov,BF_ratio_noMaxCov_log,"k-",label="BF. noMaxCov")
		ax2.plot(self.cov,G_ratio_noMaxCov_log,"k--",label='G. noMaxCov')
		x1,x2,y1,y2 = ax2.axis()
		ax2.axis((x1,max(self.cov),y1,y2))

		#Plot the red lines (if we're using a stop filter):
		if self.maxCov!=-1:
			ax1.plot(self.cov,self.num_kmers_InGraph_maxCov,"r--",label='G. maxCov='+str(maxCov))
			G_ratio_maxCov_log = [ratio_1_log_scale if x==1 else -log(1-x) for x in self.G_ratio_maxCov]
			ax2.plot(self.cov,G_ratio_maxCov_log,"r--",label='G. maxCov='+str(maxCov))

		ax1.annotate(' %i' % (self.num_kmers_in_genome_counted), xy=(0,0), xytext=(max(self.cov),self.num_kmers_in_genome_counted))
		ax2.annotate(' %0.2f' % ratio_1_log_scale, xy=(0,0), xytext=(max(self.cov),ratio_1_log_scale))
		if self.maxCov==-1:
			vars_x1 = [self.num_kmers_InGraph_noMaxCov]
			vars_x2 = [BF_ratio_noMaxCov_log,G_ratio_noMaxCov_log]
		else:
			vars_x1 = [self.num_kmers_InGraph_noMaxCov,self.num_kmers_InGraph_maxCov]
			vars_x2 = [BF_ratio_noMaxCov_log,G_ratio_noMaxCov_log,G_ratio_maxCov_log]
		for var in vars_x1:
			ax1.annotate(' %i' % max(var), xy=(0,0), xytext=(max(self.cov),max(var)))
		for var in vars_x2:
			ax2.annotate(' %0.2f' % max(var), xy=(0,0), xytext=(max(self.cov),max(var)))


		#Set axis labels and legends:
		ax1.set_title("#k-mers / cov")
		ax2.set_title("-log(1-Ratio) / cov")
		ax1.set_xlabel("Coverage")
		ax2.set_xlabel("Coverage")
		ax1.set_ylabel("#k-mers")
		ax2.set_ylabel("-log(1-Ratio)")
		ax1.grid()
		ax2.grid()
		ax1.legend(loc='upper left')#, bbox_to_anchor=(1.2, -0.05))
		ax2.legend(loc='lower right')#, bbox_to_anchor=(1.1, 0))
		
		#Set figure title and save the figure
		if self.maxCov==-1:
			fig.suptitle(self.genomeName+". No maximum coverage")
			fig.savefig(self.outDirectory+"/"+self.genomeName+"_maxCov_"+str(maxCov)+".png", bbox_inches='tight')
		else:
			fig.suptitle(self.genomeName+". maxCov="+str(maxCov))
			fig.savefig(self.outDirectory+"/"+self.genomeName+"_maxCov_"+str(maxCov)+".png", bbox_inches='tight')

#dasdf
#Inputs:
#	k: kmerLength
#	fn: a list of fastq files containing reads
#returns all substrings from segments whose kmers we've seen 2 or more times
#(with a false positive rate of p)
#whatToRun=0: Keyra BF_counter_naive
#whatToRun=1: Keyra BF_counter
def BF_counter(fn,k,BF,G,IK,maxCov,pfn=False,printProgress=False,startAtLine=0,skipPictures=False):
	if pfn:
		#print "BF_counter",locals().keys(),"\n"
		print "BF_counter",locals()
		#print "geraAllt(fn="+str(fn)+", k="+str(k)+", BF, G, pfn="+str(pfn)+", printProgress="+str(printProgress)+", startAtLine="+str(startAtLine)+")"
	if not isinstance(fn, list):
		print "fn:",fn
		raise Exception('fn has to be a list')
	if maxCov!=-1:
		IK.setMaxCov(maxCov)
	sampleDensity = IK.sampleDensity
	covFactor = float(IK.num_kmersPerLine) / IK.num_kmers_in_genome_counted
	cov = 0
	count = 0
	kd_twice_BF = collections.defaultdict(int)		#A dict storing every kmer that we see twice according to the BF
	kd_twice = collections.defaultdict(int)			#A dict storing every kmer we actually see twice
	#												len(kd_twice_BF)<=len(kd_twice)*1.01 should be true. Otherwise the BF code is not working
	for f in fn:
		h = open(f, "rU")
		for lineNr,line in enumerate(h,start=startAtLine):
			if (lineNr%4 == 1):
				segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
				for s in segments:
					assert isinstance(s, basestring)

					#Prentum út stöðuna annað slagið:
					if printProgress and (lineNr%50000==1):
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
					if (not skipPictures) and (count%sampleDensity==0):
						#if printProgress:
						#	print "Taking a sample. LineNr="+str(lineNr)+", count="+str(count)
						IK.append(G,BF,cov)

					#Bætum öllum k-merum úr s við BF
					#Bætum öllum k-merum úr s sem við höfum séð tvisvar eða 
					#oftar við G (bætum segmenti sem inniheldur k-mera)
					#B1=True means we're going to add the kmers from s to the BF
					#if they hadn't already been added
					if maxCov==-1:
						add_moreToBF = True
					else:
						add_moreToBF = (cov<=maxCov)
					L = len(s)
					start = 0
					for i, kmer in enumerate(dbg.kmers(s,k)):
						rep_kmer = min(kmer,dbg.twin(kmer))
						if add_moreToBF:
							kd_twice[rep_kmer] += 1
						#if not BF.kmerInBF_also_AddIf_B_is_True(kmer, tkmer, B1):
						if not rep_kmer in BF:
							if add_moreToBF:
								BF.add(rep_kmer)
							if i-start>0:
								goodSequence = s[start:i+k-1]
								G.addSegmentToGraph(goodSequence)
							start = i+1
						else:
							if add_moreToBF:
								assert rep_kmer in BF
								kd_twice_BF[rep_kmer] += 1
					#If we reach the end of segment we add the current sequence
					if L-start>=k:
						assert len(s[start:])>=k
						G.addSegmentToGraph(s[start:])
				count += 1
				cov = count * covFactor
		if not skipPictures:
			#Geymum einnig niðurstöðurnar í lokin
			IK.append(G,BF,cov)
	if not skipPictures:
		#Prentum niðurstöðurnar í skrá sem við getum notað síðar til að búa til mynd
		IK.printResults()
	if add_moreToBF:
		print "\nThere are "+str(len(kd_twice_BF))+" kmers we saw >=2 times according to the BF"
		print "We raise assertion if there aren't equally many kmers in G"
		assert len(kd_twice_BF)==G.numKmerPairs()
		print "There are "+str(len(kd_twice))+" kmers actually in the reads"
		temp = [x for x in kd_twice if kd_twice[x] <= 1]
		for x in temp:
			del kd_twice[x]
		print "There are "+str(len(kd_twice))+" kmers we actually saw >=2 times"
		for km in kd_twice:
			assert km in BF, "Every kmer we actually saw twice should be in the BF (otherwise we have a False Negative)"
		print "len(kd_twice_BF)/len(kd_twice)="+str(float(len(kd_twice_BF))/len(kd_twice))+". It should be less than p="+str(p)
		print "We raise assertion if len(kd)<=len(kd_twice)*1.01 isn't true"
		#assert len(kd_twice_BF)<=len(kd_twice)*1.01, "len(kd_twice_BF)/len(kd_twice)="+str(float(len(kd_twice_BF))/len(kd_twice))+" but should be less than p="+str(p)

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
	#python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 $maxCov True False Input/t/t.fa Output/t
	start = time.time()
	k = 31
	maxCovs = [-1, 5, 10, 15, 20, 30, 40, 50]
	fn = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
	genomeFile = "Input/Staphylococcus_aureus/genome.fasta"
	numberOfKmers = 100000000
	fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
	genomeFile = "Input/t/t.fa"
	numberOfKmers = 50000000
					#var áður 6000000
	genomeName = os.path.dirname(genomeFile).split("/")[-1]
	outDirectory = "Output/"+os.path.dirname(genomeFile).split("/")[-1]
	IK = infoKeeper(fn,k,outDirectory,genomeFile)
	pfn = True
	printProgress = True
	start_i = time.time()
	p = 0.01
	for maxCov in maxCovs:
		start_i = time.time()
		print "\n-------------------------------------------------------"
		print "Starting on maxCov="+str(maxCov)
		BF = Bloom.Bloom(p,numberOfKmers,pfn=True)
		G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=False)
		assert len(G)==0
		assert len(G.kmers)==0
		BF_counter(fn,k,BF,G,IK,maxCov,pfn,printProgress,startAtLine=0,skipPictures=False)
		end_i = time.time()
		print "Finished maxCov="+str(maxCov)
		print helpers.returnTime(int(end_i-start_i))
		print "-------------------------------------------------------\n"
	print "Finished everything:"
	print helpers.returnTime(int(end_i-start_i))
    
