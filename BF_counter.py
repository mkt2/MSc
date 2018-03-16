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

alphabet = ["A","C","G","T","N"]

class infoKeeper:
	def __init__(self,outDirectory,preprocFiles,kmersInGenome=-1):
		print "Initializing the infoKeeper"
		#Stuff we need to keep track of for the pictures:
		self.cov = []
		self.num_kmers_InGraph_noMaxCov = []
		self.num_kmers_InBF_noMaxCov = []
		self.G_ratio_noMaxCov = []
		self.BF_ratio_noMaxCov = []
		
		#Other stuff we need to initialize:
		self.outDirectory = outDirectory
		#self.noMaxCovFile = "figureData_maxCov_-1.csv"
		#self.noMaxCovFilePath = self.outDirectory+"/"+self.noMaxCovFile
		#self.genomeName = os.path.dirname(genomeFile).split("/")[-1]
		self.maxCov = -1	#max coverage for adding to bloom filter
		self.maxCov2 = -1	#max coverage for splitting

		#Stuff we read directly from the genome info file
		genomeInfo = helpers.readGenomeInfoFromFile(preprocFiles[0])
		self.genomeName = genomeInfo[0]
		self.k = int(genomeInfo[1])
		self.numKmersInGenome = int(genomeInfo[2])
		self.numReadsPerFile = int(genomeInfo[3])
		self.numKmersPerRead = int(genomeInfo[4])
		self.numKmersInReads = int(genomeInfo[5])
		self.numKmersInReads_twice = int(genomeInfo[6])
		self.numKmersInReads_twice_BF = int(genomeInfo[7])

		#A dict storing every k-mer in the genome:
		if kmersInGenome==-1:
			self.kmersInGenome = helpers.readKmersFromFileToDict(preprocFiles[1])
		else:
			self.kmersInGenome = kmersInGenome

		#Stuff we can calculate from the info initialized above:
		self.num_kmers_in_genome_counted = self.numKmersInGenome
		numberOfMeasurements = 10
		self.sampleDensity = self.numReadsPerFile/(numberOfMeasurements-1)
		
		#Create the output directory if it didn't already exist
		if not os.path.exists(self.outDirectory):
			os.makedirs(self.outDirectory)
	
	#When we have gathered all information about the current maxCov value
	#we call this function to update the maxCov and get ready to find all
	#information about the new maxCov value
	def setMaxCov(self,newMaxCov):
		self.maxCov = newMaxCov
		self.num_kmers_InGraph_maxCov = []
		self.G_ratio_maxCov = []

	#Add the next values about our current maxCov to the lists of 
	#things we need to store for the figures
	def append(self,G,BF,cov):
		g_fraction = helpers.fractionFromGenomeInObject(self.kmersInGenome,self.numKmersInGenome,Object=G,isGraph=1)
		#g,bf = self.ratioInGenome(G,BF)
		if self.maxCov==-1:
			bf_fraction = helpers.fractionFromGenomeInObject(self.kmersInGenome,self.numKmersInGenome,Object=BF,isGraph=0)
			self.cov.append(cov)
			self.num_kmers_InGraph_noMaxCov.append(G.numKmerPairs())
			self.num_kmers_InBF_noMaxCov.append(len(BF))
			self.G_ratio_noMaxCov.append(g_fraction)
			self.BF_ratio_noMaxCov.append(bf_fraction)
		else:
			self.num_kmers_InGraph_maxCov.append(G.numKmerPairs())
			self.G_ratio_maxCov.append(g_fraction)

	def createFigure(self):
		def toLogScale(myList,ratio_1_log_scale):
			return [ratio_1_log_scale if x==1 else -log(1-x,10) for x in myList]

		#C: Coverage
		#N: Number of k-mers in the genome
		def lw(C,N):
			out = (1-exp(-C))#*N
			assert out>=0
			assert out<=1
			return out

		#print "IK.createFigure()"
		#Initialize the plot:
		fig, (ax1,ax2) = plt.subplots(2,1)
		#fig1.subplots_adjust(hspace=.5)
		#fig.subplots_adjust(wspace=2)

		#Plot on the subplot to the left
		ax1.axhline(y=self.num_kmers_in_genome_counted,label="Genome")
		ax1.plot(self.cov,self.num_kmers_InBF_noMaxCov,"k-",label="BF. noMaxCov")
		ax1.plot(self.cov,self.num_kmers_InGraph_noMaxCov,"k--",label='G. noMaxCov')
		x1,x2,y1,y2 = ax1.axis()
		ax1.axis((x1,max(self.cov),0,max(self.num_kmers_InGraph_noMaxCov)*1.3))

		#Plot on the subplot to the right:
		ratio_1_log_scale = -log(1-0.9999,10)
		ax2.axhline(-log(1-0.99999,10),label="Fraction of 99.999%")	#XX
		ax2.axhline(y=ratio_1_log_scale,label="Fraction of 99.99%")
		ax2.axhline(-log(1-0.999,10),label="Fraction of 99.9%")	    #XX
		ax2.axhline(-log(1-0.99,10),label="Fraction of 99%")	    #XX
		BF_ratio_noMaxCov_log = toLogScale(self.BF_ratio_noMaxCov,ratio_1_log_scale)
		G_ratio_noMaxCov_log = toLogScale(self.G_ratio_noMaxCov,ratio_1_log_scale)
		ax2.plot(self.cov,BF_ratio_noMaxCov_log,"k-",label="BF. noMaxCov")
		ax2.plot(self.cov,G_ratio_noMaxCov_log,"k--",label='G. noMaxCov')
		x1,x2,y1,y2 = ax2.axis()
		ax2.axis((x1,max(self.cov),y1,y2))
		#Bætum lw við myndirnar:
		G_lw = [lw(C,self.num_kmers_in_genome_counted) for C in self.cov]
		G_lw_log = toLogScale(G_lw,ratio_1_log_scale)
		ax2.plot(self.cov,G_lw_log,"g-",label='G.lw')

		#Plot the red lines (if we're using a stop filter):
		if self.maxCov!=-1:
			ax1.plot(self.cov,self.num_kmers_InGraph_maxCov,"r--",label='G. maxCov='+str(self.maxCov))
			G_ratio_maxCov_log = [ratio_1_log_scale if x==1 else -log(1-x,10) for x in self.G_ratio_maxCov]
			ax2.plot(self.cov,G_ratio_maxCov_log,"r--",label='G. maxCov='+str(self.maxCov))

		#Setjum tölur inn á gröfin við hliðina á línunum:
		ax1.annotate(' %s' % format(self.num_kmers_in_genome_counted,','), xy=(0,0), xytext=(max(self.cov),self.num_kmers_in_genome_counted))
		ax2.annotate(' %s' % format(float('{:.2f}'.format(ratio_1_log_scale,','))), xy=(0,0), xytext=(max(self.cov),ratio_1_log_scale))
		if self.maxCov==-1:
			vars_x1 = [self.num_kmers_InGraph_noMaxCov]
			vars_x2 = [BF_ratio_noMaxCov_log,G_ratio_noMaxCov_log]
		else:
			vars_x1 = [self.num_kmers_InGraph_noMaxCov,self.num_kmers_InGraph_maxCov]
			vars_x2 = [BF_ratio_noMaxCov_log,G_ratio_noMaxCov_log,G_ratio_maxCov_log]
		for var in vars_x1:
			ax1.annotate(' %s' % format(max(var),','), xy=(0,0), xytext=(max(self.cov),max(var)))
		for var in vars_x2:
			ax2.annotate(' %s' % format(float('{:.2f}'.format(max(var),','))), xy=(0,0), xytext=(max(self.cov),max(var)))
		"""
		for var in vars_x1:
			ax1.annotate(' %i' % max(var), xy=(0,0), xytext=(max(self.cov),max(var)))
		for var in vars_x2:
			ax2.annotate(' %0.2f' % max(var), xy=(0,0), xytext=(max(self.cov),max(var)))
		"""

		#Set commas on the y-axis in the above plot:
		ax1.get_yaxis().set_major_formatter(
    	mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

		#Set axis labels and legends:
		ax1.set_title("#k-mers / cov")
		ax2.set_title("-log(1-Fraction) / cov")
		ax1.set_xlabel("Coverage")
		ax2.set_xlabel("Coverage")
		ax1.set_ylabel("#k-mers")
		ax2.set_ylabel("-log(1-Fraction)")
		ax1.grid()
		ax2.grid()
		ax1.legend(loc='upper left', prop={'size': 9})#, bbox_to_anchor=(1.2, -0.05))
		ax2.legend(loc='lower right', prop={'size': 9})#, bbox_to_anchor=(1.1, 0))
		
		plt.tight_layout()
		#Set figure title and save the figure
		if self.maxCov==-1:
			fig.suptitle(self.genomeName+". No maximum coverage", y=1.01)
			fig.savefig(self.outDirectory+"/"+self.genomeName+"_maxCov_"+str(self.maxCov)+".png", bbox_inches='tight')
		else:
			fig.suptitle(self.genomeName+". maxCov="+str(self.maxCov), y=1.01)
			fig.savefig(self.outDirectory+"/"+self.genomeName+"_maxCov_"+str(self.maxCov)+".png", bbox_inches='tight')

def BF_counter(fn,k,BF,G,IK,maxCov,pfn=False,printProgress=False,startAtLine=0,skipPictures=False):
	if pfn:
		print "BF_counter",locals()
	if not isinstance(fn, list):
		print "fn:",fn
		raise Exception('fn has to be a list')
	if maxCov!=-1:
		IK.setMaxCov(maxCov)
	sampleDensity = IK.sampleDensity
	alpha = 0.005
	covFactor = (float(IK.numKmersPerRead) / IK.num_kmers_in_genome_counted) * pow(1-alpha,k)
	cov = 0		#Hversu oft við höfum séð hvern k-mer að meðaltali
	count = 0
	for f in fn:
		h = open(f, "rU")
		for lineNr,line in enumerate(h,start=startAtLine):
			if (lineNr%4 == 1):
				segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
				for s in segments:
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
						#Check whether rep_kmer was in BF. Add it if it wasn't
						if not BF.kmerInBF_also_AddIf_B_is_True(rep_kmer,add_moreToBF):
							#rep_kmer wasn't in BF but has now been added
							if i>start:
								#Add the current goodSequence to G
								goodSequence = s[start:i+k-1]
								G.addSegmentToGraph(goodSequence)
							start = i+1
					#If we reach the end of segment we add the current sequence
					if L-start>=k:
						G.addSegmentToGraph(s[start:])
				count += 1
				cov = count * covFactor
				#cov = N * L / G
		if not skipPictures:
			#Geymum einnig niðurstöðurnar í lokin
			IK.append(G,BF,cov)
	if not skipPictures:
		#Prentum niðurstöðurnar í skrá sem við getum notað síðar til að búa til mynd
		IK.createFigure()

def BF_counter_naive(fn,BF,k,G_naive,pfn=True,printProgress=False):
	if pfn:
		print "BF_counter_naive(fn, BF, k="+str(k)+", G_naive, pfn="+str(pfn)+", printProgress="+str(printProgress)+")"
	kd = collections.defaultdict(int)
	for f in fn:
		h = open(f, "rU")
		for lineNr,line in enumerate(h):
			if (lineNr%4 == 1):
				#segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
				segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
				for s in segments:
					for km in dbg.kmers(s,k):
						rep_km = min(km,dbg.twin(km))
						if not rep_km in BF:
							BF.add(rep_km)
						else:
							kd[rep_km] += 1

	G,cs = dbg.all_contigs(kd,k)
	dbg.createGraphObject(G,cs,k,G_naive,pfn,ps=False)

def selectGenome(genomeName):
	if genomeName=="t":
		maxCovs = [-1, 5, 10, 15, 20, 30]
		fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
		numberOfKmers = 8000000
	elif genomeName=="sa":
		maxCovs = [-1, 5, 10, 15, 20, 30]#, 35, 40, 45, 50]
		fn = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
		numberOfKmers = 100000000
	else:
		raise Exception("The genomeName must be either 't' or 'sa'!")
	return maxCovs, fn, numberOfKmers

#Need to run preprocessInput.py before running this file
if __name__ == "__main__":
	genomeName = sys.argv[1]
	maxCovs, fn, numberOfKmers = selectGenome(genomeName)
	runTimes = [-1]*len(maxCovs)
	pfn = True
	printProgress = False

	#Skilgreinum breytur sem þarf að passa að séu þær sömu í þessu main falli, preproc skránum og IK
	p = 0.01
	outDirectory = "Output/"+genomeName
	assert(pathlib2.Path(outDirectory+"/genome_info.csv").is_file()), "We have to create this file using preprocessInput.py before running BF_counter.py"
	assert(pathlib2.Path(outDirectory+"/kmers_genome.txt").is_file()), "We have to create this file using preprocessInput.py before running BF_counter.py"
	assert(pathlib2.Path(outDirectory+"/kmers_reads.txt").is_file()), "We have to create this file using preprocessInput.py before running BF_counter.py"
	preprocFiles = [outDirectory+"/genome_info.csv",outDirectory+"/kmers_genome.txt",outDirectory+"/kmers_reads.txt"]
	genomeInfo = helpers.readGenomeInfoFromFile(preprocFiles[0])
	k = int(genomeInfo[1])
	kmersInGenome = helpers.readKmersFromFileToDict(preprocFiles[1])
	IK = infoKeeper(outDirectory,preprocFiles,kmersInGenome)
	BF = Bloom.Bloom(p,numberOfKmers,pfn=False)
	
	#Keyrum BF_counter fyrir sérhvert maxCov og tökum tímana
	#Búum einnig til tvær skrár fyrir hvert maxCov:
	#Sleppum skrá 2
	#	Skrá 1: Geymir alla k-mera í G fyrir núverandi maxCov
	#	Skrá 2: Geymir mismengi k-meranna úr erfðamenginu og k-meranna í G
	#Fyrir maxCov=5 verður nöfnin á skránum t.d:
	#	Skrá 1: kmersFromG_maxCov_5.txt
	#	Skrá 2: difference_Genome_G_maxCov_5.txt
	start = time.time()
	for i, maxCov in enumerate(maxCovs):
		start_i = time.time()
		print "\n-------------------------------------------------------"
		print "Starting on maxCov="+str(maxCov)
		print "Initializing an empty bloom filter and Graph"
		BF = Bloom.Bloom(p,numberOfKmers,pfn=False)
		G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=False)

		print "Running BF_counter"
		BF_counter(fn,k,BF,G,IK,maxCov,pfn,printProgress,startAtLine=0,skipPictures=False)
		timeInSeconds = int(time.time()-start_i)
		runTimes[i] = timeInSeconds
		print "G now stores all k-mers occurring twice in the reads (according to the BF and current maxCov)"

		print "Saving G to a file"
		if maxCov==-1:
			G.printToFile(outDirectory+"/G.txt")
		else:
			G.printToFile(outDirectory+"/G"+str(maxCov)+".txt")

		#helpers.printKmerdictToFile(G.kmers,outDirectory+"/kmersFromG_maxCov_"+str(maxCov)+".txt")
		print "Finished maxCov="+str(maxCov)
		print helpers.returnTime(timeInSeconds)
		print "-------------------------------------------------------\n"
	print "Finished everything:"
	print helpers.returnTime(int(time.time()-start))
	helpers.printRuntimesToFile(genomeName, maxCovs, runTimes)
    
