#coding:utf8
import collections, sys
import dbg
import compareGraphs
from math import ceil
from sets import Set
import helpers
import Path

#An object of this class represents a legal DBG graph
class Graph:
	def __init__(self,k,pfn=False,ps=False,al=True,pil=False,printInit=False):
		assert(isinstance(k,int))
		if printInit:
			print "Initializing the Graph: k="+str(k)+", pfn="+str(pfn)+", ps="+str(ps)+", al="+str(al)+", pil="+str(pil)+", printInit="+str(printInit)
		self.k = k 		#the kmer length
		self.ID = 0		#the lowest available contig ID
		
		#dictionary of the contigs in our graph with information about which contigs they connect to:
		self.contigs = collections.defaultdict(list)
		#contigID -> [contig, IN, OUT, COV]
		#IN:  List of tuples. Example: IN = [(0,True),(1,False),(99,True)]
		#OUT: List of tuples
		#COV: The total number of occurrences of the k-mers in the contig
		#	  Note that when a k-mer is added for the first time we say it has occurred twice
		#	  because it has to be seen twice to be added in the first place
		#	  Example: Say we have a contig C of length k+3 (therefore containing 3 k-mers)
		#	  Say the first k-mer has been added once, the second 10 times and the third two times
		#	  COV will therefore be 2+11+3=16
		#	  Note that we don't store the coverage of individual k-mers so when we do operations such
		#	  as splitting contigs, then COV will become an estimation

		#dictionary storing the kmers in our graph:
		self.kmers = collections.defaultdict(list)
		#kmer -> [contigID, index, B]
		#contigID is the ID of contig C with ID=contigID
		#B=False:
		#	kmer is in the twin of C
		#	index represents its location in twin(C)
		#B=True:
		#	kmer is in C
		#	index represents its location in C
		
		#Variables only used for debugging
		self.printFunctionNames = pfn
		self.printStatus = ps
		self.assertLegal = al
		self.printIsLegal = pil

		#Things for maxSplitCov:
		#self.carefulSplit = False		#When self.carefulSplit==True then we only split a contig when we have been told to do so before
		self.splitDict = collections.defaultdict(int)
		#(c_ID, i) -> count
		#c_ID is the ID of contig C
		#i is the index in the contig
		#count stores how often we have been told before that we need to split C at index i

	def __len__(self):
		return len(self.contigs)

	#Skilar hvort km sé í grafinu
	def __contains__(self,km):
		return km in self.kmers

	#--------------------------------------------------------------------------
	#--------------------------functions for printing--------------------------
	#--------------------------------------------------------------------------
	def __str__(self):
		s = "k="+str(self.k)+", ID="+str(self.ID)+", numberOfContigs="+str(len(self.contigs))+"\n"
		for cID in self.contigs:
			L = self.contigs[cID]
			s = s + "cID="+str(cID)+", Length="+str(len(L[0]))
			s = s + "\n   contig: "+str(L[0]) \
			+ "\n   twin:   "+str(dbg.twin(L[0])) \
			+ "\n   IN:     "+str(L[1]) \
			+ "\n   OUT:    "+str(L[2]) \
			+ "\n   COV:    "+str(L[3])+"\n"
		return s

	def printKmers(self,printText=-1):
		if printText!=-1:
			print printText + ":"
		for L in self.kmers:
			print L,self.kmers[L]
		print "\n"

	def printContigs(self, graphName=-1):
		print "\n"
		if graphName!=-1:
			print graphName + ":"
		print "k="+str(self.k)+", ID="+str(self.ID)+", numberOfContigs="+str(len(self.contigs))
		for cID in self.contigs:
			L = self.contigs[cID]
			print "cID="+str(cID)+", Length="+str(len(L[0]))
			print "   contig:", str(L[0])
			print "   twin:  ", str(dbg.twin(L[0]))
			print "   IN:    ", L[1]
			print "   OUT:   ", L[2]
			print "   COV:   ", L[3]
		print "\n"

	def printContigsWithRatings(self,theContigs=-1):
		if theContigs==-1:
			theContigs = self.contigs
		for c_ID in theContigs:
			[c,c_IN,c_OUT,c_COV] = theContigs[c_ID]
			c_inD, c_outD, R = self.degrees[c_ID]
			print "["+str(c_ID)+", "+c+", "+str(c_IN)+", "+str(c_OUT)+", "+str(c_COV)+"]",c_inD,c_outD,R,helpers.cLen(c,self.k)


	#--------------------------------------------------------------------------
	#--------------------functions for storing graphs in files-----------------
	#--------------------------------------------------------------------------
	def printToFile(self,fileName):
		f2 = open(fileName, 'w')
		f2.write(str(self.k)+"\n")
		f2.write(str(self.ID)+"\n")
		f2.write(str(len(self))+"\n")
		f2.write("ID;CONTIG;IN;OUT;COV\n")
		for ID, [c,IN,OUT,COV] in self.contigs.iteritems():	
			f2.write(str(ID)+";"+str(c)+";"+str(IN)+";"+str(OUT)+";"+str(COV)+"\n")
		f2.close()
	
	def createGraphFromFile(self,fileName):
		#Before:	The file is of the correct form
		#				0 k
		#				1 ID
		#				2 numberOfContigs
		#				3 ID;CONTIG;IN;OUT;COV
		#				...data of the form ID;CONTIG;IN;OUT;COV
		#After:
		#	for all lines we set self.contigs[ID] as [CONTIG,IN,OUT,COV]
		#print "createGraphFromFile(fileName="+str(fileName)+")"
		#Read header lines and assert everything is as expected:
		assert(self.isEmpty()), "This function only works on an empty Graph"
		f = open(fileName, 'r')
		k = int(f.readline())
		ID = int(f.readline())
		numberOfContigs = int(f.readline())
		h = f.readline()
		assert h=="ID;CONTIG;IN;OUT;COV\n"
		assert isinstance(k, int)
		assert isinstance(ID, int)
		assert isinstance(numberOfContigs, int)
		assert k==self.k

		#Read the contigs:
		for line in f:
			[ID,CONTIG,IN,OUT,COV] = line.strip().split(";")
			ID = int(ID)
			COV = int(COV)
			IN = eval(IN)
			OUT = eval(OUT)
			assert isinstance(ID, int)
			assert isinstance(CONTIG, str)
			assert isinstance(IN, list)
			assert isinstance(OUT, list)
			assert isinstance(COV, int)
			self.contigs[ID] = [CONTIG,IN,OUT,COV]
		f.close()
		self.addKmersFromAllContigs()
		self.ID = ID+1
		assert len(self)==numberOfContigs
		if self.assertLegal:
			if not self.isLegalDBG():
				#self.printContigs()
				raise Exception("The graph we read from the file is not legal!")
		#print "Done creating the Graph from the file "+fileName

	def saveAs_GFA_toFile(self,fileName):
		def edge(a_ID,S1,b_ID,S2):
			return "L\t"+str(a_ID)+"\t"+S1+"\t"+str(b_ID)+"\t"+S2+"\t"+str(self.k-1)+"M\n"

		f = open(fileName, 'w')

		#Skrifa haus:
		f.write("H\tVN:Z:1.0\n")

		#Skrifa contigana:
		for c_ID, [c,c_IN,c_OUT,c_COV] in self.contigs.iteritems():
			f.write("S\t%d\t%s\t*\n"%(c_ID,c))

		#Skrifa tengingarnar: (c->d)
		#edgeDict = {}
		edgeDict = collections.defaultdict(list)
		#13->twin(139)
		for c_ID, [c,c_IN,c_OUT,c_COV] in self.contigs.iteritems():
			for (d_ID,d_B) in c_OUT:
				if d_B:
					#C->D og twin(C)->twin(D)
					edgeDict[c_ID].append(edge(c_ID,"+",d_ID,"+"))
					edgeDict[d_ID].append(edge(d_ID,"-",c_ID,"-"))
				else:
					#C->twin(D) og D->twin(C)
					#Got to make sure we don't add the connection twice
					e1 = edge(c_ID,"+",d_ID,"-")
					e2 = edge(d_ID,"+",c_ID,"-")
					if not e1 in edgeDict[c_ID]:
						edgeDict[c_ID].append(e1)
						edgeDict[d_ID].append(e2)
			for (d_ID,d_B) in c_IN:
				if d_B:
					#D->C og twin(C)->twin(D)
					pass	#this has already been added
				else:
					#twin(D)->C og twin(C)->D
					e1 = edge(d_ID,"-",c_ID,"+")
					e2 = edge(c_ID,"-",d_ID,"+")
					if not e1 in edgeDict[c_ID]:
						edgeDict[d_ID].append(e1)
						edgeDict[c_ID].append(e2)

		for key in sorted(edgeDict.iterkeys()):
			for value in edgeDict[key]:
				f.write(value)
		
		f.close()

							#a->twin(b) kemur tvisvar


	#--------------------------------------------------------------------------
	#----------------functions for creating and comparing to G_naive-----------
	#--------------------------------------------------------------------------
	def createNaive(self):
		#create G_naive from the same contigs as G consists of
		if self.printFunctionNames:
			print "createNaive(G)"
		
		#Make sure we have all kmers are in the kmerdict
		self.addKmersFromAllContigs()
		G,cs = dbg.all_contigs(self.kmers,self.k)
		G_naive = Graph(k=self.k,pfn=self.printFunctionNames,ps=self.printStatus,al=self.assertLegal,pil=self.printIsLegal)
		dbg.createGraphObject(G,cs,self.k,G_naive,self.printFunctionNames,ps=False)
		return G_naive

	def equalsNaive(self,G_naive=-1):
		if self.printFunctionNames:
			print "graphEqualsNaive(G,G_naive)"
		if G_naive==-1:
			G_naive = self.createNaive()
		if self.printStatus:
			self.printContigs("G")
			G_naive.printContigs("G_naive")
		return compareGraphs.isSameGraph(self,G_naive,False,self.printFunctionNames,self.printStatus)


	#--------------------------------------------------------------------------
	#-----------------Simple functions for working with graphs-----------------
	#--------------------------------------------------------------------------
	def setIN(self,ID,IN_list):
		#Before:	ID is the ID of a contig c in the graph
		#			IN_list is a list of tuples. E.g. [(0,True)]
		#After:		the IN of c has been set as IN_list
		if self.printFunctionNames:
			print "setIN(ID="+str(ID)+", IN_list="+str(IN_list)+")"
		temp = self.contigs[ID]
		temp[1] = IN_list

	def setOUT(self,ID,OUT_list):
		if self.printFunctionNames:
			print "setOUT(ID="+str(ID)+", OUT_list="+str(OUT_list)+")"
		temp = self.contigs[ID]
		temp[2] = OUT_list

	def addIN(self,ID,IN_tuple):
		#Before:	ID is the ID of a contig c in the graph
		#			IN_tuple is a tuple. E.g. (0,True)
		#After:		IN_tuple has been added to c's IN
		if self.printFunctionNames:
			print "(addIN(ID="+str(ID)+", IN_tuple="+str(IN_tuple)+")"
		temp = self.contigs[ID]
		temp[1].append(IN_tuple)

	def addOUT(self,ID,OUT_tuple):
		if self.printFunctionNames:
			print "(addOUT(ID="+str(ID)+", OUT_tuple="+str(OUT_tuple)+")"
		temp = self.contigs[ID]
		temp[2].append(OUT_tuple)

	def ID_has_IN(self,ID,IN_tuple):
		return IN_tuple in self.contigs[ID][1]

	def ID_has_OUT(self,ID,IN_tuple):
		return IN_tuple in self.contigs[ID][2]

	def getOtherIN(self,c_ID,c_IN):
		#Returns any IN connection except a connection from C->C
		assert(len(c_IN)>0),"c_IN can't be empty"
		for (x_ID,x_B) in c_IN:
			if x_ID!=c_ID:
				return (x_ID,x_B)
		assert(False),"We should never get here"

	def getOtherOUT(self,c_ID,c_OUT):
		assert(len(c_OUT)>0),"c_OUT can't be empty"
		for (x_ID,x_B) in c_OUT:
			if x_ID!=c_ID:
				return (x_ID,x_B)
		assert(False),"We should never get here"

	def addKmersFromAllContigs(self):
		for ID, values in self.contigs.iteritems():
			self.addKmersFromContig(ID,values[0])

	def addKmersFromContig(self,ID,c="-1",c_twin="-1"):
		#Before: 	c is a contig with ID ID
		#After:		all kmers from c and their twins have been added to the kmerDict
		if c=="-1":
			c = self.contigs[ID][0]
		if c_twin=="-1":
			c_twin = dbg.twin(c)
		for i, km in enumerate(dbg.kmers(c,self.k)):
			self.kmers[km] = [ID,i,True]
		for i, km in enumerate(dbg.kmers(c_twin,self.k)):
			self.kmers[km]= [ID,i,False]

	def getID(self):
		#returns the lowest available ID in the graph. Maintains self.ID as
		#the lowest available ID
		ID = self.ID
		self.ID += 1
		return ID

	def setID(self,newID):
		self.ID = newID

	#XX spurning með að fækka þessum eitthvað
	#Þetta er allt útgáfur af __len__
	def isEmpty(self):
		#returns true if the graph is empty. False otherwise
		return len(self.contigs)==0

	def hasNoKmers(self):
		return len(self.kmers)==0

	def numKmerPairs(self):
		return len(self.kmers)/2

	def contigSum(self):
		S = 0
		for c_ID,[c,c_IN,c_OUT,c_COV] in self.contigs.iteritems():
			S += helpers.cLen(c,self.k)
		return S

	def num_bps(self):
		S = 0
		for values in self.contigs.values():
			S += len(values[0])
		return S

	#--------------------------------------------------------------------------
	#-------------------Functions for working with coverage--------------------
	#--------------------------------------------------------------------------
	def getAverageKmerCoverageOfContig(self,cID):
		#Before: cID is the ID of a contig C in the graph
		#After:  Returns the average coverage of the k-mers in C
		[c,IN,OUT,COV] = self.contigs[cID]
		return COV / helpers.cLen(c,self.k)

	def getCoverageOfKmer(self,km):
		#Returns the estimated coverage of km
		#by returning the average coverage of the
		#contig km occurs in
		[cID, i, B] = self.kmers[km]
		return self.getAverageKmerCoverageOfContig(cID)


	#--------------------------------------------------------------------------
	#----------Functions for checking whether the graph is a legal DBG---------
	#--------------------------------------------------------------------------
	def isLegalDBG(self, skipKmers=False):
		#assert(False), "We should never call this"
		#if self.printIsLegal:
		#	print "isLegalDBG(skipKmers="+str(skipKmers)+")"
		#We start by making sure that we have a legal Graph
		#	-Make sure the dict of contigs has a list of 4 values for each contig
		#	and all of them of the correct type
		#	-Make sure the dict of kmers has a list of 3 values for each kmer and all
		#	of them of the correct type
		#   -Make sure all the kmers are in the correct place in the kmerdict and occur exactly once
		#	-Make sure the highest ID in the Graph is lower then self.ID
		#	-Make sure connections go both ways
		#	-Make sure connections from a contig to itself are of the correct form
		if not self.isLegalGraph():
			if self.printIsLegal:
				print "The graph is not legal"
			return False

		#Now we can make sure we have a legal DBG:

		#Make sure that for every connected pairs of contigs that they cannot
		#merge together into one contig. For example:
		#	ID0.OUT = [(ID1,True)]	and  ID1.IN = [(ID0,True)]
		#	can merge into one contig
		if not self.isLegal_merge():
			if self.printStatus:
				print "For every pairs of connected contigs, we're not supposed to be able to merge them"
			return False

		return True

	def isLegalGraph(self, skipKmers=False):
		#assert(False), "We should never call this"
		#if self.printIsLegal:
		#	print "isLegalGraph()"
		for c_ID in self.contigs:
			#Make sure the dict has a list of 4 values for each contig
			if not len(self.contigs[c_ID])==4:
				if self.printIsLegal:
					print "self.contigs[x] must have 4 values for every x"
				return False

			#Make sure each list is of the form:
			#	ID->[c,IN,OUT,COV]
			#	c is a string, IN and OUT are lists and COV is a number
			[c,c_IN,c_OUT,c_COV] = self.contigs[c_ID]
			if not isinstance(c, basestring):
				if self.printIsLegal:
					print "The contig must be a string"
				return False
			if not isinstance(c_IN, list):
				if self.printIsLegal:
					print "c_IN must be a list"
				return False
			if not isinstance(c_OUT, list):
				if self.printIsLegal:
					print "c_OUT must be a list"
				return False
			if not isinstance(c_COV, int):
				if self.printIsLegal:
					print "c_COV must be a number"
				return False

		#Make sure all kmers in the dict of the correct form
		for km in self.kmers:
			#Make sure the dict has a list of 4 values for each contig
			if not len(self.kmers[km])==3:
				if self.printIsLegal:
					print "self.kmers[x] must have 3 values for every x"
				return False

			#Make sure each list is of the form:
			#	#kmer -> [contigID, index, B]
			#	c is a string, IN and OUT are lists and COV is a number
			[contigID, index, B] = self.kmers[km]
			if not isinstance(contigID, int):
				if self.printIsLegal:
					print "The contigID must be an integer"
				return False
			if not isinstance(index, int):
				if self.printIsLegal:
					print "the index must be an integer"
				return False
			if not isinstance(B, bool):
				if self.printIsLegal:
					print "B must be a boolean"
				return False

		#Make sure that the kmerDict stores all the kmers from the contigs
		#and that every km in the kmerDict is in the correct place in the contigs
		if not skipKmers:
			if not self.isLegal_kmerDict():
				if self.printIsLegal:
					print "Illegal kmerDict"
				return False

		#make sure self.ID has a value higher than those in the dict
		if not self.isLegal_ID():
			if self.printIsLegal:
				print "Illegal ID"
			return False

		#Make sure that connections go both ways
		if not self.isLegal_IN_and_OUT():
			if self.printIsLegal:
				"Illegal IN and/or OUT"
			return False

		#For each contig in the graph with a connection to itself
		#make sure it's of the correct form [(0,False),(0,False)] etc
		if not self.isLegal_connectionsToSelf():
			if self.printIsLegal:
				print "Every contig must have legal connections to itself"
			return False

		return True

	def isLegal_kmerDict(self):
		#assert(False), "We should never call this"
		#if self.printIsLegal:
		#	print "isLegal_kmerDict()"
		#make sure all kmers in correct place
		for c_ID in self.contigs:
			[c,c_IN,c_OUT,c_COV] = self.contigs[c_ID]

			for i, km in enumerate(dbg.kmers(c,self.k)):
				try:
					[ID,index,B] = self.kmers[km]
				except ValueError:
					if self.printIsLegal:
						print "kmer table is broken - ValueError"
					return False
				if not ID==c_ID:
					if self.printIsLegal:
						print "kmer table is broken - incorrect ID"
						print i,km,ID,c_ID
						print c_ID, self.contigs[c_ID]
						print ID, self.contigs[ID]
						print c==self.contigs[ID][0]
					return False
				if not index==i:
					if self.printIsLegal:
						print "kmer table is broken - incorrect index"
					return False
				if not B==True:
					if self.printIsLegal:
						print "kmer table is broken - incorrect B"
					return False

			for i, km in enumerate(dbg.kmers(dbg.twin(c),self.k)):
				try:
					[ID,index,B] = self.kmers[km]
				except ValueError:
					if self.printIsLegal:
						print "kmer table is broken - ValueError - twin"
					return False
				if not ID==c_ID:
					if self.printIsLegal:
						print "kmer table is broken - incorrect ID - twin"
					return False
				if not index==i:
					if self.printIsLegal:
						print "kmer table is broken - incorrect index - twin"
					return False
				if not B==False:
					if self.printIsLegal:
						print "kmer table is broken - incorrect B - twin"
					return False

		#make sure all kmers in the kmerdict indeed occur in the place they're supposed to
		for km, [ID,index,B] in self.kmers.iteritems():
			try:
				[c,IN,OUT,COV] = self.contigs[ID]
			except ValueError:
				if self.printIsLegal:
					print "kmer table is broken - ID is not in self.contigs"
				return False
			if B:
				if not c[index:index+self.k]==km:
					if self.printIsLegal:
						print "kmer table is broken - kmer is not where it's supposed to be in the contig"
					return False
			elif not B:
				if not dbg.twin(c)[index:index+self.k]==km:
					if self.printIsLegal:
						print "kmer table is broken - twin(kmer) is not where it's supposed to be in the contig"
					return False

		return True

	def isLegal_ID(self):
		#assert(False), "We should never call this"
		#if self.printIsLegal:
		#	print "isLegal_ID()"
		ID_max = -1
		for ID in self.contigs:
			if ID>ID_max:
				ID_max = ID

		if ID_max >= self.ID:
			if self.printIsLegal:
				assert isinstance(ID_max, int)
				assert isinstance(self.ID, int)
				print "ID_max >= self.ID"
				print ID_max, self.ID
			return False
		else:
			return True

	def isLegal_IN_and_OUT(self):
		#if self.printIsLegal:
		#	print "isLegal_IN_and_OUT()"
		#make sure connections go both ways:
		for c_ID in self.contigs:
			[c,c_IN,c_OUT,c_COV] = self.contigs[c_ID]
			#Make sure both c_IN and c_OUT have max 4 connectioins
			assert(len(c_IN)<=4)
			assert(len(c_OUT)<=4)

			#IN:
			for (x_ID,x_B) in c_IN:
				legal = False
				try:
					[x,x_IN,x_OUT,x_COV] = self.contigs[x_ID]
				except ValueError:
					if self.printIsLegal:
						print "there is no contig in the graph with ID: "+str(x_ID)
					return False

				#print "x_ID:",x_ID, "x_B:",x_B
				if x_B:
					#x->c
					for (y_ID,y_B) in x_OUT:
						#print "y_ID,y_B:",y_ID,y_B
						if c_ID==y_ID and y_B==True:
							legal = True
							break
				elif not x_B:
					#twin(x)->c
					#twin(c)->x
					for (y_ID,y_B) in x_IN:
						if c_ID==y_ID and y_B==False:
							legal = True
							break
				if legal==False:
					if self.printIsLegal:
						print "x connects to y must imply that y has x connecting to it. For example if x.OUT = (y,True), then y.IN = (x,True)"
					return False

			#print "Done checking IN"
			#OUT:
			for (x_ID,x_B) in c_OUT:
				legal = False
				try:
					[x,x_IN,x_OUT,x_COV] = self.contigs[x_ID]
				except ValueError:
					if self.printIsLegal:
						print "there is no contig in the graph with ID: "+str(x_ID)
					return False
				if x_B:
					#c->x
					for (y_ID,y_B) in x_IN:
						if c_ID==y_ID and y_B==True:
							legal = True
							break
				elif not x_B:
					#c->twin(x)
					#x->twin(c)
					for (y_ID,y_B) in x_OUT:
						if c_ID==y_ID and y_B==False:
							legal = True
							break
				if legal==False:
					if self.printIsLegal:
						print "x connects to y must imply that y has x connecting to it. For example if x.OUT = (y,True), then y.IN = (x,True)"
					return False

		return True

	def isLegal_connectionsToSelf(self):
		#assert(False), "We should never call this"
		#Before:	self is a Graph object
		#After:		if some contig contains a connection to itself
		#			return False if it's not of the correct form. [(0,False),(0,False)] instead of
		#			[(0,False)] for example
		#if self.printIsLegal:
		#	print "isLegal_connectionsToSelf()"
		#if self.printStatus:
		#	self.printContigs("Inside isLegal_connectionsToSelf")
		#Iterate through all contigs in the graph
		for c_ID, [c,IN,OUT,COV] in self.contigs.iteritems():
			#make sure all IN connections are legal
			for i, (i_ID,i_B) in enumerate(IN):
				if (i_ID==c_ID) and (not i_B):
					#We have found a tuple representing a connection from twin(c)->c
					#Make sure we have exactly two occurrences of (c_ID,False)
					count = 0
					for j, (j_ID,j_B) in enumerate(IN):
						if (i_ID,False)==(j_ID,j_B):
							count += 1
					if not count==2:
						if self.printIsLegal:
							print "IN connection ("+str(i_ID)+","+str(i_B)+") from contig " + str(c_ID) +" occurs illegal number of times:", count
						return False


				if (i_ID==c_ID) and i_B:
					#We have found a tuple representing a connection from c->c
					#Make sure (c_ID,True) occurs exactly once in OUT
					count = 0
					for j, (j_ID,j_B) in enumerate(OUT):
						if (i_ID,True)==(j_ID,j_B):
							count += 1
					if not count==1:
						if self.printIsLegal:
							print "IN connection ("+str(i_ID)+","+str(i_B)+") from contig " + str(c_ID) +" occurs illegal number of times:", count
						return False

			#make sure all OUT connections are legal
			for i, (i_ID,i_B) in enumerate(OUT):
				if (i_ID==c_ID) and (not i_B):
					#We have found a tuple representing a connection from c->twin(c)
					#Make sure we have exactly two occurrences of (c_ID,False)
					count = 0
					for j, (j_ID,j_B) in enumerate(OUT):
						if (i_ID,False)==(j_ID,j_B):
							count += 1
					if not count==2:
						if self.printIsLegal:
							print "OUT connection ("+str(i_ID)+","+str(i_B)+") from contig " + str(c_ID) +" occurs illegal number of times:", count
						return False

				if (i_ID==c_ID) and i_B:
					#We have found a tuple representing a connection from c->c
					#Make sure (c_ID,True) occurs exactly once in IN
					count = 0
					for j, (j_ID,j_B) in enumerate(IN):
						if (i_ID,True)==(j_ID,j_B):
							count += 1
					if not count==1:
						if self.printIsLegal:
							print "OUT connection ("+str(i_ID)+","+str(i_B)+") from contig " + str(c_ID) +" occurs illegal number of times:", count
						return False
		return True

	def isLegal_merge(self):
		#assert(False), "We should never call this"
		#if self.printIsLegal:
		#	print "isLegal_merge()"
		for ID, [c,IN,OUT,COV] in self.contigs.iteritems():
			canMerge1 = self.canMergeOUT(ID)
			canMerge2 = self.canMergeIN(ID)
			if canMerge1 or canMerge2:
				if self.printIsLegal:
					print "isLegal_merge returning False"
					print "contig ID: "+str(ID)+" can merge"
				return False
		return True


	#--------------------------------------------------------------------------
	#--------------More complex functions for working with graphs--------------
	#--------------------------------------------------------------------------
	def connect_a_to_b(self,aID,bID,A,B):
		#Before:	aID and bID are the IDs of contigs a and b
		#			a and b are contigs already added to the graph
		#			A and B are Booleans
		#Note:		It is possible that aID=bID. I.e. we may be connecting a to itself
		#After:		a and b have been connected together
		#			A and B represent how we're connecting a and b:
		#			| A=T, B=T |   A=T, B=F  |  A=F, B=T  |     A=F, B=F     |
		#			|   a->b   |  a->twin(b) | twin(a)->b | twin(a)->twin(b) |
		if self.printFunctionNames:
			print "connect_a_to_b(aID=" +str(aID)+ ", bID=" + str(bID) + ", A=" + str(A) + ", B=" + str(B)+")"
		if A and B:
			self.contigs[aID][2].append((bID,True))
			self.contigs[bID][1].append((aID,True))
		elif A and not B:
			self.contigs[aID][2].append((bID,False))
			self.contigs[bID][2].append((aID,False))
		elif not A and B:
			self.contigs[aID][1].append((bID,False))
			self.contigs[bID][1].append((aID,False))
		elif not A and not B:
			self.contigs[bID][2].append((aID,True))
			self.contigs[aID][1].append((bID,True))

	def connect_a_to_b_ifNotAlreadyConnected(self,aID,bID,A,B):
		#Function that connects a and b if and only if they are not already connected
		if self.printFunctionNames:
			print "connect_a_to_b_ifNotAlreadyConnected(aID=" +str(aID)+ ", bID=" + str(bID) + ", A=" + str(A) + ", B=" + str(B)+")"
		if A and B:
			if (not self.ID_has_OUT(aID,(bID,True))) and (not self.ID_has_IN(bID,(aID,True))):
				self.connect_a_to_b(aID,bID,A,B)
		elif A and not B:
			if (not self.ID_has_OUT(aID,(bID,False))) and (not self.ID_has_OUT(bID,(aID,False))):
				self.connect_a_to_b(aID,bID,A,B)
		elif not A and B:
			if (not self.ID_has_IN(aID,(bID,False))) and (not self.ID_has_IN(bID,(aID,False))):
				self.connect_a_to_b(aID,bID,A,B)
		elif not A and not B:
			if (not self.ID_has_OUT(bID,(aID,True))) and (not self.ID_has_IN(aID,(bID,True))):
				self.connect_a_to_b(aID,bID,A,B)

	def flipContig(self,c_ID):
		#Before:	The graph is a legal DBG and
		#			self.contigs[ID] = [c,IN,OUT,Cov]
		#After:		The graph is still a legal DBG and
		#			self.contigs[ID] = [twin(c),IN*,OUT*,COV]
		#			The INs/OUTs of c have been changed so that they connect
		#			to c* instead of c where c* is c after flipping
		if self.printFunctionNames:
			print "flipContig(c_ID="+str(c_ID)+")"
		if self.assertLegal:
			assert self.isLegalGraph()

		[c,IN,OUT,COV] = self.contigs[c_ID]
		new_c = dbg.twin(c)
		new_IN = []		#The IN of c after we have flipped it
		new_OUT = []	#The OUT of c after we have flipped it

		#Booleans used to make sure we only flip each connection from c to c once
		alreadyFlipped_1 = False	#c->c
		alreadyFlipped_2 = False	#c->twin(c)
		alreadyFlipped_3 = False	#twin(c)->c

		#In the for loop below we:
		#	change all INs of c so they instead connect to c*
		#	add tuples to new_IN and new_OUT
		for (i_ID,i_B) in IN:
			[i,i_IN,i_OUT,i_COV] = self.contigs[i_ID]
			if i_B:
				#case 1
				if i_ID==c_ID:	#c->c
					if not alreadyFlipped_1:
						new_IN.append((c_ID,True))
						new_OUT.append((c_ID,True))
						alreadyFlipped_1 = True
					continue
				for index, (j_ID,j_B) in enumerate(i_OUT):
					if (j_ID,j_B)==(c_ID,True):
						new_OUT.append((i_ID,False))
						i_OUT[index] = (c_ID,False)
						self.setOUT(i_ID,i_OUT)
						break
			elif not i_B:
				#case 2
				if i_ID==c_ID:	#twin(c)->c
					if not alreadyFlipped_3:
						new_OUT.append((c_ID,False))
						new_OUT.append((c_ID,False))
						alreadyFlipped_3 = True
					continue
				for index, (j_ID,j_B) in enumerate(i_IN):
					if (j_ID,j_B)==(c_ID,False):
						new_OUT.append((i_ID,True))
						i_IN[index] = (c_ID,True)
						self.setIN(i_ID,i_IN)
						break

		#Now we do the same for OUT:
		for (i_ID,i_B) in OUT:
			[i,i_IN,i_OUT,i_COV] = self.contigs[i_ID]
			if i_B:
				#case 3
				if i_ID==c_ID:	#c->c
					if not alreadyFlipped_1:
						new_IN.append((c_ID,True))
						new_OUT.append((c_ID,True))
						alreadyFlipped_1 = True
					continue
				for index, (j_ID,j_B) in enumerate(i_IN):
					if (j_ID,j_B)==(c_ID,True):
						new_IN.append((i_ID,False))
						i_IN[index] = (c_ID,False)
						self.setIN(i_ID,i_IN)
						break
			elif not i_B:
				#case 4
				if i_ID==c_ID:	#c->twin(c)
					if not alreadyFlipped_2:
						new_IN.append((c_ID,False))
						new_IN.append((c_ID,False))
						alreadyFlipped_2 = True
					continue
				for index, (j_ID,j_B) in enumerate(i_OUT):
					if (j_ID,j_B)==(c_ID,False):
						new_IN.append((i_ID,True))
						i_OUT[index] = (c_ID,True)
						self.setOUT(i_ID,i_OUT)
						break

		#Finally we update c
		self.contigs[c_ID] = [new_c,new_IN,new_OUT,COV]
		self.addKmersFromContig(c_ID,new_c)

		if self.assertLegal:
			assert self.isLegalGraph()

	#--------------------------------------------------------------------------
	#-------------------changeID_FromTo and helper functions-------------------
	#--------------------------------------------------------------------------
	def changeID_FromTo(self,ID,ID_new):
		#Before:	ID is the id of contig c in the legal DBG
		#After:		contig c now has the id ID_new and all appropriate changes have been made
		#				the INs/OUTs of c now connect to ID_new instead of ID
		#				the graph is a legal DBG
		if self.printFunctionNames:
			print "changeID_FromTo(ID="+str(ID)+", ID_new="+str(ID_new)+")"
		if ID==ID_new:
			return

		assert ID_new >= self.ID, "ID_new must be a free ID. ID="+str(self.ID)+", ID_new="+str(ID_new)

		#Find all occurrences of ID in the IN and OUT connections of c and change them into ID_new
		self.changeID_ofConnections(ID,ID_new)

		#initialize self.contigs[ID_new] with the same values as self.contigs[ID]:
		#[c,IN,OUT,COV] = self.contigs[ID]
		#self.contigs[ID_new] = [c,IN,OUT,COV]
		temp = self.contigs[ID]
		self.contigs[ID_new] = temp

		#delete self.contigs[ID] from the graph
		del self.contigs[ID]

		#add the kmers from c again (now with ID_new)
		self.addKmersFromContig(ID_new,temp[0])

		self.setID(max(ID_new+1,self.ID))

		#now self.contigs[ID_new] has all the same attributes as self.contigs[ID] had before we called the function and the graph is a legal DBG

	def changeID_ofConnections(self,ID,ID_new):
		if self.printFunctionNames:
			print "changeID_ofConnections(ID="+str(ID)+", ID_new="+str(ID_new)+")"
		#Create a list of all occurences of ID in the Graph (plus self.contigs[ID] which is not in this list by default)
		conns = self.changeID_findConnections(ID)

		#change all occurences of ID into ID_new:
		self.changeID_changeConnections(ID,ID_new,conns,ID_old=-1)

	def changeID_findConnections(self,ID):
		if self.printFunctionNames:
			print "changeID_findConnections(ID="+str(ID)+")"
		temp = self.contigs[ID]
		conns = []	#a list of the IDs of all the contigs that c connects to
		for (i_ID,i_B) in temp[1]:
			conns.append(i_ID)
		for (i_ID,i_B) in temp[2]:
			conns.append(i_ID)
		return conns

	def changeID_changeConnections(self,ID,ID_new,conns,ID_old=-1):
		if self.printFunctionNames:
			print "changeID_changeConnections(ID="+str(ID)+", ID_new="+str(ID_new)+", conns="+str(conns)+", ID_old="+str(ID_old)+")"
		if ID_old==-1:
			ID_old = ID
		for i_ID in conns:
			[i,i_IN,i_OUT,i_COV] = self.contigs[i_ID]
			i_IN_new = list(i_IN)
			for index, (x_ID,x_B) in enumerate(i_IN):
				if x_ID==ID_old:
					i_IN_new[index] = (ID_new,x_B)
			self.setIN(i_ID,i_IN_new)

			i_OUT_new = list(i_OUT)
			for index, (x_ID,x_B) in enumerate(i_OUT):
				if x_ID==ID_old:
					i_OUT_new[index] = (ID_new,x_B)
			self.setOUT(i_ID,i_OUT_new)



	#--------------------------------------------------------------------------
	#--------------------addSegmentToGraph and subfunctions--------------------
	#--------------------------------------------------------------------------
	def addSegmentToGraph(self,s,CS=False):
		#print "addSegmentToGraph(s="+str(s)+")"
		#Before:	segment is a DNA string
		#			this Graph is a legal DBG graph
		#After:		segment has been added to the graph
		#			the kmers from segment have been added to self.kmers
		#			the graph is a legal DBG graph
		if self.printFunctionNames:
			print "addSegmentToGraph(s="+str(s)+")"
		if self.printStatus:
			self.printContigs("inside addSegmentToGraph")
		if self.assertLegal:
			assert self.isLegalDBG()
		L = len(s)
		if L < self.k:
			return
		#if CS==True:
		#	self.carefulSplit = True
		#	print "We're setting carefulSplit=True"
		#	print "s:", s
		#assert(self.carefulSplit==False),"Haven't implemented maxSplitCov yet"

		#Fyrir hvert s í segment sem við höfum aldrei séð áður:
		#	1. Splittum s á tengingum við sjálft sig og skilum einu s2 í einu
		#	2. Splittum s2 á tengingum við þá sem í grafinu og skilum einu s3 í einu
		#	3. Bætum s3 við grafið
		#for s,s_COV in SL.iteritems():
		s_twin = dbg.twin(s)
		for i,j in self.splitSegmentOnSeenKmers(s,s_twin,L):
			assert(i>=0)
			assert(j<=L)
			assert(j-i>=self.k)
			if self.assertLegal:
				s_test = s[i:j]
				assert(len(s_test)>=self.k)
				for temp in dbg.kmers(s_test,self.k):
					assert(temp not in self)
					assert(dbg.twin(temp) not in self)
			for i2,j2 in helpers.splitOnConnToSelf(s,s_twin,L,self.k,i,j,i):
				assert(i2>=i)
				assert(j2<=j)
				assert(j2-i2>=self.k)
				if self.assertLegal:
					s2 = s[i2:j2]
					assert(len(s2)>=self.k)
					for temp in dbg.kmers(s2,self.k):
						assert(temp in s_test), "\n"+str(s)+"\n"+str(s_test)+"\n"+str(s2)
						assert(temp not in self)
						assert(dbg.twin(temp) not in self)
				for i3,j3 in self.splitSegmentOnConnToOthers(s,L,i2,j2,i2):
					#s3 does not need to be split any further
					assert(i3>=i2)
					assert(j3<=j2)
					assert(j3-i3>=self.k)
					s3 = s[i3:j3]
					if self.assertLegal:
						assert(len(s3)>=self.k)
						for temp in dbg.kmers(s3,self.k):
							assert(temp not in self)
							assert(dbg.twin(temp) not in self)
					if CS:
						#Either add s3 as normal or mark that we have been told to split
						skipDueTo_CS = self.splitOthers(s3,CS=True)
						if skipDueTo_CS:
							continue
					else:
						self.splitOthers(s3,CS=False)		#<-- Eina fallið sem breytist þegar nota maxSplitCov
					s3_ID = self.getID()
					self.contigs[s3_ID] = [s3,[],[],helpers.cLen(s3,self.k)*2]
					self.addKmersFromContig(s3_ID,s3)
					self.connectSegment(s3_ID,s3)
					s3_ID = self.mergeSegment(s3_ID)				#kalla tvisvar á merge því við getum lent í að þurfa tvö merge þegar bætum við nýrri nóðu í grafið
					s3_ID = self.mergeSegment(s3_ID)
					if self.assertLegal:
						assert self.isLegalDBG()

	def splitSegmentOnSeenKmers(self,s,s_twin,L):
		#print "splitSegmentOnSeenKmers(s="+str(s)+")"
		#Pre:	s is a segment
		#	yields prefixes of s which don't contain any seen k-mers.
		#	Increments the cov of seen k-mers by 1
		if L<self.k:
			return
		elif L==self.k:
			if s in self:
				contigID = self.kmers[s][0]
				self.contigs[contigID][3] += 1
				return
			else:
				#print "yielding 2, s:           ", s
				yield 0, L
				return
		kmersIn_s = set()
		start = 0
		for i, km in enumerate(dbg.kmers(s,self.k)):
			#print i,km,start
			km_twin = s_twin[L-self.k-i:L-i]
			km = min(km,km_twin)
			B1 = km in self
			B2 = km in kmersIn_s
			if not (B1 or B2):		#We have not seen km before
				kmersIn_s.add(km)
			else:					#We have seen km before
				if B1:				#increment cov by one
					contigID = self.kmers[km][0]
					self.contigs[contigID][3] += 1
				if i-start!=0:
					#				         start
					#split s into:   | seen |  s0  |          | where
					#                |  	   | km |         |
					#                           i    i+k
					if start<=L-self.k:
						#yield s[start:i-1+self.k]
						#print start, i-1+self.k
						yield start, i-1+self.k
				start = i+1
		if start<=L-self.k:
			#print "yielding 4, s[start:]:   ", s[start:]
			#yield s[start:]
			#print start, L
			yield start, L
			if self.assertLegal:
				assert(len(s[start:])>=self.k)

	def splitSegmentOnConnToOthers(self,s,L,i,j,start):
		if j-start < self.k:
			return
		elif j-start == self.k:
			yield start, j
			return
		#for i, km in enumerate(dbg.kmers(s,self.k)):
		km = s[i:i+self.k]
		if i>start:
			for y in dbg.bw(km):
				if y in self: #a non first kmer in s has a connection from a node in the graph
					#s0, s = helpers.splitString(s,i+self.k-1,self.k)	#split in front of km
					#yield s0
					yield start, i+self.k-1
					for temp1,temp2 in self.splitSegmentOnConnToOthers(s,L,i,j,i):
						yield temp1, temp2	#hugsa að þetta gæti orðið óendanleg lúppa
					return
		if i!=j-self.k:
			for y in dbg.fw(km):
				if y in self:	#a non last kmer in s connects to a node in the graph
					#s0, s = helpers.splitString(s,i+self.k,self.k)	#split behind km
					#yield s0
					yield start, i+self.k
					for temp1,temp2 in self.splitSegmentOnConnToOthers(s,L,i+1,j,i+1):
						yield temp1, temp2
					return
		#If a is the last kmer we return s[start:]
		if i+self.k==j:
			yield start, j
			return
		#Recursively call to get the next a
		#print "case 9"
		for temp1,temp2 in self.splitSegmentOnConnToOthers(s,L,i+1,j,start):
			yield temp1,temp2

	def splitOthers(self,s,CS=False):
		#Helper function for addSegmentAlreadySplit
		#Pre:		The graph is a legal DBG
		# 			s does not need to be split further before being added to the graph
		#			We may need to split some of the contigs in the graph before adding s
		#Post:		s is unchanged (no recursive calls)
		#			No more contigs need to be split in order to add s to the graph
		#			The graph is not necessarily a legal DBG
		if self.printFunctionNames:
			print "splitOthers(s="+str(s)+")"
		if self.printStatus:
			self.printContigs("inside splitOthers")
		first = s[0:self.k]			#first is the first km in s
		first_twin = dbg.twin(first)
		last = s[len(s)-self.k:]	#last is the last km in s
		last_twin = dbg.twin(last)
		splitSomething = True
		while splitSomething:
			splitSomething = False
			for km in dbg.fw(last):
				if km in self.kmers:		#s->km
					[c_ID, c_i, c_B] = self.kmers[km]
					if c_B and (c_i!=0):
						if CS:
							count = self.splitDict[c_ID,c_i+self.k-1]
							if count==0:
								self.splitDict[c_ID,c_i+self.k-1] += 1
								return True
						self.splitContig(c_ID,c_i+self.k-1)		#split before km
						splitSomething = True
			for km in dbg.fw(first_twin):
				if km in self.kmers:		#twin(s)->km
					[c_ID, c_i, c_B] = self.kmers[km]
					if c_B and (c_i!=0):
						if CS:
							count = self.splitDict[c_ID,c_i+self.k-1]
							if count==0:
								self.splitDict[c_ID,c_i+self.k-1] += 1
								return True
						self.splitContig(c_ID,c_i+self.k-1)		#split before km
						splitSomething = True
			for km in dbg.bw(first):
				if km in self.kmers:		#km->s
					[c_ID, c_i, c_B] = self.kmers[km]
					c_L = len(self.contigs[c_ID][0])
					if c_B and (c_i!=(c_L-self.k)):
						if CS:
							count = self.splitDict[c_ID,c_i+self.k]
							if count==0:
								self.splitDict[c_ID,c_i+self.k] += 1
								return True
						self.splitContig(c_ID,c_i+self.k)		#split after km
						splitSomething = True
			for km in dbg.bw(last_twin):	#km->twin(s)
				if km in self.kmers:
					[c_ID, c_i, c_B] = self.kmers[km]
					c_L = len(self.contigs[c_ID][0])
					if c_B and (c_i!=(c_L-self.k)):
						if CS:
							count = self.splitDict[c_ID,c_i+self.k]
							if count==0:
								self.splitDict[c_ID,c_i+self.k] += 1
								return True
						self.splitContig(c_ID,c_i+self.k)		#split after km
						splitSomething = True
		assert(splitSomething==False)
		return False

	def splitContig(self,ID,i):
		#helper function for splitOthers
		#Before: 	ID is the ID of contig c
		#			i is an integer
		#			IF c has a connection to itself it is from end to end
		#			(because we started with a legal DBG and c was part of it. Otherwise c would have been split up earlier)
		#			SHOULD I MAKE IT NECESSARY THAT self.k <= i < L (so c0 and c1 both have at least 1 kmer)?!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		#After:
		#	c has been deleted from the graph
		#		c has been deleted from self.contigs
		#		The kmers of c have been deleted from self.kmers
		#	c0 and c1 have been added to the graph
		#		and their kmers added to self.kmers
		#	c0 has been split up into
		#	   0	 i               L
		#	c0[0:i-1]  ->  s[i-k:L-1]
		#	  c0           c1
		#	where
		#		c0.IN = c.IN (and the contigs in c.IN get c0 as OUT instead of c)
		#		c0.OUT = c1
		#		c1.IN = c0
		#		c1.OUT = c.OUT (and the contigs in c.OUT get c1 as IN instead of c)
		#Need to test this for special cases. Might need to handle connections from c->c in
		#a different way
		if self.printFunctionNames:
			print "splitContig(ID="+str(ID)+", i="+str(i)+")"
		if self.printStatus:
			print self
		if self.assertLegal:
			assert self.isLegalGraph()
		[c,IN,OUT,COV] = self.contigs[ID]

		c0, c1 = helpers.splitString(c,i,self.k)
		if self.printStatus:
			print c0,c1
			print dbg.twin(c0),dbg.twin(c1)

		#Delete c
		del self.contigs[ID]

		if self.printStatus:
			self.printContigs("After deleting ID")

		#Create c0 and c1 and add correct connections
		ID0 = self.getID()
		ID1 = self.getID()
		l0 = helpers.cLen(c0,self.k)
		l1 = helpers.cLen(c1,self.k)
		if COV<=0:
			raise Exception("COV must be greater than zero")
		#if COV==1 or COV==2:
		#	COV0 = 1
		#	COV1 = 1
		else:
			COV0 = int(ceil(COV * ( float(l0) / (l0+l1) )))
			COV1 = COV - COV0
		if self.assertLegal:
			assert(COV0+COV1==COV), "The total COV must be the same before and after the split"
			assert(COV0>0 and COV1>0), "COV must be greater than zero"
		self.contigs[ID0] = [c0, IN, [(ID1,True)], COV0]
		self.contigs[ID1] = [c1, [(ID0,True)], OUT, COV1]
		self.addKmersFromContig(ID0,c0)
		self.addKmersFromContig(ID1,c1)

		#If c had a connection from c->twin(c) i.e. c.OUT=[...(ID,False),(ID,False)...]
		#we change:
		#	(ID,False),(ID,False) in ID1.OUT into (ID1,False),(ID1,False)
		#Faster code:
		for i,y in enumerate(OUT):
			if y==(ID,False):
				OUT[i] = (ID1,False)
			elif y==(ID,True):
				OUT[i] = (ID0,True)
		for i,y in enumerate(IN):
			if y==(ID,False):
				IN[i] = (ID0,False)
			elif y==(ID,True):
				IN[i] = (ID1,True)

		if self.printStatus:
			print ID,ID0,ID1
			self.printContigs("After creating c0 and c1. Before changing neighbour connections")

		#Now everything that's left is to update the connections to c0 and
		#from c1 so that they indeed connect to c0/c1 instead of c

		for (x_ID,x_B) in IN:
			[x,x_IN,x_OUT,x_COV] = self.contigs[x_ID]
			if x_B:
				#x_OUT[:] = [y if y != (ID,True) else (ID0,True) for y in x_OUT]
				for i,y in enumerate(x_OUT):
					if y==(ID,True):
						x_OUT[i] = (ID0,True)
			else:
				#x_IN[:] = [y if y != (ID,False) else (ID0,False) for y in x_IN]
				for i,y in enumerate(x_IN):
					if y==(ID,False):
						x_IN[i] = (ID0,False)

		for (x_ID,x_B) in OUT:
			[x,x_IN,x_OUT,x_COV] = self.contigs[x_ID]
			if x_B:
				#x_IN[:] = [y if y != (ID,True) else (ID1,True) for y in x_IN]
				for i,y in enumerate(x_IN):
					if y==(ID,True):
						x_IN[i] = (ID1,True)
			else:
				#x_OUT[:] = [y if y != (ID,False) else (ID1,False) for y in x_OUT]
				for i,y in enumerate(x_OUT):
					if y==(ID,False):
						x_OUT[i] = (ID1,False)

		if self.printStatus:
			print ID,ID0,ID1
			self.printContigs("After changing neighbour connections")

		if self.assertLegal:
			assert self.isLegalGraph()

		return ID0, ID1, c0, c1

	def connectSegment(self,ID,segment=""):
		#Helper function for addSegmentAldreadySplit
		#Before:	ID is the ID of segment
		#			self.contigs[ID] = [segment,[],[],0]
		#			The kmers from segment have been added to self.kmers
		#			The graph is a legal DBG except for one thing:
		#				segment has not been connected to the other contigs
		#After:		segment has been connected to the other contigs in the graph
		if segment=="":
			[segment,IN,OUT,COV] = self.contigs[ID]
		if self.printFunctionNames:
			print "connectSegment(ID="+str(ID)+", segment="+str(segment)+")"
		if self.printStatus:
			self.printContigs("inside connectSegment")
		if self.assertLegal:
			assert self.isLegalGraph()

		#Add connections of the form:
		#	segment->N
		#	segment->twin(N)
		L = len(segment)
		first = segment[0:self.k]
		last = segment[L-self.k:]
		for y in dbg.fw(last):
			if y in self.kmers:
				[N_ID, N_i, N_B] = self.kmers[y]
				if N_i==0:
					self.connect_a_to_b_ifNotAlreadyConnected(ID,N_ID,True,N_B)

		#Add connections of the form:
		#	N->segment
		#	twin(N)->segment
		for y in dbg.bw(first):
			if y in self.kmers:
				[N_ID, N_i, N_B] = self.kmers[y]
				N = self.contigs[N_ID][0]
				if N_i==len(N)-self.k:
					self.connect_a_to_b_ifNotAlreadyConnected(N_ID,ID,N_B,True)

		if self.assertLegal:
			assert self.isLegalGraph()

	def mergeSegment(self, c_ID):
		#Helper function for addSegmentAldreadySplit
		#Before:	c_ID is the c_ID of C
		#			The graph is a legal DBG except we might need to merge C
		#			with it's connections
		#After:		C has been merged with it's neighbours wherever possible.
		#			The graph is a legal DBG
		#			Note: We need to make sure that the INs and OUTs of new_c have
		#			new_ID instead of ID/N_ID in their INs / OUTs
		#Note:
		#	It is possible we can merge with an IN connection as well
		#	if we flipped the contig when merge-ing it will become an
		#	OUT connection. Therefore we don't want to throw an exception
		#	in this case (the error discovered 22.03.2017)
		assert c_ID in self.contigs, "c_ID must be the ID of a contig in the graph"
		if self.printFunctionNames:
			print "mergeSegment(c_ID="+str(c_ID)+")"
		if self.printStatus:
			self.printContigs("inside mergeSegment")
		if self.assertLegal:
			assert self.isLegalGraph()

		CMO = self.canMergeOUT(c_ID)
		if CMO:
			[n_ID,n_B] = CMO
			if n_B:	#C->N
				c_ID = self.merge(c_ID,n_ID,True,True,"")
			else:	#C->twin(N)
				c_ID = self.merge(c_ID,n_ID,True,False,"A")
			return c_ID
		CMI = self.canMergeIN(c_ID)
		if CMI:
			[n_ID,n_B] = CMI
			if n_B:	#N->C
				c_ID = self.merge(n_ID,c_ID,True,True,"")
			else:	#twin(N)->C
				c_ID = self.merge(n_ID,c_ID,False,True,"B")
			return c_ID

		if self.assertLegal:
			if not self.isLegalDBG():
				self.printContigs("Not a legal DBG")
				print "The ID of the contig failing to merge:",c_ID
				raise Exception("The graph is supposed to be a legal DBG")

		return c_ID		#We also return c_ID if we were unable to merge

	def canMergeOUT(self,c_ID):
		assert c_ID in self.contigs, "ID must be the ID of a contig in the graph"
		c_OUT = self.contigs[c_ID][2]
		if len(c_OUT)==1:
			#C has exactly one OUT: N. Find all the info we have about N:
			[n_ID,n_B] = c_OUT[0]
			[n, n_IN, n_OUT, n_COV] = self.contigs[n_ID]
			if c_ID==n_ID:
				return []

			#Check whether we can merge C->N:
			if n_B and len(n_IN)==1:
				if self.printStatus:
					print "C->N. C has exactly one OUT and N has exactly one IN. Can merge"
				return [n_ID,True]

			#Check whether we can merge C->twin(N):
			elif not n_B and len(n_OUT)==1:
				if self.printStatus:
					print "C->twin(N) and N->twin(C). both have exactly one OUT. Can merge"
				return [n_ID,False]
		return []

	def canMergeIN(self,c_ID):
		assert c_ID in self.contigs, "ID must be the ID of a contig in the graph"
		c_IN = self.contigs[c_ID][1]
		if len(c_IN)==1:
			#C has exactly one IN: N. Find all the info we have about N:
			[n_ID,n_B] = c_IN[0]
			[n, n_IN, n_OUT, n_COV] = self.contigs[n_ID]
			if c_ID==n_ID:
				return []

			#Check whether we can merge N->C:
			if n_B and len(n_OUT)==1:
				if self.printStatus:
					print "N->C. C has exactly one IN and N has exactly one OUT. Can merge"
				return [n_ID,True]

			#Check whether we can merge twin(N)->C:
			elif not n_B and len(n_IN)==1:
				if self.printStatus:
					print "twin(N)->C and twin(C)->N. both have exactly one IN. Can merge"
				return [n_ID,False]
		return []

	#flipSelect:
	#	"" flip either one
	#	"A" flip A
	#	"B" flip B
	def merge(self,a_ID,b_ID,A,B,flipSelect=""):
		#Helper function for mergeSegment
		#A=True, B=True:   a->b
		#A=True, B=False:  a->twin(b)
		#A=False, B=True:  twin(a)->b
		#A=False, B=False: twin(a)->twin(b) = b->a
		if self.printFunctionNames:
			print "merge(a_ID="+str(a_ID)+", b_ID="+str(b_ID)+", A="+str(A)+", B="+str(B)+")"
		assert a_ID in self.contigs, "a_ID must be the ID of a contig in the graph"
		assert b_ID in self.contigs, "b_ID must be the ID of a contig in the graph"
		assert isinstance(A,bool)
		assert isinstance(B,bool)
		if A and B:
			assert(flipSelect=="")
			#merging a->b into new_c (by deleting b, adding it to a and updating a and it's connections)
			#self.changeID_OUTs(b_ID,a_ID,conns)		#make the OUTs of b instead connect to a
			self.changeID_ofConnections(b_ID,a_ID)
			[a, a_IN, a_OUT, a_COV] = self.contigs[a_ID]
			[b, b_IN, b_OUT, b_COV] = self.contigs[b_ID]
			new_c = a + b[self.k-1:]
			self.contigs[a_ID] = [new_c, a_IN, b_OUT, a_COV+b_COV]
			del self.contigs[b_ID]
			#add kmers from new_c to kmerDict
			self.addKmersFromContig(a_ID, new_c)
			return a_ID
		elif A and (not B):
			#merging a->twin(b) into new_c
			if flipSelect in ["","A"]:
				self.flipContig(a_ID)
				return self.merge(b_ID,a_ID,True,True,"")
			elif flipSelect=="B":
				self.flipContig(b_ID)
				return self.merge(a_ID,b_ID,True,True,"")
		elif (not A) and B:
			#merging twin(a)->b into new_c (by deleting a)
			if flipSelect in ["","A"]:
				self.flipContig(a_ID)
				return self.merge(a_ID,b_ID,True,True,"")
			elif flipSelect=="B":
				self.flipContig(b_ID)	#twin(a)->twin(b)==b->a
				return self.merge(b_ID,a_ID,True,True,"")
		elif (not A) and (not B):
			return self.merge(b_ID,a_ID,True,True,"")

	
	#--------------------------------------------------------------------------
	#----------functions for checking for isolated, tips and bubbles-----------
	#--------------------------------------------------------------------------
	#---------------------------------
	#markAllIsolated and it's helpers:
	#---------------------------------
	def createAllPathsInContigsCollection(self, MPL, c_ID, skipNodes=set()):
		assert(isinstance(skipNodes,set)), "skipNodes must be a set"
		#Create all paths with length < MPL which include C
		allPathsIncludingC = self.createAllPathsIncludingContig(MPL, c_ID, skipNodes)
		if allPathsIncludingC.isEmpty():	
			return Path.PathList()

		#Take a copy of allPathsIncludingC
		allPaths = Path.PathList()		#This will store every path in C's collection of contigs

		#Now we know that all the paths C is a part of (stored in allPathsIncludingC) each have a length shorter than MPL.
		#Below we make sure that none of the contigs in any of the paths is part of another path which has length >= MPL
		#We add all the paths to allPaths
		for p in allPathsIncludingC.getPaths():
			for (x_ID,x_B) in p.contigs:
				xPaths = self.createAllPathsIncludingContig(MPL, x_ID, skipNodes)
				if xPaths.isEmpty():
					#X is part of a path with length >= MPL
					return Path.PathList()
				else:
					for p in xPaths.getPaths():
						if not p in allPaths:
							allPaths.append(p)
		return allPaths

	def createAllPathsIncludingContig(self, MPL, c_ID, skipNodes=set(),ignoringLengthOfLast=False):
		#Returns an empty PathList if any of the paths has a total lengh > MPL
		#print "createAllPathsIncludingContig(MPL="+str(MPL)+", c_ID="+str(c_ID)+", skipNodes="+str(skipNodes)+", ignoringLengthOfLast="+str(ignoringLengthOfLast)+")"
		assert(isinstance(skipNodes,set)), "skipNodes must be a set"
		[c, c_IN, c_OUT, c_COV] = self.contigs[c_ID]
		withinMaxLength = helpers.cLen(c,self.k) < MPL
		if not withinMaxLength:
			return Path.PathList()
		visited = {c_ID}		#a set storing all contigs we have visited
		temp = Path.Path()
		temp.append(c_ID,True,c,self.k)
		pathsInProgress = Path.PathList([temp])
		allPaths = Path.PathList()

		#Create all possible paths which include C
		while not pathsInProgress.isEmpty() and withinMaxLength:
			nodeHasChildAncestorNeither = 2		#0 means child, 1 ancestor, 2 neither
			p = pathsInProgress.pop()
			
			#Get all child nodes of the last contig in the path
			(e_ID,e_B) = p.getLastContigTuple()
			DCorA = self.getDirectChildren(e_ID,e_B)
			DCorA = helpers.DCorA_toAdd(DCorA,visited,skipNodes)
			if DCorA:
				nodeHasChildAncestorNeither = 0
			else:
				#If there were no child nodes we instead get all ancestor nodes of the first contig in the path
				(e_ID,e_B) = p.getFirstContigTuple()
				DCorA = self.getDirectAncestors(e_ID,e_B)
				DCorA = helpers.DCorA_toAdd(DCorA,visited,skipNodes)
				if DCorA:
					nodeHasChildAncestorNeither = 1
			#print DCorA

			#Remove nodes which had already been marked
			newDCorA = []
			for (x_ID,x_B) in DCorA:
				if (self.degrees[x_ID][2] not in [1,2,-3]):
					newDCorA.append((x_ID,x_B))
			DCorA = newDCorA

			#If there were no ancestor nodes then the path we just popped from pathsInProgress neither had an ancestor nor a child.
			#We therefore add it to allPaths
			if not DCorA:
				allPaths.append(p)
			else:
				for (x_ID,x_B) in DCorA:
					assert(not x_ID in skipNodes), "made sure of this above"
					assert(not x_ID in visited), "made sure of this above"
					visited.add(x_ID)
					assert(self.degrees[x_ID][2] not in [1,2,-3]), "We had already removed nodes which had already been marked"
					[x, x_IN, x_OUT, x_COV] = self.contigs[x_ID]
					if nodeHasChildAncestorNeither==0:	#child
						p.append(x_ID,x_B,x,self.k)
					elif nodeHasChildAncestorNeither==1:	#ancestor
						p.prepend(x_ID,x_B,x,self.k)
					if ignoringLengthOfLast:
						lenP = p.lengthOfAllButLast()
					else:
						lenP = len(p)
					if lenP<MPL:
						pathsInProgress.append(p)
					else:
						withinMaxLength = False
						allPaths.append(p)
						return Path.PathList()
		return allPaths

	def isPartOfIsolated(self, MPL, c_ID, skipNodes=set()):
		assert(isinstance(skipNodes,set)), "skipNodes must be a set"
		#Create all paths with length < MPL in C's collection of contigs
		allPaths = self.createAllPathsInContigsCollection(MPL, c_ID, skipNodes)
		if allPaths.isEmpty():	
			return False
		else:
			return True

	def markAllIsolated(self, MPL):
		for c_ID in self.contigs:
			if not (self.degrees[c_ID][2] in [-3,1]):	#no need to mark the same collection twice
				#Create all paths with length < MPL in C's collection of contigs
				allPaths = self.createAllPathsInContigsCollection(MPL, c_ID, skipNodes=set())
				if allPaths.isEmpty():	
					isIso = False
				else:
					isIso = True
				if isIso:
					assert(isinstance(allPaths,Path.PathList)), "allPaths must be a PathList object"
					theContigs = allPaths.getContigIDs()
					assert(len(theContigs)>0),"theContigs can't be empty"
					self.markAllContigsInSetWithNew_R(theContigs,new_R=1)
				else:
					theContigs2 = self.findAllConnectedContigs(c_ID)
					assert(len(theContigs2)>0),"theContigs2 can't be empty"
					self.markAllContigsInSetWithNew_R(theContigs2,new_R=-3)

		#self.printContigsWithRatings()
		#Now we have marked all contigs which are not part of an isolated collection of contigs with a -3
		#The rest is marked with a -1. Those are all isolated
		#Therefore we change -1 into 1
		#         and change -3 into -1
		for c_ID in self.contigs:
			R = self.degrees[c_ID][2]
			assert(R in [1,-3]), "We have marked all non isolated with -3 and the rest should be marked with a -1. R="+str(R)
		for c_ID in self.contigs:
			R = self.degrees[c_ID][2]
			if R==-3:
				self.degrees[c_ID][2] = -1
		for c_ID in self.contigs:
			R = self.degrees[c_ID][2]
			assert(R in [-1,1]), "We have marked all isolated with a 1 and the rest with a -1. R="+str(R)

	#---------------------------------
	#--isPartOfTip and it's helpers:--
	#---------------------------------
	def existsAltPath_fw(self,MPL,a_ID,b_ID):
		#Returns True if there exists a path fw from a_ID which doesn't go through b_ID
		#and contains a path which breaks the MPL limit
		#print "existsAltPath_fw(MPL="+str(MPL)+", a_ID="+str(a_ID)+", b_ID="+str(b_ID)+")"
		#Create all paths which we can reach from a_ID, skipping B and all incoming connections of A
		[a,a_IN,a_OUT,a_COV] = self.contigs[a_ID]
		skipNodes = {a_ID,b_ID}
		for (x_ID,x_B) in a_IN:
			skipNodes.add(x_ID)
		DC = self.getDirectChildren(a_ID,True,skipNodes)
		for (x_ID,x_B) in DC:
			allPaths = self.createAllPathsInContigsCollection(MPL, x_ID, skipNodes)
			if allPaths.isEmpty():
				return True
		return False

	def existsAltPath_bw(self,MPL,a_ID,b_ID):
		[a,a_IN,a_OUT,a_COV] = self.contigs[a_ID]
		skipNodes = {a_ID,b_ID}
		for (x_ID,x_B) in a_OUT:
			skipNodes.add(x_ID)
		DC = self.getDirectAncestors(a_ID,True,skipNodes)
		for (x_ID,x_B) in DC:
			allPaths = self.createAllPathsInContigsCollection(MPL, x_ID, skipNodes)
			if allPaths.isEmpty():
				return True
		return False

	def isPartOfTip(self,c_ID,MPL):
		#Returns True if C is part of a tip. Also marks all contigs in the tip with R=2
		#False otherwise (possibility of a false negative but not false positive)
		#Don't mark anything when returning False
		#print "isPartOfTip(c_ID="+str(c_ID)+", MPL="+str(MPL)+")"
		assert(c_ID in self.contigs), "c_ID must be in the Graph"
		[c, c_IN, c_OUT, c_COV] = self.contigs[c_ID]
		c_inD, c_outD, R = self.degrees[c_ID]
		assert(R in [-1,2]), "We have marked all isolated with a 1, tips with a 2 and the rest with a -1. We never call this function if C had been marked as isolated"
		if R==2:	return True
		if helpers.cLen(c,self.k)>=MPL: return False
		if c_inD==1:
			#passa að skila ekki false hér því að ég gæti haft c_inD==1 og c_outD==1 og þá vil ég einnig skoða tilfellið fyrir neðan
			#L = helpers.cLen(c,self.k)+self.longest_fw(c_ID,MPL)
			#C er hluti af tip ef C er hluti af isolated þegar skippa X og (til önnur leið áfram með lengd >= MPL eða isPartOfTip(x))
			(x_ID,x_B) = self.getOtherIN(c_ID,c_IN)
			#print c_ID, x_ID
			if self.isPartOfIsolated(MPL,c_ID,skipNodes={x_ID}):
				#print "Skipping x_ID="+str(x_ID)+" makes c_ID="+str(c_ID)+" isolated and therefore potentially a tip"
				existsAltPath = False
				if x_B and self.existsAltPath_fw(MPL,x_ID,c_ID):
					existsAltPath = True
				elif (not x_B) and self.existsAltPath_bw(MPL,x_ID,c_ID):
					existsAltPath = True
				#print "existsAltPath", existsAltPath
				#if self.existsAltPath_fw(MPL,x_ID,c_ID):
				if existsAltPath:
					#print "there exists an alt path fw from x_ID="+str(x_ID)+" skipping c_ID="+str(c_ID)
					#Mark the entire tip before returning True
					#print "found a tip and about to mark it"
					allPaths = self.createAllPathsInContigsCollection(MPL,c_ID,skipNodes={x_ID})
					theContigs = allPaths.getContigIDs()
					self.markAllContigsInSetWithNew_R(theContigs,new_R=2)
					return True
				#if calledFrom!=x_ID:							#To prevent infinite recursions
				#	if self.isPartOfTip(x_ID,MPL,calledFrom=c_ID):	return True	#Endurkvæma kallið sér um að merkja collectionið sem tip
		if c_outD==1:
			(x_ID,x_B) = self.getOtherOUT(c_ID,c_OUT)
			if self.isPartOfIsolated(MPL,c_ID,skipNodes={x_ID}):
				existsAltPath = False
				if x_B and self.existsAltPath_bw(MPL,x_ID,c_ID):
					existsAltPath = True
				elif (not x_B) and self.existsAltPath_fw(MPL,x_ID,c_ID):
					existsAltPath = True
				#print "Skipping X makes C isolated and therefore potentially a tip"
				if existsAltPath:
					#print "found a tip and about to mark it"
					allPaths = self.createAllPathsInContigsCollection(MPL,c_ID,skipNodes={x_ID})
					theContigs = allPaths.getContigIDs()
					self.markAllContigsInSetWithNew_R(theContigs,new_R=2)
					return True
				#if calledFrom!=x_ID:							#To prevent infinite recursions
				#	if self.isPartOfTip(x_ID,MPL,calledFrom=c_ID):	return True	#Endurkvæma kallið sér um að merkja collectionið sem tip
		return False

	#---------------------------------
	#-markAllBubbles and it's helpers:
	#---------------------------------
	def existsPathFrom_n1_to_n2_skipping_intNodes(self,c_ID,end_ID,intNodes1,MPL):
		#Returns intNodes2 or an empty Path
		#c_ID is an integer
		#end_ID is an integer
		#intNodes1 is a Path
		skipNodes = intNodes1.getContigIDs()

		[c, c_IN, c_OUT, c_COV] = self.contigs[c_ID]
		visited = {c_ID}		#a set storing all contigs we have visited
		temp = Path.Path()
		temp.append(c_ID,True,c,self.k)
		#lenFront = helpers.cLen(c,self.k)
		pathsInProgress = Path.PathList([temp])

		#Create all possible paths which include C
		while not pathsInProgress.isEmpty():
			p = pathsInProgress.pop()
			
			#Get all child nodes of the last contig in the path
			(e_ID,e_B) = p.getLastContigTuple()
			DC = self.getDirectChildren(e_ID,e_B)
			DC = helpers.DCorA_toAdd(DC,visited,skipNodes)

			for (x_ID,x_B) in DC:
				if x_ID==end_ID:
					#We have found an alternate path from c_ID to end_ID skipping intNodes1
					#p stores c_ID and intNodes2. if len(p)==1 intNodes2 is empty
					if p.numContigs()>=2:
						return p.allButFirst()
				assert(not x_ID in skipNodes), "made sure of this above"
				assert(not x_ID in visited), "made sure of this above"
				visited.add(x_ID)
				[x, x_IN, x_OUT, x_COV] = self.contigs[x_ID]
				p.append(x_ID,x_B,x,self.k)
				if p.lenOfAllButFirst()<MPL:
					pathsInProgress.append(p)
		return Path.Path()

	def isFront(self, c_ID, MPL):
		#Checks whether C is the front node in a bubble
		#Returns [end_ID,primary,secondary] or an empty list
		#	1. Búa til öll paths sem byrja í C (en ekki hafa C með í þeim)
		#	   og hafa lengd < MPL þegar síðasti contig-inn er ekki talinn með
		#      þ.e. við pælum bara í lengdinni á int. Dæmi:
		#		              front    int    end
		#		C->A->B->D      C->   A->B->   D
		#		C->E->F         C->      E->   F
		#		C->E->G	        C->      E->   G
		#	2. Athuga hvort til alternate path frá front til end
			#print "isFront(c_ID="+str(c_ID)+", MPL="+str(MPL)+")"
		[c, c_IN, c_OUT, c_COV] = self.contigs[c_ID]
		c_inD, c_outD, R = self.degrees[c_ID]
		assert(R in [-1,1,2,3,4]), "The legal options. R="+str(R)
		if R in [1,2]:
			return []
		if c_outD<1:
			return []

		#Set skipNodes as the set containing C and all of it's incoming connections
		skipNodes = {c_ID}
		for (x_ID,x_B) in c_IN:
			skipNodes.add(x_ID)

		firstInts = set()	#Would be {A,E} in the example above
		for (x_ID,x_B) in c_OUT:
			firstInts.add(x_ID)

		allPaths = Path.PathList()		#Will be {[A->B->D, E->F, E->G]|} in the example above
		for x_ID in firstInts:
			xPaths = self.createAllPathsIncludingContig(MPL, x_ID, skipNodes, ignoringLengthOfLast=True)
			allPaths += xPaths

		#Now we have found all potential bubbles. Need to convert them from paths into sets of
		#front_ID, I1, end_ID. Currently each Path contains I1 and end_ID
		
		allPaths = allPaths.filterOutPathsOfLengthLessThan2()	#Make sure each Path has a minimum length of 2 (otherwise I1 will be empty)
		allPaths = allPaths.all2()	#Create all possible lengths from each path (i.e. all possibilities of shortening each path)
		potentialBubbles = allPaths.toListOf_frontIntEnd(c_ID)

		#Now check for alternate paths (i.e. intNodes2)
		for front_ID, I1, end_ID in potentialBubbles:
			I2 = self.existsPathFrom_n1_to_n2_skipping_intNodes(front_ID,end_ID,I1,MPL)
			if not I2.isEmpty():
				#Now we return [end_ID,primary,secondary]
				COV1 = self.pCov(pathAsSetOfIDs=I1.getContigIDs())
				COV2 = self.pCov(pathAsSetOfIDs=I2.getContigIDs())
				if COV1>COV2:
					return [end_ID,I1,I2]
				else:
					return [end_ID,I2,I1]
		return []

	def markBubble(self,front_ID,end_ID,primary,secondary):
		#Helper function which marks a bubble we have found
		#Front_ID and end_ID are the IDs of the front and end respectively
		#primary and secondary are the two paths in the bubble. primary is the path with the higher pCov of the two paths
		#use: self.markBubble(front_ID,end_ID,primary, secondary)
		#print "markBubble(front_ID="+str(front_ID)+", end_ID="+str(end_ID)+", path1="+str(path1)+", path2="+str(path2)+")"
		assert isinstance(front_ID, int), "front_ID must be an integer"
		assert isinstance(end_ID, int), "end_ID must be an integer"
		assert isinstance(primary, Path.Path), "primary must be a Path object"
		assert isinstance(secondary, Path.Path), "secondary must be a Path object"
		assert (not primary.isEmpty()), "primary can't be empty"
		assert (not secondary.isEmpty()), "secondary can't be empty"
		self.degrees[front_ID][2] = 3
		self.degrees[end_ID][2] = 3
		set1 = primary.getContigIDs()
		set2 = secondary.getContigIDs()
		self.markAllContigsInSetWithNew_R(set1,3,listOfTuples=False)
		self.markAllContigsInSetWithNew_R(set2,4,listOfTuples=False)

	def markAllBubbles(self,MPL):
		#After running this function all the contigs in the set should
		#have been marked with 0,1,2,3 or 4 depending on isolated,tip,bubble.
		for c_ID in self.contigs:
			c_inD, c_outD, R = self.degrees[c_ID]
			assert(R in {-1,1,2}), "All contigs should be unanalyzed, isolated or tips when we call this function. R="+str(R)

		#Finnum allar nóður C sem eru Front í bubble. Þegar við finnum slíkt
		#C þá merkjum við C og allar nóðurnar í búbblunni.
		for c_ID in self.contigs:
			c_inD, c_outD, R = self.degrees[c_ID]
			if R==-1:
				self.isPartOfBubble(c_ID,MPL)

		#Mark the contigs which are still unmarked and those we marked as not part of a bubble with a 0 instead of -1 and -2
		for c_ID in self.contigs:
			if self.degrees[c_ID][2] in [-1,-2]:
				self.degrees[c_ID][2] = 0
	
	def isPartOfBubble(self,c_ID,MPL):
		#Returns True if C is part of a bubble. False otherwise
		#Marks the bubble with 3, 4 when returning True. 3 primary. 4 secondary.
		#Marks C with a -2 when returning False to denote that we've ran this function on C
		c_inD, c_outD, R = self.degrees[c_ID]
		assert(R==-1), "R is unmarked. R="+str(R)

		#Finnum allar nóður C sem eru Front í bubble. Þegar við finnum slíkt
		#C þá merkjum við C og allar nóðurnar í búbblunni.
		end_primary_secondary = self.isFront(c_ID,MPL)
		if end_primary_secondary:	#c_ID is the front of a bubble
			[end_ID,primary,secondary] = end_primary_secondary
			self.markBubble(c_ID,end_ID,primary,secondary)
			return True
		else:
			#Look for other potential front nodes
			if c_inD>0:
				PotentialFronts = self.getAllPotentialFrontsFromC(c_ID,MPL)
				for front_ID in PotentialFronts:
					end_primary_secondary = self.isFront(front_ID,MPL)
					if end_primary_secondary:	#There is a bubble starting at front_ID
						[end_ID,primary,secondary] = end_primary_secondary
						self.markBubble(front_ID,end_ID,primary,secondary)

						#If c_ID is part of the bubble: return True
						if ((c_ID == end_ID) or (c_ID in primary) or (c_ID in secondary) ):
							return True

		#Mark that we have ran this function on C before returning False
		self.degrees[c_ID][2] = -2
		return False

	def getAllPotentialFrontsFromC(self,c_ID,MPL):
		#Finds and returns all nodes which are close enough to C to be candidates
		#as first nodes in a bubble which C is a part of
		#A->B->C
		#Here we'd return {(B,True),(A,True)}
		#if len(a)+len(b)+len(c) < MPL
		[c, c_IN, c_OUT, c_COV] = self.contigs[c_ID]
		c_inD, c_outD, R = self.degrees[c_ID]
		if c_inD<1:
			return []
		
		#Find all possible paths with length < MPL ending at C
		potentialFronts = set()
		visited = {c_ID}
		L = helpers.cLen(c,self.k)
		pathsInProgress = [( [(c_ID,True)] , L )]
		underMaxLength = L < MPL
		while pathsInProgress and underMaxLength:
			contigs,L = pathsInProgress.pop()
			(e_ID,e_B) = contigs[0]
			DA = self.getDirectAncestors(e_ID,e_B)
			if DA:
				for (x_ID,x_B) in DA:
					if not (x_ID in visited):
						visited.add(x_ID)
						if (self.degrees[x_ID][2] in [1,2,-2,-3]):
							continue
					else:
						continue
					[x, x_IN, x_OUT, x_COV] = self.contigs[x_ID]
					new_contigs = list(contigs)			#create a copy of contigs
					new_contigs.insert(0 , (x_ID,x_B))	#add x_ID to the front of the new path of contigs
					new_L = L+helpers.cLen(x,self.k)
					if new_L<MPL:
						pathsInProgress.append( (new_contigs,new_L) )
						potentialFronts.add(x_ID)	#sleppi fyrsta og síðasta (þ.e. sleppi front og end því við geymum þá sér)
		return list(potentialFronts)

	#-----breytingar--------------------------
	def markAllInGenome(self,kmersInGenome):
		#markAllContigsInSetWithNew_R(self,contigSet,new_R,listOfTuples=False)
		#Uppfærum ratings á öllum contigs sem voru rate-aðir með 0 eða 3
		#Breytum ekki ratingi á contigs sem voru með 1, 2 eða 4
		#Ný rating:
		#	0: complex
		#	5: genomic
		#	6: partial
		for c_ID, [c,c_IN,c_OUT,c_COV] in self.contigs.iteritems():
			R = self.degrees[c_ID][2]
			if R in [0,3]:
				genomic = True
				partial = False
				#L = helpers.cLen(c,self.k)	#number of k-mers in c
				#inGenomeCounter = 0
				for km in dbg.kmers(c,self.k):
					km = min(km, dbg.twin(km))
					if not (km in kmersInGenome):
						genomic = False
					else:
						partial = True
				if genomic:
					#assert(L==inGenomeCounter)
					self.degrees[c_ID][2] = 5
				elif partial:
					self.degrees[c_ID][2] = 6
				else:
					self.degrees[c_ID][2] = 0

	#----------------------------------------------------
	#---analyzeAllContigsInCollection and it's helpers:--
	#----------------------------------------------------
	def createDegreeTable(self):
		def getDegreesOfContig(c_ID,c_IN,c_OUT):
			inD = 0
			for (x_ID,x_B) in c_IN:
				if x_ID!=c_ID:
					inD += 1
			outD = 0
			for (x_ID,x_B) in c_OUT:
				if x_ID!=c_ID:
					outD += 1
			return inD,outD

		self.degrees = collections.defaultdict(list)
		for c_ID,[c,c_IN,c_OUT,c_COV] in self.contigs.iteritems():
			c_inD,c_outD = getDegreesOfContig(c_ID,c_IN,c_OUT)
			self.degrees[c_ID] = [c_inD,c_outD, -1]

	#kmersInGenome is a set of all rep_km of k-mers in the genome
	def analyzeAllContigsInCollection(self,kmersInGenome):
		#	self.degrees[c_ID] -> [c_inD,c_outD, R]
		#	R stores info about whether C is part of a collection of contigs which is isolated, a tip or a bubble
		#	R = -1:	We haven't analyzed C yet
		#	R = 1: C is a part of an isolated collection of contigs
		#	R = 2: C is a part of a tip
		#	R = 3: C is a primary path in a bubble
		#	R = 4: C is a secondary path in a bubble
		#	R = 0: We were unable to identify C as any of the above

		#Create and initialize self.degrees for all contigs in the graph:
		self.createDegreeTable()
		for c_ID in self.contigs:
			assert(self.degrees[c_ID][2]==-1), "All contigs should initially be marked with a -1"

		#Set the maximum length for a path to be part of an isolated collection
		#or a tip. This is also the maximum length through a bubble
		MPL = self.k+1

		#Find and mark all isolated contigs in the graph:
		isoMPL = MPL
		#print "We're using MPL="+str(isoMPL)+" for markAllIsolated instead of MPL="+str(MPL)
		self.markAllIsolated(MPL=isoMPL)

		#Make sure the ratings are as expected after marking all isolated contigs
		for c_ID in self.contigs:
			assert(self.degrees[c_ID][2] in [-1,1]), "All contigs should either be marked as isolated, or still with a -1"

		#Find and mark all tips in the graph:
		for c_ID in self.contigs:
			R = self.degrees[c_ID][2]
			if R==-1:
				self.isPartOfTip(c_ID,MPL)	#Merkir í leiðinni C með 2 ef C er tip
			else:
				assert(R in [1,2]), "All contigs should be marked with a -1, 1 or 2. R="+str(R)

		#Make sure the ratings are as expected after marking all tips
		for c_ID in self.contigs:
			R = self.degrees[c_ID][2]
			assert(R in [-1,1,2]), "All contigs should either be marked as isolated, tips or with a -1"

		#Find and mark all bubbles:
		bubMPL = MPL
		#print "We're using MPL="+str(bubMPL)+" for markAllBubbles instead of MPL="+str(MPL)
		self.markAllBubbles(MPL=bubMPL)
		
		#Make sure the ratings are as expected after marking all bubbles
		for c_ID in self.contigs:
			R = self.degrees[c_ID][2]
			assert(R in [0,1,2,3,4]), "R="+str(R)+". We should have marked all contigs as part of isolated, tip, bubble or non of the three"

		self.markAllInGenome(kmersInGenome)
		#Make sure the ratings are as expected
		for c_ID in self.contigs:
			R = self.degrees[c_ID][2]
			assert(R in [0,1,2,4,5,6]), "R="+str(R)

	#---------------------------------------------------------------
	#---------------------------------------------------------------
	#Helpful functions for working with isolated tips and bubbles:
	#---------------------------------------------------------------
	#---------------------------------------------------------------
	def getDirectAncestors(self,c_ID,c_B,skipNodes=set()):
		A = set()
		[c, c_IN, c_OUT, c_COV] = self.contigs[c_ID]
		if c_B:
			for (x_ID,x_B) in c_IN:
				if not x_ID in skipNodes:
					A.add((x_ID,c_B&x_B))
		elif not c_B:
			for (x_ID,x_B) in c_OUT:
				if not x_ID in skipNodes:
					A.add((x_ID,c_B&x_B))
		return A

	def getDirectChildren(self,c_ID,c_B,skipNodes=set()):
		A = set()
		[c, c_IN, c_OUT, c_COV] = self.contigs[c_ID]
		if c_B:
			for (x_ID,x_B) in c_OUT:
				if not x_ID in skipNodes:
					A.add((x_ID,c_B&x_B))
		elif not c_B:
			for (x_ID,x_B) in c_IN:
				if not x_ID in skipNodes:
					A.add((x_ID,c_B&x_B))
		return A

	def markAllContigsInSetWithNew_R(self,contigSet,new_R,listOfTuples=False):
		#print "markAllContigsInSetWithNew_R(contigSet="+str(contigSet)+", new_R="+str(new_R)+", listOfTuples="+str(listOfTuples)+")"
		if listOfTuples:
			for (c_ID,c_B) in contigSet:
				assert(c_ID in self.contigs), "c_ID must be in the graph. c_ID="+str(c_ID)
				assert(c_ID in self.degrees), "c_ID must have been assigned a degree. c_ID="+str(c_ID)
				self.degrees[c_ID][2] = new_R
		else:
			for c_ID in contigSet:
				assert(c_ID in self.contigs), "c_ID must be in the graph. c_ID="+str(c_ID)
				assert(c_ID in self.degrees), "c_ID must have been assigned a degree. c_ID="+str(c_ID)
				self.degrees[c_ID][2] = new_R

	def findAllConnectedContigs(self,c_ID,skipNode=-1):
		#returns all contigs which can be reached from C
		#print "findAllConnectedContigs(c_ID="+str(c_ID),", skipNode="+str(skipNode)+")"
		visited = set()
		toExplore = {c_ID}
		if skipNode!=-1:
			visited.add(skipNode)
		while toExplore:
			e_ID = toExplore.pop()
			[e, e_IN, e_OUT, e_COV] = self.contigs[e_ID]
			visited.add(e_ID)
			for (x_ID,x_B) in e_IN:
				if not x_ID in visited:
					toExplore.add(x_ID)
			for (x_ID,x_B) in e_OUT:
				if not x_ID in visited:
					toExplore.add(x_ID)
		if skipNode!=-1:
			visited.remove(skipNode)
			assert(not skipNode in visited), "skipNode should not be in the set of visited nodes"
		assert(c_ID in visited), "bla"
		return visited

	def pCov(self,pathAsSetOfIDs):
		#We define pCov(path P) as the sum of cCov(C) for all contigs C in P.
		assert(isinstance(pathAsSetOfIDs,set))
		sum = 0
		for c_ID in pathAsSetOfIDs:
			[c,c_IN,c_OUT,c_COV] = self.contigs[c_ID]
			sum+=c_COV
		return sum

if __name__ == "__main__":
	k = 5
	MPL = 2*k
	G = Graph(k)
	G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
	G.contigs[G.getID()] = ["CCATTT",[(0,True)],[],0]       #2 k-mers
	G.addKmersFromAllContigs()
	G.createDegreeTable()
	P0_1 = G.createAllPathsIncludingContig(MPL, 0, {1})
	#P1_0 = G.createAllPathsIncludingContig(MPL, 1, [0])
	print P0_1
	#print P1_0
	"""
	G = Graph(5)
	G.contigs[G.getID()] = ["CAGGT",[],[(1,True)],0]
	G.contigs[G.getID()] = ["AGTTA",[(0,True)],[],0]
	G.addKmersFromAllContigs()
	MPL = 2*G.k 
	G.degrees = collections.defaultdict(list)
	for c_ID,[c,c_IN,c_OUT,c_COV] in G.contigs.iteritems():
		c_inD,c_outD = G.getDegreesOfContig(c_ID,c_IN,c_OUT)
		G.degrees[c_ID] = [c_inD,c_outD, -1]

	print G.isPartOfTip(0,MPL)
	#G = Graph(5)
	#G.addSegmentToGraph("AAACGGGTTGGCGATTTG")
	#G.addSegmentToGraph("CCCCC")
	#print G.findAllConnectedContigs(30)
	#print G.findAllConnectedContigs(31)
	#print G.isPartOfIsolated("CCCCC")
	#print G.isPartOfIsolated("TTGAA")
	#print G.LP_fw(0, 2*5, 0)
	#print G.isPartOfIsolated("AAACG")
	#print G.isPartOfIsolated("CCCCC")
	"""


	'''
	G = Graph(7)
	#Just so the bubble isn't isolated
	G.addSegmentToGraph("AAGCGTTGGTAAGTCAGTGGATC")	
	#The bubble itself
	G.addSegmentToGraph("TTGGATC")
	G.addSegmentToGraph( "TGGATCTTAGACT")
	G.addSegmentToGraph( "TGGATCATAGACT")
	G.addSegmentToGraph(        "TAGACTT")
	G.printContigs()
	G.analyzeAllContigsInCollection()
	G.printContigsWithRatings()
	'''

	
