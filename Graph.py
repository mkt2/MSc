#coding:utf8
import collections, sys
import dbg
import compareGraphs
#an object of this class represents a legal DBG graph
class Graph:
	def __init__(self,k,pfn=False,ps=False,al=True,pil=False,printInit=False):
		if printInit:
			print "Initializing the Graph: k="+str(k)+", pfn="+str(pfn)+", ps="+str(ps)+", al="+str(al)+", pil="+str(pil)+", printInit="+str(printInit)
		self.k = k 		#the kmer length
		self.ID = 0		#the lowest available contig ID
		self.halfAdded = -99	#a number representing that a kmer has
								#not been fully added to the kmerDict

		#contigID -> [contig, IN, OUT, coverage]
		self.contigs = collections.defaultdict(list)
		#dictionary of the contigs in our graph with information
		#over which contigs they connect to
		#IN: List of tuples. Example: IN = [(0,True),(1,False),(99,True)]
		#OUT: List of tuples

		#kmer -> [contigID, index, B]
		#B=False:
		#	kmer is in the twin of the contig with contigID
		#	index represents its location in the twin
		#B=True:
		#	kmer is in the contig with contigID
		#	index represents its location in the contig
		self.kmers = collections.defaultdict(list)
		#dictionary of the kmers in our graph, and where each of them occurs
		#in the contigs

		self.printFunctionNames = pfn
		self.printStatus = ps
		self.assertLegal = al
		self.printIsLegal = pil
		self.smartSplit = False		#When self.smartSplit==True then we only split a contig when we have been told to do so before

	def __len__(self):
		return len(self.contigs)

	def __contains__(self,contig):
		#this is a comment
		return contig in self.contigs

	def numKmerPairs(self):
		return len(self.kmers)/2
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

	def printToFile(self,f):
		f.write(str(self.k)+"\n")
		f.write(str(self.ID)+"\n")
		f.write(str(len(self))+"\n")
		f.write("ID;CONTIG;IN;OUT;COV\n")
		for ID, [c,IN,OUT,COV] in self.contigs.iteritems():	
			f.write(str(ID)+";"+str(c)+";"+str(IN)+";"+str(OUT)+";"+str(COV)+"\n")

	#Before:	The file is of the correct form
	#				0 k
	#				1 ID
	#				2 numberOfContigs
	#				3 ID;CONTIG;IN;OUT;COV
	#				...data of the form ID;CONTIG;IN;OUT;COV
	#After:
	#	for all lines we set self.contigs[ID] as [CONTIG,IN,OUT,COV]
	def createGraphFromFile(self,fileName):
		print "createGraphFromFile(fileName="+str(fileName)+")"
		f = open("Output/"+fileName, 'r')
		k = int(f.readline())
		ID = int(f.readline())
		numberOfContigs = int(f.readline())
		h = f.readline()
		assert h=="ID;CONTIG;IN;OUT;COV\n"
		assert isinstance(k, int)
		assert isinstance(ID, int)
		assert isinstance(numberOfContigs, int)
		assert k==self.k
		for i, line in enumerate(f):
			#print line
			[ID,CONTIG,IN,OUT,COV] = line.strip().split(";")
			#print ID,CONTIG,IN,OUT,COV
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
		if not self.isLegalDBG():
			self.printContigs()
			raise Exception("The graph we read from the file is not legal!")
		print "Done creating the Graph from the file"

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

	#create G_naive from the same contigs as G consists of
	#We can use self.kmers to get all the kmers from the contigs
	def createNaive(self):
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
		return compareGraphs.isSameGraph(self,G_naive,False,self.printFunctionNames,self.printStatus,relaxAssertions=False)

	#--------------------------------------------------------------------------
	#-----------------Simple functions for working with graphs-----------------
	#--------------------------------------------------------------------------

	#Before:	ID is the ID of a contig c in the graph
	#			IN_list is a list of tuples. E.g. [(0,True)]
	#After:		the IN of c has been set as IN_list
	def setIN(self,ID,IN_list):
		if self.printFunctionNames:
			print "setIN(ID="+str(ID)+", IN_list="+str(IN_list)+")"
		temp = self.contigs[ID]
		temp[1] = IN_list

	def setOUT(self,ID,OUT_list):
		if self.printFunctionNames:
			print "setOUT(ID="+str(ID)+", OUT_list="+str(OUT_list)+")"
		temp = self.contigs[ID]
		temp[2] = OUT_list

	def increaseCOV_by(self,ID,x):
		if self.printFunctionNames:
			print "increaseCOV_by(ID="+str(ID)+", x="+str(x)+")"
		temp = self.contigs[ID]
		self.contigs[ID] = [temp[0],temp[1],temp[2],temp[3]+x]

	#Before:	ID is the ID of a contig c in the graph
	#			IN_tuple is a tuple. E.g. (0,True)
	#After:		IN_tuple has been added to c's IN
	def addIN(self,ID,IN_tuple):
		if self.printFunctionNames:
			print "(addIN(ID="+str(ID)+", IN_tuple="+str(IN_tuple)+")"
		temp = self.contigs[ID]
		temp[1].append(IN_tuple)

	def addOUT(self,ID,OUT_tuple):
		if self.printFunctionNames:
			print "(addOUT(ID="+str(ID)+", OUT_tuple="+str(OUT_tuple)+")"
		temp = self.contigs[ID]
		temp[2].append(OUT_tuple)

	"""
	#Before:	ID is the ID of a contig c in the graph
	#			IN_tuple is a tuple. E.g. (0,True)
	#After:		IN_tuple has been deleted from c's IN (if it is one of c's IN)
	def deleteIN(self,ID,IN_tuple):
		if self.printFunctionNames:
			print "deleteIN(ID+"+str(ID)+", IN_tuple="+str(IN_tuple)+")"
		[c,IN,OUT,COV] = self.contigs[ID]
		for t in IN:
			if t==IN_tuple:
				IN.remove(IN_tuple)
				break
		self.contigs[ID] = [c,IN,OUT,COV]

	def deleteOUT(self,ID,OUT_tuple):
		if self.printFunctionNames:
			print "deleteOUT(ID+"+str(ID)+", OUT_tuple="+str(OUT_tuple)+")"
		[c,IN,OUT,COV] = self.contigs[ID]
		for t in OUT:
			if t==OUT_tuple:
				OUT.remove(OUT_tuple)
				break
		self.contigs[ID] = [c,IN,OUT,COV]
	"""

	def ID_has_IN(self,ID,IN_tuple):
		return IN_tuple in self.contigs[ID][1]

	def ID_has_OUT(self,ID,IN_tuple):
		return IN_tuple in self.contigs[ID][2]


	def addKmersFromAllContigs(self):
		for ID, values in self.contigs.iteritems():
			self.addKmersFromContig(ID,values[0])

	#Before: 	c is a contig with ID ID
	#After:		all kmers from c and their twins have been added to the kmerDict
	def addKmersFromContig(self,ID,c="-1"):
		if c=="-1":
			c = self.contigs[ID][0]
		for i, km in enumerate(dbg.kmers(c,self.k)):
			self.kmers[km] = [ID,i,True]
		for i, km in enumerate(dbg.kmers(dbg.twin(c),self.k)):
			self.kmers[km]= [ID,i,False]

	#Before: 	c is a contig with ID ID
	#After:		all kmers from c and their twins have been deleted from self.kmers
	def deleteKmersFromContig(self,ID,c="-1"):
		if c=="-1":
			c = self.contigs[ID][0]
		for km in dbg.kmers(c,self.k):
			if km in self.kmers:
				del self.kmers[km]
				del self.kmers[dbg.twin(km)]

	#returns the lowest available ID in the graph. Maintains self.ID as
	#the lowest available ID
	def getID(self):
		ID = self.ID
		self.ID += 1
		return ID

	def setID(self,newID):
		self.ID = newID

	#returns true if the graph is empty. False otherwise
	def isEmpty(self):
		return len(self.contigs)==0

	def hasNoKmers(self):
		return len(self.kmers)==0



	def isLegalDBG(self, skipKmers=False):
		if self.printIsLegal:
			print "isLegalDBG(skipKmers="+str(skipKmers)+")"

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
		if self.printIsLegal:
			print "isLegalGraph()"
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
		if self.printIsLegal:
			print "isLegal_kmerDict()"
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
		if self.printIsLegal:
			print "isLegal_ID()"
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
		if self.printIsLegal:
			print "isLegal_IN_and_OUT()"
		#make sure connections go both ways:
		for c_ID in self.contigs:
			[c,c_IN,c_OUT,c_COV] = self.contigs[c_ID]
			#print "c_ID:",c_ID
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


	#Before:	self is a Graph object
	#After:		if some contig contains a connection to itself
	#			return False if it's not of the correct form. [(0,False),(0,False)] instead of
	#			[(0,False)] for example
	def isLegal_connectionsToSelf(self):
		if self.printIsLegal:
			print "isLegal_connectionsToSelf()"
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
		if self.printIsLegal:
			print "isLegal_merge()"
		for ID, [c,IN,OUT,COV] in self.contigs.iteritems():
			[canMerge1,a_ID,b_ID,A,B] = self.canMergeOUT(ID)
			[canMerge2,a_ID,b_ID,A,B] = self.canMergeIN(ID)
			if canMerge1 or canMerge2:
				if self.printIsLegal:
					print "isLegal_merge returning False"
					print "contig ID: "+str(ID)+" can merge"
				return False
		return True



	#Before:	aID and bID are the IDs of contigs a and b
	#			a and b are contigs already added to the graph
	#			A and B are Booleans
	#Note:		It is possible that aID=bID. I.e. we may be connecting a to itself
	#After:		a and b have been connected together
	#			A and B represent how we're connecting a and b:
	#			| A=T, B=T |   A=T, B=F  |  A=F, B=T  |     A=F, B=F     |
	#			|   a->b   |  a->twin(b) | twin(a)->b | twin(a)->twin(b) |
	def connect_a_to_b(self,aID,bID,A,B):
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
		#update COV:
		#[a,a_IN,a_OUT,a_COV] = self.contigs[aID]
		#[b,b_IN,b_OUT,b_COV] = self.contigs[bID]
		#self.increaseCOV_by(aID, int(((self.k-1)/float(len(b)))*b_COV))
		#self.increaseCOV_by(bID, int(((self.k-1)/float(len(a)))*a_COV))




	#Function that connects a and b if and only if they are not already connected
	def connect_a_to_b_ifNotAlreadyConnected(self,aID,bID,A,B):
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





	#--------------------------------------------------------------------------
	#--------------More complex functions for working with graphs--------------
	#--------------------------------------------------------------------------

	#--------------------changeID_FromTo and helper functions--------------------

	#Before:	ID is the id of contig c in the legal DBG
	#After:		contig c now has the id ID_new and all appropriate changes have been made
	#				the INs/OUTs of c now connect to ID_new instead of ID
	#				the graph is a legal DBG
	def changeID_FromTo(self,ID,ID_new):
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

	#--------------------changeID_FromTo functions finished--------------------

	#Before:	The graph is a legal DBG and
	#			self.contigs[ID] = [c,IN,OUT,Cov]
	#After:		The graph is still a legal DBG and
	#			self.contigs[ID] = [twin(c),IN*,OUT*,COV]
	#			The INs/OUTs of c have been changed so that they connect
	#			to c* instead of c where c* is c after flipping
	def flipContig(self,c_ID):
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
	#--------------------addSegmentToGraph and subfunctions--------------------
	#--------------------------------------------------------------------------

	#Before:	segment is a DNA string
	#			this Graph is a legal DBG graph
	#After:		segment has been added to the graph
	#			the kmers from segment have been added to self.kmers
	#			the graph is a legal DBG graph
	def addSegmentToGraph(self,segment):
		#print "addSegmentToGraph(segment="+str(segment)+")"
		#self.printContigs("bla")
		if self.printFunctionNames:
			print "addSegmentToGraph(segment="+str(segment)+")"
		if self.printStatus:
			self.printContigs("inside addSegmentToGraph")
			#self.printKmers("Kmers")
		if self.assertLegal:
			assert self.isLegalDBG()
		if len(segment) < self.k:
			return

		start = 0
		tmp = collections.defaultdict(list)
		for i, km in enumerate(dbg.kmers(segment,self.k)):
			B1 = km in self.kmers
			B2 = km in tmp
			if not (B1 or B2):
				#We have not seen this kmer before
				tmp[km] = [None,None,None]
				tmp[dbg.twin(km)] = [None,None,None]
			else:
				#We have seen this kmer before
				if B1:
					[contigID, index, B] = self.kmers[km]
					self.contigs[contigID][3] += 1
				if i-start>=1:
					self.addSegmentWithNoSeenKmers(segment[start:i-1+self.k])
				start = i+1
			
		self.addSegmentWithNoSeenKmers(segment[start:])
		if self.assertLegal:
			assert self.isLegalDBG()



	#Helper function for addSegmentToGraph
	#Before:	segment is a DNA string
	#			this Graph is a legal DBG graph
	#			The kmers from segment (or its twin) are neither already in
	#			the graph or in a different place in segment
	#			No kmer from segment is in self.kmers
	#			Note:
	#				the fw/bw of some kmer in segment may be in the graph
	#				the fw/bw of some kmer in segment may be another kmer in segment
	#					(i.e. segment might connect to itself)
	#After:		Same as addSegmentToGraph
	#			(addSegmentToGraph(S) finds which parts of S we haven't seen before and
	#			calls addSegmentWithNoSeenKmers to actually add those parts)
	def addSegmentWithNoSeenKmers(self,segment):
		#print "addSegmentWithNoSeenKmers(segment="+str(segment)+")"
		#self.printContigs("blee")
		if self.printFunctionNames:
			print "addSegmentWithNoSeenKmers(segment="+str(segment)+")"
		if self.printStatus:
			self.printContigs("inside addSegmentWithNoSeenKmers")
		if self.assertLegal:
			assert self.isLegalDBG()
		L = len(segment)
		if L < self.k:
			return
		elif L == self.k:
			self.addSegmentAlreadySplit(segment)
			return
		else:
			twinSegment = dbg.twin(segment)
			for i, km in enumerate(dbg.kmers(segment,self.k)):
				if self.printStatus:
					print i, km
				for y in dbg.fw(km):
					if y in segment:
						if self.printStatus:
							print "km:"+str(km), "fw(km) in segment. y:",str(y)
						if y==km:
							#km connects to itself
							if i==0:
								s0, s1 = splitString(segment,self.k,self.k)
								self.addSegmentAlreadySplit(s0)
								self.addSegmentWithNoSeenKmers(s1)
								return
							else:
								s0, s1 = splitString(segment,i+self.k-1,self.k)
								self.addSegmentWithNoSeenKmers(s0)	#Má ég segja addSegmentAlreadySplitOnSelf(s0) hér?
								self.addSegmentWithNoSeenKmers(s1)
								return
						else:
							#km connects to another km in segment
							if i==L-self.k:
								pass
								#raise Exception("We should have done this split as a bw of a previous kmer earlier in the code")
							elif not (y==segment[i+1:i+1+self.k]):	#viljum ekki splitta ef y er næsti kmer í segment
								s0, s1 = splitString(segment,i+self.k,self.k)
								self.addSegmentWithNoSeenKmers(s0)
								self.addSegmentWithNoSeenKmers(s1)
								return

					if y in twinSegment:
						if self.printStatus:
							print "km:"+str(km), "fw(km) in twin(segment). y:",str(y)
						if y==dbg.twin(km):
							#km connects to the twin of itself
							if i==0:
								s0, s1 = splitString(segment,self.k,self.k)
								self.addSegmentAlreadySplit(s0)
								self.addSegmentWithNoSeenKmers(s1)
								return
							elif i<L-self.k:			#i=L-self.k er indexid a sidasta kmernum
								s0, s1 = splitString(segment,i+self.k,self.k)
								self.addSegmentWithNoSeenKmers(s0)
								self.addSegmentWithNoSeenKmers(s1)
								return
							else:
								pass					#don't need to split if the last kmer connects to it's twin
						else:
							#km connects to the twin of another km in segment
							if i==L-self.k:
								raise Exception("km->twinSegment. This can't happen. We should have split earlier. km="+str(km)+", y="+str(y))
							else:
								s0, s1 = splitString(segment,i+self.k,self.k)
								self.addSegmentWithNoSeenKmers(s0)
								self.addSegmentWithNoSeenKmers(s1)
								return

				for y in dbg.bw(km):
					if y in segment:
						if self.printStatus:
							print "km:"+str(km), "bw(km) in segment. y:",str(y)
						if y==km:
							#km connects to itself
							pass			#already covered in fw
						else:
							#another km in segment connects to km
							if i==0:
								pass		#don't need to split
							elif not (y==segment[i-1:i-1+self.k]):
								if i<L-self.k:
									s0, s1 = splitString(segment,i+self.k-1,self.k)
									self.addSegmentWithNoSeenKmers(s0)
									self.addSegmentWithNoSeenKmers(s1)
									return
								elif i==L-self.k:
									raise Exception("segment->km. This can't happen. We should have found this split with a fw(km) earlier on. km="+str(km)+", y="+str(y))

					if y in twinSegment:
						if self.printStatus:
							print "km:"+str(km), "bw(km) in twin(segment). y:",str(y)+", i:"+str(i)
						if y==dbg.twin(km):
							#twin(km) connects to km
							if i==0:
								pass
							else:
								s0, s1 = splitString(segment,i+self.k-1,self.k)
								self.addSegmentWithNoSeenKmers(s0)
								self.addSegmentWithNoSeenKmers(s1)
								return
						else:
							#twin(another kmer) connects to km
							if i==0:
								if self.printStatus:
									print "twin(another km) connects to the first km in segment. We'll make bw of that later kmer do the split"
								pass		#bw of a later km will catch this (then we'll need to split)
							else:
								s0, s1 = splitString(segment,i+self.k-1,self.k)
								self.addSegmentWithNoSeenKmers(s0)
								self.addSegmentWithNoSeenKmers(s1)
								return

		#If we didn't need to split anywhere we know that segment
		#doesn't need to be split further due to connections to itself
		if self.printStatus:
			print "Adding due to no split"				
		self.addSegmentAlreadySplitOnSelf(segment)
		if self.assertLegal:
			assert self.isLegalDBG()


	#Before:	segment does not need to be split further because of connections to itself
	#			We might need to split segment further due to connections to other segments
	#			in the graph
	def addSegmentAlreadySplitOnSelf(self,segment):
		if self.printFunctionNames:
			print "addSegmentAlreadySplitOnSelf(segment="+str(segment)+")"
		if self.printStatus:
			self.printContigs("inside addSegmentAlreadySplitOnSelf")
			#self.printKmers("Kmers")
		if self.assertLegal:
			assert self.isLegalDBG()
		L = len(segment)
		if L < self.k:
			return
		if L == self.k:
			self.addSegmentAlreadySplit(segment)
			return
		else:
			twinSegment = dbg.twin(segment)	
			for i, km in enumerate(dbg.kmers(segment,self.k)):
				if self.printStatus:
					print i, km
				for y in dbg.fw(km):
					#km i in segment connects to another segment
					if (y in self.kmers) and (not i==L-self.k):	#we split at i+k unless km is the last kmer in segment, then we skip
						if self.printStatus:
							print "fw, y in kmers. y =", y
						#return i+self.k
						s0, s1 = splitString(segment,i+self.k,self.k)
						self.addSegmentAlreadySplitOnSelf(s0)
						self.addSegmentAlreadySplitOnSelf(s1)
						return
				for y in dbg.bw(km):
					#another segment connects to km i in segment
					if (y in self.kmers) and (not i==0): #We split at i+k-1 unless km is the first km in segment, then we skip
						if self.printStatus:
							print "bw, y in kmers. y =", y
						#return i+self.k-1
						s0, s1 = splitString(segment,i+self.k-1,self.k)
						self.addSegmentAlreadySplitOnSelf(s0)
						self.addSegmentAlreadySplitOnSelf(s1)
						return

		#If we didn't need to split anywhere we know that segment
		#doesn't need to be split further
		if self.printStatus:
			print "Adding due to no split"				
		self.addSegmentAlreadySplit(segment)
		if self.assertLegal:
			assert self.isLegalDBG()


	#helper function for addSegmentWithNoSeenKmers
	#Before:	segment is a DNA string
	#			The graph is a legal DBG graph
	#			the kmers in segment have not been seen before, neither in previous
	#			contigs or segment itself
	#			Segment does not need to be split up further
	#				(did that in addSegmentWithNoSeenKmers)
	#				We may have to split up contigs which the ends of segment can connect to
	#After:
	#			Segment has been added to the graph
	#			the kmers from segment have been added to self.kmers
	#			All contigs in this graph (including segment) have been split up as much as is needed
	#			Note:
	#				The graph is NOT necessary a legal DBG Graph (we might need to merge segment with it's connections)
	def addSegmentAlreadySplit(self,segment):
		if self.printFunctionNames:
			print "addSegmentAlreadySplit(segment="+str(segment)+")"
		if self.printStatus:
			self.printContigs("inside addSegmentAlreadySplit")
			#self.printKmers("Kmers")
		if self.assertLegal:
			assert self.isLegalDBG()
		if len(segment) < self.k:
			if self.printStatus:
				print "We don't add segment because it's too short"
			return

		self.splitOthers(segment)

		#give segment an ID and add it to self.contigs
		ID = self.getID()
		self.contigs[ID] = [segment,[],[],len(segment)-self.k+1]

		#add the kmers from segment to self.kmers
		self.addKmersFromContig(ID,segment)

		if self.printStatus:
			self.printContigs("After calling splitOthers(segment) inside addSegmentAlreadySplit and adding segment to the graph")
			self.printKmers()

		#connect segment to other contigs in the graph
		self.connectSegment(ID,segment)

		if self.assertLegal:
			assert self.isLegalGraph()	#The graph is not a legal DBG because we might be able to merge some segments

		#merge with connecting contigs where possible
		self.mergeSegment(ID)

		#Now segment has been fully added to the graph.
		#Note: The graph is now a legal DBG so we can now think of segment as a contig

		if self.assertLegal:
			assert self.isLegalDBG()

	#Helper function for addSegmentAlreadySplit
	#Before:	segment has NOT been added to the graph (or connected to any contigs)
	#			We will not have to split segment further in order to add it to the graph
	#			NO kmer from segment is already in self.kmers
	#			some kmers in the fw/bw of segment might already be in some contigs already
	#			in the graph (in those cases we have to split the previous contigs)
	#After:		segment is unchanged (no recursive calls)
	#			We have split all other contigs in the graph wherever segment can connect to a kmer in them.
	#				Note: We don't split contigs where segment can connect to their first/last kmers
	#			i.e. NO more contigs need to be split in order to add segment to the graph
	def splitOthers(self,segment,relaxAssertions=False):
		if self.printFunctionNames:
			print "splitOthers(segment="+str(segment)+", relaxAssertions="+str(relaxAssertions)+")"
		if self.printStatus:
			self.printContigs("inside splitOthers")
			#self.printKmers("kmers")
		if self.assertLegal and not relaxAssertions:
			assert self.isLegalDBG()
		elif self.assertLegal and relaxAssertions:
			assert self.isLegalGraph()
		splitSomething = False
		#Case 1 and 2
		for y in dbg.fw(segment[len(segment)-self.k:]):
			if y in self.kmers:
				if self.printStatus:
					print "y in fw(last). y:", y
				#Note: Since segment hasn't been added to graph we know that this is returning true because y is in some other contig than segment itself.
				#	Therefore we don't need to worry about splitting segment itself
				#kmer -> [contigID, index, B]

				#Get information about the contig c where we've seen y before
				[c_ID, c_i, c_B] = self.kmers[y]
				if self.printStatus:
					print c_ID, c_i, c_B
				[c, c_IN, c_OUT, c_COV] = self.contigs[c_ID]
				c_L = len(c)

				#Find where to split c
				if c_i == 0:
					continue	#Don't need to split because c_i is the location of the first kmer in c (or twin(c)).

				if c_B:
					splitIndex = c_i + self.k - 1
					#splitIndex = c_i + self.k
				else:
					c_i = c_L - self.k - c_i 	#let c_i be the index of y in c instead of in twin(c)
					splitIndex = c_i + self.k

				#Split c into c[0:splitIndex] -> c[splitIndex-k:L]
				self.splitContig(c_ID, splitIndex, relaxAssertions)
				splitSomething = True

		#If we split in the first loop we end the function and call it again
		#instead of going into the second loop
		if splitSomething:
			if self.printStatus:
				print "We call splitOthers(segment) again because splitOthers(segment) did at least one split"
			self.splitOthers(segment,relaxAssertions=True)
			if self.assertLegal:
				assert self.isLegalGraph()
			return
		else:
			if self.assertLegal and not relaxAssertions:
				assert self.isLegalDBG()


		#Case 3 and 4
		for y in dbg.bw(segment[:self.k]):
			if y in self.kmers:
				if self.printStatus:
					print "y in bw(first). y:", y
				#Get information about the contig c where we've seen y before
				[c_ID, c_i, c_B] = self.kmers[y]
				if self.printStatus:
					print c_ID, c_i, c_B
				[c, c_IN, c_OUT, c_COV] = self.contigs[c_ID]
				c_L = len(c)

				#Find where to split c
				if c_i == c_L-self.k:
					continue	#Don't need to split because c_i is the location of the last kmer in c (or twin(c)).

				if c_B:
					#splitIndex = c_i + self.k - 1
					splitIndex = c_i + self.k
				else:
					if self.printStatus:
						print c_L, c_i
					c_i = c_L - self.k - c_i 	#let c_i be the index of y in c instead of in twin(c)
					splitIndex = c_i + self.k - 1

				#Split c into c[0:splitIndex] -> c[splitIndex-k:L]
				self.splitContig(c_ID, splitIndex, relaxAssertions)
				splitSomething = True

		if splitSomething:
			if self.printStatus:
				print "We call splitOthers(segment) again because splitOthers(segment) did at least one split"
			self.splitOthers(segment,relaxAssertions=True)
			if self.assertLegal:
				assert self.isLegalGraph()
		else:
			if self.assertLegal and not relaxAssertions:
				assert self.isLegalDBG()

		#Note: Quite surprisingly it seems that we select splitIndex the same way in all cases. If testing confirms this
		#then we can simplify this function A LOT! :D

		#Note: There are technically 4 more cases but since a->b = twin(b)->twin(a) etc
		#they are already covered

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
	def splitContig(self,ID,i,relaxAssertions=False):
		#Need to test this for special cases. Might need to handle connections from c->c in
		#a different way
		if self.printFunctionNames:
			print "splitContig(ID="+str(ID)+", i="+str(i)+", relaxAssertions="+str(relaxAssertions)+")"
		if self.printStatus:
			print self
		if self.assertLegal:
			if relaxAssertions:
				assert self.isLegalGraph()
			else:
				assert self.isLegalDBG()
		[c,IN,OUT,COV] = self.contigs[ID]

		c0, c1 = splitString(c,i,self.k)
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
		self.contigs[ID0] = [c0, IN, [(ID1,True)], COV/2]
		self.contigs[ID1] = [c1, [(ID0,True)], OUT, COV/2]
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

		#OUT[:] = [y if y != (ID,False) else (ID1,False) for y in OUT]
		#now do the same for twin(c)->c
		#IN[:] = [y if y != (ID,False) else (ID0,False) for y in IN]
		#now do the same for c->c:
		#IN[:] = [y if y != (ID,True) else (ID1,True) for y in IN]
		#OUT[:] = [y if y != (ID,True) else (ID0,True) for y in OUT]

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


	#Helper function for addSegmentAldreadySplit
	#Before:	ID is the ID of segment
	#			self.contigs[ID] = [segment,[],[],0]
	#			The kmers from segment have been added to self.kmers
	#			The graph is a legal DBG except for one thing:
	#				segment has not been connected to the other contigs
	#				(after this function finishes we need to run mergeSegment
	#				to make sure this is a legal DBG)
	#			Note: No contig needs to be split up further in order to make the
	#				  graph a legal dBG
	#After:		segment has been connected to the other contigs in the graph using
	#			self.connect_a_to_b(aID,bID,A,B)
	#			i.e. self.contigs[ID] = [segment,IN,OUT,0]
	def connectSegment(self,ID,segment=""):
		#connect(Segment)
		#	ID=46, segment = CGAATC
		if segment == "":
			[segment,IN,OUT,COV] = self.contigs[ID]
		if self.printFunctionNames:
			print "connectSegment(ID="+str(ID)+", segment="+str(segment)+")"
		if self.printStatus:
			self.printContigs("inside connectSegment")
		if self.assertLegal:
			assert self.isLegalGraph()
		#self.printKmers("Kmers")
		#We need 4 for loops for this function.
		#fw(last)
		#bw(first)
		#and fw/bw of the twins of last and first

		#Add connections of the form:
		#	segment->N
		#	segment->twin(N)

		L = len(segment)
		for y in dbg.fw(segment[L-self.k:]):
			if y in self.kmers:
				#Find the information we need about the neighbouring contig N where we've seen y before
				[N_ID, N_i, N_B] = self.kmers[y]
				[N, N_IN, N_OUT, N_COV] = self.contigs[N_ID]
				N_L = len(N)

				#Stop if y is not the first kmer in N (or twin(N))
				if not N_i==0:
					if not N_ID==ID:
						print "about to throw an exception"
						print "segment:", segment
						print "segment[L-self.k:]:", segment[L-self.k:]
						print list(dbg.fw(segment[L-self.k:]))
						print "y [N_ID,N_i,N_B]:", y, [N_ID, N_i, N_B]
						print "[N, N_IN, N_OUT, N_COV]:", [N, N_IN, N_OUT, N_COV]
						raise Exception("We're not supposed to have to split any more contigs")
					continue

				#Now we know that either segment->N or segment->twin(N)
				self.connect_a_to_b_ifNotAlreadyConnected(ID,N_ID,True,N_B)

		#Add connections of the form:
		#	N->segment
		#	twin(N)->segment
		for y in dbg.bw(segment[0:self.k]):
			if y in self.kmers:
				#Find the information we need about the neighbouring contig N where we've seen y before
				[N_ID, N_i, N_B] = self.kmers[y]
				[N, N_IN, N_OUT, N_COV] = self.contigs[N_ID]
				N_L = len(N)

				#Stop if y is not the last kmer in N (or twin(N))
				if not N_i==N_L-self.k:
					continue

				self.connect_a_to_b_ifNotAlreadyConnected(N_ID,ID,N_B,True)

		if self.assertLegal:
			assert self.isLegalGraph()


	#Helper function for addSegmentAldreadySplit
	#Before:	ID is the ID of segment
	#			The graph is a legal DBG except we might need to merge segment
	#			with it's connections
	#After:		segment has been merged with it's neighbours wherever possible.
	#			The graph is a legal DBG
	#			Note: We need to make sure that the INs and OUTs of new_c have
	#			new_ID instead of ID/N_ID in their INs / OUTs
	#Note:
	#	It is possible we can merge with an IN connection as well
		#if we flipped the contig when merge-ing it will become an
		#OUT connection. Therefore we don't want to throw an exception
		#in this case (the error discovered 22.03.2017)
	def mergeSegment(self, ID):
		assert ID in self.contigs, "ID must be the ID of a contig in the graph"
		if self.printFunctionNames:
			print "mergeSegment(ID="+str(ID)+")"
		if self.printStatus:
			self.printContigs("inside mergeSegment")
		if self.assertLegal:
			assert self.isLegalGraph()

		count = 0
		maxCount = 1	#We can merge with OUT at most 1 time (2 times if we flip)
		[canMerge,a_ID,b_ID,A,B] = self.canMergeOUT(ID)
		while not ([canMerge,a_ID,b_ID,A,B]==[False,None,None,None,None]):
			count += 1
			if self.printStatus:
				print "About to merge with OUT for the "+str(count)+" time"
			ID, flipped = self.merge(a_ID,b_ID,A,B)
			if self.printStatus:
				self.printContigs("After merging with OUT")
			if (count==1) and flipped:
				maxCount = 2
			[canMerge,a_ID,b_ID,A,B] = self.canMergeOUT(ID)

		if count>maxCount:
			self.printContigs()
			print ID
			raise Exception("We can only merge OUT once if we dont flip and twice if we flip")

		count = 0
		maxCount = 1	#We can merge with IN at most 1 time
		[canMerge,a_ID,b_ID,A,B] = self.canMergeIN(ID)
		while not ([canMerge,a_ID,b_ID,A,B]==[False,None,None,None,None]):
			count += 1
			if self.printStatus:
				print "About to merge with IN for the "+str(count)+" time"
			ID, flipped = self.merge(a_ID,b_ID,A,B)
			if self.printStatus:
				self.printContigs("After merging with IN")
			[canMerge,a_ID,b_ID,A,B] = self.canMergeIN(ID)

		if count>maxCount:
			self.printContigs()
			print ID
			raise Exception("We can only merge IN once if we dont flip and twice if we flip")

		if self.assertLegal:
			if not self.isLegalDBG():
				self.printContigs("Not a legal DBG")
				print "The ID of the contig failing to merge:",ID
				raise Exception("The graph is supposed to be a legal DBG")



	def canMergeOUT(self,ID):
		assert ID in self.contigs, "ID must be the ID of a contig in the graph"
		if self.printFunctionNames:
			pass
			#print "canMergeOUT(ID="+str(ID)+")"
		if self.printStatus:
			pass
			#self.printContigs("inside canMergeOUT")
			#self.printKmers("Kmers")

		#Get the info about segment:
		[segment, IN, OUT, COV] = self.contigs[ID]
		L = len(segment)
		#print ID, segment, IN, OUT, COV

		#Try to merge segment with its OUT:
		if len(OUT)==1:
			if self.printStatus:
				print "segment has exactly one OUT"
			#segment has exactly one OUT: N. Find all the info we have about N:
			[N_ID,N_B] = OUT[0]
			[N, N_IN, N_OUT, N_COV] = self.contigs[N_ID]

			#check whether we can merge segment and N into 1 contig
			if N_B and len(N_IN)==1:
				if self.printStatus:
					print "We can merge segment with it's one OUT"
				if ID==N_ID:
					if self.printStatus:
						print "Nema nei! Segment er að tengjast í sjálfan sig svo við viljum ekki merge-a"
				else:
					#segment->N
					#N has exactly one IN: segment
					#Now we can merge segment and N into 1 contig new_c
					return [True,ID,N_ID,True,True]

			#check whether we can merge segment and twin(N) into 1 contig
			elif not N_B and len(N_OUT)==1:
				if self.printStatus:
					print "We can merge segment with the twin of it's one OUT"
				if ID==N_ID:
					if self.printStatus:
						print "Nema nei! Segment er að tengjast í sjálfan sig svo við viljum ekki merge-a"
				else:
					#segment->twin(N)    or N->twin(segment)
					#N has exactly one OUT: twin(segment)
					#Now we can merge segment and twin(N) into 1 contig new_c. Better: merge N and twin(segment)
					return [True,ID,N_ID,True,False]

		return [False,None,None,None,None]


	def canMergeIN(self, ID):
		assert ID in self.contigs, "ID must be the ID of a contig in the graph"
		if self.printFunctionNames:
			pass
			#print "canMergeIN(ID="+str(ID)+")"
		if self.printStatus:
			pass
			#self.printContigs("inside canMergeIN")
			#self.printKmers("Kmers")

		#Get the info about segment (might have changed if we merged with some OUT):
		[segment, IN, OUT, COV] = self.contigs[ID]
		L = len(segment)

		#Try to merge segment with its IN:
		if len(IN)==1:
			if self.printStatus:
				print "segment has exactly one IN"
			#segment has exactly one IN: N. Find all the info we have about N:
			[N_ID,N_B] = IN[0]
			[N, N_IN, N_OUT, N_COV] = self.contigs[N_ID]

			#check whether we can merge N and segment into 1 contig
			if N_B and len(N_OUT)==1:
				if self.printStatus:
					print "We can merge segment with it's one IN"
				if ID==N_ID:
					if self.printStatus:
						print "Nema nei! Segment er að tengjast í sjálfan sig svo við viljum ekki merge-a"
				else:
					#N->segment
					#N has exactly one OUT: segment
					#Now we can merge N and segment into 1 contig new_c
					return [True,N_ID,ID,True,True]

			elif not N_B and len(N_IN)==1:
				if self.printStatus:
					print "We can merge segment with the twin of it's one IN"
				if ID==N_ID:
					if self.printStatus:
						print "Nema nei! Segment er að tengjast í sjálfan sig svo við viljum ekki merge-a"
				else:
					#twin(N)->segment or twin(segment)->N
					#N has exactly one IN: segment
					#now we can merge twin(segment) and N into 1 contig new_c
					#new_ID = self.getID()
					return [True,N_ID,ID,False,True]

		return [False,None,None,None,None]


	#A=True, B=True:   a->b
	#A=True, B=False:  a->twin(b)
	#A=False, B=True:  twin(a)->b
	#A=False, B=False: twin(a)->twin(b) = b->a
	def merge(self,a_ID,b_ID,A,B,flipped=False):
		if self.printFunctionNames:
			print "merge(a_ID="+str(a_ID)+", b_ID="+str(b_ID)+", A="+str(A)+", B="+str(B)+")"
		assert a_ID in self.contigs, "a_ID must be the ID of a contig in the graph"
		assert b_ID in self.contigs, "b_ID must be the ID of a contig in the graph"
		assert isinstance(A,bool)
		assert isinstance(B,bool)
		if A and B:
			#merging a->b into new_c (by deleting b, adding it to a and updating a and it's connections)
			#self.changeID_OUTs(b_ID,a_ID,conns)		#make the OUTs of b instead connect to a
		#def changeID_ofConnections(self,ID):
			self.changeID_ofConnections(b_ID,a_ID)
			[a, a_IN, a_OUT, a_COV] = self.contigs[a_ID]
			[b, b_IN, b_OUT, b_COV] = self.contigs[b_ID]
			new_c = a + b[self.k-1:]
			self.contigs[a_ID] = [new_c, a_IN, b_OUT, a_COV+b_COV]
			del self.contigs[b_ID]
			#add kmers from new_c to kmerDict
			self.addKmersFromContig(a_ID, new_c)
			return a_ID, flipped
		elif A and (not B):
			#merging a->twin(b) into new_c
			self.flipContig(a_ID)
			return self.merge(b_ID,a_ID,True,True,True)
		elif (not A) and B:
			#merging twin(a)->b into new_c (by deleting a)
			self.flipContig(a_ID)
			return self.merge(a_ID,b_ID,True,True,True)
		elif (not A) and (not B):
			return self.merge(b_ID,a_ID,True,True)


"""
	#This is not used anywhere anymore. Could be useful later though
	#-----Helper function for splitContig-----
	#Before:	ID is the ID of a contig c
	#			The graph is NOT necessary a legal DBG
	#After:
	#		c has been removed from the graph
	#		c has been removed from self.contigs
	#		c's kmers have been removed from self.kmers
	#		the connection to c has been removed from the IN's and OUT's of c
	#Note:
	#	It is NOT guaranteed that the graph is a legal DBG
	def deleteContig(self,ID):
		if self.printFunctionNames:
			print "deleteContig(ID="+str(ID)+")"
		if self.printStatus:
			self.printContigs("inside deleteContig")
		#self.printKmers("Kmers")
		assert ID in self.contigs
		#Get notað setOUT og setIN á samtals 4 stöðum hér til að einfalda kóðann!!!!!!!!!!!!!!!!!!!!!!!1
		[c, IN, OUT, coverage] = self.contigs[ID]
		if self.printStatus:
			print c,IN,OUT,coverage

		#for every contig x connecting to c
		#	delete c from x's list of connections
		for [xID,xB] in IN:
			if self.printStatus:
				print xID,xB
			#x is a contig connecting to c
			[x, xIN, xOUT, xCov] = self.contigs[xID]
			#print x, xIN, xOUT, xCov
			if xB:
				#x->c so we remove ((ID,True)) from x.OUT
				xOUT.remove((ID,True))
				temp = self.contigs[xID]
				self.contigs[xID] = [temp[0],temp[1],xOUT,temp[3]]
			else:
				#twin(x)->c so we remove ((ID,False)) from x.IN
				xIN.remove((ID,False))
				temp = self.contigs[xID]
				self.contigs[xID] = [temp[0],xIN,temp[2],temp[3]]

		#for every contig x c connects to
		#	delete c from x's list of connections
		for [xID,xB] in OUT:
			#x is a contig c connects to
			[x, xIN, xOUT, xCov] = self.contigs[xID]
			if xB:
				xIN.remove((ID,True))
				temp = self.contigs[xID]
				self.contigs[xID] = [temp[0],xIN,temp[2],temp[3]]
			else:
				xOUT.remove((ID,False))
				temp = self.contigs[xID]
				self.contigs[xID] = [temp[0],temp[1],xOUT,temp[3]]

		#delete all kmers of c from self.kmers
		self.deleteKmersFromContig(ID,c)

		#delete c (and therefore we don't need to worry about deleting every individual IN/OUT of c)
		del self.contigs[ID]
"""


#------------------------Helper functions------------------------
#Before:	s is a string representing a segment
#				0    i    L
#			s: [    |    ]
#After:
#	s has been split up into
#	  0	    i               L
#	s[0:i-1]  ->  s[i-k:L-1]
def splitString(s,i,k,ps=False):
	if ps:
		print "splitString(s="+str(s)+", i="+str(i)+", k="+str(k)+")"
	L = len(s)
	#print "\ns: " + str(s) + ". i: " + str(i) + ". k: " + str(k) + ". L: " + str(L)
	assert(i>=k)	#i=k is the lowest i allowed
	#assert(i<=L-k)
	assert(i<=L-1)	#L-1 is the highest i allowed
	s0 = s[0:i]
	s1 = s[i-k+1:]
	if ps:
		print s0,s1
	return s0,s1

def reverseList(L):
	for i, (ID,B) in enumerate(L):
		L[i] = (ID,not B)


if __name__ == "__main__":
	G = Graph(3)
	#b: b_ID=0
	G.contigs[G.getID()] = ["AGT",[(99,True)],[(3,True),(4,True)],0]
	G.contigs[G.getID()] = ["TCA",[],[(99,True)],0] #1
	G.contigs[G.getID()] = ["ACA",[],[(99,True)],0] #2
	G.contigs[G.getID()] = ["GTT",[(0,True)],[],0]  #3
	G.contigs[G.getID()] = ["GTA",[(0,True)],[],0]  #4
	G.contigs[G.getID()] = ["TGG",[(99,False)],[],0] #5

	"""
	AGT
	ACT
	TCA
	TGA
	ACA
	TGT
	GTT x
	AAC x
	GTA
	TAC
	TGG
	CCA
	CAG
	CTG
	"""



	#TGG
	#CCA
	#a: a_ID=99
	G.setID(99)
	G.contigs[G.getID()] = ["CAG",[(1,True),(2,True),(5,False)],[(0,True)],0]
	G.addKmersFromAllContigs()
	print "Before merge:"
	print G.isLegal()
	G.printContigs("G")
	#G.mergeSegment(99)
	#print "After merge:"
	#G.printContigs("G")
	#self.assertTrue(G.isLegal())
	#self.assertTrue(G_correct.isLegal())
	#self.assertTrue(isSameGraph(G,G_correct))

	#NOTE: I need to write a test to actually check whether this result is correct. Felt creating a
	#G_correct was too much work. Checked this manually so this test is quite informal
