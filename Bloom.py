#coding:utf8
from math import log
from math import pow
from math import ceil
from random import randint
import bitarray
import dbg

#To make this actually random I have to use random seeds

#p: The desired false positive probability
#k: number of hash functions	
#m: number of bits in the array
#n: number of elements inserted
#	I will let n be the combined length of segments we're going to add to the BF (úrelt)


class Bloom:
	def __init__(self, p, n ,pfn=False):
		if pfn:
			print "Initializing the Bloom Filter: p="+str(p)+", n="+str(n)+", pfn="+str(pfn)
		self.p = p
		assert not p==0, "Changed the BF so that it can't handle p=0 anymore. Use BF=Set() if you need p=0"		
		self.k = int(ceil(-log(p)/pow(log(2),2)))
		self.m = int(ceil(-n*log(p)/pow(log(2),2)))
		self.m = int(2**(ceil(log(self.m, 2))))
		self.bitarray = bitarray.bitarray(self.m)
		self.bitarray.setall(False)
		self.randNum = randint(10, 100)		#used for myHash
		self.count = 0	#the number of kmers added

	def __len__(self):
		"""Return the number of keys stored by this bloom filter."""
		#return self.count	<--- this is the old approximation
		return int( float(-self.m)/self.k * log(1.0-self.bitarray.count()/float(self.m)) )

	def __str__(self):
		return \
		"\nPrinting the Bloom BF (except ratio):" \
		+   "  p:          "+str(self.p) \
		+"\n  k:          "+str(self.k) \
		+"\n  m:          "+str(self.m) \
		+"\nAdded "+str(len(self))+" kmers\n"
		+"\n"


	def computeRatio(self):
		cZero = 0
		cOne = self.bitarray.count()
		cZero = len(self.bitarray) - cOne
		ratio = float(cOne)/(cZero+cOne)
		if ratio >= 0.5:
			assert cOne>=cZero
		else:
			assert cOne<cZero
		return ratio,cZero,cOne

	def hasAcceptableRatio(self,ratio=None):
		if not ratio:
			ratio,cZero,cOne= self.computeRatio()
		if ratio>0.5:
			return False
		if ratio <0:
			return False
		return True

	def print_bitarray(self,ratio=None):
		if not ratio:
			ratio,cZero,cOne= self.computeRatio()
		print "\nPrinting the 0/1 ratio in the BF bitarray:" \
        +"\n  Number of ones:  "+str(cOne) \
        +"\n  Number of zeros: "+str(cZero) \
        +"\n  Ratio:           "+str(ratio)+"\n"

	def bitarray_str(self,ratio=None):
		if not ratio:
			ratio,cZero,cOne= self.computeRatio()
		return "\nPrinting the 0/1 ratio in the BF bitarray:" \
        +"\n  Number of ones:  "+str(cOne) \
        +"\n  Number of zeros: "+str(cZero) \
        +"\n  Ratio:           "+str(ratio)+"\n"


	#def myHash(self,km,i):
	#	return (self.randNum + i*hash(km)) % self.m
		

	#Add kmer to the Bloom Filter
	def add(self,kmer):
		self.count += 1
		h = hash(kmer)
		for i in xrange(0,self.k):
			#generate a hash value for hash function i
			#value = (self.randNum + i*h) % self.m
			value = (self.randNum + i*h) & (self.m-1)
			
			#Set the value at index 'value' as 1
			self.bitarray[value] = 1
		#assert kmer in self
				
					
	#Checks whether the Bloom Filter contains kmer
	#Either "Maybe it contains the kmer" or
	#		"It does not contain the kmer"
	def __contains__(self,kmer):
		#if not isinstance(kmer, basestring):
		#	print "kmer: ",kmer
		#	raise Exception('Each kmer has to be a string')
		#kmerToInt = lambda x: int(sum([pow({'A': 0, 'C': 1, 'G': 2, 'T': 3}[B],2) for B in x]))
		h = hash(kmer)
		for i in xrange(0,self.k):
			#generate a hash value for hash function i
			#value = (self.randNum + i*h) % self.m
			value = (self.randNum + i*h) & (self.m-1)
			if not (self.bitarray[value] == 1):
				return False
		return True

	#Purpose:
	#	Do the same as:
	#		if not kmer in BF and B:
	#			BF.add(kmer)
	#	...just faster...
	#Adds the kmer to the BF and returns whether or not it was there already
	#This function is supposed to run 
	def kmerInBF_also_AddIf_B_is_True(self,kmer,B,alsoAddTwin=True):
		alreadyInBF = True
		h = hash(kmer)
		for i in xrange(0,self.k):
			#generate a hash value for hash function i
			#value = (self.randNum + i*h) % self.m
			value = (self.randNum + i*h) & (self.m-1)
			if not (self.bitarray[value] == 1):
				alreadyInBF = False
				#Only add the kmer to BF if B==True
				if B:
					self.bitarray[value] = 1
		#assert (B==False) or (kmer in self), "B="+str(B)+", kmer in self="+str(kmer in self)
		if not alreadyInBF and alsoAddTwin:
			self.add(dbg.twin(kmer))
		return alreadyInBF
		

if __name__ == "__main__":
	BF = Bloom(0.01,100,True)
	BF.add("AAA")
	print len(BF)
	BF.add("AAA")
	print len(BF)
	BF.add("ATCGA")
	print len(BF)
	BF.add("TTT")
	print len(BF)
	

