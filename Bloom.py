#coding:utf8
from math import log
from math import pow
from math import ceil
from sets import Set
from random import randint
import bitarray

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
		if p==0:
			self.theSet = Set()
		else:		
			self.k = int(ceil(-log(p)/pow(log(2),2)))
			self.m = int(ceil(-n*log(p)/pow(log(2),2)))
			self.bitarray = bitarray.bitarray(self.m)
			self.bitarray.setall(False)
			self.randNum = randint(10, 100)		#used for myHash
		self.count = 0	#the number of kmers added

	def __len__(self):
		"""Return the number of keys stored by this bloom filter."""
		return self.count

	def __str__(self):
		if self.p==0:
			return "\np: "+str(self.p)+"\ntheSet: "+str(self.theSet)+"\n"
		else:
			#print "\nPrinting the Bloom BF (except ratio):"
		
			return \
			"\nPrinting the Bloom BF (except ratio):" \
			+   "  p:          "+str(self.p) \
			+"\n  k:          "+str(self.k) \
			+"\n  m:          "+str(self.m) \
			+"\nAdded "+str(self.count)+" kmers\n"
			+"\n"


	def computeRatio(self):
		if self.p==0:
			raise Exception("This function is not defined for p=0")
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
		if self.p==0:
			raise Exception("This function is not defined for p=0")
		if not ratio:
			ratio,cZero,cOne= self.computeRatio()
		if ratio>0.5:
			return False
		if ratio <0:
			return False
		return True

	def print_bitarray(self,ratio=None):
		if self.p==0:
			raise Exception("This function is not defined for p=0")
		if not ratio:
			ratio,cZero,cOne= self.computeRatio()
		print "\nPrinting the 0/1 ratio in the BF bitarray:" \
        +"\n  Number of ones:  "+str(cOne) \
        +"\n  Number of zeros: "+str(cZero) \
        +"\n  Ratio:           "+str(ratio)+"\n"

	def bitarray_str(self,ratio=None):
		if self.p==0:
			raise Exception("This function is not defined for p=0")
		if not ratio:
			ratio,cZero,cOne= self.computeRatio()
		return "\nPrinting the 0/1 ratio in the BF bitarray:" \
        +"\n  Number of ones:  "+str(cOne) \
        +"\n  Number of zeros: "+str(cZero) \
        +"\n  Ratio:           "+str(ratio)+"\n"


	def myHash(self,km,i):
		return (self.randNum + i*hash(km)) % self.m
		

	#Add kmer to the Bloom Filter
	def add(self,kmer):
		if not isinstance(kmer, basestring):
			raise Exception('Each kmer has to be a string')
		#print "add(kmer="+str(kmer)+")"
		if self.p==0:
			self.theSet.add(kmer)
		else:
			if not (kmer in self):
				#kmerToInt = lambda x: int(sum([pow({'A': 0, 'C': 1, 'G': 2, 'T': 3}[B],2) for B in x]))
				for i in xrange(0,self.k):
					#generate a hash value for hash function i
					value = int(self.myHash(kmer,i))
					
					#Set the value at index 'value' as 1
					self.bitarray[value] = 1
				assert kmer in self
				self.count += 1
					
	#Checks whether the Bloom Filter contains kmer
	#Either "Maybe it contains the kmer" or
	#		"It does not contain the kmer"
	def __contains__(self,kmer):
		if not isinstance(kmer, basestring):
			print "kmer: ",kmer
			raise Exception('Each kmer has to be a string')
		if self.p==0:
			return kmer in self.theSet
		else:
			#kmerToInt = lambda x: int(sum([pow({'A': 0, 'C': 1, 'G': 2, 'T': 3}[B],2) for B in x]))
			for i in xrange(0,self.k):
				#generate a hash value for hash function i
				value = int(self.myHash(kmer,i))
				if not (self.bitarray[value] == 1):
					return False
			return True


if __name__ == "__main__":
	pass
	

