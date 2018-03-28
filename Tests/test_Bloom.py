#coding:utf8
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import unittest
import Graph
from Bloom import *
import random
import itertools
import dbg
import collections

def computeFPR_simple(kd,BF,p,printFPR=False):
	fpos = 0
	L = len(kd)
	for i,km in enumerate(kd.keys()):
		if i<L/2:
			BF.add(km)
			assert km in BF
		else:
			if km in BF:
				fpos += 1
			else:
				assert not (km in BF)
	#FPR = (fjoldi rangra km in BF) / Fjöldi km í heild
	FPR = float(fpos)/L
	if printFPR:
		printResults(BF,fpos,FPR,p)
	return FPR

def computeFPR_complex(kd,BF,p,printFPR=False):
	fpos = 0
	L = len(kd)
	for i,km in enumerate(kd.keys()):
		if i<L/2:
			BF.kmerInBF_also_AddIf_B_is_True(km,True)
			assert km in BF
		else:
			if BF.kmerInBF_also_AddIf_B_is_True(km,False):
				fpos += 1
			else:
				assert not (km in BF)
	#FPR = (fjoldi rangra km in BF) / Fjöldi km í heild
	FPR = float(fpos)/L
	if printFPR:
		printResults(BF,fpos,FPR,p)
	return FPR

def printResults(BF,fpos,FPR,p):
	print BF
	BF.print_bitarray()
	print "Number of False Positives:",fpos
	print "Actual false positive rate:", FPR
	print "Desired false positive proability p:", p

#--------------------------------------------------------------------------
#-----------------Simple functions for working with graphs-----------------
#--------------------------------------------------------------------------
class Test_Bloom(unittest.TestCase):
	#test add and contains
	def test_01(self):
		#print "\nTest_Bloom.test_01"
		BF = Bloom(0.01,10)
		BF.add("ATC")
		self.assertTrue("ATC" in BF)
		self.assertFalse("AAA" in BF)
		ratio,cZero,cOne = BF.computeRatio()
		self.assertTrue(BF.hasAcceptableRatio(ratio))

	#check the false positive rate
	def test_02(self):
		#print "\nTest_Bloom.test_02"
		n = 5000    #the number of kmers to be added to the bloom filter
		p = 0.1
		kmerLength=5
		BF = Bloom(p,n)

		#generate n kmers and add them to BF
		kmers = generate_n_kmers_of_length_k(n*2,k=kmerLength)
		for i in xrange(0,len(kmers)/2):
			BF.add(kmers[i])

		#check whether the false positive rate holds
		fpos = 0   #number of false positives

		for i in xrange(len(kmers)/2,len(kmers)):
			km = kmers[i]
			if km in BF:
				fpos += 1

		FPR = float(fpos)/n
		#print "Printing the bloom filter"
		#print BF
		#print "Number of False Positives:",fpos
		#print "Actual false positive rate:", FPR
		#print "Desired false positive proability p:", p
		self.assertTrue(FPR<=p)
		ratio,cZero,cOne = BF.computeRatio()
		self.assertTrue(BF.hasAcceptableRatio(ratio))

	def test_03(self):
		#print "\nTest_Bloom.test_03"
		BF = Bloom(0.01,10)
		BF.add("ATC")
		ratio,cZero,cOne = BF.computeRatio()
		self.assertTrue(BF.hasAcceptableRatio(ratio))

	def test_04(self):
		#print "\nTest_Bloom.test_04"
		#n = 16000; kmerLength = 9; p=0.1   #0.582 sek
		#n = 100000; kmerLength = 9; p=0.1 #3.8 sec
		#n = 1000000; kmerLength = 9; p=0.1 #4.7 sec
		#n = 10000000; kmerLength = 9; p=0.1  #5 sec
		n = 100000000; kmerLength = 9; p=0.1  #6.537 sec
		BF = Bloom(p,n)

		#generate n kmers and add them to BF
		kmers = generate_n_kmers_of_length_k(n*2,kmerLength)
		for i in xrange(0,len(kmers)/2):
			BF.add(kmers[i])

		#check whether the false positive rate holds
		fpos = 0   #number of false positives

		for i in xrange(len(kmers)/2,len(kmers)):
			km = kmers[i]
			if km in BF:
				fpos += 1

		FPR = float(fpos)/n
		#print "Printing the bloom filter"
		#print BF
		#print "Number of False Positives:",fpos
		#print "Actual false positive rate:", FPR
		#print "Desired false positive proability p:", p
		#print "Number of kmers of given length: 4^kmerLength =", pow(4,kmerLength)
		self.assertTrue(FPR<=p)
		ratio,cZero,cOne = BF.computeRatio()
		self.assertTrue(BF.hasAcceptableRatio(ratio))
		#print BF

	def test_05(self):
		#print "\nTest_Bloom.test_05"
		BF = Bloom(0.01,2,pfn=False)
		BF.add("ATCAAAAAACCCTCTCGGGTGGGGAAATCAT")
		BF.add("GATTATAGACACACGATATGATGGTTCCCCC")
		ratio,cZero,cOne = BF.computeRatio()
		self.assertTrue(BF.hasAcceptableRatio(ratio))

	def test_06(self):
		#print "\nTest_Bloom.test_06"
		BF = Bloom(0.01,2,pfn=False)
		BF.add("ATCAAAAAACCCTCTCGGGTGGGGAAATCAT")
		BF.add("GATTATAGACACACGATATGATGGTTCCCCC")
		self.assertTrue(BF.hasAcceptableRatio())

	def test_07(self):
		#print "\nTest_Bloom.test_07"
		BF = Bloom(0.01,2,pfn=False)
		BF.add("ATCAAAAAACCCTCTCGGGTGGGGAAATCAT")
		BF.add("GATTATAGACACACGATATGATGGTTCCCCC")
		BF.add("ATATATCCGGCGATAGTAGCGCGTTTTAGGC")
		self.assertFalse(BF.hasAcceptableRatio())

	def test_08(self):
		#print "\nTest_Bloom.test_08"
		alphabet = ["A","C","G","T","N"]
		k = 31
		fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
		numberOfKmers = 6000000
		p = 0.01
		BF = Bloom(p,numberOfKmers,False)
		kd = collections.defaultdict(int)	#stores all kmers in the reads
		for f in fn:
			h = open(f, "rU")
			for lineNr,line in enumerate(h,start=0):
				if (lineNr%4 == 1):
					segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
					for s in segments:
						assert isinstance(s, basestring)
						for i, kmer in enumerate(dbg.kmers(s,k)):
							rep_kmer = min(kmer,dbg.twin(kmer))
							kd[rep_kmer] += 1
		FPR = computeFPR_simple(kd,BF,p,printFPR=False)
		self.assertTrue(FPR<=p)
		self.assertTrue(BF.hasAcceptableRatio())

	def test_09(self):
		#print "\nTest_Bloom.test_09"
		alphabet = ["A","C","G","T","N"]
		k = 31
		fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
		numberOfKmers = 6000000
		p = 0.01
		BF = Bloom(p,numberOfKmers,False)
		kd = collections.defaultdict(int)	#stores all kmers in the reads
		for f in fn:
			h = open(f, "rU")
			for lineNr,line in enumerate(h,start=0):
				if (lineNr%4 == 1):
					segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
					for s in segments:
						assert isinstance(s, basestring)
						for i, kmer in enumerate(dbg.kmers(s,k)):
							rep_kmer = min(kmer,dbg.twin(kmer))
							kd[rep_kmer] += 1
		FPR = computeFPR_complex(kd,BF,p,printFPR=False)
		self.assertTrue(FPR<=p)
		self.assertTrue(BF.hasAcceptableRatio())

	def test_10(self):
		#print "\nTest_Bloom.test_10"
		alphabet = ["A","C","G","T","N"]
		k = 31
		fn = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
		numberOfKmers = 100000000
		p = 0.01
		BF = Bloom(p,numberOfKmers,False)
		kd = collections.defaultdict(int)	#stores all kmers in the reads
		for f in fn:
			h = open(f, "rU")
			for lineNr,line in enumerate(h,start=0):
				if (lineNr%4 == 1):
					segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
					for s in segments:
						assert isinstance(s, basestring)
						for i, kmer in enumerate(dbg.kmers(s,k)):
							rep_kmer = min(kmer,dbg.twin(kmer))
							kd[rep_kmer] += 1
		FPR = computeFPR_simple(kd,BF,p,printFPR=True)
		self.assertTrue(FPR<=p)
		self.assertTrue(BF.hasAcceptableRatio())

	def test_11(self):
		#print "\nTest_Bloom.test_11"
		alphabet = ["A","C","G","T","N"]
		k = 31
		fn = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
		numberOfKmers = 100000000
		p = 0.01
		BF = Bloom(p,numberOfKmers,False)
		kd = collections.defaultdict(int)	#stores all kmers in the reads
		for f in fn:
			h = open(f, "rU")
			for lineNr,line in enumerate(h,start=0):
				if (lineNr%4 == 1):
					segments = filter(lambda x: len(x)>=k and all([c in alphabet for c in x]),line.strip().split("N"))
					for s in segments:
						assert isinstance(s, basestring)
						for i, kmer in enumerate(dbg.kmers(s,k)):
							rep_kmer = min(kmer,dbg.twin(kmer))
							kd[rep_kmer] += 1
		FPR = computeFPR_complex(kd,BF,p,printFPR=True)
		self.assertTrue(FPR<=p)
		self.assertTrue(BF.hasAcceptableRatio())

def generate_n_kmers_of_length_k(n,k):
	return [''.join(i) for i in itertools.product("ACGT",repeat=k)][0:n]

if __name__ == '__main__':
	unittest.main()