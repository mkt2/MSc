#coding:utf8
import unittest
import Graph
from Bloom import *
import random
import itertools

#--------------------------------------------------------------------------
#-----------------Simple functions for working with graphs-----------------
#--------------------------------------------------------------------------
class Test_Bloom(unittest.TestCase):
    #test add and contains
    def test_1(self):
        print "Test_Bloom.test_1"
        BF = Bloom(0.01,10)
        BF.add("ATC")
        self.assertTrue("ATC" in BF)
        self.assertFalse("AAA" in BF)
        ratio,cZero,cOne = BF.computeRatio()
        self.assertTrue(BF.hasAcceptableRatio(ratio))

    #check the false positive rate
    def test_2(self):
        print "\nTest_Bloom.test_2"
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

    #try a false positive rate of 0
    def test_3(self):
        print "\nTest_Bloom.test_3"
        #n = 1000000, kmerLength=13     #17 sek
        #n = 10000000, kmerLength=13    #26 sek
        n = 1000    #the number of kmers to be added to the bloom filter
        kmerLength = 13
        p = 0
        BF = Bloom(p,n,pfn=False)

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
        self.assertTrue(FPR<=p)


    def test_4(self):
        print "\nTest_Bloom.test_4"
        BF = Bloom(0.01,10)
        BF.add("ATC")
        ratio,cZero,cOne = BF.computeRatio()
        self.assertTrue(BF.hasAcceptableRatio(ratio))

    def test_5(self):
        print "\nTest_Bloom.test_5"
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

    def test_6(self):
        print "\nTest_Bloom.test_6"
        BF = Bloom(0.01,2,pfn=False)
        BF.add("ATCAAAAAACCCTCTCGGGTGGGGAAATCAT")
        BF.add("GATTATAGACACACGATATGATGGTTCCCCC")
        ratio,cZero,cOne = BF.computeRatio()
        self.assertTrue(BF.hasAcceptableRatio(ratio))

    def test_7(self):
        print "\nTest_Bloom.test_7"
        BF = Bloom(0.01,2,pfn=False)
        BF.add("ATCAAAAAACCCTCTCGGGTGGGGAAATCAT")
        BF.add("GATTATAGACACACGATATGATGGTTCCCCC")
        self.assertTrue(BF.hasAcceptableRatio())

    def test_8(self):
        print "\nTest_Bloom.test_8"
        BF = Bloom(0.01,2,pfn=False)
        BF.add("ATCAAAAAACCCTCTCGGGTGGGGAAATCAT")
        BF.add("GATTATAGACACACGATATGATGGTTCCCCC")
        BF.add("ATATATCCGGCGATAGTAGCGCGTTTTAGGC")
        self.assertFalse(BF.hasAcceptableRatio())
        

def generate_n_kmers_of_length_k(n,k):
    return [''.join(i) for i in itertools.product("ACGT",repeat=k)][0:n]


if __name__ == '__main__':
    unittest.main()