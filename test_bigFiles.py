#coding:utf8
#import sys
#sys.path.insert(0, '/home/markus/Project-hack')
import unittest
import Graph
from compareGraphs import *
from BF_counter import *
import Bloom
from pybloom import BloomFilter

#--------------------------------------------------------------------------
#------------------------Tests for addSegmentToGraph-----------------------
#--------------------------------------------------------------------------
#Since we don't call generateSegmentsWithOnlySeenKmers we can use BF_counter and 
#BF_counter_naive directly. Reasons:
#   BF_counter:        Adds one segment at a time to the Graph using addSegmentToGraph
#   BF_counter_naive:  Creates a dict of all kmers from all segments and creates a Graph 
#                      from the results using dbg.py

#Test addSegmentToGraph by comparing the results of BF_counter with the results from
#BF_counter_naive (aka dbg.py)
class Test_BF_counter(unittest.TestCase):
    def test_1(self):
        print "test_bigFiles.py Test_BF_counter.test_1"
        k = 5
        fn = ["Input/read1.fq"]
        BF = Bloom.Bloom(p=0.01,n=1000)
        G = Graph.Graph(k,pfn=False,ps=False,al=True)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=True)
        BF_counter_and_naive(fn,BF,k,G,G_naive,printProgress=False,pfn=False)
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_2(self):
        print "test_bigFiles.py Test_BF_counter.test_2"
        fn = ["Input/read4.fq"]
        k = 5
        BF = Bloom.Bloom(p=0.01,n=1000)
        G = Graph.Graph(k,pfn=False,ps=False,al=True)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=True)
        BF_counter_and_naive(fn,BF,k,G,G_naive,printProgress=False,pfn=False)
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_3(self):
        print "test_bigFiles.py Test_BF_counter.test_3"
        fn = ["Input/read4.fq"]
        k = 5
        BF = Bloom.Bloom(p=0.01,n=1000)
        G = Graph.Graph(k,pfn=False,ps=False,al=True)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=True)
        BF_counter_and_naive(fn,BF,k,G,G_naive,printProgress=False,pfn=False)
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_4(self):
        print "test_bigFiles.py Test_BF_counter.test_4"
        fn = ["Input/read4.fq"]
        k = 5
        BF = Bloom.Bloom(p=0.01,n=1000)
        G = Graph.Graph(k,pfn=False,ps=False,al=True)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=True)
        BF_counter_and_naive(fn,BF,k,G,G_naive,printProgress=False,pfn=False)
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_5(self):
        print "test_bigFiles.py Test_BF_counter.test_5"
        k = 5
        G = Graph.Graph(k,pfn=False,ps=False)
        G.addSegmentToGraph("TGAAA")
        G.addSegmentToGraph("GCTGA")
        G.addSegmentToGraph("CTGAA")
        G_naive = G.createNaive()
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))

    def test_6(self):
        print "test_bigFiles.py Test_BF_counter.test_6"
        k = 5
        G = Graph.Graph(k,pfn=False,ps=False)
        G.addSegmentToGraph("GCTGAAA")
        G_naive = G.createNaive()
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))

    def test_7(self):
        print "test_bigFiles.py Test_BF_counter.test_7"
        k = 5
        G = Graph.Graph(k,pfn=False,ps=False)
        G.addSegmentToGraph("GCTGAAA")
        G.addSegmentToGraph('TCAGATTG')
        G_naive = G.createNaive()
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))
        
    def test_8(self):
        print "test_bigFiles.py Test_BF_counter.test_8"
        k = 5
        G = Graph.Graph(k,pfn=False,ps=False)
        G.contigs[G.getID()] = ["GCTGAAA",[],[],3]
        G.addKmersFromAllContigs()
        G.addSegmentAlreadySplit('TCAGATTG')
        
    def test_9(self):
        print "test_bigFiles.py Test_BF_counter.test_9"
        G = Graph.Graph(5,pfn=False,ps=False)
        G.addSegmentToGraph("GCTGA")
        G.addSegmentToGraph("TTCAGA")
        G.addSegmentToGraph("TTTCAT")
        self.assertTrue(graphEqualsNaive(G,-1,False,False))
        
    def test_10(self):
        print "test_bigFiles.py Test_BF_counter.test_10"
        G = Graph.Graph(5,pfn=False,ps=False,al=True,pil=False)
        G.addSegmentToGraph("ATCTTAGATT")
        G_naive = G.createNaive()
        #G.printContigs("G")
        #G_naive.printContigs("G_naive")
        self.assertTrue(graphEqualsNaive(G,G_naive=-1,pfn=False,ps=False))
        #Sýnist þetta vera villa í firstIndexOfSelfSplit
        
    def test_11(self):
        print "test_bigFiles.py Test_BF_counter.test_11"
        fn = ["Input/read8.fq"]
        k = 5
        BF = Bloom.Bloom(p=0.01,n=1000)     #Ratio: 0.30
        G = Graph.Graph(k,pfn=False,ps=False,al=True)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=True)
        BF_counter_and_naive(fn,BF,k,G,G_naive,printProgress=False,pfn=False)
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))
        self.assertTrue(BF.hasAcceptableRatio())
        
    def test_12(self):
        print "test_bigFiles.py Test_BF_counter.test_12"
        fn = ["Input/read72lines.fq"]
        k = 5
        BF = Bloom.Bloom(p=0.001,n=1300)    #Ratio: 0.46
        G = Graph.Graph(k,pfn=False,ps=False,al=True)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=True)
        BF_counter_and_naive(fn,BF,k,G,G_naive,printProgress=False,pfn=False)
        #G.printContigs("G")
        #G_naive.printContigs("G_naive")
        #print len(G)
        #print len(G_naive)
        #print BF
        #BF.print_bitarray()
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))
        self.assertTrue(BF.hasAcceptableRatio())
        
    def test_13(self):
        print "test_bigFiles.py Test_BF_counter.test_13"
        fn = ["Input/read240lines.fq"]
        k = 5
        BF = Bloom.Bloom(p=0.01,n=2800)     #Ratio: 0.28
        G = Graph.Graph(k,pfn=False,ps=False,al=True)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=True)
        BF_counter_and_naive(fn,BF,k,G,G_naive,printProgress=False,pfn=False)
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=False))
        self.assertTrue(BF.hasAcceptableRatio())
        
    def test_14(self):
        print "test_bigFiles.py Test_BF_counter.test_14"
        fn = ["Input/read_1.fq"]
        k = 5
        BF = Bloom.Bloom(p=0.01,n=100000)   #Keyrir án villna á 352 sekúndum (sleppti að búa til naive)
                                            #Ratio: 0.01. m=958506. 1226 segment sem er bætt við Grafið
                                            #Grafið inniheldur 512 contig-a
        G = Graph.Graph(k,pfn=False,ps=False,al=True)
        BF_counter(fn,BF,k,G,printProgress=True,pfn=True,startAtLine=0)
        print "Number of contigs:",len(G)
        print BF
        BF.print_bitarray()
        self.assertTrue(BF.hasAcceptableRatio())

    def test_15(self):
        print "test_bigFiles.py Test_BF_counter.test_15"
        fn = ["Input/read_1.fq"]
        k = 31
        BF = Bloom.Bloom(p=0.01,n=100000)   #Keyrir án villna á 1776 sekúndum (sleppti að búa til naive)
                                            #Ratio: 0.35. m=958506. 1578 segment sem er bætt við Grafið
                                            #Grafið inniheldur 762 contig-a
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=True)
        BF_counter_naive(fn,BF,k,G_naive,printProgress=False,pfn=False)
        print "Number of contigs:",len(G_naive)
        print BF
        BF.print_bitarray()
        self.assertTrue(BF.hasAcceptableRatio())

class Test_BloomFilters(unittest.TestCase):
    #Prófa hvort pybloom hegði sér eðlilega með því að búa til BF
    #úr öllum segmentum í r1.fq og r2.fq
    def test_1(self):
        fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
        k = 31
        numberOfKmers = 6000000
        BF = BloomFilter(capacity=numberOfKmers, error_rate=0.01)
        for segment,lineNr in generateSegmentsFrom_fq(fn,k,BF,pfn=True,printProgress=True):
            pass
        print BF
        BF.print_bitarray()
        self.assertTrue(BF.hasAcceptableRatio())
        #55 sek. 0.20 ratio

    #Prófa hvort Bloom hegði sér eðlilega með því að búa til BF
    #úr öllum segmentum í r1.fq og r2.fq
    def test_2(self):
        fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
        k = 31
        numberOfKmers = 6000000
        BF = Bloom.Bloom(0.01,numberOfKmers,pfn=True)
        for segment,lineNr in generateSegmentsFrom_fq(fn,k,BF,pfn=True,printProgress=True):
            pass
        print BF
        BF.print_bitarray()
        self.assertTrue(BF.hasAcceptableRatio())
        #83 sek. 0.4 ratio

    #Athuga hvort generateSegmentsFrom_fq búi til eðlilegan fjölda af kmerum
    #Bera saman við töluna frá Páli
    def test_3(self):
        fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
        k = 31
        numberOfKmers = 6000000
        BF = Bloom.Bloom(0.01,numberOfKmers,pfn=True)
        kmerDict = createKmerDict(fn,BF,k,pfn=True,printProgress=True)
        print len(kmerDict)
        self.assertTrue(BF.hasAcceptableRatio())
        #Vantar að vita hvaða tala er eðlileg. Get t.d. gert fall til að athuga hvort
        #þetta sé á eðlilegu bili
        self.assertTrue(len(kmerDict)==10)
        #83 sek


        
"""
    def test_16(self):
        fn = ["Input/read8.fq"]
        k = 5
        segments = generateSegmentsFrom_fq(fn,k)
        numberOfKmers = countKmers(segments,k)
        BF = Bloom.Bloom(0.01,numberOfKmers,k)
        segments = list(generateSegmentsWithOnlySeenKmers(k,segments,BF))
        segments = list( segments[i] for i in [37])
        print segments

        G = Graph.Graph(k,pfn=True,ps=False,al=True,pil=False)
        L = len(segments)
        for i, s in enumerate(segments):
            print "Adding segment", str(i)+" out of", str(L)
            G.addSegmentToGraph(s)
            #assert graphEqualsNaive(G,G_naive=-1,pfn=False,ps=False)
            self.assertTrue(graphEqualsNaive(G,G_naive=-1,pfn=False,ps=False))
    
    def test_17(self):
        fn = ["Input/read8.fq"]
        k = 5
        segments = generateSegmentsFrom_fq(fn,k)
        numberOfKmers = countKmers(segments,k)
        BF = Bloom.Bloom(0.01,numberOfKmers,k)
        segments = list(generateSegmentsWithOnlySeenKmers(k,segments,BF))
        #segments = list( segments[i] for i in [13,36])
        segments = ["CTTTCAGATTG","CGCCGCCATGCCGACCATCCCTTTCATCCCCGTACCAGACACGCTGACCATTGCCATGTTATTCAGATTG"]
        print segments

        G = Graph.Graph(k,pfn=False,ps=False,al=True,pil=False)
        L = len(segments)
        for i, s in enumerate(segments):
            print "Adding segment", str(i)+" out of", str(L)
            G.addSegmentToGraph(s)
            #G.printContigs()
            #assert graphEqualsNaive(G,G_naive=-1,pfn=False,ps=False)
            self.assertTrue(graphEqualsNaive(G,G_naive=-1,pfn=False,ps=False))

    def test_18(self):
        fn = ["Input/read8.fq"]
        k = 5
        segments = generateSegmentsFrom_fq(fn,k)
        numberOfKmers = countKmers(segments,k)
        BF = Bloom.Bloom(0.01,numberOfKmers,k)
        segments = list(generateSegmentsWithOnlySeenKmers(k,segments,BF))
        segments = list( segments[i] for i in [13,35,36])
        print segments

        G = Graph.Graph(k,pfn=False,ps=False,al=True,pil=False)
        L = len(segments)
        for i, s in enumerate(segments):
            print "Adding segment", str(i)+" out of", str(L)
            G.addSegmentToGraph(s)
            #G.printContigs()
            #assert graphEqualsNaive(G,G_naive=-1,pfn=False,ps=False)
            self.assertTrue(graphEqualsNaive(G,G_naive=-1,pfn=False,ps=False))
"""


if __name__ == '__main__':
    unittest.main()
    x=3