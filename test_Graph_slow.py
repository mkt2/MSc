#coding:utf8
import unittest
import Graph
from compareGraphs import *
import BF_counter
import helpers, Bloom

class Test_readGraphFromFile(unittest.TestCase):
    def test_1(self):
        fn = ["Input/read_1.fq"]
        k = 31
        BF = Bloom.Bloom(p=0.01,n=100000)
        G1 = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
        BF_counter.BF_counter_naive(fn,BF,k,G1,printProgress=False,pfn=True)
        print "Number of contigs:",len(G1)
        print BF
        BF.print_bitarray()
        self.assertTrue(BF.hasAcceptableRatio())
        outFolder = "defaultOutFolder"
        helpers.printResultsToFile(BF,G1,outFolder,timeInSeconds=-1,nextLine=-1)
        G2 = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
        G2.createGraphFromFile(outFolder+"/G.txt")
        self.assertTrue(G2.isLegalDBG())
        self.assertTrue(isSameGraph(G1,G2))

    def test_2(self):
        fn = ["Input/read_1.fq"]
        k = 31
        BF = Bloom.Bloom(p=0.01,n=100000)
        G1 = Graph.Graph(k,pfn=False,ps=False,al=True,pil=False,printInit=True)
        BF_counter.BF_counter_naive(fn,BF,k,G1,printProgress=False,pfn=True)
        print "Number of contigs:",len(G1)
        print BF
        BF.print_bitarray()
        self.assertTrue(BF.hasAcceptableRatio())
        outFolder = "defaultOutFolder"
        helpers.printResultsToFile(BF,G1,outFolder,timeInSeconds=-1,nextLine=-1)
        G2 = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
        G2.createGraphFromFile(outFolder+"/G.txt")
        self.assertTrue(G2.isLegalDBG())
        self.assertTrue(isSameGraph(G1,G2))

    def test_3(self):
        fn = ["Input/read_1.fq"]
        k = 31
        BF = Bloom.Bloom(p=0.01,n=100000)
        G1 = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
        BF_counter.BF_counter(fn,BF,k,G1,printProgress=False,pfn=True)
        print "Number of contigs:",len(G1)
        print BF
        BF.print_bitarray()
        self.assertTrue(BF.hasAcceptableRatio())
        outFolder = "defaultOutFolder"
        helpers.printResultsToFile(BF,G1,outFolder,timeInSeconds=-1,nextLine=-1)
        G2 = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
        G2.createGraphFromFile(outFolder+"/G.txt")
        self.assertTrue(G2.isLegalDBG())
        self.assertTrue(isSameGraph(G1,G2))

if __name__ == '__main__':
    unittest.main()
