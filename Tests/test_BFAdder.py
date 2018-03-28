#coding:utf8
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import unittest
import Graph
import BFAdder
import dbg
import collections



"""Man ekki hvernig ég bjó til þessar skrár sem testin nota þ.a. ég kommenta þetta út til að byrja með
class Test_BF_counter_naive(unittest.TestCase):
    def test_1(self):
        fn = ["Input/Test_BF_counter_naive/test_1.fq"]
        segments = ["AAAAA","TAG","TAG"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        G_correct.contigs[G_correct.getID()] = ["TAG",[(1,False),(1,False)],[],0] #twin: CTA
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G_naive,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_2(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_9.fq"]
        segments = ["AAA","AAA"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G_naive,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_3(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_8.fq"]
        segments = ["AAAAAA"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G_naive,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_4(self):
        fn = ["Input/Test_BF_counter_naive/test_4.fq"]
        segments = ["ATC","TACG","ATC","TACG"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ATC",[(0,False),(0,False)],[],0] #twin: GAT
        G_correct.contigs[G_correct.getID()] = ["TACG",[(1,False),(1,False)],[(1,False),(1,False)],0] #twin: CGTA
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G_naive,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_5(self):
        fn = ["Input/Test_BF_counter_naive/test_5.fq"]
        segments = ["ATC","TCC","ATC","TCC"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["GGAT",[],[(0,False),(0,False)],0] #twin: ATCC
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G_naive,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_6(self):
        fn = ["Input/Test_BF_counter_naive/test_6.fq"]
        segments = ["CAT","ATA","CAT","ATA"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["CAT",[],[(0,False),(0,False),(1,True)],0] #twin: ATG
        G_correct.contigs[G_correct.getID()] = ["ATA",[(0,True),(1,False),(1,False)],[(1,False),(1,False)],0] #twin: ATG
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G_naive,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_7(self):
        fn = ["Input/Test_BF_counter_naive/test_7.fq"]
        segments = ["ATCAT","ATCAT"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ATCAT",[(0,True),(0,False),(0,False)],[(0,True),(0,False),(0,False)],0] #twin: ATGAT
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G_naive,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_8(self):
        fn = ["Input/Test_BF_counter_naive/test_8.fq"]
        segments = ["ATCC","ATCC"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ATCC",[(0,False),(0,False)],[],0] #twin: GGAT
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G_naive,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_9(self):
        fn = ["Input/Test_BF_counter_naive/test_9.fq"]
        segments = ["AAT","GAT","AAT","GAT"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAT",[],[(0,False),(0,False),(1,False)],0] #twin: GGAT
        G_correct.contigs[G_correct.getID()] = ["GAT",[],[(0,False),(1,False),(1,False)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G_naive,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

Sjá print skipanir í main fallinu
class Test_BF_counter(unittest.TestCase):
    def test_1(self):
        fn = ["Input/Test_BF_counter/test_1.fq"]
        segments = ["AAA","TAG","AAA","TAG"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        #def BF_counter(fn,k,BF,G,IK,maxCov,pfn=False,printProgress=False,startAtLine=0,skipPictures=False)
        BF_counter(fn,k,BF,G,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        G_correct.contigs[G_correct.getID()] = ["TAG",[(1,False),(1,False)],[],0] #twin: CTA
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_2(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_9.fq"]
        segments = ["AAA","AAA"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter(fn,BF,k,G,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_3(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_8.fq"]
        segments = ["AAA","AAA"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter(fn,BF,k,G,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))
        self.assertTrue(BF.hasAcceptableRatio())

    def test_4(self):
        fn = ["Input/Test_BF_counter_naive/test_4.fq"]
        segments = ["ATC","TACG","ATC","TACG"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter(fn,BF,k,G,pfn=False,printProgress=False)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ATC",[(0,False),(0,False)],[],0] #twin: GAT
        G_correct.contigs[G_correct.getID()] = ["TACG",[(1,False),(1,False)],[(1,False),(1,False)],0] #twin: CGTA
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=True,pfn=False,ps=False,relaxAssertions=False))
        self.assertTrue(BF.hasAcceptableRatio())
"""


#First we test the naive functions in order:
class Test_naive_createDict(unittest.TestCase):
    def helper(self):
        #fn = ["testData/num3_len5_r1.fq","testData/num3_len5_r2.fq"]
        fn = ["Input/testData/num3_len5_r1.fq","Input/testData/num3_len5_r2.fq"]
        k = 5
        return fn,k

    def test_1(self):
        fn,k = self.helper()        
        kd = BFAdder.naive_createDict(fn,k)
        print kd

        kmerList = ["AAAAA","TTTTT","CCCCC","GGGGG","AAAAT","ATTTT"]
        for km in kmerList:
            self.assertTrue(km in kd),"All kmers in kmerList should be in the kmerdict"
        for km in kd:
            self.assertTrue(km in kmerList),"All kmers in the kmerdict should be in the kmerList"
        print kmerList

    def test_2(self):
        pass


if __name__ == '__main__':
    unittest.main()