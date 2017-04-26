#coding:utf8
import unittest
import Graph
from compareGraphs import *
from BF_counter import *
import random
import itertools
import dbg

def createCorrectDict(segments,k):
    kmerDict_correct = collections.defaultdict(int)
    for segment in segments:
        for km in dbg.kmers(segment,k):
            if km not in kmerDict_correct:
                kmerDict_correct[km] = 0
            else:
                kmerDict_correct[km] += 1
        for km in dbg.kmers(dbg.twin(segment),k):
            if km not in kmerDict_correct:
                kmerDict_correct[km] = 0
            else:
                kmerDict_correct[km] += 1

    #Throw all kmers that occur only once:
    kmerDict_correct = {k:v for k,v in kmerDict_correct.items() if v != 0}
    return kmerDict_correct

def createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn):
    if pfn:
        print "createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn)"
    kmerDict = collections.defaultdict(int)
    for segment,lineNr in generateSegmentsFrom_fq(fn,k,BF,pfn):
        assert isinstance(segment,str)
        for km in dbg.kmers(segment,k):
            kmerDict[km] += 1
        for km in dbg.kmers(dbg.twin(segment),k):
            kmerDict[km] += 1
    if pfn:
        print "Printing the kmerDict"
        for km in kmerDict:
            print km
        print kmerDict
    return kmerDict

def checkCorrectness(dict1,dict2):
    if not len(dict1)==len(dict2):
        print len(dict1), len(dict2)
        raise Exception("The dicts must have the same length")
    for km in dict1:
        if not km in dict2:
            raise Exception("The dicts must have the same kmers")

class Test_generateSegmentsFrom_fq(unittest.TestCase):
    def test_1(self):
        fn = ["Input/read1.fq"]
        k = 3
        numberOfKmers=100
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        segments = ["AAGTTGCGCTAGGGTTAAACTCGGCTAACTCGATTAACATCAGCCGTTTGGTGGCGCAGATTTGCTACTA"]
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())
        #self.assertEqual(segments,["AAGTTGCGCTAGGGTTAAACTCGGCTAACTCGATTAACATCAGCCGTTTGGTGGCGCAGATTTGCTACTA"])

    def test_2(self):
        fn = ["Input/read4.fq"]
        k = 3
        numberOfKmers=1000
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        segments = ["AAGTTGCGCTAGGGTTAAACTCGGCTAACTCGATTAACATCAGCCGTTTGGTGGCGCAGATTTGCTACTA","CGCCGCCATGCCGACCATCCCTTTCATCCCCGTACCAGACACGCTGACCATTGCCATGTTATTCAGATTG","TGATCTTTCAGATTGTAGAGTTTCATTTAGTTTACCAGTACTCGTGCGCCCGCCGAATCCAGGCGTCAAA","ACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCTGAAATTACGTTAGTCCCGTCATTAAAA"]
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())

    def test_3(self):
        fn = ["Input/read4.fq","Input/read1.fq"]
        k = 3
        numberOfKmers=1000
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        segments = ["AAGTTGCGCTAGGGTTAAACTCGGCTAACTCGATTAACATCAGCCGTTTGGTGGCGCAGATTTGCTACTA","CGCCGCCATGCCGACCATCCCTTTCATCCCCGTACCAGACACGCTGACCATTGCCATGTTATTCAGATTG","TGATCTTTCAGATTGTAGAGTTTCATTTAGTTTACCAGTACTCGTGCGCCCGCCGAATCCAGGCGTCAAA","ACTTTCCCTGGTTCTGGTCGCTCCCATGGCAGCACAGGCTGCTGAAATTACGTTAGTCCCGTCATTAAAA"]
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())

    def test_4(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_4.fq"]
        segments=["ATC","AAAAAA","CTA","TAG"]
        k = 3
        numberOfKmers=20
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())

    def test_5(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_5.fq"]
        segments=["AAA"]
        k = 3
        numberOfKmers=20
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())

    def test_6(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_6.fq"]
        segments=["AAAA"]
        k = 3
        numberOfKmers=20
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())

    def test_7(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_7.fq"]
        segments=["AAAAA"]
        k = 3
        numberOfKmers=20
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())

    def test_8(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_8.fq"]
        segments=["AAAAAA"]
        k = 3
        numberOfKmers=20
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())

    def test_9(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_9.fq"]
        segments=["AAA","AAA"]
        k = 3
        numberOfKmers=20
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())

    def test_10(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_10.fq"]
        segments=["AAA","TTT"]
        k = 3
        numberOfKmers=20
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)
        self.assertTrue(BF.hasAcceptableRatio())

    def test_11(self):
        fn = ["Input/Test_generateSegmentsFrom_fq/test_11.fq"]
        segments=["AAACCGTACGATCA"]
        k = 3
        numberOfKmers=30
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        kmerDict_correct = createCorrectDict(segments,k)
        kmerDict = createDictUsing_generateSegmentsFromfq(fn,k,BF,pfn=False)
        checkCorrectness(kmerDict,kmerDict_correct)        
        self.assertTrue(BF.hasAcceptableRatio())

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

class Test_BF_counter(unittest.TestCase):
    def test_1(self):
        fn = ["Input/Test_BF_counter/test_1.fq"]
        segments = ["AAA","TAG","AAA","TAG"]
        k = 3
        numberOfKmers = 50
        BF = Bloom.Bloom(0.0001,numberOfKmers,pfn=False)
        G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
        BF_counter(fn,BF,k,G,pfn=False,printProgress=False)
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

if __name__ == '__main__':
    unittest.main()