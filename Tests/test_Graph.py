#coding:utf8
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import unittest
import Graph
from compareGraphs import *
import BFAdder
import helpers
import Bloom
import Path

def myHelper(k=5,assertLegal=True):
    G = Graph.Graph(k,al=assertLegal)
    G_correct = Graph.Graph(k,al=assertLegal)
    return k, G, G_correct

#--------------------------------------------------------------------------
#-----------------Simple functions for working with graphs-----------------
#--------------------------------------------------------------------------
class Test_contains(unittest.TestCase):
    def test_1(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["AAAAA",[],[],2]
        G.addKmersFromAllContigs()
        self.assertTrue("AAAAA" in G)
        self.assertTrue("TTTTT" in G)
        self.assertFalse("AAAAT" in G)

    def test_2(self):
        k = 5
        G = Graph.Graph(k)
        s = "AAAAACGTTT"
        for km in dbg.kmers(s,k):
            G.contigs[G.getID()] = [km,[],[],2]
        G.addKmersFromAllContigs()
        for km in dbg.kmers(s,k):
            self.assertTrue(km in G)
            self.assertTrue(dbg.twin(km) in G)
        self.assertFalse("AAAAT" in G)

class Test_setINandOUT(unittest.TestCase):
    #def setOUT(self,ID,OUT):
    def allIsEqual(self,G1,G2):
        for index in G1.contigs:
            [c1,IN1,OUT1,COV1] = G1.contigs[index]
            [c2,IN2,OUT2,COV2] = G2.contigs[index]

            self.assertEqual(c1,c2)
            self.assertEqual(sorted(IN1),sorted(IN2))
            self.assertEqual(sorted(OUT1),sorted(OUT2))
            self.assertEqual(COV1,COV2)

        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()

        self.assertTrue(isSameGraph(G1,G2,alsoCompareWithNaive=False,relaxAssertions=True))

    def test_1(self):
        G1 = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["AAA",[],[],0]
        G1.setIN(0,[(0,True)])
        G1.setOUT(0,[(0,True)])
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        self.allIsEqual(G1,G_correct)
        
        
    def test_2(self):
        G1 = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[],[],0]
        G1.contigs[G1.getID()] = ["TTC",[],[],0]
        G1.setIN(0,[(0,False),(0,False)])
        G1.setOUT(0,[(1,True)])
        G1.setIN(1,[(0,True)])
        G_correct.contigs[G_correct.getID()] = ["ATT",[(0,False),(0,False)],[(1,True)],0]
        G_correct.contigs[G_correct.getID()] = ["TTC", [(0,True)], [], 0]
        self.allIsEqual(G1,G_correct)

class Test_COV(unittest.TestCase):
    def helper(self,G,correct_C_COV,correct_km_COV):
        G.addKmersFromAllContigs()
        C_COV = G.getAverageKmerCoverageOfContig(0)
        km_COV = G.getCoverageOfKmer("AAAAA")
        self.assertTrue(C_COV==correct_C_COV)
        self.assertTrue(km_COV==correct_km_COV)

    def assert_C_Cov_equals_trueCov(self,G,cID,trueCOV):
        C_COV = G.getAverageKmerCoverageOfContig(cID)
        self.assertTrue(C_COV==trueCOV)

    def assert_km_cov_equals_true_km_cov(self,G,km,trueCOV):
        km_COV = G.getCoverageOfKmer(km)
        self.assertTrue(km_COV==trueCOV)

    #Check that getAverageKmerCoverageOfContig, and
    #getCoverageOfKmer return the right values for a
    #value of COV which we added manually
    def test_m1(self):
        G = Graph.Graph(5)
        G.contigs[G.getID()] = ["AAAAA",[(0,True)],[(0,True)],1]
        self.helper(G,1,1)

    def test_m2(self):
        G = Graph.Graph(5)
        G.contigs[G.getID()] = ["AAAAA",[(0,True)],[(0,True)],5]
        self.helper(G,5,5)

    #Check whether addSegmentToGraph assigns the right COV value
    #to our contigs
    def test_a1(self):
        G = Graph.Graph(5)
        G.addSegmentToGraph("AAAAA")
        self.helper(G,2,2)

    def test_a2(self):
        G = Graph.Graph(5)
        G.addSegmentToGraph("AAAAA")
        G.addSegmentToGraph("AAAAA")
        self.helper(G,3,3)

    #Test splitting a contig
    def test_s1(self):
        G = Graph.Graph(5)
        G.addSegmentToGraph("AACCGA")
        G.addKmersFromAllContigs()
        C_COV = G.getAverageKmerCoverageOfContig(0)
        km_COV1 = G.getCoverageOfKmer("AACCG")
        km_COV2 = G.getCoverageOfKmer("ACCGA")
        self.assertTrue(C_COV==2)
        self.assertTrue(km_COV1==2)
        self.assertTrue(km_COV2==2)

    #Test splitting a contig
    def test_s2(self):
        G = Graph.Graph(5)
        G.addSegmentToGraph("AACCGA")
        G.addKmersFromAllContigs()
        C_COV = G.getAverageKmerCoverageOfContig(0)
        km_COV1 = G.getCoverageOfKmer("AACCG")
        km_COV2 = G.getCoverageOfKmer("ACCGA")
        self.assertTrue(C_COV==2)
        self.assertTrue(km_COV1==2)
        self.assertTrue(km_COV2==2)
        G.splitContig(0,5)
        self.assert_C_Cov_equals_trueCov(G,cID=1,trueCOV=2)
        self.assert_C_Cov_equals_trueCov(G,cID=2,trueCOV=2)

    def test_s3(self):
        G = Graph.Graph(5)
        G.contigs[G.getID()] = ["AACCGA",[],[],2]
        G.addKmersFromAllContigs()
        G.splitContig(0,5)
        self.assert_C_Cov_equals_trueCov(G,cID=1,trueCOV=1)
        self.assert_C_Cov_equals_trueCov(G,cID=2,trueCOV=1)

    def test_s4(self):
        G = Graph.Graph(5)
        G.contigs[G.getID()] = ["AACCGA",[],[],10]
        G.addKmersFromAllContigs()
        G.splitContig(0,5)
        self.assert_C_Cov_equals_trueCov(G,cID=1,trueCOV=5)
        self.assert_C_Cov_equals_trueCov(G,cID=2,trueCOV=5)

    def test_s5(self):
        G = Graph.Graph(5)
        G.contigs[G.getID()] = ["AACCGA",[],[],9]
        G.addKmersFromAllContigs()
        G.splitContig(0,5)
        self.assert_C_Cov_equals_trueCov(G,cID=1,trueCOV=5)
        self.assert_C_Cov_equals_trueCov(G,cID=2,trueCOV=4)

    def test_s6(self):
        G = Graph.Graph(5)
        G.addSegmentToGraph("AACCGAC")
        G.addKmersFromAllContigs()
        C_COV = G.getAverageKmerCoverageOfContig(0)
        km_COV1 = G.getCoverageOfKmer("AACCG")
        km_COV2 = G.getCoverageOfKmer("ACCGA")
        km_COV3 = G.getCoverageOfKmer("CCGAC")
        self.assertTrue(C_COV==2)
        self.assertTrue(km_COV1==2)
        self.assertTrue(km_COV2==2)
        self.assertTrue(km_COV3==2)
        G.splitContig(0,5)
        self.assert_C_Cov_equals_trueCov(G,cID=1,trueCOV=2)
        self.assert_C_Cov_equals_trueCov(G,cID=2,trueCOV=2)

    def test_s7(self):
        G = Graph.Graph(5)
        G.addSegmentToGraph("AACCGAC")
        G.addKmersFromAllContigs()
        C_COV = G.getAverageKmerCoverageOfContig(0)
        km_COV1 = G.getCoverageOfKmer("AACCG")
        km_COV2 = G.getCoverageOfKmer("ACCGA")
        km_COV3 = G.getCoverageOfKmer("CCGAC")
        self.assertTrue(C_COV==2)
        self.assertTrue(km_COV1==2)
        self.assertTrue(km_COV2==2)
        self.assertTrue(km_COV3==2)
        G.splitContig(0,6)
        self.assert_C_Cov_equals_trueCov(G,cID=1,trueCOV=2)
        self.assert_C_Cov_equals_trueCov(G,cID=2,trueCOV=2)

    def test_f1(self):
        G = Graph.Graph(5)
        G.addSegmentToGraph("AACCGAC")
        G.addKmersFromAllContigs()
        C_COV = G.getAverageKmerCoverageOfContig(0)
        self.assert_C_Cov_equals_trueCov(G,0,2)
        self.assert_km_cov_equals_true_km_cov(G,"AACCG",2)
        self.assert_km_cov_equals_true_km_cov(G,"ACCGA",2)
        self.assert_km_cov_equals_true_km_cov(G,"CCGAC",2)
        G.flipContig(0)
        self.assert_C_Cov_equals_trueCov(G,0,2)
        self.assert_km_cov_equals_true_km_cov(G,"AACCG",2)
        self.assert_km_cov_equals_true_km_cov(G,"ACCGA",2)
        self.assert_km_cov_equals_true_km_cov(G,"CCGAC",2)

    def test_c1(self):
        G = Graph.Graph(5)
        G.contigs[G.getID()] = ["AACCG",[],[],1]
        G.contigs[G.getID()] = ["ACCGA",[],[],1]
        G.addKmersFromAllContigs()
        G.connectSegment(0)
        self.assert_km_cov_equals_true_km_cov(G,"AACCG",1)
        self.assert_km_cov_equals_true_km_cov(G,"ACCGA",1)

    def test_m3(self):
        G = Graph.Graph(5)
        G.contigs[G.getID()] = ["AACCG",[],[(1,True)],1]
        G.contigs[G.getID()] = ["ACCGA",[(0,True)],[],1]
        G.addKmersFromAllContigs()
        G.connectSegment(0)
        G.mergeSegment(0)
        self.assert_km_cov_equals_true_km_cov(G,"AACCG",1)
        self.assert_km_cov_equals_true_km_cov(G,"ACCGA",1)

class Test_addINandOUT(unittest.TestCase):
    def allIsEqual(self,G1,G2,skipNaive=False):
        for index in G1.contigs:
            [c1,IN1,OUT1,COV1] = G1.contigs[index]
            [c2,IN2,OUT2,COV2] = G2.contigs[index]

            self.assertEqual(c1,c2)
            self.assertEqual(sorted(IN1),sorted(IN2))
            self.assertEqual(sorted(OUT1),sorted(OUT2))
            self.assertEqual(COV1,COV2)

        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()

        self.assertTrue(isSameGraph(G1,G2,not skipNaive))

    def test_1(self):
        #1 contig with connections to self
        G1 = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["AAA",[],[],0]
        G1.addIN(0,(0,True))
        G1.addOUT(0,(0,True))
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        self.allIsEqual(G1,G_correct)

    def test_2(self):
        #1 contig with same tuples more than once
        G1 = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["AAA",[],[],0]
        G1.addIN(0,(0,False))
        G1.addOUT(0,(0,False))
        G1.addIN(0,(0,False))
        G1.addOUT(0,(0,False))
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,False),(0,False)],[(0,False),(0,False)],0]
        self.allIsEqual(G1,G_correct,skipNaive=True)

    def test_3(self):
        #2 contigs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["AAA",[],[],0]
        G1.contigs[G1.getID()] = ["AAC",[],[],0]
        G1.addIN(0,(0,True))
        G1.addOUT(0,(0,True))
        G1.addIN(1,(0,True))
        G1.addOUT(0,(1,True))
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True),(1,True)],0]
        G_correct.contigs[G_correct.getID()] = ["AAC",[(0,True)],[],0]
        self.allIsEqual(G1,G_correct)

"""
class Test_kmerFunctions(unittest.TestCase):
    def test_1(self):
        #1 contig with connections to self
        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["AAA",[],[],0]
        self.assertTrue(G1.kmers=={})
        G1.addKmersFromContig(0)
        self.assertTrue(G1.isLegalDBG())
        self.assertTrue(G1.kmers=={"AAA": [0,0,True], "TTT": [0,0,False]})
        G1.deleteKmersFromContig(0)
        self.assertFalse(G1.kmers=={"AAA": [0,0,True], "TTT": [0,0,False]})
        self.assertTrue(G1.kmers=={})

    def test_2(self):
        #1 contig with connections to self
        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["AAA",[],[],0]
        G1.contigs[G1.getID()] = ["AAT",[],[],0]
        G1.addKmersFromContig(0)
        G1.addKmersFromContig(1)
        self.assertTrue(G1.isLegalDBG())
        self.assertTrue(G1.kmers=={"AAA": [0,0,True], "TTT": [0,0,False], "AAT": [1,0,True], "ATT": [1,0,False]})
        G1.deleteKmersFromContig(1)
        self.assertFalse(G1.kmers=={"AAA": [0,0,True], "TTT": [0,0,False], "AAT": [1,0,True], "ATT": [1,0,False]})
        self.assertTrue(G1.kmers=={"AAA": [0,0,True], "TTT": [0,0,False]})

    def test_3(self):
        #longer contig
        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["AAACCC",[],[],0]
        self.assertTrue(G1.kmers=={})
        G1.addKmersFromContig(0)
        self.assertTrue(G1.isLegalDBG())
        self.assertTrue(G1.kmers=={"AAA": [0,0,True], "AAC": [0,1,True], "ACC": [0,2,True], "CCC": [0,3,True] , "GGG": [0,0,False], "GGT": [0,1,False], "GTT": [0,2,False], "TTT": [0,3,False]})
        G1.deleteKmersFromContig(0)
        self.assertFalse(G1.kmers=={"AAA": [0,0,True], "AAC": [0,1,True], "ACC": [0,2,True], "CCC": [0,3,True] , "GGG": [0,0,False], "GGT": [0,1,False], "GTT": [0,2,False], "TTT": [0,3,False]})
        self.assertTrue(G1.kmers=={})
"""

class Test_getID(unittest.TestCase):
    def test_1(self):
        G1 = Graph.Graph(3)
        i = G1.getID()
        self.assertTrue(i==0)
        i = G1.getID()
        self.assertTrue(i==1)

class Test_isEmpty(unittest.TestCase):
    def test_1(self):
        G1 = Graph.Graph(3)
        self.assertTrue(G1.isEmpty())
        G1.contigs[0] = ["AAACCC",[],[],0]
        self.assertFalse(G1.isEmpty())

class Test_isLegal(unittest.TestCase):
    def test_1(self):
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[],[(1,True)],0]
        G.contigs[G.getID()] = ["AGT",[(0,True)],[],0]
        G.addKmersFromAllContigs()
        self.assertTrue(G.isLegalGraph())
        self.assertFalse(G.isLegalDBG())

    def test_2(self):
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[],[(1,True)],0]
        G.contigs[G.getID()] = ["AGT",[(0,True)],[],0]
        self.assertFalse(G.isLegalDBG())

    def test_3(self):
        G = Graph.Graph(k=3,pil=False)
        G.contigs[G.getID()] = ["CAG",[],[(1,True)],0]
        G.contigs[G.getID()] = ["AGT",[],[],0]
        G.addKmersFromAllContigs()
        self.assertFalse(G.isLegalDBG())

    def test_4(self):
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[],[(1,True)],0]
        G.contigs[G.getID()] = ["AGT",[(0,True)],[],0]
        G.addKmersFromAllContigs()
        self.assertTrue(G.isLegalGraph())
        self.assertFalse(G.isLegalDBG())

    def test_5(self):
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAGTTT",[],[],0]
        G.addKmersFromAllContigs()
        self.assertTrue(G.isLegalDBG())

    def test_6(self):
        G1 = Graph.Graph(k=3,pil=False)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[(1,True)],0]
        self.assertFalse(G1.isLegalDBG())

    def test_7(self):
        #Two graphs with 1 contig. Not identical (twins of each other but not same IN and OUT)
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[(0,False),(0,False)],[(1,True)],0]
        G2.contigs[G2.getID()] = ["AAT",[],[(0,False),(0,False)],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()
        self.assertFalse(G1.isLegalDBG())
        self.assertTrue(G2.isLegalDBG())

class Test_isLegal_connectionsToSelf(unittest.TestCase):
    def test_1(self):
        #Empty graph and twin(c)->c
        #ATT
        #AAT
        G1 = Graph.Graph(3)
        self.assertTrue(G1.isLegal_connectionsToSelf())
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[],0]
        self.assertTrue(G1.isLegal_connectionsToSelf())

    def test_2(self):
        #twin(c)->c
        G1 = Graph.Graph(3)
        self.assertTrue(G1.isLegal_connectionsToSelf())
        G1.contigs[0] = ["ATT",[(0,False)],[],0]
        self.assertFalse(G1.isLegal_connectionsToSelf())

    def test_3(self):
        #c->twin(c)
        G1 = Graph.Graph(3)
        self.assertTrue(G1.isLegal_connectionsToSelf())
        G1.contigs[0] = ["ATT",[],[(0,False),(0,False)],0]
        self.assertTrue(G1.isLegal_connectionsToSelf())

    def test_4(self):
        #c->twin(c)
        G1 = Graph.Graph(3)
        self.assertTrue(G1.isLegal_connectionsToSelf())
        G1.contigs[0] = ["ATT",[],[(0,False)],0]
        self.assertFalse(G1.isLegal_connectionsToSelf())

    def test_5(self):
        #c->
        G1 = Graph.Graph(3)
        self.assertTrue(G1.isLegal_connectionsToSelf())
        G1.contigs[0] = ["AAA",[(0,True)],[(0,True)],0]
        self.assertTrue(G1.isLegal_connectionsToSelf())

    def test_6(self):
        #c->c
        G1 = Graph.Graph(3)
        self.assertTrue(G1.isLegal_connectionsToSelf())
        G1.contigs[0] = ["AAA",[(0,True)],[],0]
        self.assertFalse(G1.isLegal_connectionsToSelf())

    def test_7(self):
        #c->c
        G1 = Graph.Graph(3)
        self.assertTrue(G1.isLegal_connectionsToSelf())
        G1.contigs[0] = ["AAA",[],[(0,True)],0]
        self.assertFalse(G1.isLegal_connectionsToSelf())

class Test_connect_a_to_b(unittest.TestCase):
    #a->b
    def test_1(self):
        #ATT
        #AAT
        #TTC
        #GAA
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ATT", [], [(1,True)], 2]
        G_correct.contigs[G_correct.getID()] = ["TTC", [(0,True)], [], 2]
        G.contigs[G.getID()] = ["ATT",[],[],1]
        G.contigs[G.getID()] = ["TTC",[],[],1]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connect_a_to_b(0,1,True,True)
        self.assertEqual(G.contigs[0][2],[(1,True)])
        self.assertEqual(G.contigs[1][1],[(0,True)])
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    #a->a
    def test_2(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAA", [(0,True)], [(0,True)], 0]
        G.contigs[G.getID()] = ["AAA",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connect_a_to_b(0,0,True,True)
        self.assertEqual(G.contigs[0][2],[(0,True)])
        self.assertEqual(G.contigs[0][1],[(0,True)])
        self.assertTrue(isSameGraph(G,G_correct))

    #a->twin(b)
    def test_3(self):
        #CAAA
        #TTTG
        #GATT
        #AATC
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["CAAA", [], [(1,False)], 0]
        G_correct.contigs[G_correct.getID()] = ["GATT", [], [(0,False)], 0]
        G.contigs[G.getID()] = ["CAAA",[],[],0]
        G.contigs[G.getID()] = ["GATT",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connect_a_to_b(0,1,True,False)
        self.assertEqual(G.contigs[0][2],[(1,False)])
        self.assertEqual(G.contigs[1][2],[(0,False)])
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))
    #a->twin(a)
    def test_4(self):
        #AAT
        #ATT
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAT", [], [(0,False),(0,False)], 0]
        G.contigs[G.getID()] = ["AAT",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connect_a_to_b(0,0,True,False)
        self.assertEqual(G.contigs[0][2],[(0,False),(0,False)])
        self.assertTrue(isSameGraph(G,G_correct))

    #twin(a)->b
    def test_5(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["CAAA", [(1,False)], [], 0]
        G_correct.contigs[G_correct.getID()] = ["TGGG", [(0,False)], [], 0]
        G.contigs[G.getID()] = ["CAAA",[],[],0]
        G.contigs[G.getID()] = ["TGGG",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connect_a_to_b(0,1,False,True)
        self.assertEqual(G.contigs[0][1] , [(1,False)])
        self.assertEqual(G.contigs[1][1] , [(0,False)])
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))


    #twin(a)->a
    def test_6(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ATT", [(0,False),(0,False)], [], 0]
        G.contigs[G.getID()] = ["ATT",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connect_a_to_b(0,0,False,True)
        self.assertEqual(G.contigs[0][1] , [(0,False),(0,False)])
        self.assertTrue(isSameGraph(G,G_correct))

    #twin(a)->twin(b)
    def test_7(self):
        #CAAA
        #TTTG

        #CCCA
        #TGGG
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["CAAA", [(1,True)], [], 0]
        G_correct.contigs[G_correct.getID()] = ["TGGG", [], [(0,True)], 0]
        G.contigs[G.getID()] = ["CAAA",[],[],0]
        G.contigs[G.getID()] = ["TGGG",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connect_a_to_b(0,1,False,False)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    #twin(a)->twin(a)
    def test_8(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["AAA", [(0,True)], [(0,True)], 0]
        G.contigs[G.getID()] = ["AAA",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connect_a_to_b(0,0,False,False)
        self.assertEqual(G.contigs[0][2],[(0,True)])
        self.assertEqual(G.contigs[0][1],[(0,True)])
        self.assertTrue(isSameGraph(G,G_correct))  

class Test_changeID_FromTo(unittest.TestCase):
    #má ekki nota isSameGraph af því að changeID_FromTo er hluti af því!
    def allIsEqual(self,G1,G2):
        for index in G1.contigs:
            [c1,IN1,OUT1,COV1] = G1.contigs[index]
            [c2,IN2,OUT2,COV2] = G2.contigs[index]

            self.assertEqual(c1,c2)
            self.assertEqual(sorted(IN1),sorted(IN2))
            self.assertEqual(sorted(OUT1),sorted(OUT2))
            self.assertEqual(COV1,COV2)
            self.assertTrue(G1.isLegalDBG())
            self.assertTrue(G2.isLegalDBG())
            m1 = 0
            m2 = 0
            for ID in G1.contigs:
                if ID>m1:
                    m1 = ID
            for ID in G2.contigs:
                if ID>m2:
                    m2 = ID
            self.assertEqual(G1.ID,m1+1)
            self.assertEqual(G2.ID,m2+1)

    def test_1(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[0] = ["AAA",[(0,True)],[(0,True)],0]
        G_correct.contigs[10] = ["AAA",[(10,True)],[(10,True)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G_correct.setID(11)
        G.changeID_FromTo(0,10)
        self.allIsEqual(G,G_correct)
        #self.assertTrue(isSameGraph(G,G_correct))

    def test_2(self):
        #change ID of INs:
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G_correct.contigs[10] = ["ATT",[(10,False),(10,False)],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G_correct.setID(11)
        G.changeID_FromTo(0,10)
        self.allIsEqual(G,G_correct)
        #self.assertTrue(isSameGraph(G,G_correct))

    def test_3(self):
        #change ID of OUTs:
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["AAT",[],[(0,False),(0,False)],0]
        G_correct.contigs[10] = ["AAT",[],[(10,False),(10,False)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G_correct.setID(11)
        G.changeID_FromTo(0,10)
        self.allIsEqual(G,G_correct)
        #self.assertTrue(isSameGraph(G,G_correct))

    def test_4(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G_correct.contigs[3567] = ["ATT",[(3567,False),(3567,False)],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G_correct.setID(3568)
        G.changeID_FromTo(0,3567)
        self.allIsEqual(G,G_correct)
        #self.assertTrue(isSameGraph(G,G_correct))

    def test_5(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["AAT",[],[(0,False),(0,False)],0]
        G_correct.contigs[3567] = ["AAT",[],[(3567,False),(3567,False)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G_correct.setID(3568)
        G.changeID_FromTo(0,3567)
        self.allIsEqual(G,G_correct)
        #self.assertTrue(isSameGraph(G,G_correct))

    def test_6(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[0] = ["TTT",[(0,True)],[(0,True)],0]
        G_correct.contigs[10] = ["TTT",[(10,True)],[(10,True)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G_correct.setID(11)
        G.changeID_FromTo(0,10)
        self.allIsEqual(G,G_correct)
        #self.assertTrue(isSameGraph(G,G_correct))

    def test_7(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[0] = ["TTT",[(0,True)],[(0,True),(1,True)],0]
        G.contigs[1] = ["TTA",[(0,True)],[(1,False),(1,False)],0]
        #twin: TAA
        G_correct.contigs[10] = ["TTT",[(10,True)],[(10,True),(1,True)],0]
        G_correct.contigs[1] = ["TTA",[(10,True)],[(1,False),(1,False)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G_correct.setID(11)
        G.changeID_FromTo(0,10)
        self.allIsEqual(G,G_correct)
        #self.assertTrue(isSameGraph(G,G_correct))

    #More complex tests
    def test_8(self):
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["ATGAT",[(1,True),(0,True),(0,False),(0,False)],
        [(1,False),(0,True),(0,False),(0,False)],0]
        G.contigs[G.getID()] = ["AAT",[],[(1,False),(1,False),(0,True),(0,False)],0]
        G.addKmersFromAllContigs()

        G_naive = G.createNaive()

        self.allIsEqual(G,G_naive)
        self.assertTrue(isSameGraph(G,G_naive))

    def test_9(self):
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["ATT",[(0,False),(0,False),(1,True),(1,False)],[],0]
        G.contigs[G.getID()] = ["ATCAT",[(0,False),(1,True),(1,False),(1,False)],[(0,True),(1,True),(1,False),(1,False)],0]
        G.addKmersFromAllContigs()
        
        G_naive = G.createNaive()

        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False))

    def test_10(self):
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G.addKmersFromAllContigs()
        G.changeID_FromTo(0,1)
        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["ATT",[(1,False),(1,False)],[],0]
        G_correct.setID(2)
        G_correct.addKmersFromAllContigs()
        self.allIsEqual(G,G_correct)

    def test_11(self):
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["ATCAT",[(0,True),(0,False),(0,False)],[(0,True),(0,False),(0,False)],0]
        G.addKmersFromAllContigs()
        G.changeID_FromTo(0,1)
        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["ATCAT",[(1,True),(1,False),(1,False)],[(1,True),(1,False),(1,False)],0]
        G_correct.setID(2)
        G_correct.addKmersFromAllContigs()
        self.allIsEqual(G,G_correct)

    def test_12(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ATT",[(0,False),(0,False),(1,True),(1,False)],[],0]
        G.contigs[G.getID()] = ["ATCAT",[(0,False),(1,True),(1,False),(1,False)],[(0,True),(1,True),(1,False),(1,False)],0]
        G.addKmersFromAllContigs()
        G_correct.contigs[2] = ["ATT",[(2,False),(2,False),(3,True),(3,False)],[],0]
        G_correct.contigs[3] = ["ATCAT",[(2,False),(3,True),(3,False),(3,False)],[(2,True),(3,True),(3,False),(3,False)],0]
        G_correct.addKmersFromAllContigs()
        G_correct.setID(4)
        G.changeID_FromTo(0,2)
        G.changeID_FromTo(1,3)
        self.allIsEqual(G,G_correct)

    def test_13(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["TTT",[(0,True),(1,True)],[(0,True),(1,True)],0]
        G.contigs[G.getID()] = ["TTCTT",[(0,True),(1,True)],[(0,True),(1,True)],0]
        G_correct.contigs[0] = ["TTT",[(0,True),(2,True)],[(0,True),(2,True)],0]
        G_correct.contigs[2] = ["TTCTT",[(0,True),(2,True)],[(0,True),(2,True)],0]
        G_correct.setID(2)
        G.changeID_FromTo(1,2)

class Test_flipContig(unittest.TestCase):
    #Note: We can't use compareGraphs to test this function because we use 
    #flipContig as a part of fixTwins in compareGraphs.
    #we therefore need to manually verify the correctness of the tests 
    def printFlipPrint(self,G,cID,text):
        #print "-----------" + text + "-----------"
        #print "Before:"
        G.addKmersFromAllContigs()
        #G.printContigs()
        G.flipContig(cID)
        #print "\nAfter:"
        #G.printContigs()
        #print "\n\n"
        self.assertTrue(G.isLegalDBG())

    def test_1(self):
        #1 contig c->c
        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        self.printFlipPrint(G1,0,"Test 1")

        #Test for correctness:
        [c,IN,OUT,COV] = G1.contigs[0]
        self.assertEqual(c,"TTT")
        self.assertEqual(IN,[(0,True)])
        self.assertEqual(OUT,[(0,True)])

    def test_2(self):
        #1 contig twin(c)->c
        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATC",[(0,False),(0,False)],[],0]
        self.printFlipPrint(G1,0,"Test 2")

        #Test for correctness:
        [c,IN,OUT,COV] = G1.contigs[0]
        self.assertEqual(c,"GAT")
        self.assertEqual(IN,[])
        self.assertEqual(OUT,[(0,False),(0,False)])

    def test_3(self):
        #1 contig c->twin(c)
        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["GAT",[],[(0,False),(0,False)],0]
        self.printFlipPrint(G1,0,"Test 3")

        #Test for correctness:
        [c,IN,OUT,COV] = G1.contigs[0]
        self.assertEqual(c,"ATC")
        self.assertEqual(IN,[(0,False),(0,False)])
        self.assertEqual(OUT,[])

    def test_4(self):
        #2 contigs. c0->c1
        #c0: ATT
        #c1: TTC
        #ATT->TTC
        #AAT<-GAA
        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[],[(1,True),(2,True)],0]
        G1.contigs[G1.getID()] = ["TTC",[(0,True)],[],0]
        G1.contigs[G1.getID()] = ["TTG",[(0,True)],[],0]
        self.printFlipPrint(G1,0,"Test 4")

        #Test for correctness:
        [c0,IN0,OUT0,COV0] = G1.contigs[0]
        [c1,IN1,OUT1,COV1] = G1.contigs[1]
        [c2,IN2,OUT2,COV2] = G1.contigs[2]
        self.assertEqual(c0,"AAT")
        self.assertEqual(c1,"TTC")
        self.assertEqual(c2,"TTG")
        self.assertEqual(IN0,[(1,False),(2,False)])
        self.assertEqual(OUT0,[])
        self.assertEqual(IN1,[(0,False)])
        self.assertEqual(OUT1,[])
        self.assertEqual(IN2,[(0,False)])
        self.assertEqual(OUT2,[])

    def test_5(self):
        #2 contigs. c0->twin(c1)
        #c0: ATT
        #c1: GAA

        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[],[(1,False),(2,True)],0]
        G1.contigs[G1.getID()] = ["GAA",[],[(0,False)],0]
        G1.contigs[G1.getID()] = ["TTT",[(0,True),(2,True)],[(2,True)],0]
        self.printFlipPrint(G1,0,"Test 5")

        #After flipping c0:
        #   c0: AAT
        #   c1: GAA

        #Test for correctness:
        [c0,IN0,OUT0,COV0] = G1.contigs[0]
        [c1,IN1,OUT1,COV1] = G1.contigs[1]
        [c2,IN2,OUT2,COV2] = G1.contigs[2]
        self.assertEqual(c0,"AAT")
        self.assertEqual(c1,"GAA")
        self.assertEqual(c2,"TTT")
        self.assertEqual(IN0,[(1,True),(2,False)])
        self.assertEqual(OUT0,[])
        self.assertEqual(IN1,[])
        self.assertEqual(OUT1,[(0,True)])
        self.assertEqual(IN2,[(0,False),(2,True)])
        self.assertEqual(OUT2,[(2,True)])

    def test_6(self):
        #2 contigs. twin(c0)->c1
        #c0: GGC. twin GCC
        #c1: CCA. twin TGG

        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["GGC",[(1,False),(2,True)],[],0]
        G1.contigs[G1.getID()] = ["CCA",[(0,False)],[],0]
        G1.contigs[G1.getID()] = ["AGG",[],[(0,True)],0]
        self.printFlipPrint(G1,1,"Test 5")

        #After flipping c1:
        #   c0: GGC. twin GCC
        #   c1: TGG. twin CCA

        #Test for correctness:
        [c0,IN0,OUT0,COV0] = G1.contigs[0]
        [c1,IN1,OUT1,COV1] = G1.contigs[1]
        [c2,IN2,OUT2,COV2] = G1.contigs[2]
        self.assertEqual(c0,"GGC")
        self.assertEqual(c1,"TGG")
        self.assertEqual(c2,"AGG")
        self.assertEqual(IN0,[(1,True),(2,True)])
        self.assertEqual(OUT0,[])
        self.assertEqual(IN1,[])
        self.assertEqual(OUT1,[(0,True)])
        self.assertEqual(IN2,[])
        self.assertEqual(OUT2,[(0,True)])

    def test_7(self):
        #just like test 4 but with different second contig (prep for test 8)
        #c0: ATT
        #c2: TTG
        G1 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[(2,True),(3,True)],0]
        G1.contigs[2] = ["TTG",[(0,True)],[],0]
        G1.contigs[3] = ["TTC",[(0,True)],[],0]
        G1.setID(4)
        self.printFlipPrint(G1,0,"Test 7")

        #After flip:
        #c0: AAT. twin ATT
        #c1: TTC
        #c2: TTG

        #Test for correctness:
        [c0,IN0,OUT0,COV0] = G1.contigs[0]
        [c2,IN2,OUT2,COV2] = G1.contigs[2]
        [c3,IN3,OUT3,COV3] = G1.contigs[3]
        self.assertEqual(c0,"AAT")
        self.assertEqual(c2,"TTG")
        self.assertEqual(c3,"TTC")
        self.assertEqual(IN0,[(2,False),(3,False)])
        self.assertEqual(OUT0,[])
        self.assertEqual(IN2,[(0,False)])
        self.assertEqual(OUT2,[])
        self.assertEqual(IN3,[(0,False)])
        self.assertEqual(OUT3,[])

    def test_8(self):
        #3 contigs. c0->c1
        #           c0->c2
        #(we took c0 and c1 from test 4)
        #c0: ATT
        #c1: TTC
        #c2: TTG
        #ATT->TTC
        #AAT<-GAA
        G1 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[],[(1,True),(2,True)],0]
        G1.contigs[G1.getID()] = ["TTC",[(0,True)],[],0]
        G1.contigs[G1.getID()] = ["TTG",[(0,True)],[],0]
        self.printFlipPrint(G1,0,"Test 8")

        #After flip:
        #c0: AAT. twin ATT
        #c1: TTC
        #c2: TTG

        #Test for correctness:
        [c0,IN0,OUT0,COV0] = G1.contigs[0]
        [c1,IN1,OUT1,COV1] = G1.contigs[1]
        [c2,IN2,OUT2,COV2] = G1.contigs[2]
        self.assertEqual(c0,"AAT")
        self.assertEqual(c1,"TTC")
        self.assertEqual(c2,"TTG")
        self.assertEqual(IN0,[(1,False),(2,False)])
        self.assertEqual(OUT0,[])
        self.assertEqual(IN1,[(0,False)])
        self.assertEqual(OUT1,[])
        self.assertEqual(IN2,[(0,False)])
        self.assertEqual(OUT2,[])

    def test_9(self):
        #the contigs from test_8 except now we have a few more contigs
        #3 contigs. c0->c1
        #           c0->c2
        #           c3->c0
        #           twin(c4)->c0
        #           c0->twin(c5)
        #c0: ATT. twin AAT
        #c1: TTC. twin GAA
        #c2: TTG. twin CAA
        #c3: CAT. twin ATG
        #c4: ATC. twin GAT
        #c5: TAA. twin TTA
        #ATT->TTC
        #AAT<-GAA
        G1 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(3,True),(4,False)],[(1,True),(2,True),(5,False)],0]
        G1.contigs[1] = ["TTC",[(0,True)],[],0]
        G1.contigs[2] = ["TTG",[(0,True)],[],0]
        G1.contigs[3] = ["CAT",[],[(0,True)],0]
        G1.contigs[4] = ["ATC",[(0,False)],[],0]
        G1.contigs[5] = ["TAA",[],[(0,False)],0]
        G1.setID(6)
        self.printFlipPrint(G1,0,"Test 9")

        #After flip:
        #c0: AAT. twin ATT

        #Test for correctness:
        [c0,IN0,OUT0,COV0] = G1.contigs[0]
        [c1,IN1,OUT1,COV1] = G1.contigs[1]
        [c2,IN2,OUT2,COV2] = G1.contigs[2]
        [c3,IN3,OUT3,COV3] = G1.contigs[3]
        [c4,IN4,OUT4,COV4] = G1.contigs[4]
        [c5,IN5,OUT5,COV5] = G1.contigs[5]
        self.assertEqual(c0,"AAT")
        self.assertEqual(c1,"TTC")
        self.assertEqual(c2,"TTG")
        self.assertEqual(c3,"CAT")
        self.assertEqual(c4,"ATC")
        self.assertEqual(c5,"TAA")
        self.assertEqual(IN0,[(1,False),(2,False),(5,True)])
        self.assertEqual(OUT0,[(3,False),(4,True)])
        self.assertEqual(IN1,[(0,False)])
        self.assertEqual(OUT1,[])
        self.assertEqual(IN2,[(0,False)])
        self.assertEqual(OUT2,[])
        self.assertEqual(IN3,[])
        self.assertEqual(OUT3,[(0,False)])
        self.assertEqual(IN4,[(0,True)])
        self.assertEqual(OUT4,[])
        self.assertEqual(IN5,[])
        self.assertEqual(OUT5,[(0,True)])

    def test_10(self):
        #Two graphs with 1 contig. Not identical (twins of each other but not same IN and OUT)
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[0] = ["AAT",[],[(0,False)],0]
        self.assertTrue(G1.isLegal_connectionsToSelf())
        G2_legal = G2.isLegal_connectionsToSelf()
        self.assertFalse(G2_legal)
        with self.assertRaises(AssertionError):
            self.printFlipPrint(G2,0,"Test 10")
    
    def allIsEqual(self,G1,G2,index):
        [c1,IN1,OUT1,COV1] = G1.contigs[index]
        [c2,IN2,OUT2,COV2] = G2.contigs[index]

        self.assertEqual(c1,c2)
        self.assertEqual(sorted(IN1),sorted(IN2))
        self.assertEqual(sorted(OUT1),sorted(OUT2))
        self.assertEqual(COV1,COV2)

    def test_11(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CGTA",[(0,False),(0,False)],[(0,False),(0,False)],0]
        G_correct.contigs[G_correct.getID()] = ["TACG",[(0,False),(0,False)],[(0,False),(0,False)],0]
        self.printFlipPrint(G,0,"Test 11")
        self.allIsEqual(G,G_correct,0)
        G.addKmersFromAllContigs()
        self.assertTrue(G.isLegalDBG())

#--------------------------------------------------------------------------
#--------------------addSegmentToGraph and subfunctions--------------------
#--------------------------------------------------------------------------
"""
class Test_addSegmentToGraph_OLD(unittest.TestCase):
    #Afrit af testunum úr addSegmentWithNoSeenKmers nema bara breytt í kall á addSegmentToGraph
    def test_1(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],2]
        G.addKmersFromAllContigs()

        G.addSegmentToGraph("TCC")

        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["ACC",[],[(2,True)],0]
        G_correct.contigs[2] = ["CCA",[(1,True),(3,True)],[],0]
        G_correct.contigs[3] = ["TCC",[],[(2,True)],0]
        G_correct.addKmersFromAllContigs()
        G_correct.setID(4)
        self.assertTrue(isSameGraph(G,G_correct))
        

    def test_2(self):
        #create G and G_correct:
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],10]
        #twin:TGGT
        G.contigs[G.getID()] = ["TTT",[(1,True)],[(1,True)],10]

        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["ACC",[],[(2,True)],0]
        G_correct.contigs[2] = ["CCA",[(1,True),(3,True)],[],0]
        G_correct.contigs[3] = ["TCC",[],[(2,True)],0]
        G_correct.contigs[4] = ["TTT",[(4,True)],[(4,True)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G_correct.setID(5)

        G.addSegmentToGraph("TCC")

        self.assertTrue(isSameGraph(G,G_correct))

    def test_3(self):
        #create G and G_correct:
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],10]
        #twin:TGGT
        G.contigs[G.getID()] = ["TTT",[(1,True),(2,True)],[(1,True),(2,True)],10]
        G.contigs[G.getID()] = ["TTCTT",[(1,True),(2,True)],[(1,True),(2,True)],10]

        G.addKmersFromAllContigs()
        G.addSegmentToGraph("TCC")
        G_naive = G.createNaive()
        self.assertTrue(isSameGraph(G,G_naive,alsoCompareWithNaive=False))

    def test_4(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],10] #twin:TGGT
        G.addKmersFromAllContigs()

        G.addSegmentToGraph("TCCT") #twin:AGGA

        #ACC   ->   CCA
        #TCC   ->   CCT

        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["ACC",[],[(2,True),(4,True)],0]
        G_correct.contigs[2] = ["CCA",[(1,True),(3,True)],[],0]
        G_correct.contigs[3] = ["TCC",[],[(2,True),(4,True)],0]
        G_correct.contigs[4] = ["CCT",[(1,True),(3,True)],[],0]
        G_correct.setID(5)
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_5(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],0] #twin:TGGT
        G.addKmersFromAllContigs()

        G.addSegmentToGraph("TCTC") #twin:GAGA

        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ACCA",[],[],0]
        G_correct.contigs[G_correct.getID()] = ["TCTC",[(1,True)],[(1,True)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_6(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],0] #twin:TGGT
        G.contigs[G.getID()] = ["CTA",[],[(1,False),(1,False)],0] #twin:TAG
        G.addKmersFromAllContigs()

        G.addSegmentToGraph("TCTC") #twin:GAGA

        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ACCA",[],[],0]
        G_correct.contigs[G_correct.getID()] = ["CTA",[(2,True)],[(1,False),(1,False)],0]#twin:TAG
        G_correct.contigs[G_correct.getID()] = ["CTCT",[(2,True)],[(1,True),(2,True)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    #Þetta test er til að einfalda test_6 niður í seinni helming keyrslunnar
    #test_6 var flóknara en ég bjóst við
    def test_7(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],10] #twin:TGGT
        G.contigs[G.getID()] = ["TCTA",[],[(1,False),(1,False)],10] #twin:TAGA
        G.addKmersFromAllContigs()

        G.addSegmentToGraph("CTC") #twin:GAG

        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ACCA",[],[],0]
        G_correct.contigs[G_correct.getID()] = ["CTA",[(2,True)],[(1,False),(1,False)],0]#twin:TAG
        G_correct.contigs[G_correct.getID()] = ["CTCT",[(2,True)],[(1,True),(2,True)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    #a->a and a->twin(a) and twin(a)->a
    def test_8(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)

        G.addSegmentToGraph("ATCAT")    #twin:ATGAT

        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ATCAT",[(0,True),(0,False),(0,False)],[(0,True),(0,False),(0,False)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_9(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ATT",[(0,False),(0,False)],[],0] #twin:AAT
        G.addKmersFromAllContigs()

        G.addSegmentToGraph("ATCAT")    #twin:ATGAT

        G_correct = Graph.Graph(3)
        G_correct.contigs[G_correct.getID()] = ["ATT",[(0,False),(0,False),(1,True),(1,False)],[],0]
        G_correct.contigs[G_correct.getID()] = ["ATCAT",[(0,False),(1,True),(1,False),(1,False)],[(0,True),(1,True),(1,False),(1,False)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    #Test skrifuð sérstaklega til að kanna hvort addSegmentToGraph
    #hegði sér rétt þegar það fær inn kmera sem við höfum séð áður
    def test_10(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAA")
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],1]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_11(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAA")
        G.addSegmentToGraph("AAA")
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],2]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_12(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAAA")
        G.addKmersFromContig(0)
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],2]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_13(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAA")
        G.addSegmentToGraph("TTT")
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],2]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_14(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAACAAA")
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True),(1,True)],[(0,True),(1,True)],2]
        G_correct.contigs[G_correct.getID()] = ["AACAA",[(0,True),(1,True)],[(0,True),(1,True)],6]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_15(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAACTTT")
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True),(1,True),(1,False)],2]
        G_correct.contigs[G_correct.getID()] = ["AACTT",[(0,True)],[(0,False)],6]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_16(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAAAAAAAA")
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],8]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    #Flókinn contig sem ætti ekki að splittast
    def test_17(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("ATCAT")
        G_correct.contigs[G_correct.getID()] = ["ATCAT",[(0,True),(0,False),(0,False)],[(0,True),(0,False),(0,False)],0] #twin: ATGAT
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_18(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("ATCC")
        G_correct.contigs[G_correct.getID()] = ["ATCC",[(0,False),(0,False)],[],0] #twin: GGAT
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    def test_19(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AATAA")
        G_correct.contigs[3] = ["TAAT",[(1,True),(3,False),(3,False)],[(1,True),(3,False),(3,False)],0] #twin: GGAT
        G_correct.contigs[1] = ["ATA",[(3,True),(1,False),(1,False)],[(3,True),(1,False),(1,False)],0] #twin: GGAT
        G_correct.addKmersFromAllContigs()
        G_correct.setID(4)
        self.assertTrue(isSameGraph(G,G_correct))

    def test_20(self):
        G = Graph.Graph(3)
        G.addSegmentToGraph("AAATCCC")
        G_naive = G.createNaive()
        self.assertTrue(G.equalsNaive())

    def test_21(self):
        G = Graph.Graph(3)
        G.addSegmentToGraph("AAATCCCA")
        G_naive = G.createNaive()
        self.assertTrue(G.equalsNaive())

    def test_22(self):
        G = Graph.Graph(3,pfn=False,ps=False,al=False,pil=False)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAA")
        G.addSegmentToGraph("AAA")
        G.addSegmentToGraph("AAA")
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],7]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=True,pfn=False,ps=False,relaxAssertions=False))

    def test_23(self):
        G = Graph.Graph(3,pfn=False,ps=False,al=False,pil=False)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAAA")
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],7]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=True,pfn=False,ps=False,relaxAssertions=False))

    def test_24(self):
        G = Graph.Graph(3,pfn=False,ps=False,al=False,pil=False)
        G_correct = Graph.Graph(3)
        G.addSegmentToGraph("AAAAA")
        G_correct.contigs[G_correct.getID()] = ["AAA",[(0,True)],[(0,True)],7]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=True,pfn=False,ps=False,relaxAssertions=False))
"""

class Test_addSegmentToGraph(unittest.TestCase):
    def test_1(self):
        k, G, Gc = myHelper(5)
        G.addSegmentToGraph("AAAAA")
        Gc.contigs[Gc.getID()] = ["AAAAA",[(0,True)],[(0,True)],2]
        Gc.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,Gc))

    def test_2(self):
        k, G, Gc = myHelper(5)
        G.addSegmentToGraph("AAAAA")
        G.addSegmentToGraph("CCCCC")
        Gc.contigs[Gc.getID()] = ["AAAAA",[(0,True)],[(0,True)],2]
        Gc.contigs[Gc.getID()] = ["CCCCC",[(1,True)],[(1,True)],2]
        Gc.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,Gc))

    def test_3(self):
        k, G, Gc = myHelper(5)
        G.addSegmentToGraph("AAAAA")
        G.addSegmentToGraph("CCCCC")
        G.addSegmentToGraph("AAAAA")
        Gc.contigs[Gc.getID()] = ["AAAAA",[(0,True)],[(0,True)],3]
        Gc.contigs[Gc.getID()] = ["CCCCC",[(1,True)],[(1,True)],2]
        Gc.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,Gc))

    def test_4(self):
        k, G, Gc = myHelper(5)
        G.addSegmentToGraph("AAAAA")
        G.addSegmentToGraph("CCCCC")
        G.addSegmentToGraph("AAAAACC")
        Gc.contigs[Gc.getID()] = ["AAAAA",[(0,True)],[(0,True),(2,True)],3]
        Gc.contigs[Gc.getID()] = ["CCCCC",[(1,True)],[(1,True)],2]
        Gc.contigs[Gc.getID()] = ["AAAACC",[(0,True)],[],2]
        Gc.addKmersFromAllContigs()
        #G.printContigs()
        self.assertTrue(isSameGraph(G,Gc))
    
    def test_5(self):
        k, G, Gc = myHelper(5)
        G.addSegmentToGraph("AAAAACC")
        Gc.contigs[Gc.getID()] = ["AAAAA",[(0,True)],[(0,True),(1,True)],2]
        Gc.contigs[Gc.getID()] = [ "AAAACC",[(0,True)],[],2]
        Gc.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,Gc))

    def test_6(self):
        k, G, Gc = myHelper(5)
        G.addSegmentToGraph("AAAAA")
        G.addSegmentToGraph("AAAAACC")
        Gc.contigs[Gc.getID()] = ["AAAAA",[(0,True)],[(0,True),(1,True)],3]
        Gc.contigs[Gc.getID()] = [ "AAAACC",[(0,True)],[],2]
        Gc.addKmersFromAllContigs()
        #G.printContigs()
        self.assertTrue(isSameGraph(G,Gc))

    def test_7(self):
        k, G, Gc = myHelper(5)
        G.addSegmentToGraph("AAAAA")
        G.addSegmentToGraph("AAAACC")
        Gc.contigs[Gc.getID()] = ["AAAAA",[(0,True)],[(0,True),(1,True)],2]
        Gc.contigs[Gc.getID()] = [ "AAAACC",[(0,True)],[],2]
        Gc.addKmersFromAllContigs()
        #G.printContigs()
        self.assertTrue(isSameGraph(G,Gc))

#for s in self.splitSegmentOnSeenKmers(segment)
class Test_splitSegmentOnSeenKmers(unittest.TestCase):
    #Ekkert split
    def test_1(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["TCCAT",[],[],2]
        G.addKmersFromAllContigs()
        s = "ACCATT"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue("ACCATT" in seqs)
        self.assertTrue(len(seqs)==1)

    #len(s)==k
    def test_2(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["TCCAT",[],[],2]
        G.addKmersFromAllContigs()
        s = "ACCAT"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue("ACCAT" in seqs)
        self.assertTrue(len(seqs)==1)

    #len(s)==k og s in G
    def test_3(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["TCCAT",[],[],2]
        G.addKmersFromAllContigs()
        s = "TCCAT"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue(seqs==[])
        self.assertTrue(len(seqs)==0)
        self.assertTrue(G.contigs[0][3]==3)

    #split
    def test_4(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["TCCAT",[],[],2]
        G.addKmersFromAllContigs()
        s = "CTCCATC"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue("CTCCA" in seqs)
        self.assertTrue(  "CCATC" in seqs)
        self.assertTrue(len(seqs)==2)
        self.assertTrue(G.contigs[0][3]==3)

    #split because a k-mer occurs twice in s
    def test_5(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["TCCAT",[],[],2]
        G.addKmersFromAllContigs()
        s = "CAAAAATGGGAAAAACGC"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue("CAAAAATGGGAAAA" in seqs)
        self.assertTrue(          "AAAACGC" in seqs)
        self.assertTrue(len(seqs)==2)
        self.assertTrue(G.contigs[0][3]==2)

    #split because a k-mer which has been seen before
    #occurs twice in s
    def test_6(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["AAAAA",[],[],2]
        G.addKmersFromAllContigs()
        s = "CAAAAATGGGAAAAACGC"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue("CAAAA" in seqs)
        self.assertTrue( "AAAATGGGAAAA" in seqs)
        self.assertTrue(         "AAAACGC" in seqs)
        self.assertTrue(len(seqs)==3)
        self.assertTrue(G.contigs[0][3]==4)

    #Fyrsti veldur splitti
    def test_7(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["AAAAA",[],[],2]
        G.addKmersFromAllContigs()
        s = "AAAAATGGGAAAA"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue("AAAATGGGAAAA" in seqs)
        self.assertTrue(len(seqs)==1)
        self.assertTrue(G.contigs[0][3]==3)

    #Síðasti veldur splitti
    def test_8(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["AAAAA",[],[],2]
        G.addKmersFromAllContigs()
        s = "AAAATGGGAAAAA"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue("AAAATGGGAAAA" in seqs)
        self.assertTrue(len(seqs)==1)
        self.assertTrue(G.contigs[0][3]==3)
        
    def test_9(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["AAAAA",[(0,True)],[(0,True)],2]
        G.addKmersFromAllContigs()
        s = "AAAAATGGGAAAAACGC"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue("AAAATGGGAAAA" in seqs)
        self.assertTrue(        "AAAACGC" in seqs)
        self.assertTrue(len(seqs)==2)
        self.assertTrue(G.contigs[0][3]==4)

    def test_10(self):
        k =5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["AAAAA",[(0,True)],[(0,True)],2]
        G.contigs[G.getID()] = ["CCCCC",[(1,True)],[(1,True)],2]
        G.addKmersFromAllContigs()
        s = "AAAAACC"
        seqs = []
        for i,j in G.splitSegmentOnSeenKmers(s,dbg.twin(s),len(s)):
            seqs.append(s[i:j])
        self.assertTrue("AAAACC" in seqs)
        self.assertTrue(len(seqs)==1)
        self.assertTrue(G.contigs[0][3]==3)

class Test_splitSegmentOnConnToOthers(unittest.TestCase):
    #split
    def test_1(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["TCCAT",[],[],2]
        G.addKmersFromAllContigs()
        s = "ACCATT"
        seqs = []
        for i,j in G.splitSegmentOnConnToOthers(s,len(s),0,len(s),0):
            seqs.append(s[i:j])
        self.assertTrue("ACCAT" in seqs)
        self.assertTrue("CCATT" in seqs)
        self.assertTrue(len(seqs)==2)

    #ekkert split því len(s)==k
    def test_2(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["TCCAT",[],[],2]
        G.addKmersFromAllContigs()
        s = "CCATT"
        seqs = []
        for i,j in G.splitSegmentOnConnToOthers(s,len(s),0,len(s),0):
            seqs.append(s[i:j])
        self.assertTrue("CCATT" in seqs)
        self.assertTrue(len(seqs)==1)

    #Sama og 1 nema fullt af auka contigum
    def test_3(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["TCCAT",[],[],2]
        G.contigs[G.getID()] = ["AAAAA",[],[],2]
        G.contigs[G.getID()] = ["CCCCC",[],[],2]
        G.contigs[G.getID()] = ["ACCGGAAG",[],[],8]
        G.addKmersFromAllContigs()
        s = "ACCATT"
        seqs = []
        for i,j in G.splitSegmentOnConnToOthers(s,len(s),0,len(s),0):
            seqs.append(s[i:j])
        self.assertTrue("ACCAT" in seqs)
        self.assertTrue( "CCATT" in seqs)
        self.assertTrue(len(seqs)==2)

    #len(s)>k og ekkert split
    def test_4(self):
        k = 5
        G = Graph.Graph(k)
        G.contigs[G.getID()] = ["TCCAT",[],[],2]
        G.addKmersFromAllContigs()
        s = "AAACCAAC"
        seqs = []
        for i,j in G.splitSegmentOnConnToOthers(s,len(s),0,len(s),0):
            seqs.append(s[i:j])
        self.assertTrue("AAACCAAC" in seqs)
        self.assertTrue(len(seqs)==1)

"""
class Test_addSegmentAlreadySplit(unittest.TestCase):
    def test_1(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],10]
        G.addKmersFromAllContigs()

        G.addSegmentAlreadySplit("TCC")

        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["ACC",[],[(2,True)],0]
        G_correct.contigs[2] = ["CCA",[(1,True),(3,True)],[],0]
        G_correct.contigs[3] = ["TCC",[],[(2,True)],0]
        G_correct.addKmersFromAllContigs()
        G_correct.setID(4)
        self.assertTrue(isSameGraph(G,G_correct))

    def test_2(self):
        #create G and G_correct:
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],10]
        #twin:TGGT
        G.contigs[G.getID()] = ["TTT",[(1,True)],[(1,True)],10]

        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["ACC",[],[(2,True)],0]
        G_correct.contigs[2] = ["CCA",[(1,True),(3,True)],[],0]
        G_correct.contigs[3] = ["TCC",[],[(2,True)],0]
        G_correct.contigs[4] = ["TTT",[(4,True)],[(4,True)],0]
        G.addKmersFromAllContigs()
        G_correct.setID(5)
        G_correct.addKmersFromAllContigs()
        G.addSegmentAlreadySplit("TCC")
        self.assertTrue(isSameGraph(G,G_correct))

    def test_3(self):
        #create G and G_correct:
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],10]
        #twin:TGGT
        G.contigs[G.getID()] = ["TTT",[(1,True),(2,True)],[(1,True),(2,True)],10]
        G.contigs[G.getID()] = ["TTATT",[(1,True),(2,True)],[(1,True),(2,True)],10]

        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["ACC",[],[(2,True)],0]
        G_correct.contigs[2] = ["CCA",[(1,True),(3,True)],[],0]
        G_correct.contigs[3] = ["TCC",[],[(2,True)],0]
        G_correct.contigs[4] = ["TTT",[(4,True),(5,True)],[(4,True),(5,True)],0]
        G_correct.contigs[5] = ["TTATT",[(4,True),(5,True)],[(4,True),(5,True)],0]
        G_correct.setID(6)
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.addSegmentAlreadySplit("TCC")
        G_naive = G.createNaive()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False))

    def test_4(self):
        G = Graph.Graph(3)
        G.addSegmentAlreadySplit("AAA")
        G_correct = Graph.Graph(3)
        G_correct.contigs[0] = ["AAA",[(0,True)],[(0,True)],10]
        G_correct.addKmersFromAllContigs()
        G_correct.setID(1)
        self.assertTrue(isSameGraph(G,G_correct))
"""

class Test_splitOthers(unittest.TestCase):
    #1 contig already in the graph
    #segment -> c0[1:]
    def test_1(self):
        k, G, G_correct = myHelper(5)
        G.contigs[G.getID()] = ["ACCATT",[],[],4]
        G.addKmersFromAllContigs()
        G.splitOthers("TCCAT")
        G_correct = Graph.Graph(k)
        G_correct.contigs[G_correct.getID()] = ["ACCAT",[],[(1,True)],2]
        G_correct.contigs[G_correct.getID()] = ["CCATT",[(0,True)],[],2]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    #1 contig already in the graph
    #segment -> c0[0:] (making sure we don't try to split)
    def test_2(self):
        k, G, G_correct = myHelper(5)
        G.contigs[G.getID()] = ["CCAAC",[],[],2]
        G.addKmersFromAllContigs()
        G.splitOthers("TCCAA")
        G_correct.contigs[G_correct.getID()] = ["CCAAC",[],[],2]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    #1 contig already in the graph
    #c0[L-k-1:] -> segment
    def test_3(self):
        k, G, G_correct = myHelper(5)
        G.contigs[G.getID()] = ["AAACCA",[],[],2]
        G.addKmersFromAllContigs()
        G.splitOthers("AACCT")
        G_correct.contigs[G_correct.getID()] = ["AAACC",[],[(1,True)],1]
        G_correct.contigs[G_correct.getID()] = ["AACCA",[(0,True)],[],1]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    #1 contig already in the graph
    #c0[L-k:] -> segment (making sure we don't try to split)
    def test_4(self):
        k, G, G_correct = myHelper(5)
        G.contigs[G.getID()] = ["ACGTAC",[(0,False),(0,False)],[(0,False),(0,False)],2]
        #                        GTACGT
        G.addKmersFromAllContigs()
        G.splitOthers("GTAAC")
        G_correct.contigs[G_correct.getID()] = ["ACGTAC",[(0,False),(0,False)],[(0,False),(0,False)],2]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    #1 contig already in the graph
    #twin(segment) -> c0[1:]
    """
    def test_5(self): Held að þetta eigi ekki að splitta
        k, G, G_correct = myHelper(5)
        G.contigs[G.getID()] = ["AAACCA",[],[],2]
        #                        TGGTTT
        G.addKmersFromAllContigs()
        G.splitOthers("GGTTC")
        #        twin: GAACC
        G.printContigs()
        G_correct = Graph.Graph(k)
        G_correct.contigs[G_correct.getID()] = ["AAACC",[],[(1,True)],2]
        G_correct.contigs[G_correct.getID()] = ["AACCA",[(0,True)],[],2]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))
    """
    #1 contig already in the graph
    #twin(segment) -> c0[0:]    (making sure we don't try to split)
    def test_6(self):
        k, G, G_correct = myHelper(5,assertLegal=False)
        G.contigs[G.getID()] = ["AAACCA",[],[],4]
        #                        TGGTTT
        G.addKmersFromAllContigs()
        G.splitOthers("GTTTC")
        #twin:         GAAAC
        G_correct.contigs[G_correct.getID()] = ["AAACCA",[],[],4]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))
        self.assertTrue(isSameGraph(G,G_correct))
        #G1,G2,alsoCompareWithNaive=True,pfn=False,ps=False,relaxAssertions=False

    #1 contig already in the graph
    #c0[L-k-1:] -> twin(segment)
    def test_7(self):
        k, G, G_correct = myHelper(5)
        G.contigs[G.getID()] = ["AAACCA",[],[],4]
        #                        TGGTTT
        G.addKmersFromAllContigs()
        G.splitOthers("CGGTT")
        #twin:         AACCG
        #G_correct.contigs[G_correct.getID()] = ["AAACC",[],[(1,True)],2]
        #G_correct.contigs[G_correct.getID()] = ["AACCA",[(0,True)],[],2]
        #G_correct.addKmersFromAllContigs()
        #self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    #1 contig already in the graph
    #c0[L-k:] -> twin(segment)
    def test_8(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],0]
        G.addKmersFromAllContigs()
        G.splitOthers("GTA")
        #twin: TAC
        G_correct = Graph.Graph(3)
        G_correct.contigs[0] = ["ACCA",[],[],0]
        G_correct.addKmersFromAllContigs()
        G_correct.setID(1)
        self.assertTrue(isSameGraph(G,G_correct))

    #2 contigs already in the graph
    #segment -> c0[1:]
    def test_9(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],2]
        G.contigs[G.getID()] = ["GTG",[],[],2]
        G.addKmersFromAllContigs()
        G.splitOthers("TCC")   
        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["GTG",[],[],1]
        G_correct.contigs[2] = ["ACC",[],[(3,True)],2]
        G_correct.contigs[3] = ["CCA",[(2,True)],[],2]
        G_correct.addKmersFromAllContigs()
        G_correct.setID(4)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    def test_10(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["GACCAA",[],[],2]
        G.addKmersFromAllContigs()
        G.splitOthers("TCC")    #twin: GGA
        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["GACC",[],[(2,True)],1]
        G_correct.contigs[2] = ["CCAA",[(1,True)],[],1]
        G_correct.addKmersFromAllContigs()
        G_correct.setID(3)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    def test_11(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CACCAA",[],[],2]
        G.addKmersFromAllContigs()
        G.splitOthers("GGA")
        G_correct = Graph.Graph(3)
        G_correct.contigs[1] = ["CACC",[],[(2,True)],1]
        G_correct.contigs[2] = ["CCAA",[(1,True)],[],1]
        G_correct.addKmersFromAllContigs()
        G_correct.setID(3)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=True))


    #new tests with a higher k. I'm going to create these new tests as if 
    #the earlier tests didn't exist

    #fw(segment) in c
    def test_12(self):
        G = Graph.Graph(k=15,pfn=False,ps=False,al=False)
        G_correct = Graph.Graph(15)
        G.contigs[G.getID()] = ["ATGGCCTGTAGTATCGCTATGGTA",[],[],20]
        #                             CTGTAGTATCGCTATGGTA
        #                        ATGGCCTGTAGTATCGCTA
        #                            ACTGTAGTATCGCTA
        G.addKmersFromAllContigs()
        G.splitOthers("ACTGTAGTATCGCTA")
        G_correct.contigs[G_correct.getID()] = ["ATGGCCTGTAGTATCGCTA",[],[(1,True)],10]
        G_correct.contigs[G_correct.getID()] = ["CTGTAGTATCGCTATGGTA",[(0,True)],[],10]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=True))
        G_correct.contigs[4] = ["CTGTAGTATCGCTATGGTA",[(3,True),(5,True)],[],2]
        G_correct.contigs[G_correct.getID()] = ["ACTGTAGTATCGCTA",[],[(4,True)],2]
        G.addSegmentToGraph("ACTGTAGTATCGCTA")
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=True,pfn=False,ps=False,relaxAssertions=False))

    #fw(segment) in c (the last kmer)
    def test_13(self):
        G = Graph.Graph(15,pfn=False,ps=False,al=False)
        G_correct = Graph.Graph(15)
        G.contigs[G.getID()] = ["ATGGCCTGTAGTATCGCTATGGTA",[],[],10]
        #                                 AGTATCGCTATGGTA
        #                        ATGGCCTGTAGTATCGCTATGGT
        #                                GAGTATCGCTATGGT
        G.addKmersFromAllContigs()
        G.splitOthers("GAGTATCGCTATGGT")
        G_correct.contigs[G_correct.getID()] = ["ATGGCCTGTAGTATCGCTATGGT",[],[(1,True)],5]
        G_correct.contigs[G_correct.getID()] = ["AGTATCGCTATGGTA",[(0,True)],[],5]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=True))
        G_correct.contigs[4] = ["AGTATCGCTATGGTA",[(3,True),(5,True)],[],5]
        G_correct.contigs[G_correct.getID()] = ["GAGTATCGCTATGGT",[],[(4,True)],5]
        G_correct.addKmersFromAllContigs()
        G.addSegmentToGraph("GAGTATCGCTATGGT")
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=True,pfn=False,ps=False,relaxAssertions=False))


    #fw(segment) in c (the second kmer)
    def test_14(self):
        G = Graph.Graph(15)
        G_correct = Graph.Graph(15)
        G.contigs[G.getID()] = ["ATGGCCTGTAGTATCGCTATGGTA",[],[],10]
        #                        ATGGCCTGTAGTATC    
        #                         TGGCCTGTAGTATCGCTATGGTA
        #                        CTGGCCTGTAGTATC
        G.addKmersFromAllContigs()
        G.splitOthers("CTGGCCTGTAGTATC")
        G_correct.contigs[G_correct.getID()] = ["ATGGCCTGTAGTATC",[],[(1,True)],10]
        G_correct.contigs[G_correct.getID()] = ["TGGCCTGTAGTATCGCTATGGTA",[(0,True)],[],10]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=True))

    #bw(segment) in c
    def test_15(self):
        G = Graph.Graph(15)
        G_correct = Graph.Graph(15)
        G.contigs[G.getID()] = ["ATGGCCTGTAGTATCGCTATGGTA",[],[],10]
        #                        ATGGCCTGTAGTATCGCTA
        #                             CTGTAGTATCGCTATGGTA
        #                             CTGTAGTATCGCTAA
        G.addKmersFromAllContigs()
        G.splitOthers("CTGTAGTATCGCTAA")
        G_correct.contigs[G_correct.getID()] = ["ATGGCCTGTAGTATCGCTA",[],[(1,True)],10]
        G_correct.contigs[G_correct.getID()] = ["CTGTAGTATCGCTATGGTA",[(0,True)],[],10]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=True))

    #bw(segment) in c (the first kmer)
    def test_16(self):
        G = Graph.Graph(15)
        G_correct = Graph.Graph(15)
        G.contigs[G.getID()] = ["ATGGCCTGTAGTATCGCTATGGTA",[],[],10]
        #                        ATGGCCTGTAGTATC
        #                         TGGCCTGTAGTATCGCTATGGTA
        #                         TGGCCTGTAGTATCC
        G.addKmersFromAllContigs()
        G.splitOthers("TGGCCTGTAGTATCC")
        G_correct.contigs[G_correct.getID()] = ["ATGGCCTGTAGTATC",[],[(1,True)],10]
        G_correct.contigs[G_correct.getID()] = ["TGGCCTGTAGTATCGCTATGGTA",[(0,True)],[],10]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=True))

    #bw(segment) in c (the second to last kmer)
    def test_17(self):
        G = Graph.Graph(15)
        G_correct = Graph.Graph(15)
        G.contigs[G.getID()] = ["ATGGCCTGTAGTATCGCTATGGTA",[],[],10]
        #                        ATGGCCTGTAGTATCGCTATGGT
        #                                 AGTATCGCTATGGTA
        #                                 AGTATCGCTATGGTTCAAACCTCTCAAA
        G.addKmersFromAllContigs()
        G.splitOthers("AGTATCGCTATGGTTCAAACCTCTCAAA")
        G_correct.contigs[G_correct.getID()] = ["ATGGCCTGTAGTATCGCTATGGT",[],[(1,True)],10]
        G_correct.contigs[G_correct.getID()] = ["AGTATCGCTATGGTA",[(0,True)],[],10]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=True))

    #fw(segment) in twin(c)
    def test_18(self):
        G = Graph.Graph(15)
        G_correct = Graph.Graph(15)
        G.contigs[G.getID()] = [dbg.twin("ATGGCCTGTAGTATCGCTATGGTA"),[],[],10]
        #                 twin: ATGGCCTGTAGTATCGCTATGGT  
        #                             TGTAGTATCGCTATG   this km has index 6 in twin => index in c: L-k-6=24-15-6=3
        #               adding:      TTGTAGTATCGCTAT

        #                    c: TACCATAGCGATACTACAGGCCAT
        #                                         SI=3+15=18
        #                       TACCATAGCGATACTACA
        #                           ATAGCGATACTACAGGCCAT
        G.addKmersFromAllContigs()
        G.splitOthers("TTGTAGTATCGCTAT")
        G_correct.contigs[G_correct.getID()] = ["TACCATAGCGATACTACA",[],[(1,True)],10]
        G_correct.contigs[G_correct.getID()] = ["ATAGCGATACTACAGGCCAT",[(0,True)],[],10]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=True))

    #bw(segment) in twin(c)
    def test_19(self):
        G = Graph.Graph(15)
        G_correct = Graph.Graph(15)
        G.contigs[G.getID()] = [dbg.twin("ATGGCCTGTAGTATCGCTATGGTA"),[],[],10]
        #                 twin: ATGGCCTGTAGTATCGCTATGGT  
        #                           CCTGTAGTATCGCTA   this km has index 4 in twin => index in c: L-k-4=24-15-4=5
        #               adding:      CTGTAGTATCGCTAA

        #                    c: TACCATAGCGATACTACAGGCCAT
        #                       TACCATAGCGATACTACAG
        #                            TAGCGATACTACAGGCCAT
        G.addKmersFromAllContigs()
        G.splitOthers("CTGTAGTATCGCTAA")
        G_correct.contigs[G_correct.getID()] = ["TACCATAGCGATACTACAG",[],[(1,True)],10]
        G_correct.contigs[G_correct.getID()] = ["TAGCGATACTACAGGCCAT",[(0,True)],[],10]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,pfn=False,ps=False,relaxAssertions=True))

class Test_splitContig(unittest.TestCase):
    def test_1(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["AACC",[],[],10]
        #                        GGTT
        G_correct.contigs[1] = ["AAC", [], [(2,True)], 10]
        G_correct.contigs[2] = ["ACC", [(1,True)], [], 10]
        G_correct.setID(3)
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.splitContig(0,3)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    def test_2(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["AACCC",[],[],10]
        #                        GGTT
        G_correct.contigs[1] = ["AAC", [], [(2,True)], 10]
        G_correct.contigs[2] = ["ACCC", [(1,True)], [], 10]
        G.addKmersFromAllContigs()
        G_correct.setID(3)
        G_correct.addKmersFromAllContigs()
        G.splitContig(0,3)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    def test_3(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["AACCC",[],[],10]
        #                        GGTT
        G_correct.contigs[1] = ["AACC", [], [(2,True)], 10]
        G_correct.contigs[2] = ["CCC", [(1,True)], [], 10]
        G.addKmersFromAllContigs()
        G_correct.setID(3)
        G_correct.addKmersFromAllContigs()
        G.splitContig(0,4)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    #Test sett upp til að ath orsök villu í addSegmentWithNoSeenKmers
    def test_4(self):
        G = Graph.Graph(15)
        G_correct = Graph.Graph(15)
        G.contigs[G.getID()] = ["ACCAGTTGATATGAC",[],[],10] #twin:TGGT
        G.contigs[G.getID()] = ["TCTATTCCAGGGTAGG",[],[(1,False),(1,False)],10] #twin:TAGA
        G.addKmersFromAllContigs()
        G.splitContig(1,15)
        G_correct.contigs[0] = ["ACCAGTTGATATGAC",[],[],10]
        G_correct.contigs[1] = ["CTATTCCAGGGTAGG",[(2,True)],[(1,False),(1,False)],10] #twin:TAG
        G_correct.contigs[2] = ["TCTATTCCAGGGTAG",[],[(1,True)],10]
        G_correct.setID(3)
        G_correct.addKmersFromAllContigs()
        G_correct.setID(3)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    def test_5(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACCA",[],[],10] #twin:TGGT
        G.contigs[G.getID()] = ["TCTA",[],[(1,False),(1,False),(2,True)],10] #twin:TAGA
        G.contigs[G.getID()] = ["TAT",[(1,True),(2,False),(2,False)],[(2,False),(2,False)],10] #twin: ATA
        G.addKmersFromAllContigs()
        G.splitContig(1,3)
        G_correct = Graph.Graph(3)
        G_correct.contigs[0] = ["ACCA",[],[],10]
        G_correct.contigs[2] = ["TAT",[(3,True),(2,False),(2,False)],[(2,False),(2,False)],10]
        G_correct.contigs[3] = ["CTA",[(4,True)],[(2,True),(3,False),(3,False)],10] #twin:TAG
        G_correct.contigs[4] = ["TCT",[],[(3,True)],10]
        G_correct.setID(5)
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    def test_6(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["TTT",[(0,True),(1,True)],[(0,True),(1,True)],10]
        G.contigs[G.getID()] = ["TTCTT",[(0,True),(1,True)],[(0,True),(1,True)],10]
        G_correct.contigs[0] = ["TTT",[(0,True),(3,True)],[(0,True),(2,True)],10]
        G_correct.contigs[2] = ["TTCT",[(0,True),(3,True)],[(3,True)],10]
        G_correct.contigs[3] = ["CTT",[(2,True)],[(0,True),(2,True)],10]
        G_correct.setID(4)
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.splitContig(1,4)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    #New tests to thoroughly test this function
    def myAssertions(self,G,G_correct):
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))
        self.assertTrue(G.isLegalGraph())
        self.assertTrue(G_correct.isLegalGraph())
        self.assertFalse(G.isLegalDBG())
        self.assertFalse(G_correct.isLegalDBG())


    def test_7(self):
        #c->c
        G = Graph.Graph(9,pil=False)
        G_correct = Graph.Graph(k=9)
        G.contigs[G.getID()] = ["TTTTTTTTTCCCTTTTTTTT",[(0,True)],[(0,True)],10]
        G.addKmersFromAllContigs()
        G.splitContig(0,12)
        G_correct.contigs[0] = ["TTTTTTTTTCCC",[(1,True)],[(1,True)],10]
        G_correct.contigs[1] = ["TTTTTCCCTTTTTTTT",[(0,True)],[(0,True)],10]
        G_correct.setID(2)
        G_correct.addKmersFromAllContigs()
        self.myAssertions(G,G_correct)

    def test_8(self):
        #c->twin(c)
        G = Graph.Graph(9)
        G_correct = Graph.Graph(9)
        G.contigs[G.getID()] = ["GGGTATATATA",[],[(0,False),(0,False)],10]
        G.addKmersFromAllContigs()
        G.splitContig(0,10)
        G_correct.contigs[0] = ["GGGTATATAT",[],[(1,True)],10]
        G_correct.contigs[1] = ["GTATATATA",[(0,True)],[(1,False),(1,False)],10]
        G_correct.setID(2)
        G_correct.addKmersFromAllContigs()
        self.myAssertions(G,G_correct)

    def test_9(self):
        #twin(c)->c
        G = Graph.Graph(9)
        G_correct = Graph.Graph(9)
        G.contigs[G.getID()] = ["TATATATAGGG",[(0,False),(0,False)],[],10]
        G.addKmersFromAllContigs()
        G.splitContig(0,10)
        G_correct.contigs[0] = ["TATATATAGG",[(0,False),(0,False)],[(1,True)],10]
        G_correct.contigs[1] = ["TATATAGGG",[(0,True)],[],10]
        G_correct.setID(2)
        G_correct.addKmersFromAllContigs()
        self.myAssertions(G,G_correct)

class Test_connectSegment(unittest.TestCase):
    #1 contig already in the graph. k=3
    def test_1(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        #The original Graph:
        G.contigs[G.getID()] = ["AAT",[],[(0,False),(0,False)],0]

        #The contig we wish to add:
        G.contigs[G.getID()] = ["ATC",[],[],0]

        G_correct.contigs[G_correct.getID()] = ["AAT",[],[(0,False),(0,False),(1,True)],0]
        G_correct.contigs[G_correct.getID()] = ["ATC",[(1,False),(1,False),(0,True)],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connectSegment(1,"ATC")
        self.assertTrue(isSameGraph(G,G_correct))

    #2 contigs already in the graph. k=3
    def test_2(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        #The original Graph:
        G.contigs[G.getID()] = ["AAT",[],[(0,False),(0,False)],0]
        G.contigs[G.getID()] = ["TCC",[],[],0]

        #The contig we wish to add:
        G.contigs[G.getID()] = ["ATC",[],[],0]

        G_correct.contigs[G_correct.getID()] = ["AAT",[],[(0,False),(0,False),(2,True)],0]
        G_correct.contigs[G_correct.getID()] = ["TCC",[(2,True)],[],0]
        G_correct.contigs[G_correct.getID()] = ["ATC",[(2,False),(2,False),(0,True)],[(1,True)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connectSegment(2,"ATC")
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False,relaxAssertions=True))

    #empty graph and segment connects to self
    def test_3(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["TTT",[],[],0]
        G_correct.contigs[G_correct.getID()] = ["TTT",[(0,True)],[(0,True)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.connectSegment(0)
        self.assertTrue(isSameGraph(G,G_correct))

class Test_mergeSegment(unittest.TestCase):
    #Two contigs a and b in the graph with no connections to themselves
    #The only connection is a connection between a and b
    
    #case 1: a->b   (we're imagining b was already in the graph and we're in the process of adding a)
    def test_1(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[],[(1,True)],0]
        G.contigs[G.getID()] = ["AGT",[(0,True)],[],0]
        G_correct.contigs[G_correct.getID()] = ["CAGT",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct))

    #case 2: a->twin(b)     (b already in graph and we're in process of adding a) 
    def test_2(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[],[(1,False)],0]
        G.contigs[G.getID()] = ["ACT",[],[(0,False)],0]
        G_correct.contigs[G_correct.getID()] = ["ACTG",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct))

     #case 3: twin(a)->b     (b already in graph and we're in process of adding a) 
    def test_3(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[(1,False)],[],0]
        G.contigs[G.getID()] = ["TGT",[(0,False)],[],0]
        G_correct.contigs[G_correct.getID()] = ["CTGT",[],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct))

    #-----------------------------------------------------------------------------------------

    #Two contigs a and b in the graph with a connection between them
    #AND b has a connection to itself (so we're not supposed to merge)

    #case 1: a->b   (b already in graph and we're in process of adding a) 
    def test_4(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[],[(1,True)],0]
        G.contigs[G.getID()] = ["AGT",[(0,True),(1,False),(1,False)],[],0]   #AGT doesnt actually connect to itself so we cant compare with naive
        G_correct.contigs[G_correct.getID()] = ["CAG",[],[(1,True)],0]
        G_correct.contigs[G_correct.getID()] = ["AGT",[(0,True),(1,False),(1,False)],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False))

    #case 2: a->twin(b)     (b already in graph and we're in process of adding a) 
    def test_5(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[],[(1,False)],0]
        G.contigs[G.getID()] = ["ACT",[],[(0,False),(1,False),(1,False)],0]
        G_correct.contigs[G_correct.getID()] = ["CAG",[],[(1,False)],0]
        G_correct.contigs[G_correct.getID()] = ["ACT",[],[(0,False),(1,False),(1,False)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False))

    #case 3: twin(a)->b     (b already in graph and we're in process of adding a) 
    def test_3b(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[(1,False)],[],0]
        G.contigs[G.getID()] = ["TGC",[(0,False)],[(1,False),(1,False)],0]
        G_correct.contigs[G_correct.getID()] = ["CTGC",[],[(0,False),(0,False)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct))

    #case 3: twin(a)->b     (b already in graph and we're in process of adding a) 
    def test_6(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[(1,False)],[],0]
        G.contigs[G.getID()] = ["TGC",[(0,False),(1,False),(1,False)],[],0]
        G_correct.contigs[G_correct.getID()] = ["CAG",[(1,False)],[],0]
        G_correct.contigs[G_correct.getID()] = ["TGC",[(0,False),(1,False),(1,False)],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False))

    #-----------------------------------------------------------------------------------------

    #Two contigs a and b in the graph with a connection between them
    #AND a has a connection to itself (so we're not supposed to merge)

    #case 1: a->b   (b already in graph and we're in process of adding a) 
    def test_7(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[(0,True)],[(0,True),(1,True)],0] #added a connection so not supposed to merge
        G.contigs[G.getID()] = ["AGT",[(0,True)],[],0]   
        G_correct.contigs[G_correct.getID()] = ["CAG",[(0,True)],[(0,True),(1,True)],0]
        G_correct.contigs[G_correct.getID()] = ["AGT",[(0,True)],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False))

    #case 2: a->twin(b)     (b already in graph and we're in process of adding a) 
    def test_8(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[(0,True)],[(0,True),(1,False)],0]
        G.contigs[G.getID()] = ["ACT",[],[(0,False)],0]
        G_correct.contigs[G_correct.getID()] = ["CAG",[(0,True)],[(0,True),(1,False)],0]
        G_correct.contigs[G_correct.getID()] = ["ACT",[],[(0,False)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False))

    #case 3: twin(a)->b     (b already in graph and we're in process of adding a) 
    def test_9(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CAG",[(0,True),(1,False)],[(0,True)],0]
        G.contigs[G.getID()] = ["TGC",[(0,False),(1,False),(1,False)],[],0]
        G_correct.contigs[G_correct.getID()] = ["CAG",[(0,True),(1,False)],[(0,True)],0]
        G_correct.contigs[G_correct.getID()] = ["TGC",[(0,False),(1,False),(1,False)],[],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        G.mergeSegment(0)
        self.assertTrue(isSameGraph(G,G_correct,alsoCompareWithNaive=False))

    #-----------------------------------------------------------------------------------------

    #Two contigs a and b in the graph with a connection between them
    #along with some more contigs

    #case 1: a->b.  a has some IN connections and b has some OUT connections
    def test_10(self):
        k=3
        G = Graph.Graph(k)
        #b: b_ID=0
        G.contigs[G.getID()] = ["AGT",[(99,True)],[(3,True),(4,True)],0]
        G.contigs[G.getID()] = ["TCA",[],[(99,True)],0] #1
        G.contigs[G.getID()] = ["ACA",[(4,False),(3,False)],[(99,True)],0] #2
        G.contigs[G.getID()] = ["GTT",[(0,True),(2,False)],[],0]  #3
        G.contigs[G.getID()] = ["GTA",[(0,True),(2,False)],[(4,False),(4,False)],0]  #4
        G.contigs[G.getID()] = ["TGG",[(99,False)],[],0] #5

        #TGG
        #CCA
        #a: a_ID=99
        G.setID(99)
        G.contigs[G.getID()] = ["CAG",[(1,True),(2,True),(5,False)],[(0,True)],0]
        G.addKmersFromAllContigs()
        self.assertTrue(G.isLegalGraph())
        G.mergeSegment(99)
        self.assertTrue(G.isLegalDBG())
        segments = ["AGT","TCA","ACA","GTT","GTA","TGG","CAG"]
        kd = helpers.createKmerDictFromSegmentList(segments,k)
        G_correct = helpers.createNaiveFromKmerDict(kd,k)

        self.assertTrue(G_correct.isLegalDBG())
        self.assertTrue(isSameGraph(G,G_correct))


    #case 1: a->b.  a and b have a connection between them and b has 1 OUT connection
    #               a has 1 additional OUT connection (2 with b) so we're not supposed to merge
    def test_11(self):
        k=3
        G = Graph.Graph(k)
        #b: b_ID=0
        G.contigs[G.getID()] = ["AGT",[(99,True)],[(2,True)],0]             #0
        G.contigs[G.getID()] = ["TCA",[],[(99,True)],0]                     #1
        G.contigs[G.getID()] = ["GTT",[(0,True),(4,True)],[],0]             #2
        G.contigs[G.getID()] = ["AGG",[(99,True)],[],0]                     #3
        G.contigs[G.getID()] = ["CGT",[(4,False),(4,False)],[(2,True)],0]   #4
        #TGG
        #CCA
        #a: a_ID=99
        G.setID(99)
        G.contigs[G.getID()] = ["CAG",[(1,True)],[(0,True),(3,True)],0]
        G.addKmersFromAllContigs()

        self.assertTrue(G.isLegalGraph())
        G.mergeSegment(99)
        self.assertTrue(G.isLegalGraph())
        G.mergeSegment(2)
        self.assertTrue(G.isLegalDBG())
        segments = ["AGT","TCA","GTT","AGG","CGT","CAG"]
        kd = helpers.createKmerDictFromSegmentList(segments,k)
        G_correct = helpers.createNaiveFromKmerDict(kd,k)

        self.assertTrue(G_correct.isLegalDBG())
        self.assertTrue(isSameGraph(G,G_correct))

    #when we merge a contig which connects to itself we don't do anything
    def test_12(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["TTT",[(0,True)],[(0,True)],0]
        G_correct.contigs[G_correct.getID()] = ["TTT",[(0,True)],[(0,True)],0]
        G.addKmersFromAllContigs()
        G_correct.addKmersFromAllContigs()
        self.assertTrue(G.isLegalDBG())
        G.mergeSegment(0)
        self.assertTrue(G.isLegalDBG())
        self.assertTrue(G_correct.isLegalDBG())
        self.assertTrue(isSameGraph(G,G_correct))

    #a->b
    #b->twin(b)
    #ab->twin(ab)   need to make sure that the ID in the connection from ab to itself is the ID
    #               of ab instead of the ID of b (which has been deleted)
    def test_13(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["CTA",[(1,True)],[(0,False),(0,False)],0]
        G.contigs[G.getID()] = ["TCT",[],[(0,True)],0]
        G.addKmersFromAllContigs()
        G.mergeSegment(1)
        G_correct.contigs[2] = ["TCTA",[],[(2,False),(2,False)],0]
        G_correct.addKmersFromAllContigs()
        G_correct.setID(3)
        self.assertTrue(G.isLegalDBG())
        self.assertTrue(G_correct.isLegalDBG())
        self.assertTrue(isSameGraph(G,G_correct))



    #-----------------------------------------------------------------------------------------

    #Three contigs a, b, and c in the graph with a connection between them
    #We imagine we're in the process of adding b so they should merge into 1 contig

    #case 1: a->b->c
    def test_14(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACC",[],[(1,True)],0]
        G.contigs[G.getID()] = ["CCT",[(0,True)],[(2,True)],0]
        G.contigs[G.getID()] = ["CTT",[(1,True)],[],0]
        G.addKmersFromAllContigs()
        G.mergeSegment(1)
        G_correct.contigs[G_correct.getID()] = ["ACCTT",[],[],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    #case 1: a->b->c. a has some INs
    def test_15(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACC",[(3,True),(4,False)],[(1,True)],0]
        G.contigs[G.getID()] = ["CCT",[(0,True)],[(2,True)],0]
        G.contigs[G.getID()] = ["CTT",[(1,True)],[],0]

        #Some additional INs and OUTs of ACC and CTT
        G.contigs[G.getID()] = ["CAC",[],[(0,True)],0]
        G.contigs[G.getID()] = ["GTC",[(0,False)],[],0]
        #G.contigs[G.getID()] = ["",[],[],0]
        #G.contigs[G.getID()] = ["",[],[],0]
        G.addKmersFromAllContigs()


        G.mergeSegment(1)
        G_correct.contigs[G_correct.getID()] = ["ACCTT",[(1,True),(2,False)],[],0]
        G_correct.contigs[G_correct.getID()] = ["CAC",[],[(0,True)],0]
        G_correct.contigs[G_correct.getID()] = ["GTC",[(0,False)],[],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

    #case 1: a->b->c. a has some INs and c has some OUTs
    def test_16(self):
        G = Graph.Graph(3)
        G_correct = Graph.Graph(3)
        G.contigs[G.getID()] = ["ACC",[(3,True),(4,False)],[(1,True)],0]    #0
        G.contigs[G.getID()] = ["CCT",[(0,True)],[(2,True)],0]              #1
        G.contigs[G.getID()] = ["CTT",[(1,True)],[(5,True),(6,False)],0]    #2

        #Some additional INs and OUTs of ACC and CTT
        G.contigs[G.getID()] = ["CAC",[],[(0,True)],0]                      #3
        G.contigs[G.getID()] = ["GTC",[(0,False)],[],0]                     #4
        G.contigs[G.getID()] = ["TTC",[(2,True)],[],0]                      #5
        G.contigs[G.getID()] = ["CAA",[],[(2,False)],0]                     #6
        G.addKmersFromAllContigs()

        G.mergeSegment(1)
        G_correct.contigs[G_correct.getID()] = ["ACCTT",[(1,True),(2,False)],[(3,True),(4,False)],0]
        G_correct.contigs[G_correct.getID()] = ["CAC",[],[(0,True)],0]
        G_correct.contigs[G_correct.getID()] = ["GTC",[(0,False)],[],0]
        G_correct.contigs[G_correct.getID()] = ["TTC",[(0,True)],[],0]
        G_correct.contigs[G_correct.getID()] = ["CAA",[],[(0,False)],0]
        G_correct.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G,G_correct))

if __name__ == '__main__':
    unittest.main()
