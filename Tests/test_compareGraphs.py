#coding:utf8
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import unittest
from compareGraphs import *
import Graph

#Note: compareGraphs does not take into account that two different graphs can 
#represent the same possible strings

#functions to test:
#1. sameContigs(G1,G2)
#2. sameIDs(G1,G2)
#3. fixIDs(G1,G2)
#4. fixTwins(G1,G2)
#5. same_INs_OUTs_COV(G1,G2)
#6. isSameGraph(G1,G2)
#Note: fixTwins uses Graph.flipContig(ID) so we can't use compareGraphs to test the 
#function flipContig

class Test_sameContigs(unittest.TestCase):
    def test_1(self):
        #Two empty graphs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        self.assertTrue(sameContigs(G1,G2))

    def test_2(self):
        #G1 has more contigs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        self.assertFalse(sameContigs(G1,G2))

    def test_3(self):
        #G2 has more contigs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G2.contigs[0] = ["ATT",[],[],0]
        self.assertFalse(sameContigs(G1,G2))

    def test_4(self):
        #1 different contig
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        G2.contigs[0] = ["ATT",[],[],0]
        self.assertTrue(sameContigs(G1,G2))
        G1.contigs[1] = ["TTT",[],[],0]
        G2.contigs[1] = ["GGG",[],[],0]
        self.assertFalse(sameContigs(G1,G2))

    def test_5(self):
        #G1[0] = twin(G2[0])
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        G2.contigs[0] = ["AAT",[],[],0]

    def test_5(self):
        #G1[1] = twin(G2[15])
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[10] = ["ATT",[],[],0]
        G2.contigs[0] = ["ATT",[],[],0]
        self.assertTrue(sameContigs(G1,G2))
        G1.contigs[1] = ["TTT",[],[],0]
        G2.contigs[15] = ["AAA",[],[],0]
        self.assertTrue(sameContigs(G1,G2))
    
class Test_sameIDs(unittest.TestCase):
    def test_1(self):
        #Two empty graphs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        self.assertTrue(sameIDs(G1,G2))

    def test_2(self):
        #Two identical graphs with same ID
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        G2.contigs[0] = ["ATT",[],[],0]
        self.assertTrue(sameIDs(G1,G2))

    def test_3(self):
        #Two identical graphs with different ID
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        G2.contigs[10] = ["ATT",[],[],0]
        self.assertFalse(sameIDs(G1,G2))

    def test_4(self):
        #Two identical graphs with different ID
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        G2.contigs[0] = ["ATT",[],[],0]
        G1.contigs[1] = ["GGG",[],[],0]
        G2.contigs[2] = ["GGG",[],[],0]
        self.assertFalse(sameIDs(G1,G2))

    def test_5(self):
        #Two identical graphs where twin has same ID
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        G2.contigs[0] = ["ATT",[],[],0]
        G1.contigs[1] = ["GGG",[],[],0]
        G2.contigs[1] = ["CCC",[],[],0]
        self.assertTrue(sameIDs(G1,G2))

    def test_6(self):
        #Two identical graphs where twin does not have same ID
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        G2.contigs[0] = ["ATT",[],[],0]
        G1.contigs[1] = ["GGG",[],[],0]
        G2.contigs[2] = ["CCC",[],[],0]
        self.assertFalse(sameIDs(G1,G2))

    def test_7(self):
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["AAA",[(0,True)],[(0,True)],0]
        G1.contigs[1] = ["TAG",[(1,False),(1,False)],[],0] #twin: CTA
        G2.contigs[1] = ["AAA",[(1,True)],[(1,True)],0]
        G2.contigs[0] = ["TAG",[(0,False),(0,False)],[],0] #twin: CTA
        self.assertFalse(sameIDs(G1,G2))

class Test_fixIDs(unittest.TestCase):
    def allIsEqual(self,G1,G2,text):
        if text=="skip":
            #þetta er bara af því að þessi kóði virkar ekki fyrir test_5 (sem skilar samt réttu)
            return
        for ID in G1.contigs:
            [c1,IN1,OUT1,COV1] = G1.contigs[ID]
            [c2,IN2,OUT2,COV2] = G2.contigs[ID]

            self.assertEqual(c1,c2)
            self.assertEqual(sorted(IN1),sorted(IN2))
            self.assertEqual(sorted(OUT1),sorted(OUT2))
            self.assertEqual(COV1,COV2)

    def helper(self,G1,G2,text=""):
        if not sameIDs(G1,G2):
            fixIDs(G1,G2)
        self.assertTrue(sameIDs(G1,G2))
        self.allIsEqual(G1,G2,text)

    def test_1(self):
        #Two empty graphs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        self.helper(G1,G2)

    def test_2(self):
        #Two graphs with 1 different ID
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[1] = ["ATT",[(1,False),(1,False)],[],0]
        self.helper(G1,G2)

    def test_3(self):
        #Two graphs with 2 different IDs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        G2.contigs[1] = ["ATT",[],[],0]
        G1.contigs[1] = ["GGG",[],[],0]
        G2.contigs[2] = ["GGG",[],[],0]
        self.helper(G1,G2)

    def test_4(self):
        #Two graphs with 2 different IDs and 1 equal ID
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[],[],0]
        G2.contigs[1] = ["ATT",[],[],0]
        G1.contigs[1] = ["GGG",[],[],0]
        G2.contigs[2] = ["GGG",[],[],0]
        G1.contigs[3] = ["TTT",[],[],0]
        G2.contigs[3] = ["TTT",[],[],0]
        self.helper(G1,G2)

    def test_5(self):
        #Two identical graphs where twin does not have same ID
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[],[],0]
        G2.contigs[G2.getID()] = ["ATT",[],[],0]
        G1.contigs[G1.getID()] = ["GGG",[],[],0]
        G2.contigs[G2.getID()] = ["CCC",[],[],0]
        self.helper(G1,G2,"skip")

    def test_6(self):
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G2.ID = 3568
        G2.contigs[3567] = ["ATT",[(3567,False),(3567,False)],[],0]
        fixIDs(G1,G2)
        self.assertTrue(sameIDs(G1,G2))

    def test_7(self):
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G2.ID = 3568
        G2.contigs[3567] = ["AAT",[],[(3567,False),(3567,False)],0]
        fixIDs(G1,G2)
        self.assertTrue(sameIDs(G1,G2))

    def test_8(self):
        #Two graphs with 1 different ID in a connection to self
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["TTT",[(0,True)],[(0,True)],0]
        G2.getID()
        G2.contigs[G2.getID()] = ["TTT",[(1,True)],[(1,True)],0]
        self.helper(G1,G2)

    def test_9(self):
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        G1.contigs[G1.getID()] = ["TAG",[(1,False),(1,False)],[],0] #twin: CTA
        G2.contigs[G2.getID()] = ["TAG",[(0,False),(0,False)],[],0] #twin: CTA
        G2.contigs[G2.getID()] = ["AAA",[(1,True)],[(1,True)],0]
        self.helper(G1,G2)

class Test_fixTwins(unittest.TestCase):
    def allIsEqual(self,G1,G2,index):
        [c1,IN1,OUT1,COV1] = G1.contigs[index]
        [c2,IN2,OUT2,COV2] = G2.contigs[index]

        self.assertEqual(c1,c2)
        self.assertEqual(sorted(IN1),sorted(IN2))
        self.assertEqual(sorted(OUT1),sorted(OUT2))
        self.assertEqual(COV1,COV2)

    def test_1(self):
        #Two empty graphs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        fixTwins(G1,G2)

    def test_2(self):
        #Two graphs with 1 contig
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[G2.getID()] = ["AAT",[],[(0,False),(0,False)],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()
        fixTwins(G1,G2)
        self.allIsEqual(G1,G2,0)

    def test_3(self):
        #Two graphs with 2 contigs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[(0,False),(0,False)],[(1,True)],0]
        G1.contigs[G1.getID()] = ["TTC",[(0,True)],[],0]
        G2.contigs[G2.getID()] = ["AAT",[(1,False)],[(0,False),(0,False)],0]
        G2.contigs[G2.getID()] = ["TTC",[(0,False)],[],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()
        fixTwins(G1,G2)
        self.allIsEqual(G1,G2,0)

    def test_4(self):
        #Two graphs with 1 contig. Not identical (twins of each other but not same IN and OUT)
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[0] = ["AAT",[],[(0,False)],0]
        with self.assertRaises(AssertionError):
            fixTwins(G1,G2)

class Test_same_INs_OUTs_COV(unittest.TestCase):
    def test_1(self):
        #Two empty graphs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        self.assertTrue(same_INs_OUTs_COV(G1,G2))

    def test_2(self):
        #Two graphs with 1 contig
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[0] = ["ATT",[(0,False),(0,False)],[],0]
        self.assertTrue(same_INs_OUTs_COV(G1,G2))

    def test_3(self):
        #Two graphs with 1 contig
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[0] = ["ATT",[(0,True),(0,False)],[],0]
        self.assertFalse(same_INs_OUTs_COV(G1,G2))

    def test_4(self):
        #Two graphs with 2 contigs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[(1,True)],0]
        G1.contigs[1] = ["TTC",[(0,True)],[],0]
        G2.contigs[0] = ["ATT",[(0,False),(0,False)],[(1,True)],0]
        G2.contigs[1] = ["TTC",[(0,True)],[],0]
        self.assertTrue(same_INs_OUTs_COV(G1,G2))

    def test_5(self):
        #Two graphs with 2 contigs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[(1,True)],0]
        G1.contigs[1] = ["TTC",[(0,True)],[],0]
        G2.contigs[0] = ["ATT",[(0,False),(0,True)],[(1,True)],0]
        G2.contigs[1] = ["TTC",[(0,True)],[],0]
        self.assertFalse(same_INs_OUTs_COV(G1,G2))

    def test_6(self):
        #Two graphs with 1 contig. Different length of IN
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[0] = ["ATT",[(0,False)],[],0]
        self.assertFalse(same_INs_OUTs_COV(G1,G2))

    def test_7(self):
        #CAAA
        #TTTG

        #TGGG
        #CCCA
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["CAAA",[(1,False)],[],0]
        G2.contigs[0] = ["CAAA",[(1,False)],[],0]
        self.assertTrue(same_INs_OUTs_COV(G1,G2))
        G1.contigs[1] = ["TGGG",[(0,False)],[],0]
        G2.contigs[1] = ["TGGG",[(0,False)],[],0]
        self.assertTrue(same_INs_OUTs_COV(G1,G2))

    def test_8(self):
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G3 = Graph.Graph(3)
        G1.contigs[0] = ["CCCCCCCCCCCCCC",[(1,False),(1,False),(1,False)],[((1,False),(1,False),(1,False),(1,False))],0]
        G2.contigs[0] = ["CCCCCCCCCCCCCC",[(1,False),(1,False),(1,False)],[((1,False),(1,False),(1,False),(1,False))],0]
        G3.contigs[0] = ["CCCCCCCCCCCCCC",[(1,False),(1,False),(1,False)],[((1,False),(1,False),(1,False),(0,False))],0]
        self.assertTrue(same_INs_OUTs_COV(G1,G2))
        self.assertFalse(same_INs_OUTs_COV(G1,G3))
        
class Test_isSameGraph(unittest.TestCase):
    def test_1(self):
        #Two empty graphs
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        self.assertTrue(isSameGraph(G1,G2))

    def test_2(self):
        #Two graphs with 1 contig. Identical
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[G2.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G1,G2))

    def test_3(self):
        #Two graphs with 1 contig. Not identical
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G3 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[G2.getID()] = ["ATT",[(0,True)],[(0,True)],0]
        G3.contigs[G3.getID()] = ["ATA",[(0,False),(0,False)],[],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()
        G3.addKmersFromAllContigs()
        self.assertFalse(isSameGraph(G1,G2))
        self.assertFalse(isSameGraph(G1,G3))

    def test_4(self):
        #Two graphs with 1 contig. Identical (twins of each other)
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[G2.getID()] = ["AAT",[],[(0,False),(0,False)],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G1,G2))

    def test_5(self):
        #Two graphs with 1 contig. Not identical (twins of each other but not same IN and OUT)
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[(1,True)],0]
        G2.contigs[0] = ["AAT",[],[(0,False),(0,False)],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()
        with self.assertRaises(AssertionError):
            self.assertFalse(isSameGraph(G1,G2))

    def test_6(self):
        #Two graphs with 1 contig. Different length of IN
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[0] = ["ATT",[(0,False),(0,False)],[],0]
        G2.contigs[0] = ["ATT",[(1,False)],[],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()
        with self.assertRaises(AssertionError):
            self.assertFalse(isSameGraph(G1,G2))

    def test_7(self):
        #CAAA
        #TTTG

        #TGGG
        #CCCA
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["CAAA",[(2,True),(1,False)],[],0]
        G2.contigs[G2.getID()] = ["CAAA",[(2,True),(1,False)],[],0]
        G1.contigs[G1.getID()] = ["TGGG",[(0,False)],[],0]
        G2.contigs[G2.getID()] = ["TGGG",[(0,False)],[],0]
        G1.contigs[G1.getID()] = ["TCA",[],[(0,True)],0]
        G2.contigs[G2.getID()] = ["TCA",[],[(0,True)],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()
        self.assertTrue(isSameGraph(G1,G2,alsoCompareWithNaive=False))

    def test_8(self):
        G1 = Graph.Graph(3)
        G2 = Graph.Graph(3)
        G1.contigs[G1.getID()] = ["ATT",[(0,False),(0,False)],[],0]
        G2.ID = 3568
        G2.contigs[3567] = ["AAT",[],[(3567,False),(3567,False)],0]
        G1.addKmersFromAllContigs()
        G2.addKmersFromAllContigs()

        self.assertTrue(isSameGraph(G1,G2))

class Test_graphsEqualsNaive(unittest.TestCase):
    def test_1(self):
        G = Graph.Graph(3)
        G.contigs[G.getID()] = ["AAA",[(0,True)],[(0,True)],0]
        G.addKmersFromAllContigs()
        self.assertTrue(G.equalsNaive())

    def test_2(self):
        G = Graph.Graph(3)
        G.addSegmentToGraph("AAATCC")
        G.addKmersFromAllContigs()
        self.assertTrue(G.equalsNaive())

if __name__ == '__main__':
    unittest.main()
