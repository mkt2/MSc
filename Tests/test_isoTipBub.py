#coding:utf8
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import unittest
import Graph
import Path

#--------------------------------------------------------------------------
#------------------Functions for exploring missing kmers-------------------
#--------------------------------------------------------------------------
class Test_createAllPathsIncludingContig(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTTT",[(0,True)],[],0]      #3 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        p0 = G.createAllPathsIncludingContig(MPL, 0)
        p1 = G.createAllPathsIncludingContig(MPL, 1)
        self.assertTrue(p0.isEmpty())
        self.assertTrue(p1.isEmpty())

    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTT",[(0,True)],[],0]       #2 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        P0 = G.createAllPathsIncludingContig(MPL, 0)
        P1 = G.createAllPathsIncludingContig(MPL, 1)
        self.assertTrue(P0.equals(P1))
        p = Path.Path()
        p.append(0,True,"CAGGTATCCAT",k)
        p.append(1,True,"CCATTT",k)
        correctPathList = Path.PathList()
        correctPathList.append(p)
        self.assertTrue(P0.equals(correctPathList))
        self.assertTrue(P1.equals(correctPathList))
        self.assertTrue(correctPathList.equals(P0))
        self.assertTrue(correctPathList.equals(P1))

    def test_3(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTT",[(0,True)],[],0]       #2 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        P0_1 = G.createAllPathsIncludingContig(MPL, 0, {1})
        P1_0 = G.createAllPathsIncludingContig(MPL, 1, {0})
        self.assertFalse(P0_1.equals(P1_0))
        p0_1 = Path.Path()
        p1_0 = Path.Path()
        p0_1.append(0,True,"CAGGTATCCAT",k)
        p1_0.append(1,True,"CCATTT",k)
        P0_1_correct = Path.PathList()
        P1_0_correct = Path.PathList()
        P0_1_correct.append(p0_1)
        P1_0_correct.append(p1_0)
        self.assertTrue(P0_1_correct.equals(P0_1))
        self.assertTrue(P1_0_correct.equals(P1_0))

    def test_4(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]                  #7 k-mers
        G.contigs[G.getID()] = [       "CCATTT",[(0,True)],[(2,True)],0]        #2 k-mers
        G.contigs[G.getID()] = [         "ATTTT",[(1,True)],[(3,True)],0]       #1 k-mer
        G.contigs[G.getID()] = [          "TTTTCCGA",[(2,True)],[],0]           #4 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        P0 = G.createAllPathsIncludingContig(MPL, 0)
        P1 = G.createAllPathsIncludingContig(MPL, 1)
        P2 = G.createAllPathsIncludingContig(MPL, 2)
        P3 = G.createAllPathsIncludingContig(MPL, 3)
        self.assertTrue(P0.equals(P1))
        self.assertTrue(P0.equals(P2))
        self.assertTrue(P0.equals(P3))
        self.assertTrue(P0.isEmpty())
        P0_1 = G.createAllPathsIncludingContig(MPL, 0, {1})
        P1_0 = G.createAllPathsIncludingContig(MPL, 1, {0})
        self.assertFalse(P0_1.equals(P1_0))
        p0_1 = Path.Path()
        p1_0 = Path.Path()
        p0_1.append(0,True,"CAGGTATCCAT",k)
        p1_0.append(1,True,"CCATTT",k)
        p1_0.append(2,True,"ATTTT",k)
        p1_0.append(3,True,"TTTTCCGA",k)
        P0_1_correct = Path.PathList()
        P1_0_correct = Path.PathList()
        P0_1_correct.append(p0_1)
        P1_0_correct.append(p1_0)
        self.assertTrue(P0_1_correct.equals(P0_1))
        self.assertTrue(P1_0_correct.equals(P1_0))
        P2_13 = G.createAllPathsIncludingContig(MPL, 2, {1,3})
        p2_13 = Path.Path()
        p2_13.append(2,True,"ATTTT",k)
        P2_13_correct = Path.PathList()
        P2_13_correct.append(p2_13)
        self.assertTrue(P2_13_correct.equals(P2_13))

    def test_5(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["AGGTCAAGGT",[],[(1,True),(2,True)],0]       #0. 6 k-mers
        G.contigs[G.getID()] = [      "AGGTC",[(0,True)],[],0]               #1. 1 k-mers
        G.contigs[G.getID()] = [      "AGGTATGG",[(0,True)],[(3,True)],0]    #2. 4 k-mers
        G.contigs[G.getID()] = [          "ATGGGATGCA",[(2,True)],[],0]      #3. 6 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        P0 = G.createAllPathsIncludingContig(MPL, 0)
        P0_1 = G.createAllPathsIncludingContig(MPL, 0, {1})
        self.assertTrue(P0.isEmpty())
        self.assertTrue(P0_1.isEmpty())
        P0_2 = G.createAllPathsIncludingContig(MPL, 0, {2})
        P0_2_correct = Path.PathList()
        p = Path.Path()
        p.append(0,True,G.contigs[0][0],k)
        p.append(1,True,G.contigs[1][0],k)    
        P0_2_correct.append(p)
        self.assertTrue(P0_2.equals(P0_2_correct))
        P3_2 = G.createAllPathsIncludingContig(MPL, 3, {2})
        P3_2_correct = Path.PathList()
        p = Path.Path()
        p.append(3,True,G.contigs[3][0],k)
        P3_2_correct.append(p)
        self.assertTrue(P3_2.equals(P3_2_correct))
        P2_3 = G.createAllPathsIncludingContig(MPL, 2, {3})
        self.assertTrue(P2_3.isEmpty())
        P2_0 = G.createAllPathsIncludingContig(MPL, 2, {0})
        self.assertTrue(P2_0.isEmpty())
        P2_03 = G.createAllPathsIncludingContig(MPL, 2, {0,3})
        P2_03_correct = Path.PathList()
        p = Path.Path()
        p.append(2,True,G.contigs[2][0],k)
        P2_03_correct.append(p)
        self.assertTrue(P2_03.equals(P2_03_correct))
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"

    #Same as test_5 except we only look at P2_3 (which was causing problems)
    def test_6(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["AGGTCAAGGT",[],[(1,True),(2,True)],0]       #0. 6 k-mers
        G.contigs[G.getID()] = [      "AGGTC",[(0,True)],[],0]               #1. 1 k-mers
        G.contigs[G.getID()] = [      "AGGTATGG",[(0,True)],[(3,True)],0]    #2. 4 k-mers
        G.contigs[G.getID()] = [          "ATGGGATGCA",[(2,True)],[],0]      #3. 6 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        P2_3 = G.createAllPathsIncludingContig(MPL, 2, {3})
        self.assertTrue(P2_3.isEmpty())
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"

    def test_7(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["AGGTCAAGGT",[],[(1,True),(2,True)],0]       #0. 6 k-mers
        G.contigs[G.getID()] = [      "AGGTC",[(0,True)],[],0]               #1. 1 k-mers
        G.contigs[G.getID()] = [      "AGGTATGG",[(0,True)],[],0]            #2. 4 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        P2 = G.createAllPathsIncludingContig(MPL, 2)
        self.assertTrue(P2.isEmpty())
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"
        
class Test_createAllPathsInContigsCollection(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTTT",[(0,True)],[],0]      #3 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        p0 = G.createAllPathsInContigsCollection(MPL, 0, set())
        p1 = G.createAllPathsInContigsCollection(MPL, 1, set())
        self.assertTrue(p0.isEmpty())
        self.assertTrue(p1.isEmpty())

    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTT",[(0,True)],[],0]       #2 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        P0 = G.createAllPathsInContigsCollection(MPL, 0, set())
        P1 = G.createAllPathsInContigsCollection(MPL, 1, set())
        self.assertTrue(P0.equals(P1))
        p = Path.Path()
        p.append(0,True,"CAGGTATCCAT",k)
        p.append(1,True,"CCATTT",k)
        correctPathList = Path.PathList()
        correctPathList.append(p)
        self.assertTrue(P0.equals(correctPathList))
        self.assertTrue(P1.equals(correctPathList))
        self.assertTrue(correctPathList.equals(P0))
        self.assertTrue(correctPathList.equals(P1))

    def test_3(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTT",[(0,True)],[],0]       #2 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        P0_1 = G.createAllPathsInContigsCollection(MPL, 0, {1})
        P1_0 = G.createAllPathsInContigsCollection(MPL, 1, {0})
        self.assertFalse(P0_1.equals(P1_0))
        p0_1 = Path.Path()
        p1_0 = Path.Path()
        p0_1.append(0,True,"CAGGTATCCAT",k)
        p1_0.append(1,True,"CCATTT",k)
        P0_1_correct = Path.PathList()
        P1_0_correct = Path.PathList()
        P0_1_correct.append(p0_1)
        P1_0_correct.append(p1_0)
        self.assertTrue(P0_1_correct.equals(P0_1))
        self.assertTrue(P1_0_correct.equals(P1_0))

    def test_4(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["AGGTCAAGGT",[],[(1,True),(2,True)],0]       #0. 6 k-mers
        G.contigs[G.getID()] = [      "AGGTC",[(0,True)],[],0]               #1. 1 k-mers
        G.contigs[G.getID()] = [      "AGGTATGG",[(0,True)],[(3,True)],0]    #2. 4 k-mers
        G.contigs[G.getID()] = [          "ATGGGATGCA",[(2,True)],[],0]      #3. 6 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        P0 = G.createAllPathsInContigsCollection(MPL, 0)
        P0_1 = G.createAllPathsInContigsCollection(MPL, 0, {1})
        self.assertTrue(P0.isEmpty())
        self.assertTrue(P0_1.isEmpty())
        P0_2 = G.createAllPathsInContigsCollection(MPL, 0, {2})
        P0_2_correct = Path.PathList()
        p = Path.Path()
        p.append(0,True,G.contigs[0][0],k)
        p.append(1,True,G.contigs[1][0],k)    
        P0_2_correct.append(p)
        self.assertTrue(P0_2.equals(P0_2_correct))
        P3_2 = G.createAllPathsInContigsCollection(MPL, 3, {2})
        P3_2_correct = Path.PathList()
        p = Path.Path()
        p.append(3,True,G.contigs[3][0],k)
        P3_2_correct.append(p)
        self.assertTrue(P3_2.equals(P3_2_correct))
        P2_3 = G.createAllPathsInContigsCollection(MPL, 2, {3})
        self.assertTrue(P2_3.isEmpty())
        P2_0 = G.createAllPathsInContigsCollection(MPL, 2, {0})
        self.assertTrue(P2_0.isEmpty())
        P2_03 = G.createAllPathsInContigsCollection(MPL, 2, {0,3})
        P2_03_correct = Path.PathList()
        p = Path.Path()
        p.append(2,True,G.contigs[2][0],k)
        P2_03_correct.append(p)
        self.assertTrue(P2_03.equals(P2_03_correct))
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"

class Test_isPartOfIsolated(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #Make sure that two contigs which together are long are not classified as
    #isolated
    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTTT",[(0,True)],[],0]      #3 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B0 = G.isPartOfIsolated(MPL,0)
        B1 = G.isPartOfIsolated(MPL,1)
        self.assertEqual(B0 , False)
        self.assertEqual(B1 , False)

    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTT",[(0,True)],[],0]       #2 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B0 = G.isPartOfIsolated(MPL,0)
        B1 = G.isPartOfIsolated(MPL,1)
        self.assertEqual(B0 , True)
        self.assertEqual(B1 , True)

    #def isPartOfIsolated(self, MPL, c_ID, skipNodes={})
    def test_3(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = [       "CCATTT",[(0,True)],[(2,True)],0]       #2 k-mers
        G.contigs[G.getID()] = [         "ATTTT",[(1,True)],[(3,True)],0]       #1 k-mer
        G.contigs[G.getID()] = [          "TTTTCCGA",[(2,True)],[],0]       #4 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B0 = G.isPartOfIsolated(MPL,0)
        B1 = G.isPartOfIsolated(MPL,1)
        B2 = G.isPartOfIsolated(MPL,2)
        B3 = G.isPartOfIsolated(MPL,3)
        self.assertEqual(B0 , False)
        self.assertEqual(B1 , False)
        self.assertEqual(B2 , False)
        self.assertEqual(B3 , False)
        B0_1 = G.isPartOfIsolated(MPL,0,{1})
        self.assertEqual(B0_1,True)
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"

    def test_4(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]          #7 k-mers
        G.contigs[G.getID()] = [       "CCATTTT",[(0,True)],[],0]       #3 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B0_1 = G.isPartOfIsolated(MPL,0,{1})
        self.assertEqual(B0_1,True)
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"

    def test_5(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True),(3,True),(4,True)],0]     #0. 1 k-mers
        G.contigs[G.getID()] = [ "AGGTT",[(0,True)],[],0]                               #1. 1 k-mers
        G.contigs[G.getID()] = [ "AGGTA",[(0,True)],[],0]                               #2. 1 k-mers
        G.contigs[G.getID()] = [ "AGGTC",[(0,True)],[],0]                               #3. 1 k-mers
        G.contigs[G.getID()] = [ "AGGTG",[(0,True)],[(5,True)],0]                       #4. 1 k-mers
        G.contigs[G.getID()] = [  "GGTGGATGCGTT",[(4,True)],[],0]                       #5. 8 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B0 = G.isPartOfIsolated(MPL,0)
        B0_4 = G.isPartOfIsolated(MPL,0,{4})
        self.assertFalse(B0)
        self.assertTrue(B0_4)
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"

    def test_6(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True),(3,True),(4,True)],0]     #0. 1 k-mers
        G.contigs[G.getID()] = [ "AGGTT",[(0,True)],[],0]                               #1. 1 k-mers
        G.contigs[G.getID()] = [ "AGGTA",[(0,True)],[],0]                               #2. 1 k-mers
        G.contigs[G.getID()] = [ "AGGTC",[(0,True)],[],0]                               #3. 1 k-mers
        G.contigs[G.getID()] = [ "AGGTG",[(0,True)],[(5,True)],0]                       #4. 1 k-mers
        G.contigs[G.getID()] = [  "GGTGGATGCGTTCC",[(4,True)],[],0]                     #5. 10 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B0 = G.isPartOfIsolated(MPL,0)
        B0_4 = G.isPartOfIsolated(MPL,0,{4})
        B5_4 = G.isPartOfIsolated(MPL,5,{4})
        self.assertFalse(B0)
        self.assertTrue(B0_4)
        self.assertFalse(B5_4)
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"

    def test_7(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["AGGTCAAGGT",[],[(1,True),(2,True)],0]       #0. 6 k-mers
        G.contigs[G.getID()] = [      "AGGTC",[(0,True)],[],0]               #1. 1 k-mers
        G.contigs[G.getID()] = [      "AGGTATGG",[(0,True)],[(3,True)],0]    #2. 4 k-mers
        G.contigs[G.getID()] = [          "ATGGGATGCA",[(2,True)],[],0]      #3. 6 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertFalse(G.isPartOfIsolated(MPL,0))
        self.assertFalse(G.isPartOfIsolated(MPL,0,{1}))
        self.assertTrue(G.isPartOfIsolated(MPL,0,{2}))
        self.assertTrue(G.isPartOfIsolated(MPL,3,{2}))
        self.assertFalse(G.isPartOfIsolated(MPL,2,{3}))
        self.assertFalse(G.isPartOfIsolated(MPL,2,{0}))
        self.assertTrue(G.isPartOfIsolated(MPL,2,{0,3}))
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"

    def test_8(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["AGGTCAAGGT",[],[(1,True),(2,True)],0]       #0. 6 k-mers
        G.contigs[G.getID()] = [      "AGGTC",[(0,True)],[],0]               #1. 1 k-mers
        G.contigs[G.getID()] = [      "AGGTATGG",[(0,True)],[],0]            #2. 4 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertFalse(G.isPartOfIsolated(MPL,0))
        self.assertFalse(G.isPartOfIsolated(MPL,0,{1}))
        self.assertTrue(G.isPartOfIsolated(MPL,0,{2}))
        self.assertFalse(G.isPartOfIsolated(MPL,2))
        self.assertTrue(G.isPartOfIsolated(MPL,2,{0}))
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1), "We shouldn't have changed any markings"

class Test_markAllIsolated(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #Make sure short contigs not connected to anything are classified
    #as isolated
    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[],0]
        G.contigs[G.getID()] = ["AGTTA",[],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertEqual(G.degrees[0][2] , 1)
        self.assertEqual(G.degrees[1][2] , 1)

    #Make sure that a long (len > 2k) contig is not classified as isolated
    def test_3(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[],0]
        G.contigs[G.getID()] = ["AGTTA",[],[],0]
        G.contigs[G.getID()] = ["AGTTAAACCCGAGGTTAGACCCTAGGACT",[],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertEqual(G.degrees[0][2] , 1)
        self.assertEqual(G.degrees[1][2] , 1)
        self.assertEqual(G.degrees[2][2] , -1)

    #Make sure that two contigs which together are long are not classified as
    #isolated
    def test_4(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTTT",[(0,True)],[],0]      #3 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertEqual(G.degrees[0][2] , -1)
        self.assertEqual(G.degrees[1][2] , -1)

    #Make sure that two contigs which together are long are not classified as
    #isolated. Check off by 1 error
    def test_5(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTT",[(0,True)],[],0]       #2 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertEqual(G.degrees[0][2] , 1)
        self.assertEqual(G.degrees[1][2] , 1)

    #Make sure that a collection of connected contigs is still classified as isolated
    #so long as the longest path they're part of is < 2k. Even though the total length 
    #of the entire collection is greater than 2k
    def test_6(self):
        k, MPL, G = self.helper()
        #Longest path: 8 kmers
        #Total length: 11 kmers
        G.contigs[G.getID()] = ["CAGGTATCC",[],[(1,True),(2,True)],0]    #5 k-mers
        G.contigs[G.getID()] = ["ATCCCC",[(0,True)],[],0]                #3 k-mers
        G.contigs[G.getID()] = ["ATCCAA",[(0,True)],[],0]                #3 k-mers
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertEqual(G.degrees[0][2] , 1)
        self.assertEqual(G.degrees[1][2] , 1)
        self.assertEqual(G.degrees[2][2] , 1)

class Test_existsAltPath_fw(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #Make sure short contigs not connected to anything are classified
    #as isolated
    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["ACGGCGTTAGCCAGGT",[],[(1,True),(2,True)],0]            #12
        G.contigs[G.getID()] = [            "AGTTA",[(0,True)],[],0]                    #1
        G.contigs[G.getID()] = [            "AGTTAATGTTAGGTACCGT",[(0,True)],[],0]      #15
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1)
        self.assertTrue(G.existsAltPath_fw(MPL,0,1))
        self.assertFalse(G.existsAltPath_fw(MPL,0,2))

class Test_isPartOfTip(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #Make sure we correctly raise an Exception if c_ID isn't in G
    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[],0]
        G.contigs[G.getID()] = ["AGTTA",[],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertTrue("CAGGT" in G.kmers)
        self.assertTrue("CAGGA" not in G.kmers)
        self.assertTrue(not "CAGGA" in G.kmers)
        with self.assertRaises(AssertionError):
            G.isPartOfTip(2,MPL)

    #Make sure we correctly raise an Exception if c_ID is part of an
    #isolated collection of contigs
    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        with self.assertRaises(AssertionError):
            G.isPartOfTip(0,MPL)

    #Make sure we correctly identify a simple right tip
    def test_3(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]
        G.contigs[G.getID()] = [               "GCCATTA",[(0,True)],[],0]
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCATAAAAAAAA",[(0,True)],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertFalse(G.isPartOfTip(0,MPL))
        self.assertTrue(G.isPartOfTip(1,MPL))
        self.assertFalse(G.isPartOfTip(2,MPL))

    #Make sure we correctly identify a simple left tip
    def test_4(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [                "CAGGTATTTGTTATGGCCATT",[(1,True),(2,True)],[],0]   #A
        G.contigs[G.getID()] = [               "CCAGGTA",[],[(0,True)],0]                           #C
        G.contigs[G.getID()] = ["AAAAAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                           #B
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertFalse(G.isPartOfTip(0,MPL))
        self.assertTrue(G.isPartOfTip(1,MPL))
        self.assertFalse(G.isPartOfTip(2,MPL))

    #Make sure we got the length right for a simple right tip
    #A->B
    #A->C
    #C is a tip if:
    #   1. C is isolated without the connection to A and
    #   2. the number of k-mers in B >= MPL=2*k.    
    #           Number of k-mers in B = |b|-k+1
    def test_5(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]           #A
        G.contigs[G.getID()] = [               "GCCATTA",[(0,True)],[],0]                   #C
        G.contigs[G.getID()] = [               "GCCATTCCTATTTTTTTTCC",[(0,True)],[],0]      #B of length =2k
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertFalse(G.isPartOfTip(0,MPL))
        self.assertTrue(G.isPartOfTip(1,MPL))
        self.assertFalse(G.isPartOfTip(2,MPL))

    def test_6(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["AGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]           #A. 14=2k
        G.contigs[G.getID()] = [               "GCCATTA",[(0,True)],[],0]                  #C. 1
        G.contigs[G.getID()] = [               "GCCATTCCTATTTTTTTTC",[(0,True)],[],0]      #B. 13<2k
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.markAllIsolated(MPL)
        self.assertFalse(G.isPartOfTip(0,MPL))
        assert(G.degrees[1][2]==-1), "We shouldn't have changed any marking"
        assert(G.degrees[2][2]==-1), "We shouldn't have changed any marking"
        self.assertFalse(G.isPartOfTip(1,MPL))
        self.assertFalse(G.isPartOfTip(2,MPL))

    def test_6a(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["AGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]           #A. 14
        G.contigs[G.getID()] = [              "GCCATTA",[(0,True)],[],0]                   #C. 1
        G.contigs[G.getID()] = [              "GCCATTCCTATTTTTTTTC",[(0,True)],[],0]       #B. 13<2k=14
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertFalse(G.isPartOfIsolated(MPL, 0))
        self.assertFalse(G.isPartOfIsolated(MPL, 1))
        self.assertFalse(G.isPartOfIsolated(MPL, 2))
        self.assertFalse(G.isPartOfIsolated(MPL, 0 , {1,2}))
        self.assertTrue(G.isPartOfIsolated(MPL, 1, {0}))
        self.assertTrue(G.isPartOfIsolated(MPL, 2, {0}))
        G.markAllIsolated(MPL)
        for c_ID in G.contigs:
            assert(G.degrees[c_ID][2]==-1)

    #Make sure we got the length right for a simple left tip
    #C->A
    #B->A
    #C is a tip if:
    #   1. C is isolated without the connection to A and
    #   2. the number of k-mers in B >= MPL=2*k.    
    #           Number of k-mers in B = |b|-k+1
    def test_7(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [               "CAGGTATTTGTTATGGCCATAAAAAAAA",[(1,True),(2,True)],[],0]    #A. 22
        G.contigs[G.getID()] = [               "CCAGGTA",[],[(0,True)],0]                                  #C. 1
        G.contigs[G.getID()] = ["AAAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                                    #B. 14=2k
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertFalse(G.isPartOfTip(0,MPL))
        assert(G.degrees[1][2]==-1), "We shouldn't have changed any marking"
        assert(G.degrees[2][2]==-1), "We shouldn't have changed any marking"
        self.assertTrue(G.isPartOfTip(1,MPL))
        self.assertFalse(G.isPartOfTip(2,MPL))

    def test_8(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [               "CAGGTATTTGTTATGGCCATAAAAAAAA",[(1,True),(2,True)],[],0]   #A. 22
        G.contigs[G.getID()] = [               "CCAGGTA",[],[(0,True)],0]                                 #C. 1
        G.contigs[G.getID()] = ["AAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                                    #B. 13<2k
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertFalse(G.isPartOfTip(0,MPL))
        assert(G.degrees[1][2]==-1), "We shouldn't have changed any marking"
        assert(G.degrees[2][2]==-1), "We shouldn't have changed any marking"
        self.assertFalse(G.isPartOfTip(1,MPL))
        assert(G.degrees[0][2]==-1), "We shouldn't have changed any marking"
        assert(G.degrees[2][2]==-1), "We shouldn't have changed any marking"
        self.assertFalse(G.isPartOfTip(2,MPL))
        assert(G.degrees[0][2]==-1), "We shouldn't have changed any marking"
        assert(G.degrees[1][2]==-1), "We shouldn't have changed any marking"

    #Make sure we correctly identify a simple left tip when it's connected to the
    #reverse complement of a longer contig
    def test_9(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [              "CAGGTATTTGTTATGGCCATAAAAAAAA",[(1,False),(2,True)],[],0]    #A. 22
        G.contigs[G.getID()] = [              "ATACCTG",[(0,False)],[],0]                                  #C. 1
        G.contigs[G.getID()] = ["AAAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                                    #B. 14=2k
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertFalse(G.isPartOfTip(0,MPL))
        assert(G.degrees[1][2]==-1), "We shouldn't have changed any marking"
        assert(G.degrees[2][2]==-1), "We shouldn't have changed any marking"
        self.assertTrue(G.isPartOfTip(1,MPL))
        self.assertFalse(G.isPartOfTip(2,MPL))

    def test_9a(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [              "CAGGTATTTGTTATGGCCATAAAAAAAA",[(1,False),(2,True)],[],0]    #A. 22
        G.contigs[G.getID()] = [              "ATACCTG",[(0,False)],[],0]                                  #C. 1
        G.contigs[G.getID()] = ["AAAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                                    #B. 14=2k
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertTrue(G.isPartOfTip(1,MPL))

    #Check the length for a simple left tip
    def test_10(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [               "CAGGTATTTGTTATGGCCATAAAAAAAA",[(1,False),(2,True)],[],0]    #A. 22
        G.contigs[G.getID()] = [               "ATACCTG",[(0,False)],[],0]                                  #C. 1
        G.contigs[G.getID()] = ["AAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                                      #B. 13<2k
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertFalse(G.isPartOfTip(0,MPL))
        self.assertFalse(G.isPartOfTip(1,MPL))
        self.assertFalse(G.isPartOfTip(2,MPL))

    #Make sure we correctly identify a simple right tip when it's connected to the
    #reverse complement of a longer contig
    def test_11(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,False),(2,True)],0]           #A. 15
        G.contigs[G.getID()] = [               "AAATGGC",[],[(0,False)],0]                   #C. 1
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCAT",[(0,True)],[],0]       #B. 14
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertFalse(G.isPartOfTip(0,MPL))
        self.assertTrue(G.isPartOfTip(1,MPL))
        self.assertFalse(G.isPartOfTip(2,MPL))

    #Check the length for a simple right tip
    def test_12(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,False),(2,True)],0]           #A. 15
        G.contigs[G.getID()] = [               "AAATGGC",[],[(0,False)],0]                   #C. 1
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCA",[(0,True)],[],0]        #B. 13
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        self.assertFalse(G.isPartOfTip(0,MPL))
        self.assertFalse(G.isPartOfTip(1,MPL))
        self.assertFalse(G.isPartOfTip(2,MPL))

class Test_findAllConnectedContigs(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    def test_1(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]
        G.contigs[G.getID()] = [               "GCCATTA",[(0,True)],[],0]
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCATAAAAAAAA",[(0,True)],[],0]
        G.addKmersFromAllContigs()
        x0 = G.findAllConnectedContigs(0)
        self.assertEqual(x0 , {0,1,2})
        x1 = G.findAllConnectedContigs(1)
        self.assertEqual(x1 , {0,1,2})
        x2 = G.findAllConnectedContigs(2)
        self.assertEqual(x2 , {0,1,2})

    def test_2(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]
        G.contigs[G.getID()] = [               "GCCATTA",[(0,True)],[],0]
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCATAAAAAAAA",[(0,True)],[],0]
        G.addKmersFromAllContigs()
        x0 = G.findAllConnectedContigs(0,skipNode=1)
        self.assertEqual(x0 , {0,2})
        x1 = G.findAllConnectedContigs(1,skipNode=0)
        self.assertEqual(x1 , {1})
        x2 = G.findAllConnectedContigs(2,skipNode=0)
        self.assertEqual(x2 , {2})

    def test_3(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]                       #0
        G.contigs[G.getID()] = [               "GCCATTA",[(0,True)],[],0]                               #1
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCATAAAAAAAA",[(0,True)],[(3,True)],0]  #2
        G.contigs[G.getID()] = [               "GGGGG",[(2,True)],[(4,True)],0]                         #3
        G.contigs[G.getID()] = [               "AAAAA",[(3,True)],[],0]                                 #4
        G.addKmersFromAllContigs()
        x01 = G.findAllConnectedContigs(0,skipNode=1)
        self.assertEqual(x01 , {0,2,3,4})
        x10 = G.findAllConnectedContigs(1,skipNode=0)
        self.assertEqual(x10 , {1})
        x20 = G.findAllConnectedContigs(2,skipNode=0)
        self.assertEqual(x20 , {2,3,4})
        x23 = G.findAllConnectedContigs(2,skipNode=3)
        self.assertEqual(x23 , {0,1,2})
        x32 = G.findAllConnectedContigs(3,skipNode=2)
        self.assertEqual(x32 , {3,4})
        x34 = G.findAllConnectedContigs(3,skipNode=4)
        self.assertEqual(x34 , {0,1,2,3})

class Test_getDirectAncestors(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        G = Graph.Graph(k)
        return k, G

    #Make sure we correctly find the ancestors in the most simple case, i.e. c0->c1
    def test_1(self):
        k, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True)],0]
        G.contigs[G.getID()] = ["AGTTA",[(0,True)],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        A = G.getDirectAncestors(1,True)
        self.assertEqual(A,{(0,True)})

    #c0->c1
    #c2->c1
    def test_2(self):
        k, G = self.helper()
        G.contigs[G.getID()] = ["CAGTT",[],[(1,True)],0]
        G.contigs[G.getID()] = [ "AGTTA",[(0,True),(2,True)],[],0]
        G.contigs[G.getID()] = ["AAGTT",[],[(1,True)],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        A = G.getDirectAncestors(1,True)
        self.assertEqual(A,{(0,True),(2,True)})

    def test_3(self):
        k, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,False)],0]
        G.contigs[G.getID()] = ["AGTTA",[(0,False)],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        A = G.getDirectAncestors(1,True)
        self.assertEqual(A,{(0,False)})

class Test_getAllPotentialFrontsFromC(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #The most simple case: c0->c1
    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True)],0]
        G.contigs[G.getID()] = [ "AGGTA",[(0,True)],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        fronts = G.getAllPotentialFrontsFromC(1,MPL)
        self.assertEqual(fronts,[0])

    #A longer version of the case above: c0->c1->c2->c3
    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True)],0]
        G.contigs[G.getID()] = [ "AGGTA",[(0,True)],[(2,True)],0]
        G.contigs[G.getID()] = [  "GGTAA",[(1,True)],[(3,True)],0]
        G.contigs[G.getID()] = [   "GTAAC",[(2,True)],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        fronts1 = G.getAllPotentialFrontsFromC(1,MPL)
        self.assertEqual(fronts1,[0])
        fronts2 = G.getAllPotentialFrontsFromC(2,MPL)
        self.assertEqual(fronts2,[0,1])
        fronts3 = G.getAllPotentialFrontsFromC(3,MPL)
        self.assertEqual(fronts3,[0,1,2])

    #The most simple case traveling through a false positive
    def test_3(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,False)],0]
        G.contigs[G.getID()] = [ "TACCT",[(2,False)],[(0,False)],0]
        G.contigs[G.getID()] = [  "GGTAA",[(1,False)],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        fronts1 = G.getAllPotentialFrontsFromC(1,MPL)
        self.assertEqual(fronts1,[2])
        fronts2 = G.getAllPotentialFrontsFromC(2,MPL)
        self.assertEqual(fronts2,[0,1])

class Test_markBubble(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #The most simple case: c0->c1
    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]       #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]    #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[],0]  #3 end
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        front_ID = 0
        path1 = Path.Path()
        path1.append(2,True,G.contigs[2][0],k)
        path2 = Path.Path()
        path2.append(1,True,G.contigs[1][0],k)
        end_ID = 3
        G.markBubble(front_ID,end_ID,path1,path2)
        self.assertEqual(G.degrees[0][2],3)
        self.assertEqual(G.degrees[1][2],4)
        self.assertEqual(G.degrees[2][2],3)
        self.assertEqual(G.degrees[3][2],3)

    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,False),(2,True)],0]              #0 front 
        G.contigs[G.getID()] = [    "GACCT",[(3,False)],[(0,False)],0]          #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]            #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,False),(2,True)],[],0]         #3 end
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        front_ID = 0
        path1 = Path.Path()
        path1.append(2,True,G.contigs[2][0],k)
        path2 = Path.Path()
        path2.append(1,False,G.contigs[1][0],k)
        end_ID = 3
        G.markBubble(front_ID,end_ID,path1,path2)
        self.assertEqual(G.degrees[0][2],3)
        self.assertEqual(G.degrees[1][2],4)
        self.assertEqual(G.degrees[2][2],3)
        self.assertEqual(G.degrees[3][2],3)

class Test_existsPathFrom_n1_to_n2_skipping_intNodes(unittest.TestCase):
    #existsPathFrom_n1_to_n2_skipping_intNodes(self,c_ID,end_ID,intNodes1,MPL)
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #The most simple case of a bubble
    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]       #0 front. 1
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]    #1 short int. 1
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int. 4
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[],0]  #3 end. 1
        G.addKmersFromAllContigs()
        G.createDegreeTable()

        #Find the bubble using isFront and make sure it's working
        temp = G.isFront(0,MPL)
        self.assertFalse(temp==[])
        [end_ID,intNodes1,intNodes2] = temp
        int1_correct = Path.Path()
        int2_correct = Path.Path()
        int1_correct.append(2,True,G.contigs[2][0],k)
        int2_correct.append(1,True,G.contigs[1][0],k)
        self.assertTrue(intNodes1.equals(int1_correct))
        self.assertTrue(intNodes2.equals(int2_correct))
        self.assertEqual(end_ID , 3)

    #Manually find front_ID, end_ID and one of the paths and test whether existsPathFrom_n1_to_n2_skipping_intNodes
    #Finds the alternate path
    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]       #0 front. 1
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]    #1 short int. 1
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int. 4
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[],0]  #3 end. 1
        G.addKmersFromAllContigs()
        G.createDegreeTable()

        #Create front_ID, intNodes1 and end_ID
        front_ID = 0
        end_ID = 3
        intNodes1 = Path.Path()
        intNodes1.append(2,True,G.contigs[2][0],k)

        #Find intNodes2 using existsPathFrom_n1_to_n2_skipping_intNodes
        intNodes2 = G.existsPathFrom_n1_to_n2_skipping_intNodes(front_ID,end_ID,intNodes1,MPL)

        intNodes2_correct = Path.Path()
        intNodes2_correct.append(1,True,G.contigs[1][0],k)
        self.assertTrue(intNodes2.equals(intNodes2_correct))


    #Now test the length of the bubble. The length of each internal path must be <MPL
    #We ignore the length of front and end since they are not part of the bubble
    def test_3(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["GTAGTTCAGCAGGT",[],[(1,True),(2,True)],0]         #0 front. 10
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]               #1 short int. 1
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]               #2 long int. 4
        G.contigs[G.getID()] = [     "GGTCAAATTGTTCC",[(1,True),(2,True)],[],0]    #3 end. 10
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        front_ID = 0
        end_ID = 3
        intNodes1 = Path.Path()
        intNodes1.append(1,True,G.contigs[1][0],k)
        intNodes2 = G.existsPathFrom_n1_to_n2_skipping_intNodes(front_ID,end_ID,intNodes1,MPL)
        self.assertFalse(intNodes2.isEmpty())
        intNodes2_correct = Path.Path()
        intNodes2_correct.append(2,True,G.contigs[2][0],k)
        self.assertTrue(intNodes2.equals(intNodes2_correct))

class Test_isFront(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #The most simple case of a bubble
    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]       #0 front. 1
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]    #1 short int. 1
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int. 4
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[],0]  #3 end. 1
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        temp = G.isFront(0,MPL)
        self.assertFalse(temp==[])
        [end_ID,intNodes1,intNodes2] = temp
        int1_correct = Path.Path()
        int2_correct = Path.Path()
        int1_correct.append(2,True,G.contigs[2][0],k)
        int2_correct.append(1,True,G.contigs[1][0],k)
        self.assertTrue(intNodes1.equals(int1_correct))
        self.assertTrue(intNodes2.equals(int2_correct))
        self.assertEqual(end_ID , 3)

    #Same as test 1 except we added a tip
    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]               #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]            #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]            #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[(4,True)],0]  #3 end
        G.contigs[G.getID()] = [      "GTCAA",[(3,True)],[],0]                  #4 tip
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        temp = G.isFront(0,MPL)
        self.assertFalse(temp==[])
        [end_ID,intNodes1,intNodes2] = temp
        int1_correct = Path.Path()
        int2_correct = Path.Path()
        int1_correct.append(2,True,G.contigs[2][0],k)
        int2_correct.append(1,True,G.contigs[1][0],k)
        self.assertTrue(intNodes1.equals(int1_correct))
        self.assertTrue(intNodes2.equals(int2_correct))
        self.assertEqual(end_ID , 3)

    #Same as test 1 except we added a tip
    def test_2a(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]               #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]            #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]            #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[(4,True)],0]  #3 end
        G.contigs[G.getID()] = [      "GTCAA",[(3,True)],[],0]                  #4 tip
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.degrees[4][2] = 2
        temp = G.isFront(0,MPL)
        self.assertFalse(temp==[])


    #Same as test 1 again except now we're making sure isFront returns false for contigs 1,2,3
    def test_3(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]       #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]    #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[],0]  #3 end
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B1 = G.isFront(1,MPL)
        self.assertEqual(B1 , [])
        B2 = G.isFront(2,MPL)
        self.assertEqual(B2 , [])
        B3 = G.isFront(3,MPL)
        self.assertEqual(B3 , [])

    #Same as test 1 again except now we've removed the alternate path
    def test_4(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(2,True)],0]       #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[],[],0]                    #1 short int (without connections)
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(2,True)],[],0]  #3 end
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B = G.isFront(0,MPL)
        self.assertEqual(B , [])
        self.assertFalse(B)

    def test_5(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[],0]
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B = G.isFront(0,MPL)
        self.assertFalse(B)

    def test_6(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True)],0] 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[],0] 
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B = G.isFront(0,MPL)
        self.assertFalse(B)
        B1 = G.isFront(1,MPL)
        self.assertFalse(B1)

    def test_7(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True)],0] 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(2,True)],0]
        G.contigs[G.getID()] = [     "GGTCC",[(1,True)],[],0] 
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        B = G.isFront(0,MPL)
        self.assertFalse(B)
        B1 = G.isFront(1,MPL)
        self.assertFalse(B1)
        B2 = G.isFront(2,MPL)
        self.assertFalse(B2)

    #Now test the length of the bubble. The length of each internal path must be <MPL
    #We ignore the length of front and end since they are not part of the bubble
    def test_8(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["GTAGTTCAGCAGGT",[],[(1,True),(2,True)],0]         #0 front. 10
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]               #1 short int. 1
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]               #2 long int. 4
        G.contigs[G.getID()] = [     "GGTCAAATTGTTCC",[(1,True),(2,True)],[],0]    #3 end. 10
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        i = G.isFront(0,MPL)
        self.assertFalse(i==[])
        if i:
            [end_ID,intNodes1,intNodes2] = i
            intNodes2_correct = Path.Path()
            intNodes2_correct.append(1,True,G.contigs[1][0],k)
            intNodes1_correct = Path.Path()
            intNodes1_correct.append(2,True,G.contigs[2][0],k)
            self.assertTrue(intNodes1.equals(intNodes1_correct))
            self.assertTrue(intNodes2.equals(intNodes2_correct))
            self.assertEqual(end_ID , 3)

    def test_9(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["GTAGTTC",[],[(1,True),(2,True)],0]          #0 front. 3
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]         #1 short int. 1
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]         #2 long int. 4
        G.contigs[G.getID()] = [     "GGTCAA",[(1,True),(2,True)],[],0]      #3 end. 2
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        i = G.isFront(0,MPL)
        self.assertFalse(i==[])
        if i:
            [end_ID,intNodes1,intNodes2] = i
            intNodes2_correct = Path.Path()
            intNodes2_correct.append(1,True,G.contigs[1][0],k)
            intNodes1_correct = Path.Path()
            intNodes1_correct.append(2,True,G.contigs[2][0],k)
            self.assertTrue(intNodes1.equals(intNodes1_correct))
            self.assertTrue(intNodes2.equals(intNodes2_correct))
            self.assertEqual(end_ID , 3)

    def test_10(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["GTAGTTC",[],[(1,True),(2,True)],0]          #0 front. 3
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]         #1 short int. 1
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]         #2 long int. 4
        G.contigs[G.getID()] = [     "GGTCAAT",[(1,True),(2,True)],[],0]     #3 end. 3
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        i = G.isFront(0,MPL)
        self.assertFalse(i==[])
        if i:
            [end_ID,intNodes1,intNodes2] = i
            intNodes2_correct = Path.Path()
            intNodes2_correct.append(1,True,G.contigs[1][0],k)
            intNodes1_correct = Path.Path()
            intNodes1_correct.append(2,True,G.contigs[2][0],k)
            self.assertTrue(intNodes1.equals(intNodes1_correct))
            self.assertTrue(intNodes2.equals(intNodes2_correct))
            self.assertEqual(end_ID , 3)

class Test_markAllBubbles(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #The most simple case: c0->c1
    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]       #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]    #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[],0]  #3 end
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        #theContigs = {0,1,2,3}
        G.markAllBubbles(MPL)
        self.assertEqual(G.degrees[0][2],3)
        self.assertEqual(G.degrees[1][2],4)
        self.assertEqual(G.degrees[2][2],3)
        self.assertEqual(G.degrees[3][2],3)

    #Same as test 1 except we added a tip
    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]               #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]            #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]            #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[(4,True)],0]  #3 end
        G.contigs[G.getID()] = [      "GTCAA",[(3,True)],[],0]                  #4 tip
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.degrees[4][2] = 2
        #theContigs = {0,1,2,3}
        G.markAllBubbles(MPL)
        self.assertEqual(G.degrees[0][2],3)
        self.assertEqual(G.degrees[1][2],4)
        self.assertEqual(G.degrees[2][2],3)
        self.assertEqual(G.degrees[3][2],3)

    #Same as test 1 except theContigs only contains one contig in the bubble
    def test_3(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]       #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]    #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[],0]  #3 end
        G.addKmersFromAllContigs()
        for c_ID in [0,1,2,3]:
            #print "\ntesting for c_ID="+str(c_ID)
            G.createDegreeTable()
            #theContigs = {c_ID}
            G.markAllBubbles(MPL)
            self.assertEqual(G.degrees[0][2],3)
            self.assertEqual(G.degrees[1][2],4)
            self.assertEqual(G.degrees[2][2],3)
            self.assertEqual(G.degrees[3][2],3)

    def test_3a(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]       #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]    #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[],0]  #3 end
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        #theContigs = {1}
        G.markAllBubbles(MPL)
        self.assertEqual(G.degrees[0][2],3)
        self.assertEqual(G.degrees[1][2],4)
        self.assertEqual(G.degrees[2][2],3)
        self.assertEqual(G.degrees[3][2],3)

    #Same as test 1 except we have added a tip and theContigs only contains one contig in the bubble (not the tip)
    def test_4(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]               #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]            #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]            #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[(4,True)],0]  #3 end
        G.contigs[G.getID()] = [      "GTCAA",[(3,True)],[],0]                  #4 tip
        G.addKmersFromAllContigs()
        for c_ID in [0,1,2,3]:
            #print "testing for c_ID="+str(c_ID)
            G.createDegreeTable()
            G.degrees[4][2] = 2
            #theContigs = {c_ID}
            G.markAllBubbles(MPL)
            self.assertEqual(G.degrees[0][2],3)
            self.assertEqual(G.degrees[1][2],4)
            self.assertEqual(G.degrees[2][2],3)
            self.assertEqual(G.degrees[3][2],3)

class Test_analyzeAllContigsInCollection(unittest.TestCase):
    def helper(self,k=-1):
        if k==-1:
            k = 5
        MPL = 2*k
        G = Graph.Graph(k)
        return k, MPL, G

    #Make sure short contigs not connected to anything are classified
    #as isolated
    def test_1(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[],0]
        G.contigs[G.getID()] = ["AGTTA",[],[],0]
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 1)
        self.assertEqual(G.degrees[1][2] , 1)

    #Make sure that a long (len > 2k) contig is not classified as isolated
    def test_2(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[],0]
        G.contigs[G.getID()] = ["AGTTA",[],[],0]
        G.contigs[G.getID()] = ["AGTTAAACCCGAGGTTAGACCCTAGGACT",[],[],0]
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 1)
        self.assertEqual(G.degrees[1][2] , 1)
        self.assertEqual(G.degrees[2][2] , 0)

    #Make sure that two contigs which together are long are not classified as
    #isolated
    def test_3(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTTT",[(0,True)],[],0]      #3 k-mers
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 0)

    #Make sure that two contigs which together are long are not classified as
    #isolated. Check off by 1 error
    def test_4(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGTATCCAT",[],[(1,True)],0]  #7 k-mers
        G.contigs[G.getID()] = ["CCATTT",[(0,True)],[],0]       #2 k-mers
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 1)
        self.assertEqual(G.degrees[1][2] , 1)

    #Make sure that a collection of connected contigs is still classified as isolated
    #so long as the longest path they're part of is < 2k. Even though the total length 
    #of the entire collection is greater than 2k
    def test_5(self):
        k, MPL, G = self.helper()
        #Longest path: 8 kmers
        #Total length: 11 kmers
        G.contigs[G.getID()] = ["CAGGTATCC",[],[(1,True),(2,True)],0]    #5 k-mers
        G.contigs[G.getID()] = ["ATCCCC",[(0,True)],[],0]                #3 k-mers
        G.contigs[G.getID()] = ["ATCCAA",[(0,True)],[],0]                #3 k-mers
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 1)
        self.assertEqual(G.degrees[1][2] , 1)
        self.assertEqual(G.degrees[2][2] , 1)

    #-----------------------------------------------------------------------
    #---------------------The tests from the tip tests----------------------
    #-----------------------------------------------------------------------
    #Make sure we correctly identify a simple right tip
    def test_6(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]
        G.contigs[G.getID()] = [               "GCCATTA",[(0,True)],[],0]
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCATAAAAAAAA",[(0,True)],[],0]
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 2)
        self.assertEqual(G.degrees[2][2] , 0)

    #Make sure we correctly identify a simple left tip
    def test_7(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [                "CAGGTATTTGTTATGGCCATT",[(1,True),(2,True)],[],0]   #A
        G.contigs[G.getID()] = [               "CCAGGTA",[],[(0,True)],0]                           #C
        G.contigs[G.getID()] = ["AAAAAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                           #B
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 2)
        self.assertEqual(G.degrees[2][2] , 0)

    #Make sure we got the length right for a simple right tip
    #A->B
    #A->C
    #C is a tip if:
    #   1. C is isolated without the connection to A and
    #   2. the number of k-mers in B >= MPL=2*k.    
    #           Number of k-mers in B = |b|-k+1
    def test_8(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]           #A
        G.contigs[G.getID()] = [               "GCCATTA",[(0,True)],[],0]                   #C
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCAT",[(0,True)],[],0]      #B of length 3k-1=20
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 2)
        self.assertEqual(G.degrees[2][2] , 0)

    def test_9(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,True),(2,True)],0]           #A
        G.contigs[G.getID()] = [               "GCCATTA",[(0,True)],[],0]                   #C
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCA",[(0,True)],[],0]       #B of length <3k-1
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 0)
        self.assertEqual(G.degrees[2][2] , 0)

    #Make sure we got the length right for a simple left tip
    #C->A
    #B->A
    #C is a tip if:
    #   1. C is isolated without the connection to A and
    #   2. the number of k-mers in B >= MPL=2*k.    
    #           Number of k-mers in B = |b|-k+1
    def test_10(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [               "CAGGTATTTGTTATGGCCATAAAAAAAA",[(1,True),(2,True)],[],0]    #A
        G.contigs[G.getID()] = [               "CCAGGTA",[],[(0,True)],0]                                  #C
        G.contigs[G.getID()] = ["AAAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                                    #B of length 3k-1=20
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 2)
        self.assertEqual(G.degrees[2][2] , 0)

    def test_11(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [               "CAGGTATTTGTTATGGCCATAAAAAAAA",[(1,True),(2,True)],[],0]   #A
        G.contigs[G.getID()] = [               "CCAGGTA",[],[(0,True)],0]                                 #C
        G.contigs[G.getID()] = ["AAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                                    #B of length <3k-1
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 0)
        self.assertEqual(G.degrees[2][2] , 0)

    #Make sure we correctly identify a simple left tip when it's connected to the
    #reverse complement of a longer contig
    def test_12(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [               "CAGGTATTTGTTATGGCCATAAAAAAAA",[(1,False),(2,True)],[],0]    #A
        G.contigs[G.getID()] = [               "ATACCTG",[(0,False)],[],0]                                  #C
        G.contigs[G.getID()] = ["AAAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                                     #B of length 3k-1=20
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 2)
        self.assertEqual(G.degrees[2][2] , 0)

    #Check the length for a simple left tip
    def test_13(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = [               "CAGGTATTTGTTATGGCCATAAAAAAAA",[(1,False),(2,True)],[],0]    #A
        G.contigs[G.getID()] = [               "ATACCTG",[(0,False)],[],0]                                  #C
        G.contigs[G.getID()] = ["AAAAAAACCTAGGCAGGTA",[],[(0,True)],0]                                      #B of length <3k-1
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 0)
        self.assertEqual(G.degrees[2][2] , 0)

    #Make sure we correctly identify a simple right tip when it's connected to the
    #reverse complement of a longer contig
    def test_14(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,False),(2,True)],0]           #A
        G.contigs[G.getID()] = [               "AAATGGC",[],[(0,False)],0]                   #C
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCAT",[(0,True)],[],0]       #B of length 3k-1=20
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 2)
        self.assertEqual(G.degrees[2][2] , 0)

    #Check the length for a simple right tip
    def test_15(self):
        k, MPL, G = self.helper(7)
        G.contigs[G.getID()] = ["CAGGTATTTGTTATGGCCATT",[],[(1,False),(2,True)],0]           #A
        G.contigs[G.getID()] = [               "AAATGGC",[],[(0,False)],0]                   #C
        G.contigs[G.getID()] = [               "GCCATTCCTATTTATTCCA",[(0,True)],[],0]       #B of length <3k-1
        G.addKmersFromAllContigs()
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2] , 0)
        self.assertEqual(G.degrees[1][2] , 0)
        self.assertEqual(G.degrees[2][2] , 0)

    #-----------------------------------------------------------------------
    #---------------The tests from the mark all bubbles tests---------------
    #-----------------------------------------------------------------------
    #The most simple case: c0->c1
    def test_16(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]                   #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]                #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]                #2 long int
        G.contigs[G.getID()] = [     "GGTCAAAAAAAAAAA",[(1,True),(2,True)],[],0]              #3 end
        G.addKmersFromAllContigs()
        theContigs = {0,1,2,3}
        G.analyzeAllContigsInCollection()
        self.assertEqual(G.degrees[0][2],3)
        self.assertEqual(G.degrees[1][2],4)
        self.assertEqual(G.degrees[2][2],3)
        self.assertEqual(G.degrees[3][2],3)

    #Same as test 1 except we added a tip
    def test_17(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]               #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]            #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]            #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[(4,True)],0]  #3 end
        G.contigs[G.getID()] = [      "GTCAA",[(3,True)],[],0]                  #4 tip
        G.addKmersFromAllContigs()
        G.createDegreeTable()
        G.degrees[4][2] = 2
        theContigs = {0,1,2,3,4}
        G.markAllBubbles(MPL)
        self.assertEqual(G.degrees[0][2],3)
        self.assertEqual(G.degrees[1][2],4)
        self.assertEqual(G.degrees[2][2],3)
        self.assertEqual(G.degrees[3][2],3)

    #Same as test 1 except theContigs only contains one contig in the bubble
    def test_18(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]       #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]    #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]    #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[],0]  #3 end
        G.addKmersFromAllContigs()
        for c_ID in [0,1,2,3]:
            #print "\ntesting for c_ID="+str(c_ID)
            G.createDegreeTable()
            theContigs = {c_ID}
            G.markAllBubbles(MPL)
            self.assertEqual(G.degrees[0][2],3)
            self.assertEqual(G.degrees[1][2],4)
            self.assertEqual(G.degrees[2][2],3)
            self.assertEqual(G.degrees[3][2],3)

    #Same as test 1 except we have added a tip and theContigs only contains one contig in the bubble (not the tip)
    def test_19(self):
        k, MPL, G = self.helper()
        G.contigs[G.getID()] = ["CAGGT",[],[(1,True),(2,True)],0]               #0 front 
        G.contigs[G.getID()] = [    "AGGTC",[(0,True)],[(3,True)],0]            #1 short int
        G.contigs[G.getID()] = [ "AGGTGGTC",[(0,True)],[(3,True)],0]            #2 long int
        G.contigs[G.getID()] = [     "GGTCA",[(1,True),(2,True)],[(4,True)],0]  #3 end
        G.contigs[G.getID()] = [      "GTCAA",[(3,True)],[],0]                  #4 tip
        G.addKmersFromAllContigs()
        for c_ID in [0,1,2,3]:
            #print "testing for c_ID="+str(c_ID)
            G.createDegreeTable()
            G.degrees[4][2] = 2
            theContigs = {c_ID}
            G.markAllBubbles(MPL)
            self.assertEqual(G.degrees[0][2],3)
            self.assertEqual(G.degrees[1][2],4)
            self.assertEqual(G.degrees[2][2],3)
            self.assertEqual(G.degrees[3][2],3)

if __name__ == '__main__':
    unittest.main()
