#coding:utf8
import os, sys
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import unittest
import Graph
import helpers
import dbg
import collections

class Test_fix_i(unittest.TestCase):
    def doTheWork(self,c,k,t_i):
        t = dbg.twin(c)
        c_i = helpers.fix_i(t_i,len(c),k)
        c_km = c[c_i:c_i+k]
        t_km = t[t_i:t_i+k]
        self.assertTrue(len(c_km)==k)
        self.assertTrue(len(t_km)==k)
        self.assertTrue(c_km==dbg.twin(t_km))
        return c_i

    def test_1(self):
        k = 5
        c = "AAAAA"
        t_i = 0
        c_i = self.doTheWork(c,k,t_i)
        self.assertTrue(c_i==0)
        

    def test_2(self):
        k = 5
        #    01234
        #     12345
        #      23456
        #       34567
        c = "AAAAACGT"
        #    ACGTTTTT
        t_i = 3
        c_i = self.doTheWork(c,k,t_i)
        self.assertTrue(c_i==0)

        t_i = 2
        c_i = self.doTheWork(c,k,t_i)
        self.assertTrue(c_i==1)

        t_i = 1
        c_i = self.doTheWork(c,k,t_i)
        self.assertTrue(c_i==2)

        t_i = 0
        c_i = self.doTheWork(c,k,t_i)
        self.assertTrue(c_i==3)

class Test_splitString(unittest.TestCase):
    def test_1(self):
        s0,s1 = helpers.splitString("012345",3,3)
        self.assertEqual(s0,"012")
        self.assertEqual(s1,"12345")

    def test_2(self):
        s0,s1 = helpers.splitString("012345",4,3)
        self.assertEqual(s0,"0123")
        self.assertEqual(s1,"2345")

    def test_3(self):
        s0,s1 = helpers.splitString("012345",5,3)
        self.assertEqual(s0,"01234")
        self.assertEqual(s1,"345")

    def test_4(self):
        s0,s1 = helpers.splitString("012345678",8,4)
        self.assertEqual(s0,"01234567")
        self.assertEqual(s1,"5678")

    #Þarf að skoða þetta betur ef ég ætla að testa fyrir assertions
    """def test_5(self):
        with self.assertRaises(AssertionError) as cm:
            #do_something()
            print "adsfasdfasdf"
            s0,s1 = Graph.splitString("012345678",9,4)
            print s0

        the_exception = cm.exception
        print "bla", the_exception
        print cm
        #self.assertEqual(the_exception.error_code, 3)"""

class Test_splitCov(unittest.TestCase):
    def test_1(self):
        #    012345
        #    01234
        #     12345
        s = "AAAAAC"
        k = 5
        s0,s1 = helpers.splitString(s,5,k)
        self.assertTrue(s0=="AAAAA")
        self.assertTrue(s1=="AAAAC")
        COV0, COV1 = helpers.splitCov(s0,s1,k,4)
        self.assertEqual(COV0,2)
        self.assertEqual(COV1,2)

    def test_2(self):
        #    012345
        #    01234
        #     12345
        s = "AAAAAC"
        k = 5
        s0,s1 = helpers.splitString(s,5,k)
        self.assertTrue(s0=="AAAAA")
        self.assertTrue(s1=="AAAAC")
        COV0, COV1 = helpers.splitCov(s0,s1,k,6)
        self.assertEqual(COV0,3)
        self.assertEqual(COV1,3)

    def test_3(self):
        #    01234567
        #    012345
        #      234567
        s = "AAAAACGT"
        k = 5
        s0,s1 = helpers.splitString(s,6,k)
        self.assertTrue(s0=="AAAAAC")
        self.assertTrue(s1==  "AAACGT")
        COV0, COV1 = helpers.splitCov(s0,s1,k,8)
        self.assertEqual(COV0,4)
        self.assertEqual(COV1,4)

    def test_4(self):
        #    01234567
        #    01234
        #     1234567
        s = "AAAAACGT"
        k = 5
        s0,s1 = helpers.splitString(s,5,k)
        self.assertTrue(s0=="AAAAA")
        self.assertTrue(s1== "AAAACGT")
        COV0, COV1 = helpers.splitCov(s0,s1,k,8)
        self.assertEqual(COV0,2)
        self.assertEqual(COV1,6)

class Test_reverseList(unittest.TestCase):
    def test_1(self):
       L = [(0,True),(1,False)]
       L_correct = [(0,False),(1,True)]
       helpers.reverseList(L)
       L.sort()
       L_correct.sort()
       self.assertEqual(L,L_correct)

class Test_splitOnConnToSelf(unittest.TestCase):
    #splitOnConnToSelf(s,s_twin,k,i,j,start=0)
    def test_1(self):
        s = "AAAAAC"
        k = 5
        seqs = []
        for i,j in helpers.splitOnConnToSelf(s,dbg.twin(s),len(s),k,0,len(s),0):
            seqs.append(s[i:j])
        self.assertTrue("AAAAA" in seqs)
        self.assertTrue("AAAAC" in seqs)
        self.assertTrue(len(seqs)==2)

    def test_2(self):
        s = "AAAACG"
        #    CGTTTT
        k = 5
        seqs = []
        for i,j in helpers.splitOnConnToSelf(s,dbg.twin(s),len(s),k,0,len(s),0):
            seqs.append(s[i:j])
        self.assertTrue("AAAACG" in seqs)
        self.assertTrue(len(seqs)==1)

    def test_3(self):
        s = "AAAAACG"
        k = 5
        seqs = []
        for i,j in helpers.splitOnConnToSelf(s,dbg.twin(s),len(s),k,0,len(s),0):
            seqs.append(s[i:j])
        self.assertTrue("AAAAA" in seqs)
        self.assertTrue("AAAACG" in seqs)
        self.assertTrue(len(seqs)==2)

    def test_4(self):
        s = "CAAAAACTTTTTCC"
        #    GGAAAAAGTTTTTG
        k = 5
        seqs = []
        for i,j in helpers.splitOnConnToSelf(s,dbg.twin(s),len(s),k,0,len(s),0):
            seqs.append(s[i:j])
        self.assertTrue("CAAAA" in seqs)
        self.assertTrue( "AAAAA" in seqs)
        self.assertTrue(  "AAAACTTTT" in seqs)
        self.assertTrue(        "TTTTT" in seqs)
        self.assertTrue(         "TTTTCC" in seqs)
        self.assertTrue(len(seqs)==5)

    #See what happens when len(s)>k and still no split
    def test_5(self):
        s = "AAACCC"
        k = 5
        seqs = []
        for i,j in helpers.splitOnConnToSelf(s,dbg.twin(s),len(s),k,0,len(s),0):
            seqs.append(s[i:j])
        self.assertTrue("AAACCC" in seqs)
        self.assertTrue(len(seqs)==1)

class Test_canConn(unittest.TestCase):
    def test_1(self):
        #a->b
        k = 5
        a = "AAAAA"
        b = "AAAAC"
        self.assertTrue(helpers.canConn(a,b,True,True,k))
        self.assertFalse(helpers.canConn(b,a,True,True,k))

    def test_2(self):
        #a->b
        k = 5
        a = "AAAAA"
        b = "AAACC"
        self.assertFalse(helpers.canConn(a,b,True,True,k))

    def test_3(self):
        #a->twin(b)
        k = 5
        a = "AAAAA" 
        b = "GTTTT" #twin: AAAAC
        self.assertTrue(helpers.canConn(a,b,True,False,k))
        self.assertTrue(helpers.canConn(b,a,True,False,k))
        self.assertFalse(helpers.canConn(a,b,True,True,k))

    def test_4(self):
        #a->twin(b)
        k = 5
        a = "AAAAA"
        b = "GGTTT"
        self.assertFalse(helpers.canConn(a,b,True,False,k))
        self.assertFalse(helpers.canConn(b,a,True,False,k))
        self.assertFalse(helpers.canConn(a,b,True,True,k))

    def test_5(self):
        #twin(a)->b
        k = 5
        a = "TTTTT" #twin: AAAAA
        b = "AAAAC"
        self.assertTrue(helpers.canConn(a,b,False,True,k))
        self.assertTrue(helpers.canConn(b,a,False,True,k))
        self.assertFalse(helpers.canConn(a,b,True,True,k))
        self.assertFalse(helpers.canConn(b,a,True,True,k))

    def test_6(self):
        #twin(a)->twin(b)
        k = 5
        a = "AAAAA" #twin: TTTTT
        b = "CAAAA" #twin: TTTTG
        self.assertTrue(helpers.canConn(a,b,False,False,k))
        self.assertTrue(helpers.canConn(b,a,True,True,k))
        self.assertFalse(helpers.canConn(a,b,True,True,k))
        self.assertFalse(helpers.canConn(a,b,True,False,k))
        self.assertFalse(helpers.canConn(b,a,True,False,k))
        self.assertFalse(helpers.canConn(a,b,False,True,k))
        self.assertFalse(helpers.canConn(b,a,False,True,k))

    def test_self1(self):
        #a->a
        k = 5
        a = "AAAAA"
        self.assertTrue(helpers.canConn(a,a,True,True,k))
        self.assertFalse(helpers.canConn(a,a,True,False,k))
        self.assertFalse(helpers.canConn(a,a,False,True,k))

    def test_self2(self):
        #a->twin(a)
        k = 5
        a = "TATAT" #twin: ATATA
        self.assertTrue(helpers.canConn(a,a,True,False,k))
        self.assertFalse(helpers.canConn(a,a,True,True,k))
        self.assertTrue(helpers.canConn(a,a,False,True,k))
    
    def test_self3(self):
        #twin(a)->a
        k = 5
        a = "ATATA" #twin: TATAT
        self.assertTrue(helpers.canConn(a,a,False,True,k))
        self.assertFalse(helpers.canConn(a,a,True,True,k))
        self.assertTrue(helpers.canConn(a,a,True,False,k))

class Test_cannConnFromTo(unittest.TestCase):
    def test_1(self):
        #a->b
        k = 5
        a = "AAAAA"
        b = "AAAAC"
        c = "AAACC"
        self.assertTrue(helpers.canConnFromTo(a,b,k))
        self.assertTrue(helpers.canConnFromTo(a,a,k))
        self.assertFalse(helpers.canConnFromTo(b,a,k))
        self.assertFalse(helpers.canConnFromTo(a,c,k))

    def test_2(self):
        #a->twin(b)
        k = 5
        a = "AAAAA" 
        b = "GTTTT" #twin: AAAAC
        c = "GGTTT" #twin: AAACC
        self.assertTrue(helpers.canConnFromTo(a,b,k))
        self.assertTrue(helpers.canConnFromTo(b,a,k))
        self.assertFalse(helpers.canConnFromTo(a,c,k))
        self.assertTrue(helpers.canConnFromTo(c,b,k))

class Test_canConnToFrom(unittest.TestCase):
    def test_1(self):
        #a->b
        k = 5
        a = "AAAAA"
        b = "AAAAC"
        c = "AAACC"
        self.assertTrue(helpers.canConnToFrom(b,a,k))
        self.assertTrue(helpers.canConnToFrom(a,a,k))
        self.assertFalse(helpers.canConnToFrom(a,b,k))
        self.assertFalse(helpers.canConnToFrom(c,a,k))

    def test_2(self):
        #a->twin(b)
        k = 5
        a = "AAAAA" #twin: TTTTT
        b = "TTTTG" #twin: CAAAA
        c = "TTTGG" #twin: CCAAA
        self.assertTrue(helpers.canConnToFrom(b,a,k))
        self.assertTrue(helpers.canConnToFrom(a,b,k))
        self.assertFalse(helpers.canConnToFrom(c,a,k))
        self.assertTrue(helpers.canConnToFrom(c,b,k))

class Test_canConnToSelf(unittest.TestCase):
    def test_1(self):
        #a->a
        k = 5
        a = "AAAAA"
        self.assertTrue(helpers.canConnToSelf(a,k))
        self.assertTrue(helpers.canConnToSelfFront(a,k))
        self.assertTrue(helpers.canConnToSelfBack(a,k))

    def test_2(self):
        #a->twin(a) and twin(a)->a
        k = 5
        a = "ATATA" #twin: TATAT
        self.assertTrue(helpers.canConnToSelf(a,k))
        self.assertTrue(helpers.canConnToSelfFront(a,k))
        self.assertTrue(helpers.canConnToSelfBack(a,k))

    def test_3(self):
        #a->twin(a) only
        k = 5
        a = "CATAT" #twin: ATATG
        self.assertTrue(helpers.canConnToSelf(a,k))
        self.assertTrue(helpers.canConnToSelfFront(a,k))
        self.assertFalse(helpers.canConnToSelfBack(a,k))

    def test_4(self):
        #twin(a)->a only
        k = 5
        a = "ATATC" #twin: GATAT
        self.assertTrue(helpers.canConnToSelf(a,k))
        self.assertFalse(helpers.canConnToSelfFront(a,k))
        self.assertTrue(helpers.canConnToSelfBack(a,k))

    def test_5(self):
        #a->a
        k = 5
        a = "ATCGG"
        self.assertFalse(helpers.canConnToSelf(a,k))
        self.assertFalse(helpers.canConnToSelfFront(a,k))
        self.assertFalse(helpers.canConnToSelfBack(a,k))

if __name__ == '__main__':
    unittest.main()