#coding:utf8
import Graph
import BF_counter
import Bloom
import time
from bigData import returnTime
import helpers
import sys
import cProfile

def sim_BF_counter(fn,k,BF,G,pfn=False,printProgress=True,startAtLine=0):
    print sim_BF_counter.__name__
    #G.printFunctionNames = True
    #cProfile.run('BF_counter.BF_counter(fn,BF,k,G,pfn,printProgress,startAtLine)')
    BF_counter.BF_counter(fn,BF,k,G,pfn,printProgress,startAtLine)

def sim_BF_counter_naive(fn,k,BF,G_naive,pfn=False,printProgress=True):
    print sim_BF_counter_naive.__name__
    #cProfile.run("BF_counter.BF_counter_naive(fn,BF,k,G_naive,pfn,printProgress)")
    BF_counter.BF_counter_naive(fn,BF,k,G_naive,pfn,printProgress)

def initialize(fn,k,numberOfKmers,fileName=""):
    BF = Bloom.Bloom(0.01,numberOfKmers,pfn=True)
    G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
    if fileName!="":
        print "We read the graph from the file "+fileName
        G.createGraphFromFile(fileName)
    return BF,G

def readInput():
    k = int(sys.argv[1])
    read1 = sys.argv[2]
    read2 = sys.argv[3]
    numberOfKmers = int(sys.argv[4])
    whatToRun = int(sys.argv[5])        #0: BF_counter_naive. 1: BF_counter
    pfn = False
    printProgress = True
    startAtLine = 0
    fileName = ""
    if len(sys.argv)>6:
        pfn = map((lambda x: {"False":False,"True":True}[x]), [sys.argv[6]])[0]
    if len(sys.argv)>7:
        printProgress = map((lambda x: {"False":False,"True":True}[x]), [sys.argv[7]])[0]
    if len(sys.argv)>8:
        assert len(sys.argv)>9
        startAtLine = int(sys.argv[8])
        fileName = sys.argv[9]
        
    assert isinstance(k, int)
    assert isinstance(read1,str)
    assert isinstance(read2,str)
    assert isinstance(numberOfKmers,int)
    assert isinstance(whatToRun,int)
    assert isinstance(pfn,bool)
    assert isinstance(printProgress,bool)
    assert isinstance(startAtLine,int)
    assert isinstance(fileName,str)
    fn = [read1,read2]
    return k,fn,numberOfKmers,whatToRun,pfn,printProgress,startAtLine,fileName

if __name__ == '__main__':
    #command line syntax:
    #python simData.py k       read1           read2      numberOfKmers whatToRead [pfn] [printProgress] [startAtLine] [fileName]
    #python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 0
    #python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 0 False False 0 BF_counter_KeyboardInterrupt/G.txt
    start = time.time()
    k,fn,numberOfKmers,whatToRun,pfn,printProgress,startAtLine,fileName = readInput()

    BF,G = initialize(fn,k,numberOfKmers,fileName)
    if whatToRun==0:
        sim_BF_counter_naive(fn,k,BF,G,pfn,printProgress)
    elif whatToRun==1:
        sim_BF_counter(fn,k,BF,G,pfn,printProgress,startAtLine)
    else:
        raise Exception("whatToRun must be either 0 or 1")
    end = time.time()
    print "Checking whether the graph is legal"
    G.printIsLegal = True
    B = G.isLegalDBG()
    if B:
        print "The Graph is legal"
    else:
        print "The Graph is not legal"
    helpers.printResultsToFile(BF,G,outFolder="defaultOutFolder",timeInSeconds=int(end-start))
    print "We have finished writing the Graph and all info to a file"
    print "Checking whether the graph is legal"
    G.printIsLegal = True
    B = G.isLegalDBG()
    if B:
        print "The Graph is legal"
    else:
        print "The Graph is not legal"
    