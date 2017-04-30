#coding:utf8
import Graph
import BF_counter
import Bloom, pybloom
import time
from bigData import returnTime
import helpers
import sys
import cProfile

"""
def sim_BF_counter(fn,k,BF,G,pfn=False,printProgress=True,startAtLine=0,sizeOfGenome,maxCov):
    print sim_BF_counter.__name__
    #G.printFunctionNames = True
    #cProfile.run('BF_counter.BF_counter(fn,BF,k,G,pfn,printProgress,startAtLine)')
    BF_counter.BF_counter(fn,BF,k,G,pfn,printProgress,startAtLine,sizeOfGenome,maxCov)
"""

def sim_BF_counter_naive(fn,k,BF,G_naive,pfn=False,printProgress=True):
    print sim_BF_counter_naive.__name__
    #cProfile.run("BF_counter.BF_counter_naive(fn,BF,k,G_naive,pfn,printProgress)")
    BF_counter.BF_counter_naive(fn,BF,k,G_naive,pfn,printProgress)

def initialize(fn,k,numberOfKmers):
    BF = Bloom.Bloom(0.01,numberOfKmers,pfn=False)
    #BF = pybloom.BloomFilter(numberOfKmers,0.01)
    G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=False)
    #if fileName!="":
    #    print "We read the graph from the file "+fileName
    #    G.createGraphFromFile(fileName)
    return BF,G

def readInput():
    #print sys.argv
    k = int(sys.argv[1])
    read1 = sys.argv[2]
    read2 = sys.argv[3]
    numberOfKmers = int(sys.argv[4])
    whatToRun = int(sys.argv[5])        #0: BF_counter_naive. 1: BF_counter
    sizeOfGenome = int(sys.argv[6])
    maxCov = int(sys.argv[7])
    pfn = False
    printProgress = True
    #startAtLine = 0
    #fileName = ""
    genomeFile = ""
    if len(sys.argv)>8:
        pfn = map((lambda x: {"False":False,"True":True}[x]), [sys.argv[8]])[0]
    if len(sys.argv)>9:
        printProgress = map((lambda x: {"False":False,"True":True}[x]), [sys.argv[9]])[0]
    if len(sys.argv)>10:
        genomeFile  = sys.argv[10]
    if len(sys.argv)>11:
        outDirectory = sys.argv[11]
        
    assert isinstance(k, int)
    assert isinstance(read1,str)
    assert isinstance(read2,str)
    assert isinstance(numberOfKmers,int)
    assert isinstance(whatToRun,int)
    assert isinstance(sizeOfGenome,int)
    assert isinstance(maxCov,int)
    assert isinstance(pfn,bool)
    assert isinstance(printProgress,bool)
    #assert isinstance(startAtLine,int)
    assert isinstance(genomeFile,str)
    assert isinstance(outDirectory,str)
    #assert isinstance(fileName,str)
    fn = [read1,read2]
    return k,fn,numberOfKmers,whatToRun,sizeOfGenome,maxCov,pfn,printProgress,genomeFile,outDirectory

if __name__ == '__main__':
    #command line syntax:
    #python simData.py k       read1           read2      numberOfKmers whatToRun sizeOfGenome maxCov [pfn] [printProgress] [genomeFile] [outDirectory]
    #python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000        1          70000      10   True         True          t.fa       Output/t

    #Gera fyrst:
    #python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 70000 -1 True
    #Svo:
    #python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 70000 5 True
    start = time.time()
    k,fn,numberOfKmers,whatToRun,sizeOfGenome,maxCov,pfn,printProgress,genomeFile,outDirectory = readInput()

    BF,G = initialize(fn,k,numberOfKmers)
    if whatToRun==0:
        sim_BF_counter_naive(fn,k,BF,G,pfn,printProgress)
    elif whatToRun==1:
        #sim_BF_counter(fn,k,BF,G,pfn,printProgress,startAtLine)
        BF_counter.BF_counter(fn,k,BF,G,pfn,printProgress,0,sizeOfGenome,maxCov,genomeFile,outDirectory)
    else:
        raise Exception("whatToRun must be either 0 or 1")
    end = time.time()

    #Kóði til að athuga í leiðinni hvort grafið sé löglegt hér fyrir neðan sem ég kommenta út til að spara hraða:

    #print "Checking whether the graph is legal"
    #G.printIsLegal = True
    #B1 = G.isLegalDBG()
    #if B1:
    #    print "The Graph is legal"
    #else:
    #    print "The Graph is not legal"
    helpers.printResultsToFile(BF,G,outFolder="defaultOutFolder",timeInSeconds=int(end-start))
    #print "We have finished writing the Graph and all info to a file"
    #print "Checking whether the graph is legal"
    #B2 = G.isLegalDBG()
    #if B2:
    #    print "The Graph is legal"
    #else:
    #    print "The Graph is not legal"
    print 'finished simData.py in '+str(int(end-start))+' seconds'
    