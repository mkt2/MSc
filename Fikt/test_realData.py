#coding:utf8
import Graph
from compareGraphs import *
from BF_counter import *
import Bloom
from pybloom import BloomFilter


#--------------------------------------------------------------------------
#-------Various ways of running "the code" on Staphylococcus_aureus--------
#--------------------------------------------------------------------------
#Prófa hvort pybloom hegði sér eðlilega með því að búa til BF
#úr öllum segmentum í frag_1.fastq
def bigData_1(stopAfter=-1):
    print "bigData_1(stopAfter="+str(stopAfter)+")"
    """
    Printing the BF:
    p:          0.01
    m:          892888906


    Printing pybloom.bitArray:
    Number of ones:  242305870
    Number of zeros: 650583036
    Ratio:           0.271372920384
    Added 40383349 kmers

    Time in seconds: 549
    9 minutes and 9 second
    """
    fn = ["Input/Staphylococcus_aureus/frag_1.fastq"]
    k = 31
    numberOfKmers = 38814270*1.2
    BF = BloomFilter(capacity=numberOfKmers, error_rate=0.01)
    for segment,lineNr in generateSegmentsFrom_fq(fn,k,BF,pfn=True,printProgress=False):
        if stopAfter!=-1:
            if lineNr==stopAfter:
                break
    print BF
    BF.print_bitarray()

#Prófa hvort Bloom hegði sér eðlilega með því að búa til BF
#úr öllum segmentum í frag_1.fastq
def bigData_2(stopAfter=-1):
    print "bigData_2(stopAfter="+str(stopAfter)+")"
    """
    generateSegments_from_fastq(fn=['Input/Staphylococcus_aureus/frag_1.fastq', 'Input/Staphylococcus_aureus/frag_2.fastq'], k=31)

    Printing the BF (except ratio):
      p:          0.01
      k:          10
      m:          372037044


    Printing the 0/1 ratio in the BF bitarray:
      Number of ones:  260178608
      Number of zeros: 111858436
      Ratio:           0.699335219963

    Time in seconds: 1549
    25 minutes and 49 second
    """
    fn = ["Input/Staphylococcus_aureus/frag_1.fastq","Input/Staphylococcus_aureus/frag_2.fastq"]
    k = 31
    numberOfKmers = 38814270
    BF = Bloom.Bloom(0.01,numberOfKmers,pfn=True)
    count = 0
    for segment in generateSegmentsFrom_fq(fn,k,BF,pfn=True,printProgress=False):
        if stopAfter!=-1:
            count += 1
            if count==stopAfter:
                break
    print BF
    BF.print_bitarray()

#Prófa hvort BF_counter hegði sér eðlilega með því að búa til BF með Bloom.py úr öllum segmentum
#í frag_1.fq og bæta öllum segmentum sem innihalda eingöngu kmera sem koma fyrir tvisvar
#í graf
def bigData_3(outFolder,startAtLine=0):
    print "bigData_3(outFolder="+str(outFolder)+", startAtLine="+str(startAtLine)+")"
    assert int(startAtLine)==0
    fn = ["Input/Staphylococcus_aureus/frag_1.fastq","Input/Staphylococcus_aureus/frag_2.fastq"]
    k = 31
    numberOfKmers = 38814270
    BF = Bloom.Bloom(0.01,numberOfKmers,pfn=True)
    G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
    G = BF_counter(fn,BF,k,G,printProgress=False,pfn=True,startAtLine=startAtLine,outFolder=outFolder)
    print "What you have created so far will be added to the folder " \
    +'Output/'+str(outFolder)+'\n'

#Before:    We have stopped bigData_3 using ctrl+c
#           and want to continue from where we stopped
def bigData_3_continue(outFolder,startAtLine=0,inFolder="Output/bigData_3"):
    print "bigData_3b"
    fn = ["Input/Staphylococcus_aureus/frag_1.fastq","Input/Staphylococcus_aureus/frag_2.fastq"]
    k = 31
    numberOfKmers = 38814270
    BF = Bloom.Bloom(0.01,numberOfKmers,pfn=True)
    G = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
    G.createGraphFromFile("Output/bigData_3/G_from_BF_counter.txt")
    G = BF_counter(fn,BF,k,G,printProgress=False,pfn=True,startAtLine=startAtLine,outFolder=outFolder)
    print "What you have created so far will be added to the folder " \
    +'Output/'+str(outFolder)+'\n'

#Athuga hvort generateSegmentsFrom_fq búi til eðlilegan fjölda af kmerum
#Bera saman við töluna frá Páli
def bigData_4():
    print "bigData_4"
    fn = ["Input/Staphylococcus_aureus/frag_1.fastq"]
    k = 31
    numberOfKmers = 38814270*1.2
    BF = Bloom.Bloom(0.01,numberOfKmers,pfn=True)
    kmerDict = createKmerDict(fn,BF,k,pfn=True,printProgress=False)
    print len(kmerDict)

#Prófa hvort BF_counter_naive hegði sér eðlilega með því að gera það sama og fyrir
#BF_counter í bigData_3
def bigData_5():
    print "bigData_5"
    fn = ["Input/Staphylococcus_aureus/frag_1.fastq","Input/Staphylococcus_aureus/frag_2.fastq"]
    k = 31
    numberOfKmers = 38814270
    BF = Bloom.Bloom(0.01,numberOfKmers,pfn=True)
    print BF
    BF.print_bitarray()
    G_naive = BF_counter_naive(fn,BF,k,pfn=False,ps=False,printProgress=False)
    print BF
    BF.print_bitarray()


    
"""
if __name__ == '__main__':
    start = time.time()
    bigData_3()
    end = time.time()
    time = int(end-start)
    printTime(time)
"""

#Command line skipun til að runna mörg föll í einu og skila útkomu í skrá:
# bash runBigData.sh > Output/out.txt
if __name__ == '__main__':
    pass