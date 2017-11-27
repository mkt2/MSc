#coding:utf8
import helpers
import numpy as np
import matplotlib.pyplot as plt
import dbg
import collections
from math import log
import os.path
import Graph, Bloom
import time
import sys

alphabet = ["A","C","G","T","N"]

#Creates a k-merdict from a .fa or .fasta file
def createKmerDictFromGenomeFile(k,genomeFile):
    kmersInGenome = collections.defaultdict(int)
    h = open(genomeFile, "rU")
    genomeFileExtension = os.path.splitext(genomeFile)[1]
    if genomeFileExtension==".fa":
        genome = h.readline()
        genome = h.readline().rstrip('\n')
        h.close()
        for km in dbg.kmers(genome,k):
            rep_km = min(km,dbg.twin(km))
            kmersInGenome[rep_km] += 1
    elif genomeFileExtension==".fasta":
        genome = h.readline()
        for line in h:
            line = line.strip()
            line = line if all([c in alphabet for c in line]) else ""
            tl = dbg.twin(line)
            for km in dbg.kmers(line,k):
                rep_km = min(km,dbg.twin(km))
                kmersInGenome[rep_km] += 1
    numKmersInGenome = len(kmersInGenome)
    return kmersInGenome,numKmersInGenome

#Creates a k-merdict from two .fastq files
def createKmerDictFromReadFiles(k,readFiles):
    kmersInReads = collections.defaultdict(int)
    numReadsPerFile = 0
    for f in readFiles:
        h = open(f, "rU")
        for lineNr,line in enumerate(h):
            if (lineNr%4 == 1):
                numReadsPerFile += 1
                segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
                for s in segments:
                    for km in dbg.kmers(s,k):
                        rep_km = min(km,dbg.twin(km))
                        kmersInReads[rep_km] += 1
        h.close()
    numReadsPerFile = numReadsPerFile/2
    numKmersInReads = len(kmersInReads)
    return kmersInReads,numKmersInReads,numReadsPerFile

def writeGenomeInfoToFile(genomeName,k,p=0.01):
    start = time.time()
    #Helper function 1
    def doGenomeStuff(genomeName,k,genomeFile):
        #Create a dict for all the k-mers in the genome:
        kmersInGenome,numKmersInGenome = createKmerDictFromGenomeFile(k,genomeFile)
        kmer_file = "Output/"+str(genomeName)+"/kmers_genome.txt"
        f = open(kmer_file, 'w')
        for km in kmersInGenome:
            f.write(str(km)+"\n")
        f.close()
        return numKmersInGenome

    #Helper function 2
    def doReadStuff(genomeName,k,readFiles):
        #Create a dict for all the k-mers in the reads:
        kmersInReads,numKmersInReads,numReadsPerFile = createKmerDictFromReadFiles(k,readFiles)
        kmer_file = "Output/"+str(genomeName)+"/kmers_reads.txt"
        f = open(kmer_file, 'w')
        for km in kmersInReads:
            f.write(str(km)+"\n")
        f.close()
        #Count k-mers occuring twice or more often
        #using an actual counting of k-mers
        numKmersInReads_twice = 0
        for km,num in kmersInReads.iteritems():
            if num > 1:
                numKmersInReads_twice += 1
        return numReadsPerFile,numKmersInReads,numKmersInReads_twice

    #Helper function 3
    def readStuffUsingBF(genomeName,k,readFiles,p,numKmersInReads):
        #Create a dict for all the k-mers in the reads
        #occurring at least twice. Estimated using a BF:
        BF = Bloom.Bloom(p,n=(numKmersInReads*7),pfn=False)
        kmersInReads_twice_BF = collections.defaultdict(int)
        for f in readFiles:
            h = open(f, "rU")
            for lineNr,line in enumerate(h):
                if (lineNr%4 == 1):
                    segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
                    for s in segments:
                        for km in dbg.kmers(s,k):
                            rep_km = min(km,dbg.twin(km))
                            if not rep_km in BF:
                                BF.add(rep_km)
                            else:
                                kmersInReads_twice_BF[rep_km] += 1
            h.close()
        numKmersInReads_twice_BF = len(kmersInReads_twice_BF)
        BF_ratio_1_vs_0,cZero,cOne = BF.computeRatio()
        return numKmersInReads_twice_BF, BF_ratio_1_vs_0


    #----------------Main function part starts here-------:

    #Select a genome:
    if genomeName=="t":
        readFiles = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
        genomeFile = "Input/t/t.fa"
    elif genomeName=="sa":
        readFiles = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
        genomeFile = "Input/Staphylococcus_aureus/genome.fasta"

    #Let the helper functions do the hard work:
    numKmersInGenome = doGenomeStuff(genomeName,k,genomeFile)
    numReadsPerFile,numKmersInReads,numKmersInReads_twice = doReadStuff(genomeName,k,readFiles)
    numKmersInReads_twice_BF, BF_ratio_1_vs_0 = readStuffUsingBF(genomeName,k,readFiles,p,numKmersInReads)

    #Find the number of k-mers in each read:
    f = readFiles[0]                                #XX throw this line away
    numReadsPerFile2  = int(helpers.file_len(f))/4          #XX throw this line away
    assert(numReadsPerFile==numReadsPerFile2)       #XX throw this line away
    h = open(f, "rU")
    line = h.readline()
    line = h.readline().rstrip('\n')
    numKmersPerRead = len(list(dbg.kmers(line,k)))
    h.close()
    preprocessTime = helpers.returnTime(time.time()-start)
    

    #Print the results:
    #(we don't print the dicts for now)
    print "genomeName:               %9s" % genomeName
    print "k:                        %9s" % k
    print "numKmersInGenome:         %9i" % numKmersInGenome
    print "numReadsPerFile:          %9i" % numReadsPerFile
    print "numKmersPerRead:          %9i" % numKmersPerRead
    print "numKmersInReads:          %9i" % numKmersInReads
    print "numKmersInReads_twice:    %9i" % numKmersInReads_twice
    print "numKmersInReads_twice_BF: %9i" % numKmersInReads_twice_BF
    print "BF_ratio_1_vs_0           %9.2f" % BF_ratio_1_vs_0
    print "preprocessTime            %9s" % preprocessTime

    #Print the results to files:
    #All above to an info file
    infoFile = "Output/"+str(genomeName)+"/genome_info.csv"
    f = open(infoFile, 'w')
    f.write("genomeName,"+str(genomeName)+"\n")
    f.write("k,"+str(k)+"\n")
    f.write("numKmersInGenome,"+str(numKmersInGenome)+"\n")
    f.write("numReadsPerFile,"+str(numReadsPerFile)+"\n")
    f.write("numKmersPerRead,"+str(numKmersPerRead)+"\n")
    f.write("numKmersInReads,"+str(numKmersInReads)+"\n")
    f.write("numKmersInReads_twice,"+str(numKmersInReads_twice)+"\n")
    f.write("numKmersInReads_twice_BF,"+str(numKmersInReads_twice_BF)+"\n")
    f.write("BF_ratio_1_vs_0,"+str(BF_ratio_1_vs_0)+"\n")
    f.write("preprocessTime,"+str(preprocessTime)+"\n")
    f.close()

if __name__ == "__main__":
    k = int(sys.argv[1])
    genomeName = sys.argv[2]
    writeGenomeInfoToFile(genomeName,k)