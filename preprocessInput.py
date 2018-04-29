#coding:utf8
import helpers
import dbg
import collections
import Bloom
import time
import sys

#Helper function 1. Creates a dict for all the k-mers in the genome and saves them to a file
    #                   Returns the number of k-mers in the genome
def doGenomeStuff(genomeName,k,genomeFile):
    kmersInGenome = helpers.createKmerDictFromGenomeFile(k,genomeFile)
    genomeFile = "Output/"+str(genomeName)+"/kmers_genome.txt"
    #print "Printing to", genomeFile
    kmer_file = genomeFile
    f = open(kmer_file, 'w')
    for km in kmersInGenome:
        f.write(str(km)+"\n")
    f.close()
    return len(kmersInGenome)

#Helper function 2. Creates a dict for all the k-mers in the reads and saves them to a file
#                   Returns the number of reads in each file, the number of kmers in each read
#                   and the number of k-mers appearing twice or more often
#                   This function actually counts the k-mers (no estimation)
def doReadStuff(genomeName,k,readFiles):
    kmersInReads = helpers.createKmerDictFromReadFiles(k,readFiles)
    kmer_file = "Output/"+str(genomeName)+"/kmers_reads.txt"
    numReadsPerFile = int(helpers.file_len(readFiles[0]))/4
    numKmersInReads = len(kmersInReads)
    #print "bla",numReadsPerFile,numKmersInReads

    #print "Printing to", kmer_file
    f = open(kmer_file, 'w')
    for km in kmersInReads:
        f.write(str(km)+"\n")
    f.close()

    #Count k-mers occuring twice or more often:
    numKmersInReads_twice = 0
    for km,num in kmersInReads.iteritems():
        if num > 1:
            numKmersInReads_twice += 1

    #Find the number of k-mers in each line in the reads (assume they're all of equal length):
    h = open(readFiles[0], "rU")
    line = h.readline()
    line = h.readline().rstrip('\n')
    numKmersPerRead = len(list(dbg.kmers(line,k)))
    h.close()
    
    return numReadsPerFile,numKmersInReads,numKmersInReads_twice,numKmersPerRead

#Helper function 3. Similar to helper 2 except counts the k-mers using a BF (which means it's an estimation)
#                   Returns the number of k-mers which occur 2 or more times in the reads (using BF estimation) 
#                   and the ratio of 1 vs 0 in the BF
def doReadStuff_usingBF(genomeName,k,readFiles,p,numKmersInReads):
    #Create a dict for all the k-mers in the reads
    #occurring at least twice. Estimated using a BF:
    BF = Bloom.Bloom(p,n=(numKmersInReads*7),pfn=False)
    kmersInReads_twice_BF = collections.defaultdict(int)
    for s in helpers.segments(readFiles,k):
        for km in dbg.kmers(s,k):
            rep_km = min(km,dbg.twin(km))
            if not (rep_km in BF):
                BF.add(rep_km)
            else:
                kmersInReads_twice_BF[rep_km] += 1
    numKmersInReads_twice_BF = len(kmersInReads_twice_BF)
    BF_ratio_1_vs_0 = BF.computeRatio()[0]
    return numKmersInReads_twice_BF, BF_ratio_1_vs_0

def selectGenome(genomeName):
    if genomeName=="t":
        readFiles = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
        genomeFile = "Input/t/t.fa"
    elif genomeName=="sa":
        readFiles = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
        genomeFile = "Input/Staphylococcus_aureus/genome.fasta"
    else:
        raise Exception("The genomeName must be either 't' or 'sa'!")
    return readFiles, genomeFile

def writeGenomeInfoToFile(genomeName,k,p=0.01):
    start = time.time()
    readFiles, genomeFile = selectGenome(genomeName) 

    #Let the helper functions do the hard work:
    print "Starting on doGenomeStuff"
    numKmersInGenome                                                      = doGenomeStuff(genomeName,k,genomeFile)
    #print numKmersInGenome
    print "Starting on doReadStuff"
    numReadsPerFile,numKmersInReads,numKmersInReads_twice,numKmersPerRead = doReadStuff(genomeName,k,readFiles)
    #print numReadsPerFile,numKmersInReads,numKmersInReads_twice,numKmersPerRead
    print "Starting on doReadStuff_usingBF"
    #print genomeName,k,readFiles,p,numKmersInReads
    numKmersInReads_twice_BF, BF_ratio_1_vs_0                             = doReadStuff_usingBF(genomeName,k,readFiles,p,numKmersInReads)

    #Print the results to files:
    preprocessTime = helpers.returnTime(time.time()-start)
    infoFile = "Output/"+str(genomeName)+"/genome_info.csv"
    variables = [genomeName,k,numKmersInGenome,numReadsPerFile,numKmersPerRead,numKmersInReads,numKmersInReads_twice,numKmersInReads_twice_BF,BF_ratio_1_vs_0,preprocessTime]
    print "Printing to", infoFile
    helpers.writeGenomeInfoToFile(infoFile,variables)

if __name__ == "__main__":
    start = time.time()
    genomeName = sys.argv[1]
    k = int(sys.argv[2])
    writeGenomeInfoToFile(genomeName,k)
    print "preprocessing genome="+str(genomeName)+" took "+str(helpers.returnTime(int(time.time()-start)))