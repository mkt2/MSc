#coding:utf8
import Graph
import helpers
import BF_counter
import Bloom

if __name__ == "__main__":
    #Inputs
    k = 31
    outDir = "Output/t"
    inDir = "Input/t"
    fn = [inDir+"/r1.fastq",inDir+"/r2.fastq"]
    p = 0.01
    numberOfKmers = 8000000

    #Create G
    print "Starting on G"
    G = Graph.Graph(k,al=False)
    G.createGraphFromFile(outDir+"/G.txt")

    #Create G_naive
    print "Starting on G_naive"
    G_naive = Graph.Graph(k,al=False)
    BF = Bloom.Bloom(p,numberOfKmers,pfn=False)
    BF_counter.BF_counter_naive(fn,BF,k,G_naive,pfn=False,printProgress=False)
    print "len(G)      ",len(G), len(G.kmers)
    print "len(G_naive)",len(G_naive), len(G_naive.kmers)

    #Save the graphs to files
    print "saving the graphs to files"
    G.saveAs_GFA_toFile(outDir+"/G.gfa")
    G_naive.printToFile(outDir+"/G_naive.txt")
    G_naive.saveAs_GFA_toFile(outDir+"/G_naive.gfa")
    
