#coding:utf8
import Graph
import helpers
import BFAdder
import Bloom
import sys

def selectGenome(genomeName):
    if genomeName=="t":
        maxAddCovs = [5, 10, 15, 20, 30,float('inf')]
        fn = ["Input/t/r1.fastq", "Input/t/r2.fastq"]
        numberOfKmers = 8000000
    elif genomeName=="sa":
        maxAddCovs = [5, 10, 15, 20, 30,float('inf')]
        fn = ["Input/Staphylococcus_aureus/frag_1.fastq", "Input/Staphylococcus_aureus/frag_2.fastq"]
        numberOfKmers = 100000000
    else:
        raise Exception("The genomeName must be either 't' or 'sa'!")
    return maxAddCovs, fn, numberOfKmers
    
if __name__ == "__main__":
    #Inputs
    genomeName = sys.argv[1]
    k = int(sys.argv[2])
    maxAddcovs, fn , numberOfKmers = selectGenome(genomeName)
    #k = 31
    outDir = "Output/"+str(genomeName)
    inDir = "Input/"+str(genomeName)
    p = 0.01

    #Create G
    print "Starting on G"
    G = Graph.Graph(k,al=False)
    G.createGraphFromFile(outDir+"/Ginf.txt")
    G.saveAs_GFA_toFile(outDir+"/Ginf.GFA")
    print "contigsInG:", len(G)
    print "kmersInG:  ", G.numKmerPairs()

    #Create G_naive
    print "Starting on G_naive using a BF for filtering singletons:"
    BF = Bloom.Bloom(p,numberOfKmers,pfn=False)
    BFAdder.naive_GFA_nonSingleton_BF(fn,k,BF,outDir+"/Ginf_naive3.GFA",False)
    print "Starting on G_naive using a dict for filtering singletons"
    BFAdder.naive_GFA_nonSingleton(fn,k,outDir+"/Ginf_naive2.GFA")
    

    
