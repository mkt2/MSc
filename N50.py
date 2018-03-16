#coding:utf8
import Graph
import helpers
import sys

def findN50(G,k):
    #Find the sum of all the contigs in Gx
    S = G.contigSum()
    S2 = S/2    #compute half the sum to speed things up later

    #Sort the contigs in Gx by cLen(c)
    #   IDand_c will be a list of tuples of the form:
    #       (c_ID, c, cLen(c))
    #       sorted in ascending order
    IDand_c = zip(G.contigs.keys(),G.contigs.values())
    for i,t in enumerate(IDand_c):
        IDand_c[i] = (t[0],t[1][0],helpers.cLen(t[1][0],k))
    IDand_c.sort(key=lambda tup: tup[2])
    #for i in IDand_c:
    #    print i

    #Split the contigs in Gx into two groups
    #   Group B is the largest contigs in the sorted Gx
    #   Group T is the smallest contigs in the sorted Gx
    #   Sum(B)<S/2
    #   Sum(T) = S - Sum(B)

    #Find the first ID in Gx.contigs
    lastTuple = IDand_c[0]
    #print "first last_Tuple:",lastTuple

    #Find the contig N50
    currS = 0
    for c_ID,c,cL in IDand_c:
        new_currS = currS + cL
        if new_currS < S2:
            currS = new_currS
        else:
            N50_ID = lastTuple[0]
            N50 = lastTuple[2]
            N50_c = lastTuple[1]
            break
        lastTuple = (c_ID,c,cL)
        
    #print "The contig sum:", S
    #print "Half the contig sum:", S2
    #print "N50_ID:",N50_ID
    #print "N50:",N50
    #print "N50_c:",N50_c

    return N50,N50_ID,N50_c

def printN50ToFile(MCL,N50L,outDir):
    fileName = outDir+"/N50.txt"
    f = open(fileName, 'w')
    for i in range(0,len(MCL)):
        f.write(str(MCL[i])+" & "+str(N50L[i])+" & \\\\\n")
    f.close()

if __name__ == "__main__":
    #Inputs
    genomeName = sys.argv[1]
    k = int(sys.argv[2])

    #Define some stuff
    outDir = "Output/"+genomeName
    maxCovs = [5, 10, 15, 20, 30, -1]
    N50_list = [-2]*len(maxCovs)

    #Find N50 for each maxCov
    for i,maxCov in enumerate(maxCovs):
        Gx = Graph.Graph(k,al=False)
        if maxCov==-1:
            fileName = "/G.txt"
        else:
            fileName = "/G"+str(maxCov)+".txt"
        Gx.createGraphFromFile(outDir+fileName)
        N50,N50_ID,N50_c = findN50(Gx,k)
        N50_list[i] = N50

    printN50ToFile(maxCovs,N50_list,outDir)
    #for i, maxCov in enumerate(maxCovs):
    #    print maxCov,N50_list[i]