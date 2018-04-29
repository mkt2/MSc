#coding:utf8
import Graph
import helpers
import sys

def findN50(G,k):
    #Find the sum of all the contigs in Gx
    S = G.contigSum()
    S = G.num_bps()
    S2 = S/2    #compute half the sum to speed things up later

    #Sort the contigs in Gx by cLen(c)
    #   IDand_c will be a list of tuples of the form:
    #       (c_ID, c, cLen(c))
    #       sorted in ascending order
    IDand_c = zip(G.contigs.keys(),G.contigs.values())
    for i,t in enumerate(IDand_c):
        IDand_c[i] = (t[0],t[1][0],len(t[1][0]))
        #IDand_c[i] = (t[0],t[1][0],helpers.cLen(t[1][0],k))
    IDand_c.sort(key=lambda tup: tup[2],reverse=True)
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
    numInSum = 0
    for c_ID,c,cL in IDand_c:
        #print c_ID,c,cL
        #break
        new_currS = currS + cL
        if new_currS < S2:
            currS = new_currS
        else:
            #print "entering else statement"
            N50_ID = lastTuple[0]
            N50 = lastTuple[2]
            N50_c = lastTuple[1]
            break
        numInSum+=1
        lastTuple = (c_ID,c,cL)  
    #print "The contig sum:     ", S
    #print "Half the contig sum:", S2
    #print "numInSum:           ", numInSum
    #print "len(IDand_c:        ", len(IDand_c)
    #print "N50_ID:",N50_ID
    #print "N50:",N50
    #print "N50_c:",N50_c
    #print "about to return"
    return N50,N50_ID,N50_c,numInSum

def printN50ToFile(filters,N50L,L50L,outDir):
    fileName = outDir+"/N50.txt"
    f = open(fileName, 'w')
    #for i in range(0,len(MCL)):
    for i, (MAC,MSC) in enumerate(filters):
        gn = helpers.createGraphName_tex(MAC,MSC,True)
        f.write(gn+" & "+str(N50L[i])+" & "+str(L50L[i])+" \\\\\n")
        #if MCL[i]==-1:
        #    f.write("$\\infty$"+" & "+str(N50L[i])+" & "+str(L50L[i])+" \\\\\n")
        #else:
        #    f.write(str(MCL[i])+" & "+str(N50L[i])+" & "+str(L50L[i])+" \\\\\n")
    f.close()

"""
def selectGenome(genomeName):
	if genomeName=="t":
		filters = [ \
        (5,float('inf')), \
		(10,float('inf')), \
		(15,float('inf')), \
		(20,float('inf')), \
		(30,float('inf')), \
		(float('inf'),float('inf'))]
	elif genomeName=="sa":
		filters = [ \
        (15,float('inf')), \
		(20,float('inf')), \
		(30,float('inf')), \
		(float('inf'),float('inf'))]
	else:
		raise Exception("The genomeName must be either 't' or 'sa'!")
	return filters
"""

if __name__ == "__main__":
    #Inputs
    genomeName = sys.argv[1]
    k = int(sys.argv[2])

    #Define some stuff
    outDir = "Output/"+genomeName
    #filters = selectGenome(genomeName)
    covDict = helpers.readCovDictFromFile(fileName=outDir+"/covDict.txt")
    filters = covDict.keys()
    filters.sort(key=lambda tup: (tup[0],tup[1]), reverse=False)
    N50_list = [-2]*len(filters)
    L50_list = [-2]*len(filters)

    #Find N50 for each maxAddCov
    for i,(MAC,MSC) in enumerate(filters):
        Gx = Graph.Graph(k,al=False)
        fileName = helpers.createGraphName(MAC,MSC)+".txt"
        #fileName = "/G"+str(maxAddCov)+".txt"
        Gx.createGraphFromFile(outDir+"/"+fileName)
        N50,N50_ID,N50_c,L50 = findN50(Gx,k)
        N50_list[i] = N50
        L50_list[i] = L50

    printN50ToFile(filters,N50_list,L50_list,outDir)