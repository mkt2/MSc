#coding:utf8
import helpers
import sys

"""
def selectGenome(genomeName):
    if genomeName=="t":
        filters = [ \
        (5,float('inf')), \
        (10,12),(10,15),(10,float('inf')), #fæ error í merge fyrir (10,11)\
		(15,16),(15,20),(15,float('inf')), \
		(20,float('inf')), \
		(30,float('inf')), \
        (15,10),(20,15),(30,20), \
		(float('inf'),float('inf'))]
    elif genomeName=="sa":
        filters = [ \
        (15,16),(15,20),(15,float('inf')), \
		(20,21),(20,25),(20,float('inf')), \
		(30,float('inf')), \
        (15,10),(20,15),(30,20), \
		(float('inf'),float('inf'))]
    return filters
"""

if __name__ == "__main__":
    #Inputs
    genomeName = sys.argv[1]

    #Define some stuff
    outDir = "Output/"+genomeName
    #filters = selectGenome(genomeName)
    covDict = helpers.readCovDictFromFile(fileName=outDir+"/covDict.txt")
    filters = covDict.keys()
    filters.remove((float('inf'),float('inf')))
    filters.sort(key=lambda tup: (tup[0],tup[1]), reverse=False)
    #filters.sort( key=lambda x: (x.split("_")[1][1:].split("}")[0].split(",")[0] , x.split("_")[1][1:].split("}")[0].split(",")[1]) )
    #print filters
    numGraphs = len(filters)*2+1
    temp = []
    temp_tex = []
    texNames = []
    for i, (MAC,MSC) in enumerate(filters):
        gn = helpers.createGraphName(MAC,MSC)
        gn_tex = helpers.createGraphName_tex(MAC,MSC,True)
        texNames.append(gn_tex)
        filters[i] = gn
        temp.append(gn+"d")
        temp_tex.append(gn_tex+"d")
    filters.append("G_inf_inf")
    texNames.append("$G_{\\infty,\\infty}$")
    for i, t in enumerate(temp):
        filters.append(t)
        texNames.append(temp_tex[i])
    outFile = outDir+"/percTable.txt"
    #for t in texNames:
    #    print t
    """
    texNames.sort( key=lambda x: (x.split("_")[1][1:].split("}")[0].split(",")[0] , x.split("_")[1][1:].split("}")[0].split(",")[1]) )
    print "----------------------------"
    for t in texNames:
        print t
    for i,t in enumerate(texNames):
        if t[-1]=="d":
            texNames[i] = "$G-"+t[1:-1]
    print "----------------------------"
    for t in texNames:
        print t
    """



    #Read the data into our matrix
    numValues = 7
    matrix = [[]]*numGraphs
    #print matrix
    for lineNr,gn in enumerate(filters):
        totals, percs = helpers.readTotalsAndPercFromFile(outDir+"/"+gn+"_tot_perc.csv")
        for i, p in enumerate(percs):
            percs[i] = "{0:.1f}".format(p*100)+"\%"
        if texNames[lineNr][-1]=="d":
            texNames[lineNr] = "$G-"+str(texNames[lineNr][1:-1])
        matrix[lineNr] = [texNames[lineNr]] + [sum(totals)] + percs                         #<----Breyta þessari línu til að sleppa numKmers dálknum
    f = open(outFile, 'w')
    f.write(" & numKmers & isolated & tips & primary & secondary & none \\\\\n")            #<----Breyta þessari línu til að sleppa numKmers dálknum
    f.write("\\hline\n")

    for c in range(0,numGraphs):
        str1 = ' & '.join(str(e) for e in matrix[c])
        f.write(str1+" \\\\\n")
    f.close()