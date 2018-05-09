#coding:utf8
import helpers
import sys

if __name__ == "__main__":
    #Inputs
    genomeName = sys.argv[1]

    #Define some stuff
    outDir = "Output/"+genomeName
    covDict = helpers.readCovDictFromFile(fileName=outDir+"/covDict.txt")
    filters = covDict.keys()
    filters.remove((float('inf'),float('inf')))
    filters.sort(key=lambda tup: (tup[0],tup[1]), reverse=False)
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
    #texNames geymir nú nöfnin á öllum gröfunum

    #Read the data into our matrix
    #numValues = 7
    matrix = [[]]*numGraphs
    #print matrix
    for lineNr,gn in enumerate(filters):
        totals, percs = helpers.readTotalsAndPercFromFile(outDir+"/"+gn+"_tot_perc.csv")
        assert(abs(sum(percs)-1)<0.1), str(percs)
        for i, p in enumerate(percs):
            percs[i] = "{0:.1f}".format(p*100)+"\%"
        if texNames[lineNr][-1]=="d":
            texNames[lineNr] = "$G-"+str(texNames[lineNr][1:-1])
        matrix[lineNr] = [texNames[lineNr]] + [format(sum(totals),',')] + percs
    f = open(outFile, 'w')
    f.write(" & numKmers & isolated & tips & secondary & genomic & partial & complex \\\\\n")
    f.write("\\hline\n")

    for c in range(0,numGraphs):
        str1 = ' & '.join(str(e) for e in matrix[c])
        f.write(str1+" \\\\\n")
    f.close()