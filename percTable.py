#coding:utf8
import helpers
import sys

def selectGenome(genomeName):
    if genomeName=="t":
        filters = [ \
        (5,float('inf')), \
		(10,float('inf')), \
		(15,float('inf')), \
		(20,float('inf'))]
    elif genomeName=="sa":
        filters = [ \
        (15,float('inf')), \
		(20,float('inf')), \
		(30,float('inf'))]
    return filters

if __name__ == "__main__":
    #Inputs
    genomeName = sys.argv[1]

    #Define some stuff
    outDir = "Output/"+genomeName
    filters = selectGenome(genomeName)
    numGraphs = len(filters)*2+1
    temp = []
    temp_tex = []
    texNames = []
    for i, (MAC,MSC) in enumerate(filters):
        gn = helpers.createGraphName(MAC,MSC,False,False)
        gn_tex = helpers.createGraphName(MAC,MSC,False,True)
        texNames.append(gn_tex)
        filters[i] = gn
        temp.append(gn+"d")
        temp_tex.append(gn_tex+"d")
    filters.append("G_inf_inf")
    texNames.append("G$\\infty \\infty$")
    for i, t in enumerate(temp):
        filters.append(t)
        texNames.append(temp_tex[i])
    outFile = outDir+"/percTable.txt"
    print filters
    print texNames

    #Read the data into our matrix
    numValues = 7
    matrix = [[]]*numGraphs
    #matrix = [[]]*numValues
    #print matrix
    for lineNr,gn in enumerate(filters):
        totals, percs = helpers.readTotalsAndPercFromFile(outDir+"/"+gn+"_tot_perc.csv")
        for i, p in enumerate(percs):
            percs[i] = "{0:.1f}".format(p*100)+"\%"
        matrix[lineNr] = [texNames[lineNr]] + [sum(totals)] + percs

    f = open(outFile, 'w')
    f.write(" & numContigs & isolated & tips & primary & secondary & none \\\\\n")
    f.write("\\hline\n")
    list1 = ['a5','1234','c']
    str1 = ''.join(list1)
    print str1
    for c in range(0,numGraphs):
        #print matrix[c]
        str1 = ' & '.join(str(e) for e in matrix[c])
        f.write(str1+" \\\\\n")
        #print str1
        #s = ""
        #for r in range(numValues):
        #    s += str(matrix[c][r])
        #f.write(str(matrix[c][r]))
    """for r in range(0,numValues):
        if r==0:
            s = " "
        if r==1:
            s = "numContigs"
        elif r==2:
            s="isolated"
        elif r==3:
            s="tips"
        elif r==4:
            s="primary"
        elif r==5:
            s="secondary"
        elif r==6:
            s="none"
        for c in range(0,numGraphs):
            s += " & "+str(matrix[c][r])
        s += " \\\\\n"
        f.write(s)
        if r==0:
            f.write("\\hline\n")
    """
    f.close()