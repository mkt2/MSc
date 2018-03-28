#coding:utf8
import helpers
import sys

def selectGenome(genomeName):
    if genomeName=="t":
        maxAddCovs = [5, 10, 15, 20]
    elif genomeName=="sa":
        maxAddCovs = [15, 20, 30]
    return maxAddCovs

if __name__ == "__main__":
    #Inputs
    genomeName = sys.argv[1]

    #Define some stuff
    outDir = "Output/"+genomeName
    maxAddCovs = selectGenome(genomeName)
    numCols = len(maxAddCovs)*2+1
    temp = []
    for maxAddCov in maxAddCovs:
        temp.append(str(maxAddCov)+"d")
    maxAddCovs.append(float('inf'))
    for t in temp:
        maxAddCovs.append(t)
    outFile = outDir+"/percTable.txt"

    #Read the data into our matrix
    numRows = 7
    matrix = [[]]*numCols
    #print matrix
    for colNr,maxAddCov in enumerate(maxAddCovs):
        totals, percs = helpers.readTotalsAndPercFromFile(outDir+"/G"+str(maxAddCov)+"_tot_perc.csv")
        for i, p in enumerate(percs):
            percs[i] = "{0:.1f}".format(p*100)+"\%"
        matrix[colNr] = ["G"+str(maxAddCov)] + [sum(totals)] + percs

    f = open(outFile, 'w')
    for r in range(0,numRows):
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
        for c in range(0,numCols):
            s += " & "+str(matrix[c][r])
        s += " \\\\\n"
        f.write(s)
        if r==0:
            f.write("\\hline\n")
    f.close()