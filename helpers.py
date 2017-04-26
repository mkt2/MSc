#coding:utf8
import numpy as np
#import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
import dbg
import collections

def returnTime(timeInSeconds):
    s = timeInSeconds%60
    m = timeInSeconds/60
    return "Time in seconds: "+str(timeInSeconds)+"\n" \
    + str(m)+" minutes and "+str(s)+" seconds\n"
    assert m*60+s==timeInSeconds

def printResultsToFile(BF,G,outFolder="defaultOutFolder",timeInSeconds=-1,nextLine=-1):
    print "printResultsToFile(BF, G, outFolder="+str(outFolder)+", timeInSeconds="+str(timeInSeconds)+")"
    assert isinstance(outFolder, str)
    assert isinstance(timeInSeconds, int)
    f = open('Output/'+outFolder+"/G.txt", 'w')
    b = open('Output/'+outFolder+'/other_info.txt', 'w')
    G.printToFile(f)
    f.close()
    b.write(str(BF))
    b.write(BF.bitarray_str())
    if nextLine!=-1:
        b.write("\nWe were about to add the segment in line nr "+str(nextLine)+" from the input file\n")
    if timeInSeconds!=-1:
        b.write("\n"+str(returnTime(int(timeInSeconds))))
    b.close()

def printFigureFromFile(fileName):
    f = open(fileName,"r")
    lineNumber = []
    cov = []
    kmersInGraph = []
    kmersInBF = []
    for line in f:
        temp = line.split(",")
        lineNumber.append(int(temp[0]))
        cov.append(int(temp[1]))
        kmersInGraph.append(int(temp[2]))
        kmersInBF.append(int(temp[3]))
    plt.axhline(y=70000,label="The correct answer")
    plt.plot(cov,kmersInBF,"k-",label="k-mers in BF")
    plt.plot(cov,kmersInGraph,"k--",label='k-mers in G')
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,0,500000))
    plt.title("#k-mers as a function of cov")
    plt.xlabel("Coverage")
    plt.ylabel("Number of k-mers")
    plt.legend(loc='upper right')
    #plt.show()
    plt.savefig('Output/defaultOutFolder/figure_1.png', bbox_inches='tight')
    f.close()

def printFigureFromFile2(fileName):
    f = open(fileName,"r")
    lineNumber = []
    cov = []
    G_ratio = []
    BF_ratio = []
    for line in f:
        temp = line.split(",")
        lineNumber.append(int(temp[0]))
        cov.append(int(temp[1]))
        G_ratio.append(temp[2])
        BF_ratio.append(temp[3])
    plt.axhline(y=1,label="100% ratio")
    plt.plot(cov,BF_ratio,"k-",label="ratio in BF")
    plt.plot(cov,G_ratio,"k--",label='ratio in G')
    x1,x2,y1,y2 = plt.axis()
    plt.axis((x1,x2,0,1.1))
    plt.title("#Ratio of k-mers in G and BF")
    plt.xlabel("Coverage")
    plt.ylabel("Ratio of correct k-mers")
    plt.legend(loc='lower right')
    plt.savefig('Output/defaultOutFolder/figure_2.png', bbox_inches='tight')
    #plt.show()
    f.close()

def createKmerDictFrom_fa(fileName,k):
    f = open(fileName,"r")
    l = f.readline()
    l = f.readline().strip()
    kmerdict = collections.defaultdict(int)
    for km in dbg.kmers(l,k):
        kmerdict[km] = 1
    for km in dbg.kmers(dbg.twin(l),k):
        kmerdict[km] = 1
    return kmerdict

#Input:
#   kmerdict:  A dict storing all kmers actually occurring in the genome. Also stores twins
#   G:         Our Graph
#   BF:        Our Bloom filter
#Returns:
#   The ratio of kmers in the genome that occur in G
#   The ratio of kmers in the genome that occur in BF
def ratioInGenome(kmerdict,G,BF):
    G_count = 0
    BF_count = 0
    for km in kmerdict:
        if km in G.kmers:
            G_count += 1
        if km in BF:
            BF_count += 1
    if len(G)==0:
        G_ratio = 0
    else:
        G_ratio = float(G_count) / len(kmerdict)
    if len(BF)==0:
        BF_ratio = 0
    else:
        BF_ratio = float(BF_count) / len(kmerdict)
    #print "Ratios:",G_ratio,BF_ratio
    return G_ratio,BF_ratio

if __name__ == "__main__":
    printFigureFromFile("Output/defaultOutFolder/figureData.csv")
    printFigureFromFile2("Output/defaultOutFolder/figureData2.csv")
