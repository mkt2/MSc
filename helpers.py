#coding:utf8
import numpy as np
#import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
#import matplotlib.pyplot
import dbg
import collections
from math import log
import fileinput
from shutil import copyfile
import os.path

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

def printFigureFromFile(inputFile,sizeOfGenome,title="",outFolder=""):
    f = open(inputFile,"r")
    lineNumber = []
    cov = []
    kmersInGraph = []
    kmersInBF = []
    kmersInGraph_s = []
    kmersInBF_s = []
    for line in f:
        temp = line.split(",")
        L = len(temp)
        assert (L==4) or (L==6)
        lineNumber.append(int(temp[0]))
        cov.append(int(temp[1]))
        kmersInGraph.append(int(temp[2]))
        kmersInBF.append(int(temp[3]))
        if L==6:
            kmersInGraph_s.append(int(temp[4]))     #k-mers in the Graph using a BF-stopper
            kmersInBF_s.append(int(temp[4]))
    
    #for i, c in enumerate(cov):
    #    print i,c,kmersInGraph[i]

    fig1, (ax1,ax2) = plt.subplots(2,1,sharex=False)
    ax1.axhline(y=sizeOfGenome,label="#k-mers in genome")
    ax1.plot(cov,kmersInBF,"k-",label="#k-mers in BF")
    ax1.plot(cov,kmersInGraph,"k--",label='#k-mers in G')
    if L==6:
        ax1.plot(cov,kmersInGraph_s,"r--",label='#k-mers in G using BF-stopper')
    x1,x2,y1,y2 = ax1.axis()
    #ax1.axis((x1,max(cov),0,700000))
    ax1.axis((x1,max(cov),0,max(kmersInGraph)*2))
    if title=="":
        ax1.set_title("#k-mers as a function of cov (log scale below)")
    else:
        ax1.set_title(title)
    ax2.set_xlabel("Coverage")
    ax1.set_ylabel("#k-mers")
    ax2.set_ylabel("ln(#k-mers)")
    ax1.grid()

    kmersInGraph = [0 if x==0 else log(x) for x in kmersInGraph]
    kmersInBF = [0 if x==0 else log(x) for x in kmersInBF]
    ax2.axhline(y=log(sizeOfGenome),label="#k-mers in genome")
    ax2.plot(cov,kmersInBF,"k-",label="#k-mers in BF")
    ax2.plot(cov,kmersInGraph,"k--",label='#k-mers in G')
    if L==6:
        kmersInGraph_s = [0 if x==0 else log(x) for x in kmersInGraph_s]
        ax2.plot(cov,kmersInGraph_s,"r--",label='#k-mers in G using BF-stopper')
    x1,x2,y1,y2 = ax2.axis()
    ax2.axis((x1,max(cov),y1,y2))
    ax2.legend(loc='lower right')
    ax2.grid()
    if outFolder=="":
        fig1.savefig('Output/defaultOutFolder/figure_1.png', bbox_inches='tight')
    else:
        fig1.savefig(outFolder+"/figure_1.png", bbox_inches='tight')
    f.close()

def printFigureFromFile2(inputFile,title="",outFolder=""):
    f = open(inputFile,"r")
    lineNumber = []
    cov = []
    G_ratio = []
    BF_ratio = []
    G_ratio_s = []
    BF_ratio_s = []
    for line in f:
        temp = line.split(",")
        L = len(temp)
        assert (L==4) or (L==6)
        lineNumber.append(int(temp[0]))
        cov.append(int(temp[1]))
        G_ratio.append(temp[2])
        BF_ratio.append(temp[3])
        if L==6:
            G_ratio_s.append(temp[4])     #ratio of k-mers in the Graph using a BF-stopper
            BF_ratio_s.append(temp[4])
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    ax2.axhline(y=1,label="100% ratio")
    ax2.plot(cov,BF_ratio,"k-",label="ratio in BF")
    ax2.plot(cov,G_ratio,"k--",label='ratio in G')
    if L==6:
        ax2.plot(cov,G_ratio_s,"r--",label='ratio in G using BF-stopper')
    x1,x2,y1,y2 = ax2.axis()
    ax2.axis((x1,x2,0,1.1))

    if title=="":
        plt.title("#Ratio of k-mers in G and BF")
    else:
        plt.title(title)
    plt.xlabel("Coverage")
    plt.ylabel("Ratio of correct k-mers")
    ax2.legend(loc='lower right')
    ax2.grid()
    if outFolder=="":
        fig2.savefig('Output/defaultOutFolder/figure_2.png', bbox_inches='tight')
    else:
        fig2.savefig(outFolder+"/figure_2.png", bbox_inches='tight')
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

def changeFile(file,newCol1,newCol2,outFolder=""):
    if outFolder=="":
        outFolder = "Output/defaultOutFolder/BF_stopper"
    newFile = outFolder+"/"+os.path.basename(file)
    print newFile
    copyfile(file, newFile)
    for i,line in enumerate(fileinput.input(newFile, inplace=True)):
        assert line.count(",")==3
        print "%s,%s,%s" % (line.strip(),str(newCol1[i]),str(newCol2[i])+"\n"),

if __name__ == "__main__":
    #changeFile("Output/testfile.csv",[0]*21,[0]*21)
    printFigureFromFile("Output/defaultOutFolder/figureData.csv",sizeOfGenome=70000,title="#k-mers as a function of cov (log scale below)")
    printFigureFromFile2("Output/defaultOutFolder/figureData2.csv",title="#Ratio of k-mers in G and BF")
