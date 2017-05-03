#coding:utf8
import numpy as np
import matplotlib.pyplot as plt
import dbg
import collections
from math import log
import fileinput
from shutil import copyfile
import os.path
import Graph

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

def createFigureFromFile(inputFile,sizeOfGenome,titles,outFolder,maxCov,genomeName):
    print "createFigureFromFile(titles="+str(titles)+")"
    f = open(inputFile,"r")
    lineNumber = []
    cov = []
    kmersInGraph = []
    kmersInBF = []
    G_ratio = []
    BF_ratio = []
    kmersInGraph_s = []     #k-mers in the Graph using a BF-stopper
    kmersInBF_s = []
    G_ratio_s = []          #ratio of k-mers in the Graph using a BF-stopper
    BF_ratio_s = []
    for line in f:
        temp = line.strip().split(",")
        L = len(temp)
        assert (L==6) or (L==10)
        lineNumber.append(int(temp[0]))
        cov.append(int(temp[1]))
        kmersInGraph.append(int(temp[2]))
        kmersInBF.append(int(temp[3]))
        G_ratio.append(float(temp[4]))
        BF_ratio.append(float(temp[5]))
        if L==10:
            kmersInGraph_s.append(int(temp[6]))
            kmersInBF_s.append(int(temp[7]))
            G_ratio_s.append(float(temp[8]))
            BF_ratio_s.append((temp[9]))
    
    #for i, c in enumerate(cov):
    #    print i,c,kmersInGraph[i]

    fig1, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2)#,sharex='col', sharey='row')
    fig1.subplots_adjust(hspace=.5)
    fig1.subplots_adjust(wspace=.5)
    ax1.axhline(y=sizeOfGenome,label="#k-mers in genome")
    ax1.plot(cov,kmersInBF,"k-",label="#k-mers in BF")
    ax1.plot(cov,kmersInGraph,"k--",label='#k-mers in G')
    if L==10:
        ax1.plot(cov,kmersInGraph_s,"r--")
    x1,x2,y1,y2 = ax1.axis()
    #ax1.axis((x1,max(cov),0,700000))
    ax1.axis((x1,max(cov),0,max(kmersInGraph)*2))
    assert len(titles)==5
    if titles==["","","","",""]:
        fig1.suptitle('No maximum coverage')
        ax1.set_title("#k-mers / cov")
        ax2.set_title("Ratio / cov")
        ax3.set_title("ln(#k-mers) / cov")
        ax4.set_title("Ratio_log / cov")
    else:
        if titles[0]=="":
            fig1.suptitle('No maximum coverage')
        else:
            fig1.suptitle(titles[0])
        if titles[1]=="":
            ax1.set_title("#k-mers / cov")
        else:
            ax1.set_title(titles[1])
        if titles[2]=="":
            ax2.set_title("Ratio / cov")
        else:
            ax2.set_title(titles[2])
        if titles[3]=="":
            ax3.set_title("ln(#k-mers) / cov")
        else:
            ax3.set_title(titles[3])
        if titles[4]=="":
            ax4.set_title("Ratio_log / cov")
        else:
            ax4.set_title(titles[4])
    ax3.set_xlabel("Coverage")
    ax1.set_xlabel("Coverage")
    ax2.set_xlabel("Coverage")
    ax4.set_xlabel("Coverage")
    ax1.set_ylabel("#k-mers")
    ax3.set_ylabel("ln(#k-mers)")
    ax2.set_ylabel("Ratio")
    ax4.set_ylabel("Ratio on log scale")
    ax1.grid()

    kmersInGraph = [0 if x==0 else log(x) for x in kmersInGraph]
    kmersInBF = [0 if x==0 else log(x) for x in kmersInBF]
    ax3.axhline(y=log(sizeOfGenome),label="#k-mers in genome")
    ax3.plot(cov,kmersInBF,"k-",label="#k-mers in BF")
    ax3.plot(cov,kmersInGraph,"k--",label='#k-mers in G')
    if L==10:
        kmersInGraph_s = [0 if x==0 else log(x) for x in kmersInGraph_s]
        ax3.plot(cov,kmersInGraph_s,"r--",label='#With BF_stopper')
    x1,x2,y1,y2 = ax3.axis()
    ax3.axis((x1,max(cov),y1,y2))
    ax3.legend(loc='lower right', bbox_to_anchor=(1.2, -0.05))
    ax3.grid()

    ax2.axhline(y=1,label="99.9% ratio")
    ax2.plot(cov,BF_ratio,"k-",label="ratio in BF")
    ax2.plot(cov,G_ratio,"k--",label='ratio in G')
    if L==10:
        ax2.plot(cov,G_ratio_s,"r--",label='With BF-stopper')
    x1,x2,y1,y2 = ax2.axis()
    ax2.axis((x1,max(cov),y1,1.1))

    ratio_1_log_scale = -log(1-0.999)
    G_ratio = [ratio_1_log_scale if x==1 else -log(1-x) for x in G_ratio]
    BF_ratio = [ratio_1_log_scale if x==1 else -log(1-x) for x in BF_ratio]
    ax4.axhline(y=ratio_1_log_scale,label="99.9% ratio")
    ax4.plot(cov,BF_ratio,"k-",label="ratio in BF")
    ax4.plot(cov,G_ratio,"k--",label='ratio in G')
    if L==10:
        G_ratio_s = [ratio_1_log_scale if x==1 else -log(1-x) for x in G_ratio_s]
        ax4.plot(cov,G_ratio_s,"r--",label='With BF-stopper')
    x1,x2,y1,y2 = ax4.axis()
    ax4.axis((x1,max(cov),y1,y2))
    ax2.legend(loc='lower right', bbox_to_anchor=(1.1, 0))

    #if outFolder=="":
    #    print "Using default outFolder inside helpers.createFigureFromFile"
    #    outFolder = 'Output/defaultOutFolder'
    fig1.savefig(outFolder+"/"+genomeName+"_maxCov_"+str(maxCov)+".png", bbox_inches='tight')
    f.close()

"""
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
"""


def printAllInfoFromFiles(fn,k):
    kd = collections.defaultdict(int)       #Stores all kmers in the files
    kd_BF = collections.defaultdict(int)    #Stores all kmers in the files
                                            #seen for the second time according to BF
    numberOfKmers = 0
    BF = Bloom.Bloom(0.01,6000000,pfn=True)
    print "Running through the files:"
    for f in fn:
        h = open(f, "rU")
        for lineNr,line in enumerate(h):
            if (lineNr%4 == 1):
                segments = filter(lambda x: len(x)>=k,line.strip().split("N"))
                if len(segments)>1:
                    print "lineNr="+str(lineNr)+", segments:",segments, "len(segments)="+str(len(segments))
                for s in segments:
                    for km in dbg.kmers(s,k):
                        kd[km] += 1
                        numberOfKmers += 1
                        if not km in BF:
                            BF.add(km)
                        else:
                            kd_BF[km] += 1
                    for km in dbg.kmers(dbg.twin(s),k):
                        kd[km] += 1
                        numberOfKmers +=1
                        if not km in BF:
                            BF.add(km)
                        else:
                            kd_BF[km] += 1

    G_twice = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
    numberOfAtLeastTwice = 0
    for km,num in kd.iteritems():
        if num > 1:
            G_twice.addSegmentToGraph(km)
            numberOfAtLeastTwice += 1
    G_BF = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False,printInit=True)
    for km,num in kd_BF.iteritems():
        G_BF.addSegmentToGraph(km)
    
    print "Total number of kmers in the files:                  ", numberOfKmers
    print "Number of unique kmers in the files:                 ", len(kd)
    print "Number of kmers in the BF:                           ", len(BF)
    print "Number of kmers occuring at least twice in the files: ", numberOfAtLeastTwice
    print "Same number according to BF:                          ", len(kd_BF)
    print "Number of kmers in G_twice:                           ", len(G_twice.kmers)
    print "Number of contigs in G_twice:                          ", len(G_twice)
    print "Number of kmers in G_BF:                              ", len(G_BF.kmers)
    print "Number of contigs in G_BF:                            ", len(G_BF)
    print BF
    BF.print_bitarray()

def createKmerDictFromSegmentList(SL,k):
    kd = collections.defaultdict(int)
    for segment in SL:
        for km in dbg.kmers(segment,k):
            kd[km] += 1
        for km in dbg.kmers(dbg.twin(segment),k):
            kd[km] += 1
    return kd

def createNaiveFromKmerDict(kd,k):
    G,cs = dbg.all_contigs(kd,k)
    G_naive = Graph.Graph(k,pfn=False,ps=False,al=False,pil=False)
    dbg.createGraphObject(G,cs,k,G_naive,pfn=False,ps=False)
    return G_naive

if __name__ == "__main__":
    pass
    #changeFile("Output/testfile.csv",[0]*21,[0]*21)
    #createFigureFromFile("Output/defaultOutFolder/Without_BF_stopper/figureData.csv",sizeOfGenome=70000,titles=["","","","",""],outFolder="")
    #createFigureFromFile("Output/defaultOutFolder/BF_stopper/maxCov_5/figureData.csv",sizeOfGenome=70000,titles=["","","","",""],outFolder="Output/defaultOutFolder/BF_stopper/maxCov_5")
