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
            G_ratio_s.append(temp[8])
            BF_ratio_s.append(temp[9])
    
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

    G_ratio = [-log(1-x) for x in G_ratio]
    BF_ratio = [-log(1-x) for x in BF_ratio]
    ax4.axhline(y=-log(1-0.999),label="99.9% ratio")
    ax4.plot(cov,BF_ratio,"k-",label="ratio in BF")
    ax4.plot(cov,G_ratio,"k--",label='ratio in G')
    if L==10:
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

#Examples:
#outFolder = "Output/defaultOutFolder/BF_stopper"
#newFileName = "figureData_maxCov_5.csv"
def changeFile(file,newCols,outFolder,newFileName):
    assert len(newCols)==4
    if outFolder=="":
        print "Using default outFolder inside helpers.changeFile"
        outFolder = "Output/defaultOutFolder/BF_stopper"
    #newFile = outFolder+"/"+os.path.basename(file)
    newFile = outFolder+"/"+newFileName
    copyfile(file, newFile)
    for i,line in enumerate(fileinput.input(newFile, inplace=True)):
        assert line.count(",")==5
        print "%s,%s,%s,%s,%s" % (line.strip(),str(newCols[0][i]),str(newCols[1][i]),str(newCols[2][i]),str(newCols[3][i])+"\n"),

if __name__ == "__main__":
    pass
    #changeFile("Output/testfile.csv",[0]*21,[0]*21)
    #createFigureFromFile("Output/defaultOutFolder/Without_BF_stopper/figureData.csv",sizeOfGenome=70000,titles=["","","","",""],outFolder="")
    #createFigureFromFile("Output/defaultOutFolder/BF_stopper/maxCov_5/figureData.csv",sizeOfGenome=70000,titles=["","","","",""],outFolder="Output/defaultOutFolder/BF_stopper/maxCov_5")
