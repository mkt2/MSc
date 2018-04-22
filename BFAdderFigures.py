#coding:utf8
import helpers
import collections
import sys
import dbg, Bloom, Graph
import os.path
from shutil import copyfile
import fileinput
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from math import log, exp, pow
import time
import pathlib2


#for j in range(0,len(line)-1):
#   covDict[maxAddCov][j] = line[j+1]
#        0     1     2      3
#Geymir COV, lenG, fracG, lenBF og fracBF fyrir G
def createFigure(covDict,genomeLen,numReadsPerFile,numKmersPerRead,genome_dict,maxAddCov,maxSplitCov,figTitle,figName):
    def toLogScale(myList,ratio_1_log_scale):
        return [ratio_1_log_scale if x==1 else -log(1-x,10) for x in myList]

    #Væri fínt að nefna þetta fall eitthvað meira lýsandi. 
    #Núna er fallið að skila líkunum á að stakur k-mer X úr genome-inu (en ekki grafinu)
    #sé sequence-aður amk tvisvar
    def lw(C,N):
        #C:   Coverage
        #N:   Number of k-mers in the genome
        #P:   The probability of sequencing a given k-mer from the genome at least twice
        #     This should also be the fraction of k-mers sequenced at least once
        #P*N: The expected number of sequenced k-mers for the given coverage
        P = 1.0 - exp(-C) - C*exp(-C)    #P = 1-e^{-C}-Ce^{-C} skv. ritgerð
        
        if P==1:
            P = 0.99999    #Hér ætti P að vera hlutfall k-mera úr G sem við búumst við að sjá í öllum read-unum í heild. Þetta ætti þó að duga fyrir myndirnar
        assert(P>=0), "C="+str(C)+", P="+str(P)
        assert(P<1), "C="+str(C)+", P="+str(P)
        return P

    #Teiknar inn línur fyrir Gx og Gxs
    def addCurrent(key):
        assert(key!=(float('inf'),float('inf')))
        maxAddCov = key[0]
        maxSplitCov = key[1]

        #Plot Gx on top plot:
        ax1.plot(COV,covDict[(maxAddCov,float('inf'))][1],"r-",label='G'+str(maxAddCov))

        #Plot Gx on bottom plot:
        Gx_frac_log = [ratio_1_log_scale if x==1 else -log(1-x,10) for x in covDict[(maxAddCov,float('inf'))][2]]
        ax2.plot(COV,Gx_frac_log,"r-",label='G'+str(maxAddCov))

        if maxSplitCov!=float('inf'):
            #Plot Gxs on top plot:
            ax1.plot(COV,covDict[key][1],"g-",label='G'+str(maxAddCov)+"s"+str(maxSplitCov))

            #Plot Gxs on bottom plot:
            Gxs_frac_log = [ratio_1_log_scale if x==1 else -log(1-x,10) for x in covDict[key][2]]
            ax2.plot(COV,Gxs_frac_log,"g--",label='G'+str(maxAddCov)+"s"+str(maxSplitCov))

        #Setjum tölur inn á gröfin við hliðina á línunum fyrir Gx og Gxs:
        vars_x1 = [covDict[(float('inf'),float('inf'))][1],covDict[key][1]]
        if maxSplitCov!=float('inf'):
            vars_x1.append(covDict[(maxAddCov,float('inf'))][1])
        vars_x2 = [G_BF_frac_log,G_frac_log,Gx_frac_log]
        for var in vars_x1:
            ax1.annotate(' %s' % format(max(var),','), xy=(0,0), xytext=(max(COV),max(var)))
        for var in vars_x2:
            ax2.annotate(' %s' % format(float('{:.2f}'.format(max(var),','))), xy=(0,0), xytext=(max(COV),max(var)))

    #Initialize the plot:
    fig, (ax1,ax2) = plt.subplots(2,1)
    COV = covDict[(float('inf'),float('inf'))][0]

    #Plot Ginf on top plot
    ax1.axhline(genomeLen,label="Genome")
    ax1.plot(COV,covDict[(float('inf'),float('inf'))][3],"m-",label='BFinf')
    ax1.plot(COV,covDict[(float('inf'),float('inf'))][1],"k-",label='Ginf')
    x1,x2,y1,y2 = ax1.axis()
    ax1.axis((x1,max(COV),0,max(covDict[(float('inf'),float('inf'))][1])*1.3))
    #ax1.plot(COV,covDict[maxAddCov][1],"r--",label='G'+str(maxAddCov))
    
    #Plot Ginf on bottom plot:
    ratio_1_log_scale = -log(1-0.9999,10)
    ratio_1_log_scale_2 = -log(1-0.99,10)
    ax2.axhline(ratio_1_log_scale,label="Fraction of 99.99%")
    ax2.axhline(ratio_1_log_scale_2,label="Fraction of 99%")
    G_BF_frac_log = toLogScale(covDict[(float('inf'),float('inf'))][4],ratio_1_log_scale)
    G_frac_log = toLogScale(covDict[(float('inf'),float('inf'))][2],ratio_1_log_scale)
    ax2.plot(COV,G_BF_frac_log,"m-",label="BF frac")
    ax2.plot(COV,G_frac_log,"k-",label='G frac')
    x1,x2,y1,y2 = ax2.axis()
    ax2.axis((x1,max(COV),y1,y2))
    #Bætum lw við myndirnar:
    G_lw = [lw(C,genomeLen) for C in COV]
    G_lw_log = toLogScale(G_lw,ratio_1_log_scale)
    ax2.plot(COV,G_lw_log,"y-",label='G.lw')

    #Plot Gx and Gxs:
    key = (maxAddCov,maxSplitCov)
    addCurrent(key)

    #Plot Gx on bottom plot:
    #Gx_frac_log = [ratio_1_log_scale if x==1 else -log(1-x,10) for x in covDict[maxAddCov][2]]
    #ax2.plot(COV,Gx_frac_log,"r--",label='G'+str(maxAddCov))

    #Setjum tölur inn á gröfin við hliðina á línunum:
    #G:
    ax1.annotate(' %s' % format(genomeLen,','), xy=(0,0), xytext=(max(COV),genomeLen))
    ax2.annotate(' %s' % format(float('{:.2f}'.format(ratio_1_log_scale,','))), xy=(0,0), xytext=(max(COV),ratio_1_log_scale))
    ax2.annotate(' %s' % format(float('{:.2f}'.format(ratio_1_log_scale_2,','))), xy=(0,0), xytext=(max(COV),ratio_1_log_scale_2))
    #vars_x1 = [covDict[float('inf')][1],covDict[maxAddCov][1]]
    #vars_x2 = [G_BF_frac_log,G_frac_log,Gx_frac_log]
    #for var in vars_x1:
    #    ax1.annotate(' %s' % format(max(var),','), xy=(0,0), xytext=(max(COV),max(var)))
    #for var in vars_x2:
    #    ax2.annotate(' %s' % format(float('{:.2f}'.format(max(var),','))), xy=(0,0), xytext=(max(COV),max(var)))

    #Set commas on the y-axis in the above plot:
    ax1.get_yaxis().set_major_formatter(
    mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    #Set axis labels and legends:
    ax1.set_title("#k-mers / cov")
    ax2.set_title("-log(1-Fraction) / cov")
    ax1.set_xlabel("Coverage")
    ax2.set_xlabel("Coverage")
    ax1.set_ylabel("#k-mers")
    ax2.set_ylabel("-log(1-Fraction)")
    ax1.grid()
    ax2.grid()
    ax1.legend(loc='upper left', prop={'size': 9})#, bbox_to_anchor=(1.2, -0.05))
    ax2.legend(loc='lower right', prop={'size': 9})#, bbox_to_anchor=(1.1, 0))
    
    plt.tight_layout()

    fig.suptitle(figTitle, y=1.01)
    fig.savefig(figName, bbox_inches='tight')
    plt.close(fig)

if __name__ == "__main__":
    genomeName = sys.argv[1]
    
    #Lesum úr preprocess skránum:
    p = 0.01
    outDir = "Output/"+genomeName
    assert(pathlib2.Path(outDir+"/genome_info.csv").is_file()), "We have to create this file using preprocessInput.py before running BF_counter.py"
    assert(pathlib2.Path(outDir+"/kmers_genome.txt").is_file()), "We have to create this file using preprocessInput.py before running BF_counter.py"
    assert(pathlib2.Path(outDir+"/kmers_reads.txt").is_file()), "We have to create this file using preprocessInput.py before running BF_counter.py"
    preprocFiles = [outDir+"/genome_info.csv",outDir+"/kmers_genome.txt",outDir+"/kmers_reads.txt"]
    genomeInfo = helpers.readGenomeInfoFromFile(preprocFiles[0])
    k = int(genomeInfo[1])
    genomeLen = int(genomeInfo[2])   #Number of k-mers in the genome
    numReadsPerFile = int(genomeInfo[3])
    numKmersPerRead = int(genomeInfo[4])
    genome_dict = helpers.readKmersFromFileToDict(preprocFiles[1])  #A dictionary storing all k-mers from the genome

    
    #Búum til myndirnar:
    covDict = helpers.readCovDictFromFile(fileName=outDir+"/covDict.txt")
    keys = covDict.keys()
    print keys
    maxAddCovs = [x[0] for x in keys]
    maxSplitCovs = [x[1] for x in keys]
    #print keys
    #print maxAddCovs
    #print maxSplitCovs
    maxAddCovs.remove(float('inf'))
    print "Creating figures for every Gx"
    for maxAddCov in maxAddCovs:
        maxSplitCov = float('inf')
        figName = outDir+"/"+genomeName+"_"+str(maxAddCov)+"_inf"+".png"
        figTitle = genomeName+". MAC="+str(maxAddCov)
        createFigure(covDict,genomeLen,numReadsPerFile,numKmersPerRead,genome_dict,maxAddCov,maxSplitCov,figTitle,figName)

    print "Creating figures for every Gxs"
    for (MAC,MSC) in covDict.keys():
        if MAC==float('inf'):
            continue
        if MSC==float('inf'):
            continue
        figName = outDir+"/"+genomeName+"_"+str(MAC)+"_"+str(MSC)+".png"
        #figName = outDir+"/"+genomeName+"_maxAddCov"+str(maxAddCov)+"_maxSplitCov"+str(maxSplitCov)+".png"
        figTitle = genomeName+". MAC="+str(MAC)+", MSC="+str(MSC)
        createFigure(covDict,genomeLen,numReadsPerFile,numKmersPerRead,genome_dict,MAC,MSC,figTitle,figName)