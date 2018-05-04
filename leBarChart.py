#coding:utf8
import helpers
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def createBarChart(outDir,genomeName,MAC,MSC):
    gn = helpers.createGraphName(MAC,MSC)
    gn_tex = helpers.createGraphName_tex(MAC,MSC,False)
    #print gn
    #print gn_tex
    #Read the data from the files:
    G_tot, G_per = helpers.readTotalsAndPercFromFile(outDir+"/G_inf_inf_tot_perc.csv")
    Gx_tot, Gx_per = helpers.readTotalsAndPercFromFile(outDir+"/"+gn+"_tot_perc.csv")
    Gd_tot, Gd_per = helpers.readTotalsAndPercFromFile(outDir+"/"+gn+"d_tot_perc.csv")

    #Define values for the rows in the bar chart where row1 is the bottom row, row2 is the one one top, etc
    #[tot_iso,tot_tip,tot_bub4,tot_genomic,tot_complex]
    row1 = (G_tot[0], Gx_tot[0])    #iso
    row2 = (G_tot[1], Gx_tot[1])    #tip
    row3 = (G_tot[2], Gx_tot[2])    #secondary
    row4 = (G_tot[3], Gx_tot[3])    #genomic
    row5 = (G_tot[4], Gx_tot[4])    #complex

    #Initialize the plot. We start with the left subplot:
    fig, (ax1,ax2) = plt.subplots(1,2)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    #plt.tight_layout(h_pad=5.0)
    #plt.subplots_adjust(bottom=0.1, right=0.8, top=0.9)
    if genomeName=="t":
        spaceBetween = 0.3
    elif genomeName=="sa":
        spaceBetween = 0.4  #0.5 var fínt en aðeins óþarflega mikið
    else:
        raise Exception('illegal genome name')
    plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=spaceBetween, hspace=None)

    #Set commas on the y-axis in the left plot:
    ax1.get_yaxis().set_major_formatter(
    mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
    #Set commas on the y-axis in the right plot:
    ax2.get_yaxis().set_major_formatter(
    mpl.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

    #Define colors:
    color_iso =         '#5A9BD3'
    color_tip =         '#F15A61'
    color_secondary =   '#7AC36A'
    color_genomic =     '#9F66AB'
    color_complex =     '#CF6F57'

    #Create the bars:
    ind = [0.1, 0.6]    #the x locations for each bar #[0 1]
    width = 0.3         #the width of the bars
    ax1.bar(ind, row1, width, color=color_iso,   bottom=0,label="Isolated")
    ax1.bar(ind, row2, width, color=color_tip,  bottom=row1,label="Tips")
    ax1.bar(ind, row3, width, color=color_secondary, bottom=[sum(x) for x in zip(row1, row2)],label="Secondary")
    ax1.bar(ind, row4, width, color=color_genomic,    bottom=[sum(x) for x in zip(row1, row2, row3)],label="Genomic")
    ax1.bar(ind, row5, width, color=color_complex, bottom=[sum(x) for x in zip(row1, row2, row3, row4)],label="Complex")
    
    #Now the right subplot:

    #Create the bars:
    ind = [0.1]         #the x locations for each bar #[0 1]
    width = 0.3         #the width of the bars
    ax2.bar(ind, Gd_tot[0], width, color=color_iso,  bottom=0,label="Isolated")
    ax2.bar(ind, Gd_tot[1], width, color=color_tip,  bottom=Gd_tot[0],label="Tips")
    ax2.bar(ind, Gd_tot[2], width, color=color_secondary, bottom=sum(Gd_tot[0:2]),label="Secondary")
    ax2.bar(ind, Gd_tot[3], width, color=color_genomic,    bottom=sum(Gd_tot[0:3]),label="Genomic")
    ax2.bar(ind, Gd_tot[4], width, color=color_complex, bottom=sum(Gd_tot[0:4]),label="Complex")
    
    #Labels and legend
    ax1.set_xlim(0,1,1)
    ax1.xaxis.set_ticks([0.25,0.75])
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels[0] = '$G$'
    labels[1] = gn_tex
    ax1.set_xticklabels(labels)
    
    #ax1.set_xlim(0,0.8,0.8)
    ax1.set_ylim(0,sum(G_tot)*1.1)
    ax1.set_ylabel('Number of k-mers')
    ax1.set_title('$G$ and '+gn_tex)

    ax2.set_xlim(0,1,1)
    ax2.xaxis.set_ticks([0.25])
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    labels[0] = '$G-'+gn_tex[1:]
    ax2.set_xticklabels(labels)
    
    ax2.set_ylim(0,sum(Gd_tot)*1.1)
    ax2.set_title('$G-'+gn_tex[1:])
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[::-1], labels[::-1], loc='lower right', prop={'size': 13})
  
    #Save and show the figure
    fig.suptitle("A visual representation of $G$, "+gn_tex+", and "+"$G-"+gn_tex[1:]+" for genome="+genomeName)
    fig.savefig(outDir+"/"+genomeName+"_wim_"+gn+".png", bbox_inches='tight')
    #plt.show()

if __name__ == "__main__":
    outDir = "Output/t"
    genomeName="t"
    for MAC in [5, 10, 15, 20, 30]:
        createBarChart(outDir,genomeName,MAC,MSC=float('inf'))
    