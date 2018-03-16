import helpers
import numpy as np
import matplotlib.pyplot as plt

def createBarChart(maxCov,outDir,genomeName):
    #Read the data from the files:
    G_tot, G_per = helpers.readTotalsAndPercFromFile(outDir+"/G_tot_perc.csv")
    Gx_tot, Gx_per = helpers.readTotalsAndPercFromFile(outDir+"/G"+str(maxCov)+"_tot_perc.csv")
    Gd_tot, Gd_per = helpers.readTotalsAndPercFromFile(outDir+"/G"+str(maxCov)+"d_tot_perc.csv")

    #Define values for the rows in the bar chart where row1 is the bottom row, row2 is the one one top, etc
    row1 = (G_tot[0], Gx_tot[0])    #iso
    row2 = (G_tot[1], Gx_tot[1])    #tip
    row3 = (G_tot[2], Gx_tot[2])    #primary
    row4 = (G_tot[3], Gx_tot[3])    #secondary
    row5 = (G_tot[4], Gx_tot[4])    #none

    #Initialize the plot. We start with the left subplot:
    fig, (ax1,ax2) = plt.subplots(1,2)

    #Create the bars:
    ind = [0.1, 0.5]    #the x locations for each bar #[0 1]
    width = 0.2         #the width of the bars
    ax1.bar(ind, row1, width, color='green',   bottom=0,label="Isolated")
    ax1.bar(ind, row2, width, color='black',  bottom=row1,label="Tips")
    ax1.bar(ind, row3, width, color='blue', bottom=[sum(x) for x in zip(row1, row2)],label="Primary")
    ax1.bar(ind, row4, width, color='red',    bottom=[sum(x) for x in zip(row1, row2, row3)],label="Secondary")
    ax1.bar(ind, row5, width, color='orange', bottom=[sum(x) for x in zip(row1, row2, row3, row4)],label="Unidentified")
    
    #Now the right subplot:

    #Create the bars:
    ind = [0.1]    #the x locations for each bar #[0 1]
    width = 0.2         #the width of the bars
    ax2.bar(ind, Gd_tot[0], width, color='green',  bottom=0,label="Isolated")
    ax2.bar(ind, Gd_tot[1], width, color='black',  bottom=Gd_tot[0],label="Tips")
    ax2.bar(ind, Gd_tot[2], width, color='blue', bottom=sum(Gd_tot[0:2]),label="Primary")
    ax2.bar(ind, Gd_tot[3], width, color='red',    bottom=sum(Gd_tot[0:3]),label="Secondary")
    ax2.bar(ind, Gd_tot[4], width, color='orange', bottom=sum(Gd_tot[0:4]),label="Unidentified")
    
    #Labels and legend
    ax1.set_xlim(0,1,1)
    ax1.xaxis.set_ticks([0.2,0.6])
    labels = [item.get_text() for item in ax1.get_xticklabels()]
    labels[0] = 'G'
    labels[1] = 'G'+str(maxCov)
    ax1.set_xticklabels(labels)
    
    ax1.set_ylim(0,sum(G_tot)*1.01)
    ax1.set_ylabel('Number of contigs')
    ax1.set_title('The contigs in G and G'+str(maxCov))
    handles, labels = ax2.get_legend_handles_labels()
    ax1.legend(handles[::-1], labels[::-1], loc='upper right', prop={'size': 13})

    ax2.set_xlim(0,1,1)
    ax2.xaxis.set_ticks([0.2])
    labels = [item.get_text() for item in ax2.get_xticklabels()]
    labels[0] = '(G-G'+str(maxCov)+")*"
    ax2.set_xticklabels(labels)
    
    ax2.set_ylim(0,sum(Gd_tot)*1.01)
    ax2.set_title('The contigs in (G-G'+str(maxCov)+")*")
    handles, labels = ax2.get_legend_handles_labels()
    ax2.legend(handles[::-1], labels[::-1], loc='upper right', prop={'size': 13})
  
    #Save and show the figure
    fig.suptitle(genomeName+". maxCov="+str(maxCov))
    fig.savefig(outDir+"/wim"+str(maxCov)+".png", bbox_inches='tight')
    #plt.show()

if __name__ == "__main__":
    pass
    #outDir = "Output/t"
    #genomeName="t"
    #for maxCov in [5, 10, 15, 20, 30]:
    #    createBarChart(maxCov,outDir,genomeName)
    