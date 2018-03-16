#coding:utf8
import dbg
import collections
import helpers
import Graph, Bloom
import time
import sys
import leBarChart

def getTotals(G,theContigs=-1):
	if theContigs==-1:
		theContigs = G.contigs.keys()
	tot_iso = 0
	tot_tip = 0
	tot_bub3 = 0
	tot_bub4 = 0
	tot_non = 0
	for c_ID in theContigs:
		R = G.degrees[c_ID][2]
		if R==0:
			tot_non += 1
		elif R==1:
			tot_iso += 1
		elif R==2:
			tot_tip += 1
		elif R==3:
			tot_bub3 += 1
		elif R==4:
			tot_bub4 += 1
		else:
			raise Exception("Illegal value for R")
	return [tot_iso,tot_tip,tot_bub3,tot_bub4,tot_non]

def selectGenome(genomeName):
	if genomeName=="t":
		maxCovs = [-1, 5, 10, 15, 20, 30]
		outDir = "Output/t"
	elif genomeName=="sa":
		maxCovs = [-1, 5, 10, 15, 20, 30, 35, 40, 45]
		#maxCovs = [-1,20,30]	#Get bara notað þessi gildi tímabundið af því að mig vantar
								#kmersFromG_maxCov_.txt skrár fyrir hin gildin á maxCov
		outDir = "Output/sa"
	else:
		raise Exception("The genomeName must be either 't' or 'sa'!")
	return maxCovs, outDir

"""það er óþarfi að búa til grafið sjálft!
def create_Gdx(G,Gx):
	print "create_Gdx(G,Gx)"
	Gd = Graph.Graph(k)

	Gd_kd = set(G.kmers.keys()) - set(Gx.kmers.keys())
	assert(len(Gd_kd)==len(G.kmers)-len(Gx.kmers))
	print "Done creating Gd_kd=G.kmers-Gx.kmers"

	#Add all the k-mers from Gd_kd to Gd (they'll all get COV=2)
	count = 0
	L_Gd_kd = len(Gd_kd)
	for km in Gd_kd:
		count+=1
		if count%100==0:
			print "We have added "+str(count)+" out of "+str(L_Gd_kd)+" k-mers to Gd"
		#print km
		assert(isinstance(km,str)), "km has to be a string"
		assert(len(km)==k)
		Gd.addSegmentToGraph(km)
	print "Done adding all k-mers from G-Gx to Gdx"
	
	#Set the COV of each contig in Gd to what it is in G
	for c_ID, [c,c_IN,c_OUT,c_COV] in Gd.contigs.iteritems():		
		#update c_COV so that it will be the sum of the approximate COV of all k-mers in c in G
		#def getAverageKmerCoverageOfContig(self,cID)
		new_c_COV = 0
		for km in dbg.kmers(c,k):
			assert(km in G.kmers), "All k-mers in Gd should also be in G"
			[x_ID,x_i,x_B] = G.kmers[km]
			new_c_COV += G.getAverageKmerCoverageOfContig(x_ID)
		Gd.contigs[c_ID][3] = new_c_COV
	
	return Gd
"""

if __name__ == "__main__":
	genomeName = sys.argv[1]
	k = int(sys.argv[2])
	maxCovs, outDir = selectGenome(genomeName)

	#Read G from G.txt:
	G = Graph.Graph(k)
	G.createGraphFromFile(outDir+"/G.txt")

	#Analyze the contigs in G:
	G.analyzeAllContigsInCollection()
	totals_G = getTotals(G,-1)
	perc_G = [x/float(len(G)) for x in totals_G]
	helpers.writeTotalsAndPercToFile(outDir+"/G_tot_perc.csv",totals_G,perc_G)

	#For every maxCov x
	#	Read Gx from Gx.txt
	#	Create I_kmers as intersection of k-mers in G and Gx
	#	Create I_ids as the set containing all contigs the k-mers in I_kmers occur in
	for maxCov in maxCovs[1:]:
		print "\nA) Starting on maxCov="+str(maxCov)+":"
		#Create Gx, I_kmers and I_ids
		Gx = Graph.Graph(k)
		Gx.createGraphFromFile(outDir+"/G"+str(maxCov)+".txt")
		I_kmers = set(G.kmers.keys()) - set(Gx.kmers.keys())
		assert(len(I_kmers)==len(G.kmers)-len(Gx.kmers))
		I_ids = helpers.getIDsFromSetOfKmers(G,I_kmers)
		assert(len(G) > len(Gx)+len(I_ids))
		print "B) Done recreating G"+str(maxCov)+" I_kmers and I_ids for genome="+str(genomeName)
		
		#Check for errors:
		for km in I_kmers:
			assert(km in G.kmers), "I_kmers is supposed to be storing the k-mers which occured in G and NOT in Gx"
			assert(km not in Gx),  "I_kmers is supposed to be storing the k-mers which occured in G and NOT in Gx"
		for km in Gx.kmers:
			assert(km not in I_kmers), "I_kmers only stores k-mers which didn't occur in Gx"
		for km in G.kmers:
			assert((km in I_kmers) or (km in Gx.kmers)), "All k-mers from G should either occur in I_kmers or Gx"
		print "C) Done checking for errors for genome="+str(genomeName)+" and maxCov="+str(maxCov)

		#Create I_kmerPairs so we only store either km or twin(km) for each km in I_kmers
		I_kmerPairs = collections.defaultdict(int)
		for km in I_kmers:
			km = min(km , dbg.twin(km))
			I_kmerPairs[km] = 1
		assert(len(I_kmerPairs)*2==len(I_kmers)), "I_kmerPairs has exactly half as many values as I_kmers"

		#Now we'll give iso, tip, bub ratings to every contig in Gx
		Gx.analyzeAllContigsInCollection()
		print "D) Done assigning a rating to every contig in Gx for genome="+str(genomeName)+" and maxCov="+str(maxCov)
		
		#Compute percentages of iso, tip and bubbles for Gx and I_ids
		totals_Gx = getTotals(Gx,-1)
		totals_I = getTotals(G,I_ids)
		perc_Gx = [x/float(len(Gx)) for x in totals_Gx]
		perc_I = [x/float(len(I_ids)) for x in totals_I]

		helpers.writeTotalsAndPercToFile(outDir+"/G"+str(maxCov)+"_tot_perc.csv",totals_Gx,perc_Gx)
		helpers.writeTotalsAndPercToFile(outDir+"/G"+str(maxCov)+"d_tot_perc.csv",totals_I,perc_I)
		print "E) Done printing the results from analyzeAllContigsInCollection to a file for genome="+str(genomeName)+" and maxCov="+str(maxCov)

        leBarChart.createBarChart(maxCov,outDir,genomeName)
		
	print "Finished running whatIsMissing.py for genome="+genomeName
