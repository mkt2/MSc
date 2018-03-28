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
		maxAddCovs = [5, 10, 15, 20, 30,float('inf')]
	elif genomeName=="sa":
		maxAddCovs = [5, 10, 15, 20, 30,float('inf')]
	else:
		raise Exception("The genomeName must be either 't' or 'sa'!")
	return maxAddCovs

if __name__ == "__main__":
	genomeName = sys.argv[1]
	k = int(sys.argv[2])
	skipPrint = helpers.readSkipPrint(4)
	maxAddCovs = selectGenome(genomeName)
	outDir = "Output/"+str(genomeName)

	if not skipPrint:
		print "Reading G from Ginf.txt:"
	G = Graph.Graph(k,al=False)
	G.createGraphFromFile(outDir+"/G"+str(maxAddCovs[-1])+".txt")

	if not skipPrint:
		print "Analyzing the contigs in G:"
	G.analyzeAllContigsInCollection()
	totals_G = getTotals(G,-1)
	perc_G = [x/float(len(G)) for x in totals_G]
	helpers.writeTotalsAndPercToFile(outDir+"/G"+str(maxAddCovs[-1])+"_tot_perc.csv",totals_G,perc_G)

	#For every maxAddCov x
	#	Read Gx from Gx.txt
	#	Create I_kmers as intersection of k-mers in G and Gx
	#	Create I_ids as the set containing all contigs the k-mers in I_kmers occur in
	for maxAddCov in maxAddCovs[0:-1]:
		if not skipPrint:
			print "\nA) Starting on maxAddCov="+str(maxAddCov)+":"
		#Create Gx, I_kmers and I_ids
		Gx = Graph.Graph(k,al=False)
		Gx.createGraphFromFile(outDir+"/G"+str(maxAddCov)+".txt")
		I_kmers = set(G.kmers.keys()) - set(Gx.kmers.keys())
		assert(len(I_kmers)==len(G.kmers)-len(Gx.kmers))
		I_ids = helpers.getIDsFromSetOfKmers(G,I_kmers)
		if not skipPrint:
			print "B) Done recreating G"+str(maxAddCov)+" I_kmers and I_ids for genome="+str(genomeName)
		
		#Check for errors:
		for km in I_kmers:
			assert(km in G.kmers), "I_kmers is supposed to be storing the k-mers which occured in G and NOT in Gx"
			assert(km not in Gx),  "I_kmers is supposed to be storing the k-mers which occured in G and NOT in Gx"
		for km in Gx.kmers:
			assert(km not in I_kmers), "I_kmers only stores k-mers which didn't occur in Gx"
		for km in G.kmers:
			assert((km in I_kmers) or (km in Gx.kmers)), "All k-mers from G should either occur in I_kmers or Gx"
		if not skipPrint:
			print "C) Done checking for errors for genome="+str(genomeName)+" and maxAddCov="+str(maxAddCov)

		#Create I_kmerPairs so we only store either km or twin(km) for each km in I_kmers
		I_kmerPairs = collections.defaultdict(int)
		for km in I_kmers:
			km = min(km , dbg.twin(km))
			I_kmerPairs[km] = 1
		assert(len(I_kmerPairs)*2==len(I_kmers)), "I_kmerPairs has exactly half as many values as I_kmers"

		#Now we'll give iso, tip, bub ratings to every contig in Gx
		Gx.analyzeAllContigsInCollection()
		if not skipPrint:
			print "D) Done assigning a rating to every contig in Gx for genome="+str(genomeName)+" and maxAddCov="+str(maxAddCov)
		
		#Compute percentages of iso, tip and bubbles for Gx and I_ids
		totals_Gx = getTotals(Gx,-1)
		totals_I = getTotals(G,I_ids)
		perc_Gx = [x/float(len(Gx)) for x in totals_Gx]
		perc_I = [x/float(len(I_ids)) for x in totals_I]

		helpers.writeTotalsAndPercToFile(outDir+"/G"+str(maxAddCov)+"_tot_perc.csv",totals_Gx,perc_Gx)
		helpers.writeTotalsAndPercToFile(outDir+"/G"+str(maxAddCov)+"d_tot_perc.csv",totals_I,perc_I)
		if not skipPrint:
			print "E) Done printing the results from analyzeAllContigsInCollection to a file for genome="+str(genomeName)+" and maxAddCov="+str(maxAddCov)
		leBarChart.createBarChart(maxAddCov,outDir,genomeName)
		if not skipPrint:
			print "F Done creating the bar Chart for maxAddCov="+str(maxAddCov)
	if not skipPrint:	
		print "Finished running whatIsMissing.py for genome="+genomeName
