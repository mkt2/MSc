#coding:utf8
import dbg
import collections
import helpers
import Graph, Bloom
import time
import sys
import leBarChart

#Fram að 24. apríl var þetta að skila fjölda contiga
#sem eru iso, tip, bub3, bub4 eða ekkert af þessu

#Frá og með 24. apríl skilar þetta heildarfjölda k-mera
#í umræddum contigum
def getTotals(G,k,theContigs=-1):
	if theContigs==-1:
		theContigs = G.contigs.keys()
	tot_iso = 0
	tot_tip = 0
	#tot_bub3 = 0
	tot_bub4 = 0
	#tot_non = 0
	tot_genomic = 0
	tot_complex = 0
	for c_ID in theContigs:
		[c,c_IN,c_OUT,c_COV] = G.contigs[c_ID]
		R = G.degrees[c_ID][2]
		if R==0:
			tot_complex += helpers.cLen(c,k)
		elif R==1:
			tot_iso += helpers.cLen(c,k)
		elif R==2:
			tot_tip += helpers.cLen(c,k)
		elif R==3:
			raise Exception('R cant be 3')
			#tot_bub3 += helpers.cLen(c,k)
		elif R==4:
			tot_bub4 += helpers.cLen(c,k)
		elif R==5:
			tot_genomic += helpers.cLen(c,k)
		else:
			raise Exception("Illegal value for R")
	#return [tot_iso,tot_tip,tot_bub3,tot_bub4,tot_non]
	return [tot_iso,tot_tip,tot_bub4,tot_genomic,tot_complex]

def totalToPerc(Graph,totals,numKmers=-1):
	if numKmers==-1:
		return [x/float(Graph.numKmerPairs()) for x in totals]
	else:
		return [x/float(numKmers) for x in totals]

def createNewTable(outDir,total_G,total_Gxy,MAC,MSC):
	gn_tex = helpers.createGraphName_tex(MAC,MSC,inTable=True)
	tex_newLine = " \\\\\n"
	
	[tot_iso,tot_tip,tot_bub4,tot_genomic,tot_complex] = total_G
	trash_G = tot_iso+tot_tip+tot_bub4
	genomic_G = tot_genomic
	complex_G = tot_complex

	[tot_iso,tot_tip,tot_bub4,tot_genomic,tot_complex] = total_Gxy
	trash_Gxy = tot_iso+tot_tip+tot_bub4
	genomic_Gxy = tot_genomic
	complex_Gxy = tot_complex
	
	outFile = outDir+"/newTable_"+str(MAC)+"_"+str(MSC)+".txt"
	f = open(outFile, 'w')
	f.write(" & $G$ & "+str(gn_tex)+tex_newLine)
	f.write("\\hline\n")
	f.write("Genomic & "+str(genomic_G)+" & "+str(genomic_Gxy)+tex_newLine)
	f.write("Trash & "+str(trash_G)+" & "+str(trash_Gxy)+tex_newLine)
	f.write("Complex & "+str(complex_G)+" & "+str(complex_Gxy)+tex_newLine)
	f.close()


if __name__ == "__main__":
	genomeName = sys.argv[1]
	k = int(sys.argv[2])
	skipPrint = helpers.readSkipPrint(4)
	outDir = "Output/"+str(genomeName)
	#filters = selectGenome(genomeName)
	covDict = helpers.readCovDictFromFile(fileName=outDir+"/covDict.txt")
	filters = covDict.keys()
	filters.sort(key=lambda tup: (tup[0],tup[1]), reverse=False)
	kmersInGenome = helpers.readKmersFromFileToSet(kmerFile=outDir+'/kmers_genome.txt')

	if not skipPrint:
		print "Running whatIsMissing.py for genome="+str(genomeName)+":"

	if not skipPrint:
		print "Reading G from Ginf.txt:"
	G = Graph.Graph(k,al=False)
	G.createGraphFromFile(outDir+"/G_inf_inf.txt")

	if not skipPrint:
		print "Analyzing the contigs in G:"
	G.analyzeAllContigsInCollection(kmersInGenome)
	totals_G = getTotals(G,k,-1)
	#perc_G = [x/float(G.numKmerPairs()) for x in totals_G]
	perc_G = totalToPerc(G,totals_G)
	assert(filters[-1]==(float('inf'),float('inf'))), str(filters)
	helpers.writeTotalsAndPercToFile(outDir+"/G_inf_inf_tot_perc.csv",totals_G,perc_G)

	#For every maxAddCov x
	#	Read Gx from Gx.txt
	#	Create I_kmers as intersection of k-mers in G and Gx
	#	Create I_ids as the set containing all contigs the k-mers in I_kmers occur in
	for MAC, MSC in filters[0:-1]:
		gn = helpers.createGraphName(MAC,MSC)
		#print gn

		if not skipPrint:
			print "\nA) Starting on MAC="+str(MAC)+" and MSC="+str(MSC)+":"
		#Create Gx, I_kmers and I_ids
		Gx = Graph.Graph(k,al=False)
		if MSC==float('inf'):
			Gx.createGraphFromFile(outDir+"/"+gn+".txt")
		else:
			Gx.createGraphFromFile(outDir+"/"+gn+".txt")
		I_kmers = set(G.kmers.keys()) - set(Gx.kmers.keys())
		assert(len(I_kmers)==len(G.kmers)-len(Gx.kmers))
		I_ids = helpers.getIDsFromSetOfKmers(G,I_kmers)
		if not skipPrint:
			print "B) Done recreating "+gn+", I_kmers and I_ids"
		
		#Check for errors:
		for km in I_kmers:
			assert(km in G.kmers), "I_kmers is supposed to be storing the k-mers which occured in G and NOT in Gx"
			assert(km not in Gx),  "I_kmers is supposed to be storing the k-mers which occured in G and NOT in Gx"
		for km in Gx.kmers:
			assert(km not in I_kmers), "I_kmers only stores k-mers which didn't occur in Gx"
		for km in G.kmers:
			assert((km in I_kmers) or (km in Gx.kmers)), "All k-mers from G should either occur in I_kmers or Gx"
		if not skipPrint:
			print "C) Done checking for errors"

		#Create I_kmerPairs so we only store either km or twin(km) for each km in I_kmers
		I_kmerPairs = collections.defaultdict(int)
		for km in I_kmers:
			km = min(km , dbg.twin(km))
			I_kmerPairs[km] = 1
		assert(len(I_kmerPairs)*2==len(I_kmers)), "I_kmerPairs has exactly half as many values as I_kmers"

		#Now we'll give iso, tip, bub ratings to every contig in Gx
		Gx.analyzeAllContigsInCollection(kmersInGenome)
		if not skipPrint:
			print "D) Done assigning a rating to every contig in Gx"
		
		#Compute percentages of iso, tip and bubbles for Gx and I_ids
		totals_Gx = getTotals(Gx,k,-1)
		totals_I = getTotals(G,k,I_ids)
		#perc_Gx = [x/float(len(Gx)) for x in totals_Gx]
		perc_Gx = totalToPerc(Gx,totals_Gx)
		#perc_I = [x/float(len(I_ids)) for x in totals_I]
		numKmersInI_ids = 0
		for c_ID in I_ids:
			[c,c_IN,c_OUT,c_COV] = G.contigs[c_ID]
			numKmersInI_ids += helpers.cLen(c,k)
		perc_I = totalToPerc(I_ids,totals_I,numKmersInI_ids)

		helpers.writeTotalsAndPercToFile(outDir+"/"+gn+"_tot_perc.csv",totals_Gx,perc_Gx)
		helpers.writeTotalsAndPercToFile(outDir+"/"+gn+"d_tot_perc.csv",totals_I,perc_I)
		if not skipPrint:
			print "E) Done printing the results from analyzeAllContigsInCollection to a file"
		leBarChart.createBarChart(outDir,genomeName,MAC,MSC)
		createNewTable(outDir,totals_G,totals_Gx,MAC,MSC)
		if not skipPrint:
			print "F Done creating the bar Chart"
	if not skipPrint:	
		print "Finished running whatIsMissing.py for genome="+genomeName
