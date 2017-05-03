#coding:utf8

import Graph, dbg
import sys, collections
from BF_counter import BF_counter_naive


#Before:	G1 and G2 represent two legal DBGs
#After:		returns True if G1 and G2 are the same DBG. Otherwise returns False
#Note:		G1 and G2 can be different DBGs but still represent all the same strings.
#			This function ignores that possibility. We will use another function to check whether two different graph indeed represent the same strings

#Method 1:	1. create list L1 for all possible strings from G1. Same for G2
#			2. compare L1 and L2
#			3. if the contain the same strings. Return True
#			   otherwise: Return False

#Can't see how to implement step 1 of Method 1. Therefore try to come up with Method 2
#Method 2:
#			1. 	Check whether G1 and G2 have equally many contigs
#				1a. 	Check whether they have the same contigs (or twins)
#				1b. 	if we find contigs with different IDs, try to change the IDs so that contig c has the same ID in G1 and G2
#				1c. 	Now make sure that G1.contigs[ID] == G2.contigs[ID]
#			2. 	Check whether each contig in G1 has equally many INs/OUTs as
#				it's complementary in G2
#			2b. Check wheather the INs/OUTs of the contigs in G1 and G2 are the same
def isSameGraph(G1,G2,alsoCompareWithNaive=True,pfn=False,ps=False,relaxAssertions=False):
	if pfn:
		print "isSameGraph(G1,G2,alsoCompareWithNaive="+str(alsoCompareWithNaive)+")"
	if ps:
		print "Beginning of isSameGraph:"
		G1.printContigs("G1")
		G2.printContigs("G2")

	if not relaxAssertions:
		assert G1.isLegalDBG()==True, "G1 has to be legal"
		assert G2.isLegalDBG()==True, "G2 has to be legal"
	else:
		assert G1.isLegalGraph()==True, "G1 has to be legal"
		assert G2.isLegalGraph()==True, "G2 has to be legal"

	#compare G1 and G2
	if not len(G1)==len(G2):
		if ps:
			print "The graphs don't have the same number of contigs"
		return False

	if not sameContigs(G1,G2,pfn,ps):
		if ps:
			print "The graphs don't have the same contigs"
		return False

	if ps:
		print "Before fixing IDs"
		G1.printContigs("G1")
		G2.printContigs("G2")
	
	if not sameIDs(G1,G2,pfn,ps):
		if ps:
			print "The graphs don't have the same IDs"
		fixIDs(G1,G2,pfn,ps)

	if ps:
		print "Done fixing IDs"
		G1.printContigs("G1")
		G2.printContigs("G2")

	#Now G1 and G2 have the same IDs for each contig
	#Contigs in G1 may be twins of contigs with same ID in G2
	#All we need to do now is:
	#	flip twins so that G1 and G2 will have exactly the same contigs (and should therefore have the same IN/OUT if everything is correct)
	#	make sure they have the same IN, OUT and COV for each contig

	fixTwins(G1,G2,pfn,ps)
	if ps:
		print "Done fixing twins"
		G1.printContigs("G1")
		G2.printContigs("G2")

	#print "Return True if G1 and G2 have the same INs, OUTs and COV. False otherwise"
	if not same_INs_OUTs_COV(G1,G2,pfn,ps):
		if ps:
			print "The graphs don't have the same INs OUTs and COV"
		return False

	#Compare with G_naive
	if alsoCompareWithNaive:
		if ps:
			G1.printContigs("G1")
			G2.printContigs("G2")
		#B1 = graphEqualsNaive(G1,-1,pfn,ps)
		#B2 = graphEqualsNaive(G2,-1,pfn,ps)
		B1 = G1.equalsNaive()
		B2 = G2.equalsNaive()
		if not B1:
			print "G1 does not equal G_naive"
		elif not B2:
			print "G2 does not equal G_naive"
		return B1 and B2
	else:
		return True


#Before:	G1 and G2 are DBG graphs
#After:		return True if G1 and G2 have the same contigs (or twins). return False otherwise
def sameContigs(G1,G2,pfn=True,ps=True):
	if pfn:
		print "sameContigs(G1,G2)"
	if not len(G1.contigs)==len(G2.contigs):
		return False

	#idea:
	#for each contig c_1 in G1:
	#	make sure it's also in G2
	for ID_1 in G1.contigs:
		[c_1,IN_1,OUT_1,COV_1] = G1.contigs[ID_1]
		found_c_1 = False
		for ID_2 in G2.contigs:
			[c_2,IN_2,OUT_2,COV_2] = G2.contigs[ID_2]
			if c_1==c_2 or c_1==dbg.twin(c_2):
				found_c_1 = True
				break
		if not found_c_1:
			return False

	#if we get here it means that every c_1 in G1 is also in G2. len(G1)==len(G2) => G1 and G2 have the same contigs
	return True


#Before:	G1 and G2 are DBG graphs
#			G1 and G2 have the same contigs
#After:		return True if the contigs in G1 and G2 have the same IDs (or twins). return False otherwise
def sameIDs(G1,G2,pfn=True,ps=True):
	if pfn:
		print "sameIDs(G1,G2)"
	#idea:
	#for each ID ID_1 in G1:
	#	make sure G2[ID_1] contains the same contig (I could also check for IN, OUT and COV here)
	for ID_1 in G1.contigs:
		[c_1,IN_1,OUT_1,COV_1] = G1.contigs[ID_1]
		if ID_1 in G2.contigs:
			[c_2,IN_2,OUT_2,COV_2] = G2.contigs[ID_1]
		else:
			return False
		if not (c_1==c_2 or c_1==dbg.twin(c_2)):
			return False
	return True


#Before:	G1 and G2 are DBG graphs
#			G1 and G2 have the same contigs 
#			Some contig c does not have the same ID in G1 and G2
#After:		all contigs c have the same ID in G1 and G2
def fixIDs(G1,G2,pfn=True,ps=True):
	if pfn:
		print "fixIDs(G1,G2)"
	#Step 1: Create a list of all changes we wish to make for G1 and G2 respectively
	new_G1_IDs = []		#stores 1 tuple for each contig from G1. Each tuple is of the form (oldID,newID)
						#where oldID is the ID of the contig and newID is the ID we wish to change it to
	new_G2_IDs = []
	lowestFree_ID = max(G1.ID,G2.ID)
	ID_count = 0
	for ID_1 in G1.contigs:
		[c_1,IN_1,OUT_1,COV_1] = G1.contigs[ID_1]
		found_c_1 = False
		for ID_2 in G2.contigs:
			[c_2,IN_2,OUT_2,COV_2] = G2.contigs[ID_2]
			if c_1==c_2 or c_1==dbg.twin(c_2):
				found_c_1 = True
				break
		if not found_c_1:
			raise Exception('G1 and G2 are supposed to have the same contigs before this function is called')
		elif found_c_1:
			if ID_1==ID_2:
				#c_1 has the same ID in G1 and G2
				pass
				#new_G1_IDs.append((ID_1,ID_1))
				#new_G2_IDs.append((ID_1,ID_1))
			else:
				#c_1 and c_2 need to have the same ID so we give them the lowest ID available in both G1 and G2
				#ID_new = max(G1.getID(),G2.getID())
				ID_new = lowestFree_ID + ID_count
				new_G1_IDs.append((ID_1,ID_new))
				new_G2_IDs.append((ID_2,ID_new))
				ID_count += 1
		
	#Step 2: Make all the changes we added to the lists
	if ps:
		print "changing G1"
		print new_G1_IDs
		print new_G2_IDs
	for (oldID,newID) in new_G1_IDs:
		G1.changeID_FromTo(oldID,newID)

	if ps:
		print "changing G2"
		print new_G1_IDs
		print new_G2_IDs
	for (oldID,newID) in new_G2_IDs:
		G2.changeID_FromTo(oldID,newID)

	#newID = max(G1.ID,G2.ID)#+ID_count
	#G1.setID(newID)
	#G2.setID(newID)


#Before:	G1 and G2 are DBG graphs
#			G1 and G2 have the same contigs 
#			All contigs c have the same ID in G1 and G2
#			Contigs in G1 may be twins of contigs with same ID in G2 (and vice versa)
#After:		All contigs c1 with ID ID in G1 are the EXACT same contig as the contig with that ID in G2
def fixTwins(G1,G2,pfn=True,ps=True):
	if pfn:
		print "fixTwins(G1,G2)"
	if ps:
		G1.printContigs("G1")
		G2.printContigs("G2")

	for ID_1 in G1.contigs:
		[c_1,IN_1,OUT_1,COV_1] = G1.contigs[ID_1]
		for ID_2 in G2.contigs:
			[c_2,IN_2,OUT_2,COV_2] = G2.contigs[ID_2]
			if c_1==dbg.twin(c_2):
				G2.flipContig(ID_2)
				break



#Before:	G1 and G2 are DBG graphs
#			G1 and G2 have the same contigs 
#			All contigs c have the same ID in G1 and G2
#			All contigs c1 with ID ID in G1 are the EXACT same contig as the contig with that ID in G2
#After:		True if all contigs c have the same IN, OUT and COV in G1 and G2
def same_INs_OUTs_COV(G1,G2,pfn=True,ps=True):
	if pfn:
		print "same_INs_OUTs_COV(G1,G2)"
	if G1.isEmpty() and G2.isEmpty():
		return True
	found = False
	for ID_1 in G1.contigs:
		[c_1,IN_1,OUT_1,COV_1] = G1.contigs[ID_1]
		for ID_2 in G2.contigs:
			[c_2,IN_2,OUT_2,COV_2] = G2.contigs[ID_2]
			if c_1==c_2:
				found = True
				i = sorted(IN_1)==sorted(IN_2)
				o = sorted(OUT_1)==sorted(OUT_2)

				if not (i and o):
					if ps:
						print "ID:", ID_1
						print sorted(IN_1)
						print sorted(IN_2)
						print sorted(OUT_1)
						print sorted(OUT_2)
					return False

		if not found:
			#We didn't find c_1 in G2
			raise Exception("G1 and G2 don't have the same contigs!")
		found = False
	return True

"""
def graphEqualsNaive(G,G_naive=-1,pfn=True,ps=True):
	if pfn:
		print "graphEqualsNaive(G,G_naive)"
	if G_naive==-1:
		G_naive = G.createNaive()
	if ps:
		G.printContigs("G")
		G_naive.printContigs("G_naive")
	return isSameGraph(G,G_naive,False,pfn,ps,relaxAssertions=False)
"""

if __name__ == "__main__":
    k = 3
    fileReads = []
    for f in sys.argv[2:]:
        fileReads.append(dbg.myParse(f))    #index i stores a list of all reads from file i
    d = dbg.build(fileReads,k,1)
    G,cs = dbg.all_contigs(d,k)
    #dbg.print_GFA(G,cs,k)
    G1 = dbg.createGraphObject(G,cs,k)
    G1.printContigs()
    print sameContigs(G1,G1)
    print sameIDs(G1,G1)
    print same_INs_OUTs_COV(G1,G1)
    print isSameGraph(G1,G1)