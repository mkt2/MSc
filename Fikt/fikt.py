#coding:utf8
import Graph
from compareGraphs import *
from BF_counter import *

G = Graph.Graph(k=5,pfn=False,ps=False,al=True)

#G.addSegmentToGraph("GCCGA")	#x->twin(c0)
#			   twin: TCGGC 
#G.addSegmentToGraph("ACTCG")	#c1->x
#assert G.isLegalDBG()
#G.printContigs("Before adding the problem contig")
#G.addSegmentToGraph("CTCGG")	#x

#Þetta er testið þar sem ég fann að ég get þurft að merge-a OUT tvisvar
G.contigs[0] = ["GCCGA",[],[(2,False)],0]
G.contigs[1] = ["ACTCG",[],[(2,True)],0]
G.contigs[2] = ["CTCGG",[(1,True)],[(0,False)],0]
G.setID(3)
G.addKmersFromAllContigs()
G.printContigs("After adding the problem contig")
G.printFunctionNames=True
G.printStatus=True
G.mergeSegment(2)
