#coding:utf8
import bigData
import time
import sys

#Command line skipun til að runna mörg föll í einu og skila útkomu í skrá:
# bash runBigData.sh > Output/out.txt


def runFunction(functionName,funcInput1=-1,funcInput2=-1):
	if funcInput1==-1:
		getattr(bigData,functionName)()
	else:
		getattr(bigData,functionName)(funcInput1,funcInput2)

if __name__ == '__main__':
	start = time.time()
	funcName = sys.argv[1]
	outFolder = sys.argv[2]
	startAtLine = int(sys.argv[3])
	runFunction(funcName,outFolder,startAtLine)
	end = time.time()
	b = open('Output/'+outFolder+'/BF_and_other_info_from_BF_counter.txt', 'a')
	b.write("\n"+str(bigData.returnTime(int(end-start))))
	b.close()
	#print "\n-----------------------------------------------------------------------\n"




	
