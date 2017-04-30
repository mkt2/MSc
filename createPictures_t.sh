#!/bin/bash
#coding:utf8

#command line syntax:
#python simData.py k       read1           read2      numberOfKmers whatToRun sizeOfGenome maxCov [pfn] [printProgress] [genomeFile] [outDirectory]
#python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000        1          70000      10   True        True      Input/t/t.fa    Output/t


#Create a file with values without stop-filter
for maxCov in -1 5 10 15 20 30
do
	echo "t with maxCov=$maxCov"
	python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 70000 $maxCov True False Input/t/t.fa Output/t
	echo " "
done

#Create folders with pictures for various maxCov stop-filters
#													            maxCov=			
#python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 70000 1 True False
#python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 70000 3 True False
#python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 70000 5 True False
#python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 70000 6 True False
#python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 70000 10 True False
#echo "t with maxCov=20"
#python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000 1 70000 20 False False