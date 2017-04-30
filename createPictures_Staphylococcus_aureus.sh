#!/bin/bash
#coding:utf8

#command line syntax:
#python simData.py k       read1           read2      numberOfKmers whatToRun sizeOfGenome maxCov [pfn] [printProgress]             [genomeFile]                      [outDirectory]
#python simData.py 31 Input/t/r1.fastq Input/t/r2.fastq 6000000        1          70000      10   True       True        Input/Staphylococcus_aureus/genome.fa Output/Staphylococcus_aureus


#Create a file with values without stop-filter
#for maxCov in -1 5 10 15 20 30

#Note:  I'm only guessing the size of the genome
#Note2: Some parts of the code don't work without being able to actually create a dict storing all 
#		k-mers from the correct genome
for maxCov in -1
do
	echo "Staphylococcus aureus with maxCov=$maxCov"
	python simData.py 31 Input/Staphylococcus_aureus/frag_1.fastq Input/Staphylococcus_aureus/frag_2.fastq 100000000 1 2903081 $maxCov False False Input/Staphylococcus_aureus/genome.fa Output/Staphylococcus_aureus
	echo " "
done
