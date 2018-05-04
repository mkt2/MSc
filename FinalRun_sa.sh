#!/bin/bash
#coding:utf8
python preprocessInput.py sa 31
printf "\n---BFAdder:---\n"
python BFAdder.py sa
printf "\n---BFAdderFigures:---\n"
python BFAdderFigures.py sa
printf "\n---whatIsMissing:---\n"
python whatIsMissing.py sa 31
printf "\nN50:\n"
python N50.py sa 31
printf "\npercTable:\n"
python percTable.py sa
#printf "\n---exploreWithBandage:---\n"
#python exploreWithBandage.py sa 31             #<---naive sprengir minnid
#---exploreWithBandage:---
#Starting on G
#contigsInG: 817316
#kmersInG:   4185369
#Starting on G_naive using a BF for filtering singletons:
#print_GFA_to_file(fileName=Output/sa/G_inf_naive3.GFA)
#Starting on G_naive using a dict for filtering singletons
#FinalRun_sa.sh: line 15: 11019 Killed                  python exploreWithBandage.py sa 31