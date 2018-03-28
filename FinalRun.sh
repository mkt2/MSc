#!/bin/bash
#coding:utf8
#printf "Starting on t:\n"
#python preprocessInput.py t 31
printf "\nBFAdder:\n"
python BFAdder.py t
printf "\nBFAdderFigures:\n"
python BFAdderFigures.py t
printf "\nwhatIsMissing:\n"
python whatIsMissing.py t 31 False
printf "\nN50:\n"
python N50.py t 31
printf "\npercTable:\n"
python percTable.py t
printf "\nexploreWithBandage:\n"
python exploreWithBandage.py t 31

#printf "\nStarting on sa:\n"
#python preprocessInput.py sa 31
#python BFAdder.py sa
#python BFAdderFigures.py sa
#python whatIsMissing.py sa 31
#printf "\nN50:\n"
#python N50.py sa 31
#printf "\npercTable:\n"
#python percTable.py sa
#printf "\nexploreWithBandage:\n"
#python exploreWithBandage.py sa 31