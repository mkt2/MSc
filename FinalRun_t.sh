#!/bin/bash
#coding:utf8
python preprocessInput.py t 31
printf "\n---BFAdder:---\n"
python BFAdder.py t
printf "\n---BFAdderFigures:---\n"
python BFAdderFigures.py t
printf "\n---whatIsMissing:---\n"
python whatIsMissing.py t 31 False
printf "\n---N50:---\n"
python N50.py t 31
printf "\n---percTable:---\n"
python percTable.py t
printf "\n---exploreWithBandage:---\n"
python exploreWithBandage.py t 31