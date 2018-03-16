#!/bin/bash
#coding:utf8
python preprocessInput.py t 31
python BF_counter.py t
python whatIsMissing.py t 31

python preprocessInput.py sa 31
python BF_counter.py sa
python whatIsMissing.py sa 31