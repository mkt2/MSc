#!/bin/bash
echo "Enter your k: "
read k
echo "Enter your first read: "
read fn1
echo "Enter your second read: "
read fn2
#python bloom.py
python BFCounter_naive.py "$k" "$fn1"
python BFCounter_naive.py "$k" "$fn2"

#Check if the two reads give the same output.
#read4.fq and read8.fq should do that

#Now compare the two dbg's
gfa1=${fn1/.*/_dbg.gfa}
gfa2=${fn2/.*/_dbg.gfa}
echo $gfa1
echo $gfa2

python compare_gfa.py $gfa1 $gfa2


#I'm currently running this with:
#k=3
#fn1=read4.fq
#fn2=read8.fq