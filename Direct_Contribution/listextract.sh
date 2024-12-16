#!/bin/bash

#FILENAME=$1

#CDSlist=$1

#List of all my Smp IDs
#awk '{print$9}' $CDSlist | awk -F'=' '{print$3}' > temp2.txt

file=$1
lines=$(cat $file)

for line in $lines
do
        #echo ${line}
	gff="${line}.gff"
	txt="${line}.txt"
	cp ./anno/${gff} ./filteredfiles
	cp ./ints/${txt} ./filteredfiles
		
done

