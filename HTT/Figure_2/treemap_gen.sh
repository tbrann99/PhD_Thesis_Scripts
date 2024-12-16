#!/bin/bash

file=$1
lines=$(cat $file)
for line in $lines
do
	echo $line
	grep $line new_lib_RMout.gff | wc -l
	grep $line new_lib_RMout.gff | awk '{print$12-$11}' | awk -F: '{total+=$1} END{print total}'

done
