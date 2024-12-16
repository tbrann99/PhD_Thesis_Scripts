#!/bin/bash

#What features will be examining

file=$1
lines=$(cat $file)


for line in $lines
do
	cp sublists/$line .	

	sbatch -p icelake  --time=12:00:00 --mail-type=FAIL --cpus-per-task 8 --wrap "./batch_blast.sh $line"

done


