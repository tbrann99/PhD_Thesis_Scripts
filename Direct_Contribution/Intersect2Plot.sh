#!/bin/bash

#What features will be examining
ticker=1
ref=$1

mkdir results

print_usage()
{
	printf "Please provide appropriate RM output with -r" 
}


for file in ./input/*
do
	
	#Input (e.g. 2012CDS.gff)
	input=${file##*/}
	

	#Result name (e.g. 2012CDS)
	result=${input%.gff}
	
	#Move files to cwd
	cp ./input/${input} .	

	#echo $input
	#echo $result

	#Bedtools Intersect
	bedtools intersect -wa -wb -a $ref -b $input  > $result.txt
	
	bedo="${result}.txt"
	unique="unique"
	uniqueresults="$results$unique$txt"	
	#echo $uniqureresults

	#Processing Step
	#Numbers per gene output
	grep -v SmtRNA ${result}.txt | awk '{print$21}' | awk -F ';' '{print$2}'|  sort -u > Uni${result}.txt
	
	#Output
	temp= grep -v SmtRNA ${result}.txt | awk '{print$21}' | awk -F';' '{print$2}' | sort -u | wc -l 
	
	temp=${temp//$'\n'/}
	echo "$temp$Uni${result}.txt"
	#"results${ticker}"=$temp
	#echo"$result${ticker}"
	#echo$temp

	#clean up
	rm ${input}
	mv Uni*.txt results
			


done

