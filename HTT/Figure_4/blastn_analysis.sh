#!/bin/bash

#What features will be examining

file=$1
cat $file | while read line
do
	LHS=`echo $line | tr '\t' ' ' |cut -d' ' -f1`
	RHS=`echo $line | tr '\t' ' ' |cut -d' ' -f2`

	#echo "${LHS}_${RHS}_genomic.fna.gz_blastout.tbl"

	blast_out=`echo "${LHS}_${RHS}_genomic.fna.gz_blastout.tbl"`

	Perere=`awk '{if($1=="Perere-3#LINE_RTE"){print$0}}' $blast_out | sort -k4n | tail -1 | awk '{OFS=","}{print$3,$4,$11}'`
	Pereren=`awk '{if($1=="Perere-3#LINE_RTE"){print$0}}' $blast_out | awk '{if($4>1000){print$0}}'    | wc -l`
	#echo $Pereren

	Sr=`awk '{if($1=="Sr3#LINE_RTE"){print$0}}' $blast_out | sort -k4n | tail -1 | awk '{OFS=","}{print$3,$4,$11}'`
	Srn=`awk '{if($1=="Sr3#LINE_RTE"){print$0}}' $blast_out | awk '{if($4>1000){print$0}}'    | wc -l`	
	#echo $Srn

	echo "${line}    ${Perere},${Pereren}    ${Sr},${Srn}"

done

