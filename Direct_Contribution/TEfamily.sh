#!/bin/bash

for file in ./input/*
do

	#Input (e.g. 2012CDS.gff)
        input=${file##*/}


        #Result name (e.g. 2012CDS)
        result=${input%.gff}

	cp ./input/${input} .

	awk '{print$10}' ${result} | awk -F':'  '{print$2}' | awk -F'"' '{print$1}' | sort | uniq -c | awk -v var="${result}" '{print$0, var }' > TE${result}


done	
