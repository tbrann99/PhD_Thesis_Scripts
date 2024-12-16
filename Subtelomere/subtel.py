#!/usr/bin/env python3

###Script takes paired loci and returns the list with a tag corresponding to if it is subtelomeric or chromosome body (using bounds from subtelomere.bed), can be adjusted for other purposes
###Like bedtools intersect but returning a tagged list (which can then be differentially labeled for visualisation)
#ARG1 = file of regions to investigate / tag. In this case, paired loci from biser
#ARG2 = bedfile of which we compare, if something in ARG1 is in ARG2, tagged as "Subtelomere"

#Import modules
import sys
import csv
import statistics

count=0
checker=0



#Pairs
input = sys.argv[1]
count = len(open(input).readlines(  ))
lol = list(csv.reader(open(input, 'rt'), delimiter='\t'))

#Subtelomere Bounds
input2 = sys.argv[2]
count2 = len(open(input2).readlines(  ))
lol2 = list(csv.reader(open(input2, 'rt'), delimiter='\t'))

def ship_sub(i):
    for j in range(0,(len(lol[i]))):
        if(j+1==len(lol[i])):            
            print(lol[i][j], end="\t")
            print("Subtelomere")
        else:
            print(lol[i][j],end = "\t")

def ship(i):
	for j in range(0,(len(lol[i]))):
		if(j+1==len(lol[i])):
			print(lol[i][j],end = "\t")
			print("No")            
		else:
			print(lol[i][j],end = "\t")

def subtels(chrom,LHS,RHS):
    for j in range(0,count2):
        
		##if match to struc
        str_chr = lol2[j][0]
        str_LHS = int(lol2[j][1])
        str_RHS = int(lol2[j][2])
        #print(str_chr,chrom)
        if(str_chr==chrom and LHS<=str_RHS and RHS>=str_LHS):
            #print("hit")
            ##Overlap
            return(1)
        #else:
            #return(0)


for i in range(0,count):
    ##Iterate through pairs
    chrom=lol[i][0]
    LHS=int(lol[i][1])
    RHS=int(lol[i][2])

    func = subtels(chrom,LHS,RHS)

    #print(func)

    if(func==1):
        checker+=1

    chrom2=lol[i][3]
    LHS2=int(lol[i][4])
    RHS2=int(lol[i][5])

    func = subtels(chrom2,LHS2,RHS2)

    #These counters check both loci, script can be very easily changed to single loci by removing the second half of counters
    if(func==1):
        checker+=1

    #If either loci was subtelomeric, we then print subtelomere
    if(checker>=1):
        ship_sub(i)
    else:
        #print("no sub")
        ship(i)

    checker=0




