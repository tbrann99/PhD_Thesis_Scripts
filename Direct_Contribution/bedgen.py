#!/usr/bin/env python3

#Specifically for S. mansoni
#Generate windows of X, across chromosoems of the genome

#Import modules
import sys
import csv

#>SM_V9_1
#87984036
#>SM_V9_2
#45716228
#>SM_V9_3
#49793703
#>SM_V9_4
#46471565
#>SM_V9_5
#24128855
#>SM_V9_6
#24696071
#>SM_V9_7
#19942475
#>SM_V9_PAR1
#10680109
#>SM_V9_PAR2
#42949100
#>SM_V9_ZSR
#33063208
#>SM_V9_WSR
#5969971

#Initialise output
sys.stdout = open("newbed.bed", "w")


def ship(chrom,LHS,RHS):
	print(chrom, end="\t")
	print(LHS, end="\t")
	print(RHS)
        #print(".", end="\t")
        #print(".", end="\t")
       
chr=["SM_V9_1","SM_V9_2","SM_V9_3","SM_V9_4","SM_V9_5","SM_V9_6","SM_V9_7","SM_V9_PAR1","SM_V9_PAR2","SM_V9_ZSR","SM_V9_WSR"]
chrlengths=[87984036,45716228,49793703,46471565,24128855,24696071,19942475,10680109,42949100,33063208,5969971]
chrticker=0
for i in chrlengths:
	for j in range(0,i,100):
		LHS = int(j)
		RHS = LHS+100

		if(RHS > i):
			RHS=int(i)
		chrom=chr[chrticker]
		ship(chrom,LHS,RHS)	
	#print(chr[chrticker])
	#print(i)
	chrticker+=1
	

