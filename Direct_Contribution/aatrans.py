#!/usr/bin/env python3

#Fairly simple transformation of coordinates
#Usage:
#argv1= whatever gff, in my case, this was my domain annotations


#Import modules
import sys
import pandas
import csv
import statistics


count=0
cumoverlap=0
cumlength=0
avlength=0
correction=0
temp2=0

#Import and parse the file
#This is just TE annotation
GFFinput = sys.argv[1]
count = len(open(GFFinput).readlines(  ))
lol = list(csv.reader(open(GFFinput, 'rt'), delimiter='\t'))

#Initialise the output
#placeholder=lol2[0][17]
#parent=placeholder.split("=")
#parentID=parent[2]
#sys.stdout = open("unicat.gff", "w")


for i in range(0,count):

	#aa boundaries
	LHS=lol[i][3]
	RHS=lol[i][4]	
	LHS=int(LHS)
	RHS=int(RHS)

	#nt boundaries
	LHSv2=(LHS*3)-2
	RHSv2=(RHS*3)

	smp=lol[i][0]
	source=lol[i][1]
	strand=lol[i][6]	

	#Wrangling with $6	
	meta=str(lol[i][8])
	metav2=meta.replace(" ", "_")
	metav3=metav2.split(";")
	length=len(metav3)

	#Print output for $1-$5
	print(smp, end="\t")
	print(source, end="\t")
	print(LHSv2, end="\t")
	print(RHSv2, end="\t")
	print(strand, end="\t")

	found=0
	end=0
	for i in range(0,length):
		if(i>=3):
			print(metav3[i], end="")
			print(";", end="")
		#if(metav2[i]=="N" and metav2[i+1]=="a" and metav2[i+2]=="m" and metav2[i+3]=="e" and metav2[i+4]=="="):
			#print("Found")
			#found+=1
		#if(found==1 and end==0):
			#if(i-1==length):
				#print(metav2[i])
			#else:
				#print(metav2[i], end="")
			#print(metav2[i],end="")
	print("")
