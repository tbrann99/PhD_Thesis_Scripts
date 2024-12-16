#!/usr/bin/env python3

#Import modules
from __future__ import print_function
import random
import sys
import csv


#Take sys.arg's
#sim = int(sys.argv[1])
output = sys.argv[2]


#Ticker for multiple outputs
reps = 1
repsc = 1

#Do we request multiple?
try:
    reps = int(sys.argv[3])
except Exception:
    pass

#Initialise output
srepsc = str(repsc)
sys.stdout = open(output+srepsc+".gff", "w")

#Number of insertions is number of lines in the input gff
input = sys.argv[1]
sim  = len(open(input).readlines(  ))

#Parse the gff
lol = list(csv.reader(open(input, 'rt'), delimiter='\t'))

#Tickers to initialise
#Variable for the chromosome number
chr=0

#overlap is to decide whether we need to go back and redo the insertion
overlap=0

#Count the number of overlaps so we can keep track
overlapc=0

#tcount is the current line
#sim is the total number of lines required
tcount=0

#Array to populate (up to 3*sim, just incase)
vim=sim*3
chrarr=[0] * vim
lociarr=[0] * vim
loci2arr=[0] * vim

#print(reps)
#print(tcount)
#print("test")


#Function to check for overlap, uses current row, chromosome, both loci
def overlapcheck(tcount,chr,loci,loci2):
	overlap=0
	for j in range(0,tcount):
		#Check only the same chromosome
		if(chrarr[j]==chr):
			#Check if it's an overlap
			if(loci>lociarr[j] and loci<loci2arr[j]):
				#print("overlap loci", end="\t")
				#print(loci)
				#Then overlap on left side
				overlap+=1
				#break
			elif(loci2>lociarr[j] and loci2<loci2arr[j]):
				#print("overlap loci2", end="\t")
				#print(loci2)
				#Then overlap on right side
				overlap+=1
				#break
			#Check if it's engulfed (won't be flagged otherwise)
			elif(loci<lociarr[j] and loci2>loci2arr[j]):
				#print(loci)
				#print(lociarr[j])
				#print(loci2)
				#print(loci2arr[j])
				#print("engulfed", end="\t")
				overlap+=1
				#break
				#print("Chr match")
			elif(lociarr[j]<loci and loci2arr[j]>loci):
				#print("engulfed2", end="\t")
				overlap+=1
			#else:
				#overlap=0
	
#		if(overlap>=0):
#			print(" tcount " + str(tcount), end="\t")
#			print(" j " + str(j), end="\t")	
#			print(" lociarr[j] "+str(lociarr[j]), end="\t")
#			print(" loci2arr[j] "+str(loci2arr[j]), end="\t")
#			print(" loci "+ str(loci), end="\t")
#			print(" loci2 "+ str(loci2))
	return overlap


def ship(chr,loci,loci2,length):
	#print("test", end="\t")
	#print(tcount, end="\t")
	#print(sim, end="\t")
	#print chromosome here too
	print(chr, end="\t")
	print("GenSim", end="\t")
	print("Insertion", end="\t")
	print(loci, end="\t")
	print(loci2, end="\t")
	print(length, end="\t")
	print(".", end="\t")
	strand = random.randint(1, 2)
	if strand == 1:
	    print("+", end="\t")
	else:
	    print("-", end="\t")
	print("0", end="\t")
	print("na")


#if(reps>0):
#print(tcount)
for i in range(repsc, reps+1):
	#print(tcount)
	tcount=0
	while(tcount<sim):
		#print(tcount, end="\t")
		#print(sim, end="\t")
		srepsc = str(repsc)
		sys.stdout = open(output+srepsc+".gff", "w")
		for i in range(tcount, sim):
			loci = random.randint(1, 263000000)
			#loci= random.randint(1,1000000)
			
			#Account for failed overlaps
			#i=i-overlapc

			#Generate length of simulated insertion
			left=lol[tcount][3]
			right=lol[tcount][4]
			#print(left)
			#print(right)
			length=int(right)-int(left)
			#print(length)


			#Identify chromosome, transform loci value            
			if loci < 66567961:
				chr="Schisto_mansoni.Chr_1"
				#print("Schisto_mansoni.Chr_1", end="\t")
			elif loci < 101606849:
				chr="Schisto_mansoni.Chr_2"
				#print("Schisto_mansoni.Chr_2", end="\t")
				loci = loci-66567960
			elif loci < 130038448:
				chr="Schisto_mansoni.Chr_3"
				#print("Schisto_mansoni.Chr_3", end="\t")
				loci = loci-101606848
			elif loci < 162689081:
				chr="Schisto_mansoni.Chr_4"
				#print("Schisto_mansoni.Chr_4", end="\t")
				loci = loci-130038447
			elif loci < 172231772:
				chr="Schisto_mansoni.Chr_5"
				#print("Schisto_mansoni.Chr_5", end="\t")
				loci = loci-162689080
			elif loci < 192607189:
				chr="Schisto_mansoni.Chr_6"
				#print("Schisto_mansoni.Chr_6", end="\t")
				loci = loci-172231772
			elif loci < 202513131:
				chr="Schisto_mansoni.Chr_7"
				#print("Schisto_mansoni.Chr_7", end="\t")
				loci = loci-192607189
			elif loci < 263013206:
				chr="Schisto_mansoni.Chr_ZW"
				#print("Schisto_mansoni.Chr_ZW", end="\t")
				loci = loci-202513131


			#Generate the second location, check that it's not a negative
			#If it is negative we set coordinates to 1 and 1+length
			loci2=loci-length
			if(loci2<1):
				loci=1
				loci2=loci+length
				#loci=loci2
				#loci2=temp
			if(loci>loci2):
				temp=loci2
				loci2=loci
				loci=temp

			#Having identified chromosome and transformed loci, check for overlap
			#print(i)
			#chrarr[i]=chr
			#lociarr[i]=loci
			#loci2arr[i]=loci2
			
			temp=overlapcheck(tcount,chr,loci,loci2)		
			#print(temp)
			#print(overlap)

			#Break the loop if we've made it
			#if(tcount==sim):
				#break

			if(temp==0):
				chrarr[tcount]=chr
				lociarr[tcount]=loci
				loci2arr[tcount]=loci2
				tcount+=1
			else:
				#If there is, underwrite sim
				print("redo!")
				#tcount-=1
				temp=0
				#break
				#overlapc+=1
				#sim+=1
				#print(i)

		for i in range(0,tcount):
			chr=chrarr[i]
			loci=lociarr[i]
			loci2=loci2arr[i]
			length=loci2-loci
			ship(chr,loci,loci2,length)	
	repsc += 1
