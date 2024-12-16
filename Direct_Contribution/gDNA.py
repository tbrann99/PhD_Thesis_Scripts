#!/usr/bin/env python3

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
CDSinput = sys.argv[1]
count = len(open(CDSinput).readlines(  ))
lol = list(csv.reader(open(CDSinput, 'rt'), delimiter='\t'))

#Parse TEs
#This is intersect output
TEinput = sys.argv[2]
count2 = len(open(TEinput).readlines(  ))
lol2 = list(csv.reader(open(TEinput, 'rt'), delimiter='\t'))

lengtharr = [0]*count
CDSarr=[0]*count*7
counter=1
CDSn=1
intron_space=0
ilength=0
CDSarr[0]=0
temp=0

def transform(i,LHS,RHS,intron_space):
	counter=i+1
	#print("Counter is:", counter)
	temp=(counter-1)*6
	if(counter==1):
		intron_space=LHS-1
		#LHSv2 SHOULD be 1
		LHSv2=LHS-intron_space
		RHSv2=LHSv2+(RHS-LHS)
	elif(counter):
		ilength=LHS-CDSarr[temp-3]
		intron_space=intron_space+ilength-1
		#print(CDSarr[temp-3])
		LHSv2=CDSarr[temp-1]+1
		RHSv2=LHSv2+(RHS-LHS)
	
	return(intron_space,LHSv2,RHSv2)

def populate(i,LHS,RHS,LHSv2,RHSv2,intron_space):
	CDSn=i+1
	counter=CDSn
	if(CDSn==1):
		#print("CDSn==1")
		CDSarr[counter]=CDSn
		CDSarr[counter+1]=LHS
		CDSarr[counter+2]=RHS
		CDSarr[counter+3]=LHSv2
		CDSarr[counter+4]=RHSv2
		CDSarr[counter+5]=intron_space
	else:
		counter=(counter-1)*6
		CDSarr[counter+1]=CDSn
		CDSarr[counter+2]=LHS
		CDSarr[counter+3]=RHS
		CDSarr[counter+4]=LHSv2
		CDSarr[counter+5]=RHSv2
		CDSarr[counter+6]=intron_space



for i in range(0,count):

	#CDS boundaries
	LHS=lol[i][3]
	RHS=lol[i][4]
	
	LHS=int(LHS)
	RHS=int(RHS)	

	newvals=transform(i,LHS,RHS,intron_space)
	intron_space,LHSv2,RHSv2=newvals

	#print(i+1,LHS,RHS,LHSv2,RHSv2,intron_space)
	populate(i,LHS,RHS, LHSv2, RHSv2, intron_space)
	
#print(CDSarr)

#Initialise the output
placeholder=lol2[0][17]
parent=placeholder.split("=")
parentID=parent[2]
sys.stdout = open(parentID+".gff", "w")

for j in range(0,count2):
	
	#TE overlap boundaries
	LHSTE=lol2[j][3]
	RHSTE=lol2[j][4]

	LHSTE=int(LHSTE)
	RHSTE=int(RHSTE)
	#print(LHSTE, RHSTE)
	
	#Which CDS do they fall in?
	LHSCDS=lol2[j][12]
	RHSCDS=lol2[j][13]	

	#print(LHSCDS)
	#print(RHSCDS)
	LHSCDS=int(LHSCDS)
	RHSCDS=int(RHSCDS)

	#Identify exon location and transform
	for i in range(0,count*7):
		#Iterate through CDS looking for a match
		#print(CDSarr[i], end="\n")
		if(LHSCDS==CDSarr[i]):
			if(RHSCDS==CDSarr[i+1]):
				#print("spotted our matching CDS")
				#This means we're looking at the right CDS, extract intron space
				correction=CDSarr[i+4]
				#print(CDSarr[i+1])
				#print(CDSarr[i+2])
				#print(CDSarr[i+3])


	LHSTEv2=LHSTE-correction
	RHSTEv2=RHSTE-correction

	#print(correction)
	#print("LHSTE",LHSTE)
	#print("RHSTE",RHSTE)
	#print("LHSTEv2",LHSTEv2)
	#print("RHSTEv2",RHSTEv2)


	for y in range(0,18):
		if(y==3):
			print(LHSTEv2, end="\t")
		elif(y==4):
			print(RHSTEv2, end="\t")
		elif(y==17):
			print(lol2[j][y], end="\n")	
		else:
			print(lol2[j][y], end="\t")
