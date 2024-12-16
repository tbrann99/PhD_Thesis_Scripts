#!/usr/bin/env python3

#Generates 0 indexed genomic coordinates, in which everything is relative to the start of it's own genomic annotation
#Usage
#argv1=TE annotations of the genome (RM output)
#argv2=Intersect output
#We take our genomic annotations, find them and generate a value I've called 'intron_space', which is essentially the correction factor required for the element to be subtracted.


#Import modules
import sys
import pandas
import csv
import statistics
import os

#Set all files in the directory as a list
directory_list = os.listdir()

#Find number of txt files in directory
#for file in directory_list:
#	if(file.endswith("txt")):
#		ticker+=1

#Populate an array with these 0s
#This needs to be n(files)
txtarr=[0]*1

ticker=0

#Populate array with these files
for file in directory_list:
	if(file.endswith("txt")):
	        txtarr[ticker]=file
        	ticker+=1

#print(txtarr)
#i=0
#stripped=txtarr[i].strip("txt")
#gff=stripped+"gff"
#CDSinput=gff
#TEinput=txtarr[i]

#print(CDSinput)
#print(TEinput)


for i in txtarr:
	count=0
	cumoverlap=0
	cumlength=0
	avlength=0
	correction=0
	temp2=0

	#Manipulate txtarr
	stripped=i.strip("txt")
	gff=stripped+"gff"
	CDSinput=gff
	TEinput=i

	#Import and parse the file
	#This is just TE annotation
	#CDSinput = sys.argv[1]
	count = len(open(CDSinput).readlines(  ))
	lol = list(csv.reader(open(CDSinput, 'rt'), delimiter='\t'))

	#Parse TEs
	#This is intersect output
	#TEinput = sys.argv[2]
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

	
	def reverse(LHSTEv2,RHSTEv2,finLHSv2,finRHSv2):
		#Need to find the very end
		#CDSn is number of CDS
		#If CDSn is n(CDS), we can translate this to a pos in CDSarr
		#RHS
		length=RHSTEv2-LHSTEv2
		newRHS=finRHSv2-RHSTEv2
		newLHS=finRHSv2-(RHSTEv2-length)

		return(newRHS,newLHS)


	for i in range(0,count):

		#CDS boundaries
		LHS=lol[i][3]
		RHS=lol[i][4]
		
		LHS=int(LHS)
		RHS=int(RHS)	

		newvals=transform(i,LHS,RHS,intron_space)
		intron_space,LHSv2,RHSv2=newvals

		orientation=lol[i][6]
		if(orientation=="-"):
			#Gene is in the opposite orientation so must count relative to the 3' not 5'
			finLHSv2=LHSv2
			finRHSv2=RHSv2
		
		#print(i+1,LHS,RHS,LHSv2,RHSv2,intron_space)
		populate(i,LHS,RHS, LHSv2, RHSv2, intron_space)
		
	#print(CDSarr)

	#Initialise the output
	placeholder=lol2[0][17]
	parent=placeholder.split("=")
	parentID=parent[2]
	sys.stdout = open(parentID+"_norm.gff", "w")

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

		#Identify the actual overlap coords
		if(LHSTE<LHSCDS):
			LHSTE=LHSCDS
		if(RHSTE>RHSCDS):
			RHSTE=RHSCDS



		#Identify exon location and transform
		for i in range(0,count*7):
			#Iterate through CDS looking for a match
			#print(CDSarr[i], end="\n")
			if(LHSCDS==CDSarr[i]):
				if(RHSCDS==CDSarr[i+1]):
					#print("spotted our matching CDS")
					#This means we're looking at the right CDS, extract intron space
					correction=CDSarr[i+4]
					oldRHS=CDSarr[i+1]
					newRHS=CDSarr[i+3]
					#print(oldRHS, newRHS)
					#print(CDSarr[i+1])
					#print(CDSarr[i+2])
					#print(CDSarr[i+3])


		LHSTEv2=LHSTE-correction
		RHSTEv2=RHSTE-correction

		if(orientation=="-"):		
			#Finv2s are the final CDS numbers (transformed)
			#Take these, with the newly corrected TE values
			newTE=reverse(LHSTEv2,RHSTEv2,finLHSv2,finRHSv2)
			LHSTEv3,RHSTEv3=newTE
			#print(LHSTEv3)
			#print(RHSTEv3)
			LHSTEv2=LHSTEv3
			RHSTEv2=RHSTEv3



		if(RHSTEv2<LHSTEv2):
			tempTE=RHSTEv2
			RHSTEv2=LHSTEv2
			LHSTEv2=tempTE

			




		#print(correction)
		#print("LHSTE",LHSTE)
		#print("RHSTE",RHSTE)
		#print("LHSTEv2",LHSTEv2)
		#print("RHSTEv2",RHSTEv2)


		for y in range(0,18):
			if(y==3):
				#if(LHSTEv2<1):
					#print("1", end="\t")
				#else:
					#print(LHSTEv2, end="\t")
				print(LHSTEv2, end="\t")
			elif(y==4):
				#if(RHSTEv2>newRHS):
					#print(newRHS, end="\t")
				#else:
					#print(RHSTEv2, end="\t")
				print(RHSTEv2, end="\t")
			elif(y==17):
				print(lol2[j][y], end="\n")	
			else:
				print(lol2[j][y], end="\t")
