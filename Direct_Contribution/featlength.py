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

#Import and parse the file
input = sys.argv[1]
count = len(open(input).readlines(  ))

lol = list(csv.reader(open(input, 'rt'), delimiter='\t'))

lengtharr = [0]*count

def findoverlap(LTE,RTE,LFeat,RFeat):
	if(LFeat<RTE and LFeat>LTE):
		if(RFeat<RTE and RFeat>LTE):
			#Means feature engulfed by element
			#Length is just length of feature
			length=int(RFeat)-int(LFeat)
		else:
			#Overlap on right side of TE
			length=int(RTE)-int(LFeat)
	elif(RFeat<RTE and RFeat>LTE):
		length=int(RFeat)-int(LTE)
	elif(LTE>LFeat and RTE<RFeat):
		#Element engulfed by feature
		#Length is just length of element
		length=int(RTE)-int(LTE)
	else:
		length=0
	return(length)		

for i in range(0,count):

	#TE locations
	LTE=lol[i][3]
	RTE=lol[i][4]
	TElength=int(RTE)-int(LTE)
	#print(LTE)
	#print(RTE)	

	#Feature locations
	LFeat=lol[i][12]
	RFeat=lol[i][13]
	#Featlength=abs(int(RFeat)-int(LFeat))

	#cumFeatLength+=FeatLength

	temp=findoverlap(LTE,RTE,LFeat,RFeat)
	#print(temp)
	lengtharr[i]=temp
	cumlength=cumlength+temp


	#Sanity Checks
	#print(FeatLength)
	#print(IntronID)
	#print(Parent)
	#print(LFeat)
	#print(RFeat)

	#Overlap work
	#print(LTE)
	#print(RTE)
	#print(LFeat)
	#print(RFeat)
	
	print(temp)

#Present results
#avlength=cumlength/count
#median = statistics.median(lengtharr)
#formavlength = "{:.2f}".format(avlength)

#print("Total length is", cumlength)
#print("Average length is", formavlength)
#print("Median length is", median)






