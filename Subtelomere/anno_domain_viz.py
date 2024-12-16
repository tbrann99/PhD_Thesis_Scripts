#!/usr/bin/env python3

##Processing the TSV generated from Interproscan
##Due to the finite number of protein domains sampled, domains, colours, shapes were all predetermined. 
##Generates output that is used as an "additional dataset" on ITOL for visualisation

#Import modules
import sys
import csv
import statistics

#Import and parse the file
#This is just IPR TSV annotation
CDSinput = sys.argv[1]
count = len(open(CDSinput).readlines(  ))
lol = list(csv.reader(open(CDSinput, 'rt'), delimiter='\t'))


##Total length of the protein
lens = sys.argv[2]
count2= len(open(lens).readlines(  ))
lol2=list(csv.reader(open(lens, 'rt'), delimiter='\t'))

##domain_array

col="col"
shape="shape"


if(count>0):
	for i in range(0,count2):
		#print(lol[0][0])
		#print(lol2[i][0])	
		if(lol[0][0]==lol2[i][0]):
			print(lol2[i][0], end=",")
			print(lol2[i][1], end=",")

	for i in range(0,count):
		###Domain ID
		if(lol[i][4]=="TMhelix"):
			##Red Rectangle
			col="#FF0000"
			shape="RE"
		if(lol[i][4]=="SignalP-noTM"):
			##Blue Vertical Rectangle
			col="#3954BC"
			shape="HH"
		if(lol[i][4]=="SignalP-TM"):
			##Darker Blue
                        col="#241571"
                        shape="HH"
		if(lol[i][4]=="CYTOPLASMIC_DOMAIN"):
			shape="EL"
			col="#FF9933"
		if(lol[i][4]=="NON_CYTOPLASMIC_DOMAIN"):
			shape="EL"
			col="#AFB83B"
		if(lol[i][4]=="TRANSMEMBRANE"):
			shape="DI"
			col="#FF8D85"
		if(lol[i][4]=="Coil"):
			shape="OC"
			col="#884DFF"
		if(lol[i][4]=="mobidb-lite"):
                        shape="PD"
                        col="#B8B8B8"
		
			

		###Print Entry

		###Shape
		if(shape=="shape"):
			print(lol[i][4])
		print(shape,end="|")		

		###Loc
		print(lol[i][6],end="|")
		print(lol[i][7],end="|")

		###Col
		if(count-1==int(i)):
			print(col)
		else:
			print(col,end=",")
			



