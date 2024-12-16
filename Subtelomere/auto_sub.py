#!/usr/bin/env python3

#Import modules
import random
import sys
import csv

###Script to take pvalues per chromosome and to identify the end of enrichment at the chromosome termini
#ARG1 = pvals from wilcoxon sliding windows across chromosomes
#Output to stdout
#In a five window sliding range, with a threshold of 0.1 for enrichment



#Take sys.arg's
pvals = sys.argv[1]
count  = len(open(pvals).readlines(  ))
lol = list(csv.reader(open(pvals, 'rt'), delimiter='\t'))

chroms=[]

##Pull unique chromosomes
for i in range(0,count):
	if(lol[i][0] not in chroms):
		chroms.append(lol[i][0])

#Function to find enrichment at end
def find_subtel(pval_arr):
	on=0
	for i in range(0,len(pval_arr)):
		if(float(pval_arr[i]) < 0.1):
			on=1
		if(float(pval_arr[i]) > 0.1 and on==1):
			return(i-1)
		if(on==1 and i==len(pval_arr)-1):
			return(i-1)

##For each chromosome, pull out the first 40 entries and last 40, corresponds to the last and first 5Mbp
for j in chroms:
	counter=0
	forward=[]
	backward=[]
	all=[]
	locs=[]
	for i in range(0,count):
		if(lol[i][0]==j):
			#print(lol[i])
			counter+=1
			if(counter<=40):
				forward.append(lol[i][5])
			all.append(lol[i][5])
			locs.append(float(lol[i][2]))
	##When we're finished, chop the last 40 out of the all to form backward
	for y in reversed(range(0,len(all))):
		if(len(all)-y<=40):	
			backward.append(all[y])

	#print(locs)
	border = find_subtel(forward)
	#print(border)
	if(border==None):
		pass
	else:
		print(j,end="\t")
		print(int(locs[border]))
	
	##Backwards needs to be treated flipped
	border = find_subtel(backward)
	if(border==None):
		pass
	else:
		print(j,end="\t")
		print(int(locs[len(locs)-border]-100000))

