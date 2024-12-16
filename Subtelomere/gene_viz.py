#!/usr/bin/env python3

#Import modules
import sys
import csv


##Take the gene gff and loci and find the numbers for CDS relative to 0, the start of the first coding sequence, irrespective of orientation.
##Used to produce input for gggenes package

#Import and parse the alignment
###Input is the CDS annotation for a respective smp
input = sys.argv[1]
count = len(open(input).readlines(  ))
lol = list(csv.reader(open(input, 'rt'), delimiter='\t'))

##Allows for combining of genes and use of a wrapper script
smp=sys.argv[2]

revarr=[]

rev=0

for i in range(0,count):
    ##Is reverse?
    if(lol[i][6]=="+"):
        ispace=int(lol[0][3])
        for j in range(0, len(lol[i])):
            if(j==len(lol[i])-1):
                print(lol[i][j])
            elif(j==0):
                print(smp,end="\t")
                print("exon_"+str(i),end="\t")
            elif(j==3 or j==4):
                print(int(lol[i][j])-ispace+1, end="\t")
            else:
                print(lol[i][j],end="\t")
    else:
        ##Reverse
        rev=1
        ispace=int(lol[count-1][4])
        revarr.append(abs((int(lol[i][3])-ispace)-1))
        revarr.append(abs((int(lol[i][4])-ispace)-1))
        ##Add adjusted coords to an array

y=len(revarr)-1

##Flip a - orientation gene model
if(rev==1):
    for i in range(0,count):
        for j in range(0, len(lol[i])):
            if(j==len(lol[i])-1):
                print(lol[i][j])
            elif(j==0):
                print(smp,end="\t")
                print("exon_"+str(i),end="\t")
            elif(j==3 or j==4):
                print(revarr[y], end="\t")
                y=y-1
            else:
                print(lol[i][j],end="\t")

