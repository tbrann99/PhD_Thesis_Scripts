#!/usr/bin/env python3

#Import modules
import sys
import csv
import statistics

#Import and parse the domains
input = sys.argv[1]
count = len(open(input).readlines(  ))
lol = list(csv.reader(open(input, 'rt'), delimiter='\t'))

#Import and parse all TEs
input2 = sys.argv[2]
count2 = len(open(input2).readlines(  ))
lol2 = list(csv.reader(open(input2, 'rt'), delimiter='\t'))

for i in range(0,count2):
    found=0
    size=0
    ##Iterate through TEs
    chrom_1 = lol2[i][0]
    LHS_1 = int(lol2[i][3])
    RHS_1 = int(lol2[i][4])

    ##Does this TE have a domain?
    for j in range(0,count):
        chrom_2 = lol[j][0]
        LHS_2 = int(lol[j][3])
        RHS_2 = int(lol[j][4])

        if(chrom_1==chrom_2 and LHS_1==LHS_2 and RHS_1==RHS_2):
            ##Found TE
            found=1
            doms=lol[j][9]
            num=len(doms.split(","))
            doms=doms.split(",")
            size=0
            for y in range(0,num,2):
                if(float(doms[y+1]))<0.8:
                    size+=1
            
            if(size==0 and num>2):
                print(lol[j][0],lol[j][1],lol[j][2],lol[j][3],lol[j][4],lol[j][5],lol[j][6],lol[j][7],lol[j][8],"Complete",sep="\t")
            #print(chrom_1,LHS_1,RHS_1,chrom_2,LHS_2,RHS_2,doms)


    if(found==0 or size>0):
        print(lol2[i][0],lol2[i][1],lol2[i][2],lol2[i][3],lol2[i][4],lol2[i][5],lol2[i][6],lol2[i][7],lol2[i][8],"Incomplete",sep="\t")



