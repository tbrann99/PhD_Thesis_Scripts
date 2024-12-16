#!/usr/bin/env python3

#Import modules
import sys
import csv
import statistics

ticker=0

#Import and parse the RepeatMasker
#This is just TE annotation globally
input = sys.argv[1]
count = len(open(input).readlines(  ))
lol = list(csv.reader(open(input, 'rt'), delimiter='\t'))

for i in range(0,count,3):
	print(lol[i][0],end="\t")
	print(lol[i+1][0],end="\t")
	print(lol[i+2][0])
