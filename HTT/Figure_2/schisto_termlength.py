#!/usr/bin/env python3

#Import modules
import sys
import csv
import statistics

#Import and parse the file
input = sys.argv[1]

active=0
tbl=0


cooldown=0

def end(coord):
	for y in range(0,14):
		if(newick[coord+y]=="," or newick[coord+y]==")"):
			for x in range(1,y):
				print(newick[coord+x],end="")
			#tbl=newick[coord,coord:coord+(y-1)]
	print("")

def term(i):
	#print(newick[i],"this is i")
	found=0
	for j in range(0,25):
		if(newick[i+j]==":"):
			length_site=i+j
			if(found==0):
				end(length_site)
				found=1
			


with open(input, "r") as a_file:
	for line in a_file:
		newick=line.strip()
		##Stripped line is my chars to iterate through
		for i in range(0,len(newick)):
			if(newick[i]=="_"):
				if(cooldown==0):
					term(i)
					cooldown=6
			if(cooldown>0):
				cooldown-=1

