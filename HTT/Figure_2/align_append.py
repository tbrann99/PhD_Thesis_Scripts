#!/usr/bin/env python3

#Import modules
import sys
import csv
import statistics


count=0
cumoverlap=0
cumlength=0
avlength=0

#Import and parse the file
input = sys.argv[1]
count = len(open(input).readlines(  ))
lol = list(csv.reader(open(input, 'rt'), delimiter=' '))

#2084 17.87 10.36 0.22 SM_V10_1 37811 38225 (87945823) C Sr3#LINE_RTE (2387) 822 366 m_b1s001i78 65
#Kimura (with divCpGMod) = 20.33
#['7821', '6.93', '0.00', '0.10', 'SM_V10_Z', '86690544', '86691554', '(874)', 'Perere-3#LINE_RTE', '2193', '3202', '(5)', 'm_b3012s001i66', '278754']
#['Kimura', '(with', 'divCpGMod)', '=', '6.81']

input2 = sys.argv[2]
count2 = len(open(input2).readlines(  ))
lol2 = list(csv.reader(open(input2, 'rt'), delimiter='\t'))

for i in range(0,count,2):
	hit=0
	chrom=lol[i][4]
	LHS=int(lol[i][5])
	RHS=int(lol[i][6])
	kimura=lol[i+1][4]

	for j in range(0,count2):
		chrom_2=lol2[j][0]
		LHS_2=int(lol2[j][3])
		RHS_2=int(lol2[j][4])
		if(chrom==chrom_2 and LHS==LHS_2 and RHS==RHS_2):
			hit=1
			for y in range(0,len(lol2[j])):
				print(lol2[j][y],end="\t")
	if(hit==1):
		print(kimura)
		hit=0




