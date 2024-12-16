#!/usr/bin/env python3

#Import modules
import sys
import csv
import statistics

#Import and parse the file
input = sys.argv[1]
count = len(open(input).readlines(  ))
lol = list(csv.reader(open(input, 'rt'), delimiter='\t'))

def answer(result_arr):
    #print(result_arr)
    if(len(result_arr)>2):
        #return("yes")
        numover=int(result_arr[3])
        maxlen=int(result_arr[1])
        if(numover>5 and maxlen>2000):
            return("Yes2")
        elif(numover>1 and maxlen>1000):
            return("Yes1")
        else:
            return("No")
    else:
        return("No")

def send(name,P3,S3):
    ##Yes 2 is 1, Yes 1 is 2, No is 3
    print(name.replace(" ", "_"),end=",")

    if(P3=="No" and S3=="No"):
        print("0,0,0,0,0")
    elif(P3=="No"):
        print("0,0,100",end=",")
    elif(P3=="Yes2"):
        print("100,0,0",end=",")
    elif(P3=="Yes1"):
        print("0,100,0",end=",")

    if(S3=="Yes2"):
        print("100,0")
    elif(S3=="Yes1"):
        print("0,100")
    elif(S3=="No" and P3!="No"):
        print("0,0")

def push(TE,i):
    ##Selected best TE
    #print("")
    ##Extract names
    print(lol[i][3].split(" ")[0],end="")
    if(len(lol[i][3].split(" "))>1):
        print("_",end="")
        print(lol[i][3].split(" ")[1],end=",")
    else:
        print(",", end="")

    if(int(TE[1])>2000 and int(TE[3])>4):
        print("100,0")
    else:
        print("0,100")
"""     print(TE[0],end="\t")
    print(TE[1],end="\t")
    print(TE[2],end="\t")
    print(TE[3],end="\t")
    print(lol[i][0],end="\t")
    print(lol[i][1],end="\t")
    print(lol[i][2],end="\t")
    print(lol[i][3]) """


for i in range(0,count):
    name=lol[i][3]
    P3=lol[i][4].split(",")
    S3=lol[i][5].split(",")

    #S3_out=answer(S3.split(","))
    #P3_out=answer(P3.split(","))
    if(len(P3)>2 or len(S3)>2):
        if(int(P3[1])>1500 or int(S3[1])>1500):
            if(len(P3)>2 and len(S3)>2):
                if(int(P3[1])>int(S3[1])):
                    push(P3,i)
                elif(int(P3[1])==int(S3[1])):
                    if(float(P3[0])>float(S3[0])):
                        #print(P3)
                        push(P3,i)
                    elif(float(P3[0])==float(S3[0])):
                        #print(P3,S3)
                        push(P3,i)
                    else:
                        push(S3,i)
                        #print(S3)
                else:
                    push(S3,i)  
                    #print(P3,S3)


        #send(name,P3_out,S3_out)
