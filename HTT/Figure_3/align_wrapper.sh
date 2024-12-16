#!/bin/bash

#What features will be examining

inp=$1

count=$(wc -l < "$inp")

#echo $count

for (( c=1; c<=$count; c++ ))
do
sed "${c}q;d" $inp > ${c}.row.list
./align_sequences_v2.sh ${c}.row.list

rm ${c}.row.list

done




