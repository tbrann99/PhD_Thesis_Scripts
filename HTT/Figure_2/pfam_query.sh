#!/bin/bash

sequence=$1
coords=$2
min_size=$3

awk -v min=$min_size '{OFS="\t"}{if($12-$11>min){print$0}}' $coords > temp.bed 

bedtools getfasta -fi $sequence -bed temp.bed > temp.fa

getorf -sequence temp.fa -outseq temp.orf -minsize 300

pfam_scan.pl -fasta temp.orf -dir /Users/tobybrann/Documents/HGT/polished/new_library/2D/w_dom/libraries > pfam.results

##If we want to filter the results, we have to do it here
#grep Domain pfam.results | awk '{print$1}' | sort -u | awk -F"_" '{print$1}' | sort | uniq -c | awk -F"-" '{print$1,$2}' | awk '{OFS="\t"}{print$2,$3,$1}' > domain.matrix
grep Domain pfam.results | awk '{OFS="\t"}{print$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' > domain.matrix

./pfam_parse.py domain.matrix $coords > ${coords}_wdom

#awk '{if ($6~/^PF/) {print $1}}' < pfam.results |  sed 's/_/\//2;s/_/ /2' | awk '{print $1}' | sort > pf.domains.count

