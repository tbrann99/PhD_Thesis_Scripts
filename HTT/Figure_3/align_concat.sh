#!/bin/bash

#What features will be examining

tag=$1
file=$2
#cat *.aln > concat.fa

#./fasta_multi_to_one_line.pl concat.fa > concat_1L.fa

grep $tag -A1 $file | grep -v $tag | tr -d '\n'  > ${tag}1L_aln.fa

sed "s/XXXX/$tag/g" fasta.header | cat - ${tag}1L_aln.fa > temp

mv temp ${tag}1L_aln.fa




