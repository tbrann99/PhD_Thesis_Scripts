#!/bin/bash

#What features will be examining

file=$1
lines=$(cat $file)

link=`awk '{print$3}' $file`

#echo $link

tail_fa="_genomic.fna.gz"


for line in $link
do
	ftp=$line

	acc=`echo $ftp | cut -d '/' -f 10`

	genfile="${acc}${tail_fa}"

	genfile_td="${ftp}"/"${genfile}"

	wget $genfile_td

	genfile_uz=${genfile//.gz}

	gunzip $genfile

	blastn -query Sr3_Perere-3.fa -subject $genfile_uz -outfmt 6 -evalue 0.00002 -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 10 -dust yes -gapopen 30 -gapextend 6 > ${genfile}_blastout.tbl

	rm $genfile_uz

done


