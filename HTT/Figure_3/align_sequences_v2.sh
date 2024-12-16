#!/bin/bash

#What features will be examining

inp=$1

bov1="iso_1L_schistosoma_bovis.TD2_PRJEB44434.WBPS18.protein.fa"
cur2="iso_1L_schistosoma_curassoni.PRJEB44434.WBPS18.protein.fa"
gui3="iso_1L_schistosoma_guineensis.PRJEB44434.WBPS18.protein.fa"
hae4="iso_1L_schistosoma_haematobium.TD2_PRJEB44434.WBPS18.protein.fa"
int5="iso_1L_schistosoma_intercalatum.TD1_PRJEB44434.WBPS18.protein.fa"
jap6="iso_1L_schistosoma_japonicum.PRJNA520774.WBPS18.protein.fa"
man7="iso_1L_schistosoma_mansoni.PRJEA36577.WBPS18.protein.fa"
mar8="iso_1L_schistosoma_margrebowiei.PRJEB44434.WBPS18.protein.fa"
mat9="iso_1L_schistosoma_mattheei.PRJEB44434.WBPS18.protein.fa"
rod10="iso_1L_schistosoma_rodhaini.TD2_PRJEB44434.WBPS18.protein.fa"
spi11="iso_1L_schistosoma_spindale.PRJEB44434.WBPS18.protein.fa"
tur12="iso_1L_schistosoma_turkestanicum.PRJEB44434.WBPS18.protein.fa"
reg13="iso_1L_trichobilharzia_regenti.PRJEB44434.WBPS18.protein.fa"

awk '{print$1}' $inp | grep -A1 -f - proteins_fasta/${bov1} > ${inp}_1
awk '{print$2}' $inp | grep -A1 -f - proteins_fasta/${cur2} > ${inp}_2
awk '{print$3}' $inp | grep -A1 -f - proteins_fasta/${gui3} > ${inp}_3
awk '{print$4}' $inp | grep -A1 -f - proteins_fasta/${hae4} > ${inp}_4
awk '{print$5}' $inp | grep -A1 -f - proteins_fasta/${int5} > ${inp}_5
awk '{print$6}' $inp | grep -A1 -f - proteins_fasta/${jap6} > ${inp}_6
awk '{print$7}' $inp | grep -A1 -f - proteins_fasta/${man7} > ${inp}_7
awk '{print$8}' $inp | grep -A1 -f - proteins_fasta/${mar8} > ${inp}_8
awk '{print$9}' $inp | grep -A1 -f - proteins_fasta/${mat9} > ${inp}_9
awk '{print$10}' $inp | grep -A1 -f - proteins_fasta/${rod10} > ${inp}_10
awk '{print$11}' $inp | grep -A1 -f - proteins_fasta/${spi11} > ${inp}_11
awk '{print$12}' $inp | grep -A1 -f - proteins_fasta/${tur12} > ${inp}_12
awk '{print$13}' $inp | grep -A1 -f - proteins_fasta/${reg13} > ${inp}_13

cat ${inp}_1 ${inp}_2 ${inp}_3 ${inp}_4 ${inp}_5 ${inp}_6 ${inp}_7 ${inp}_8 ${inp}_9 ${inp}_10 ${inp}_11 ${inp}_12 ${inp}_13 > ${inp}.fa

rm ${inp}_1 ${inp}_2 ${inp}_3 ${inp}_4 ${inp}_5 ${inp}_6 ${inp}_7 ${inp}_8 ${inp}_9 ${inp}_10 ${inp}_11 ${inp}_12 ${inp}_13

mafft --auto --quiet ${inp}.fa > ${inp}.aln

mv ${inp}.aln ./alignments/${inp}.aln
mv ${inp}.fa ./fasta/${inp}.fa





