##Workbook for Brann et al. 2024 Horizontal Gene (Transposon) Transfer

##
##Figure 1
##

#Relaxed BlastN from Galbraith et al. 2020
Galbraith, J.D., Ludington, A.J., Suh, A., Sanders, K.L. and Adelson, D.L., 2020. New environment, new invaders—repeated horizontal transfer of LINEs to sea snakes. Genome Biology and Evolution, 12(12), pp.2370-2383.

BLASTN -reward 3 -penalty -4 -xdrop_ungap 80 -xdrop_gap 130 -xdrop_gap_final 150 -word_size 14 -dust yes -gapopen 30 -gapextend 6 -outfmt 6 -evalue 1e-5 -query Smansoni_Perere-3_Sr3.fa -subject other genomes

#Other genomes tested were:
Schistosoma mansoni GCA_000237925.5
Trichobilharzia regenti GCA_944472135.2
Heterobilharzia americana GCA_944470555.2
Fasciola hepatica GCA_948099385.1
Clonorchis sinensis GCA_003604175.2
Bulinus truncatus GCA_021962125.1
Biomphalaria straminea GCA_021533235.1
Lymnaea stagnalis GCA_900036025.1
Homo sapien GCF_000001405.40

##Results were stratified into Perere-3, Sr3 and Other
grep Perere-3 $file | awk '{print$4}' | awk 'a+=$1; END{print a}' | tail -1 | awk '{print$1}'
grep Perere-3 $file | awk '{print$4}' | wc -l | awk '{print$1}'

grep Sr3 $file | awk '{print$4}' | awk 'a+=$1; END{print a}' | tail -1 | awk '{print$1}'
grep Sr3 $file | awk '{print$4}' | wc -l | awk '{print$1}'

grep -v Perere-3 $file | grep -v Sr3 | awk '{print$4}' | awk 'a+=$1; END{print a}' | tail -1 | awk '{print$1}'
grep -v Perere-3 $file | grep -v Sr3 | awk '{print$4}' | wc -l | awk '{print$1}'

##
##Figure 2
##

#Colours
# #156082 - Perere-3
# #EB7131 - Sr3

RepeatMasker -no_is -gff -s -a -lib TE library for S. mansoni 

##Domain vizualisation of each Perere-3 and Sr3 - Panel A
getorf -minsize 500 $i ${i}.orf; done
#This was then used with InterProScan5 web portal (https://www.ebi.ac.uk/interpro/search/sequence/) to identify  L1-EN - cd09076, RT POL - PS50878 and DUF-6451 - PF20049	coordinates
#Coordinates were then used with a bedtools getfasta to extract the specific regions corresponding to the relevant domain
#Pairs of domain sequences originating from Perere-3 and Sr3 were aligned and percentage identity calculated using Jalview's Pairwise Alignment
mafft --reorder --adjustdirection

##Tree of full length copies - Panel B 
grep "Perere-3\|Sr3" RMaskerOutput | awk '{OFS="\t"}{if($12-$11>3200){print$0}}' > FL_copies
bedtools getfasta -bed FL_copies -fi mansoni_genome
##Addition of outgroup 
cat FL_copies Ancohuma > FL_copies_woutgroup
mafft --reorder --adjustdirection
iqtree --seqtype DNA -T AUTO --alrt 1000 --root-seq FILENAME,SR2#LINE_RTE

##Tree map - Panel C
#TE_list.txt = list of TEs in the genome (the unique names)
#Also have to specify genome within treemap_gen.sh
awk '{print$11}' RMaskerOutput | sort -u | awk -F":" '{print$2}' | awk -F"\"" '{print$1}' > TE_list.txt
./treemap_gen.sh TE_list.txt | awk '{print$1}' > treemap_out.txt
./processing.py treemap_out.txt > treemap_toviz.txt
##Where columns represents the name of the element,number of elements, number of bases of the element
##These are then visualised with David Wilkin's Treemapify in R
#https://cran.r-project.org/web/packages/treemapify/vignettes/introduction-to-treemapify.html

##Complete / Incomplete - Panel D
./pfam_query.sh 
#to generate a domain.matrix, which is then analysed by pfam_parse.py to combine domains with respective coordinates
./dom_parse.py Perere-3_Sr3.gff_wdom Perere-3_Sr3.gff > Perere-3_Sr3_wdom.gff
#takes the output of the above and the original annotation file to generate intact gff with an appended "complete" or "incomplete"
##Using the .align file generated from the -a flag of the above RM run
grep "m\_\|imura" ../../../V10_RM/TE_RM/SmansoniV10.fa.align | grep -A1 "Perere-3\|Sr3" | grep -v -- "^--$" > align.reform.txt
./align_append.py align.reform.txt Perere-3_Sr3_wdom.gff

##Terminal Branch Length - Panel E
#For each transposable element (where TE = $line)
grep $line SmansoniV10.fa.out.gff | awk '{if($5-$4>500){print$0}}' | bedtools getfasta -fi schistosoma_mansoni.PRJEA36577.WBPS18.genomic.fa -bed - > ${line}.fa
mafft --auto ${line}.fa > ${line}.aln
trimal -in ${line}.aln -gt 0.01 -keepheader > ${line}.trim.aln
sed 's/:/_/g' ${line}.trim.aln > ${line}.trim.aln.sub
FastTree -nt -gamma < ${line}.trim.aln.sub > ${line}.tree
./schisto_termlength.py ${line}.tree > ${line}.tbl
#We then appended things such as number of 0s, and family which can be found by:
awk '{if($1="0.0"){print$2}' | awk -F".tbl" '{print$1}' | sort " uniq -c

##
##Figure 3
##

##Manual curation was completed inline with Goubert et al. Field Guide to Manual Curation of Transposable elements

##Schistosoma / trematode  phylogeny and terminal branch lengths
#reciprologs was ran on the cambridge computing cluster with the longest isoforms of every gene
#Note, we had to run these on the himem cluster (more GB of RAM than usual) to prevent errors
sbatch -p icelake-himem  --time=12:00:00 --mail-type=FAIL --ntasks 76 --cpus-per-task 1 --wrap "./reciprologs.py -p 304 iso_1L_schistosoma_bovis.TD2_PRJEB44434.WBPS18.protein.fa iso_1L_schistosoma_haematobium.TD2_PRJEB44434.WBPS18.protein.fa iso_1L_schistosoma_mansoni.PRJEA36577.WBPS18.protein.fa iso_1L_schistosoma_rodhaini.TD2_PRJEB44434.WBPS18.protein.fa iso_1L_schistosoma_curassoni.PRJEB44434.WBPS18.protein.fa iso_1L_schistosoma_intercalatum.TD1_PRJEB44434.WBPS18.protein.fa iso_1L_schistosoma_margrebowiei.PRJEB44434.WBPS18.protein.fa iso_1L_schistosoma_spindale.PRJEB44434.WBPS18.protein.fa iso_1L_schistosoma_guineensis.PRJEB44434.WBPS18.protein.fa iso_1L_schistosoma_japonicum.PRJNA520774.WBPS18.protein.fa iso_1L_schistosoma_mattheei.PRJEB44434.WBPS18.protein.fa iso_1L_schistosoma_turkestanicum.PRJEB44434.WBPS18.protein.fa iso_1L_trichobilharzia_regenti.PRJEB44434.WBPS18.protein.fa blastp -o combi_out"

##Output of this first bit is just lists of reciprocal orthologs
##Combined downstream analysis - provides you with a list of reciprocal genes
./align_wrapper.sh calls ./align_sequences_v2.sh per row.list (which is the list of reciprologs)
./align_sequences_v2.sh aligns all reciprologs (in a given row.list)
./cmd_concat.sh, which calls ./align_concat.sh with every genome sampled individually (using their genomic / genic identfier), this results in concatenated amino acid alignments per species
#We then combine the alignments into a single fasta, which is each genome and all reciprologs within it concatenated into a single sequence


##Terminal Branch Lengths
Generated with the same process as figure one, using schisto_schisto_termlength.py to pull out the respective lengths
awk '{if($12-$11>2900){print$0}}' for nums per species


##Analysis of Snails - per species
#n(>2,900bp) was calculated 
awk '{if($12-$11>2900){print$0}}' RMaskerOutput | wc -l 
awk '{print$12-$11}' RMaskerOutput | sort -k1n | tail -1

##Distance matrix
cat schisto_TEs.fa snail_TEs.fa > all_TEs.fa
mafft --auto --adjustdirection all_TEs.fa > all_TEs.aln
distmat -nucmethod 1 all_TEs.aln
#This is in then used in R to produce a heatmap
#Annotations were added after assembly of the figure, as was slicing the heatmap in half


##
##Figure 4
##

##RepBase Comparative Analysis -- Panel A
##RTE sequences downloaded from 
https://www.girinst.org/repbase/update/browse.php?type=RTE&format=FASTA&autonomous=on&nonautonomous=on&simple=on&division=all > repbase_RTE.fa
#Add in our new TEs
cat ../../gastropoda/w_schisto/all_TEs.fa repbase_RTE.fa > repbase_RTE_wcur.fa
mafft --auto --adjustdirection repbase_RTE_wcur.fa > repbase_RTE_wcur.aln
trimal -in repbase_RTE_wcur.aln -out repbase_RTE_wcur_trim.aln -gt 0.25
#Tree built with http://iqtree.cibiv.univie.ac.at/ (iqtree webserver)
#Annotations were done by sequence ID, inline with ITOL documentation
grep AviRTE repbase_RTE_wcur.fa | awk '{print$1}' | awk -F">" '{print$2",branch,node,#COLOUR,2,normal"}'
grep BovB repbase_RTE_wcur.fa | awk '{print$1}' | awk -F">" '{print$2",branch,node,#COLOUR,2,normal"}'


##Metazoan Blast -- Panel B
esearch -db assembly -query '"Metazoa"[Organism] AND ("latest genbank"[filter] AND ("representative genome"[filter] OR "reference genome"[filter]) AND all[filter] NOT anomalous[filter] AND ("1000000"[ScaffoldN50] : "5000000000"[ScaffoldN50])' | esummary | xtract -pattern DocumentSummary -def "NA" -element AssemblyAccession,AssemblyName,FtpPath_GenBank,ScaffoldN50,Organism

#Generate species list
awk -F"\t" '{OFS="\t"}{print$5}' genome_list.tbl | sed 's/_/ /g' > list_of_species.txt
#Input species list into timetree.org and extract tree of the successful entries
#Some species will be dropped
#For these species, simply just grep -v the species from our blastfile
#Schistosoma mansoni (no substitute found) --> sed 's/ (no substitute found)//g' | awk -F"\t" '{print"grep -v "$1" blastfile"}'  
#Many will be renamed, these can be found in the dialog box
#e.g. Schistosoma mansoni (replaced with Schistosoma japonicum)
##Whilst this is generally fine, it's important to keep a linked list of these, to align annotations later with genomes

#Blast the esearch output
#Was ran on the Cambridge High Performance Computer with the icelake architecture
#we split the blastfile into batches of fifty lines, to be ran independently, partly to reduce the simulataneous storage requirements and length of the run
#We called a wrapper to submit all 68 files concurrently, which submits ./batch_blast.sh for each respective submit on the splitfile
./batch_submit_wrapper.sh list_of_split_files
#./batch_blast.sh iterates through each file and downloads the genome, blasts Perere-3 and Sr3 against it, saves the blast table and deletes the genome file
./blastn_analysis.sh blasttable 
##Generates numbers for each blast table, for Perere-3 and Sr3 (seperately)
##Finds the largest hit and outputs the: percentage ID, length, evalue 
##Also finds the number of hits over 1,000bp (which is later superseded by multibar_plotter.py)
./multibar_plotter.py blastn_summary
##Multibar plotter translates the csvs into a Yes1 or Yes2 (strong or weak hit), which in turn gets translated into an ITOL suitable annotation 

##
##Figure 5
##

mafft --auto --adjustdirection

iqtree -bb 1000 -wbt -alrt 1000

