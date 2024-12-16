#SessionInfo()
#R version 4.2.1 (2024-09-24)
#Platform: x86_64-apple-darwin17.0 (64-bit)
#Running under: macOS Ventura 13.0

library(treemapify)
#treemapify_2.5.6
library(ggplot2)
#ggplot2_3.4.3 
library(viridis)
#viridis_0.6.4     
#viridisLite_0.4.2 
library(dplyr)
#dplyr_1.1.3 
library(stringr)
#stringr_1.5.0
library(pheatmap)
#pheatmap_1.0.12 
library(dunn.test)
#dunn.test_1.3.6
library(MetBrewer)
#MetBrewer_0.2.0 

##Figure 2C
setwd("/Users/tobybrann/Documents/HGT/polished/new_library/2C/")
treemap <- read.delim(file="treemap_toviz.txt",sep="\t",header=FALSE)

tiff('gen_cov.tiff', units="in", width=10, height=8, res=600, compression = 'lzw')

ggplot(treemap, aes(area = (V3/391395382)*100, fill = (V3/391395382)*100, label=V1)) +
  geom_treemap(color = "black", size=2) +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 25,
                    grow=FALSE,
                    fontface="bold")  +
  labs(colour="Percentage of Genome") +
  scale_fill_viridis(option="viridis")

dev.off()

##Figure 2D
setwd("/Users/tobybrann/Documents/HGT/polished/new_library/2D/")
activity <- read.delim(file="sman_Pr3_Sr3_comps_kim.gff",sep="\t",header=FALSE)

temp<-t(data.frame(strsplit(activity$V9," ")))
temp<-as.data.frame(temp)

activity<-cbind(activity,temp$V2,temp$V3,temp$V4)
colnames(activity)[12] <- "V12"
colnames(activity)[13] <- "V13"
colnames(activity)[14] <- "V14"

activity$V13<-as.numeric(activity$V13)
activity$V14<-as.numeric(activity$V14)

activity$V10 = factor(activity$V10,levels=c("Incomplete","Complete"))

##Lets add our new changes
reg_colours <- c("#156082","#EB7131")

tiff('activity.tiff', units="in", width=12, height=4, res=600, compression = 'lzw')

ggplot(activity, aes(x=V11,fill=V12)) + 
  geom_histogram(bins=75,col=("black")) + 
  facet_wrap(~ factor(V10, levels=c("Incomplete", "Complete"))) + 
  ylab("Count") + xlab("Kimura (with divCpGMod)") +
  scale_fill_manual(values=c("#156082","#EB7131"))

dev.off()

##Figure 2E
setwd("/Users/tobybrann/Documents/HGT/polished/new_library/2E")
TEs.tbl <- read.delim(file="combined_all.tbl",sep="\t",header=FALSE)

ggplot(TEs.tbl, aes(x=V2,y=V1+0.0001)) + 
  geom_hline(yintercept=0.0001,linetype="dashed",col="red") + 
  geom_boxplot() + 
  scale_y_continuous(trans="log10") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

Zero <-TEs.tbl %>% filter(V1==0) 

medians <- TEs.tbl %>% group_by(V2) %>%
  summarise(median = median(V1, na.rm = TRUE))

copy_num <- read.delim(file="TE_copy.nums",sep="\t",header=FALSE)
family <- read.delim(file="TE_family.nums",sep="\t",header=FALSE)

temp <- medians %>% inner_join(copy_num,
                               by=c("V2"))

colnames(family)[1]<- "V2"
colnames(family)[2]<- "V1"

combi_table<-temp %>% inner_join(family,
                                 by=c("V2"))

colnames(combi_table)[1] <- "TE"
colnames(combi_table)[2] <- "Median"
colnames(combi_table)[3] <- "Copy_Number"
colnames(combi_table)[4] <- "Family"

combi_2 <- combi_table %>% filter(Family=="LINE_RTE"|Family=="LTR_Gypsy"|Family=="LINE_Jockey"|Family=="LTR_Pao"|Family=="PLE")

library(MetBrewer)

sum_0 <- TEs.tbl %>% group_by(V2) %>%
  mutate(nn=sum(V1==0))
sum_0 <- cbind(sum_0$V2,sum_0$nn)
sum_0<-as.data.frame(unique(sum_0))
sum_0$V2 <- as.numeric(sum_0$V2)
colnames(sum_0)[1] <- "TE"
combi3 <- combi_table %>% inner_join(sum_0,
                                     by=c("TE"))
colnames(combi3)[5] <- "num0"

combi_3.5 <- combi3 %>% filter(Family=="LINE_RTE"|Family=="LTR_Gypsy"|Family=="LINE_Jockey"|Family=="LTR_Pao"|Family=="PLE"|Family=="LINE_R2")

combi_3.5[10,4]<-"Perere-3"
combi_3.5[27,4]<-"Sr3"

combi_3.5$Family<-factor(combi_3.5$Family,levels=c("Perere-3","Sr3","LINE_RTE","LINE_Jockey","LINE_R2",
                                                   "LTR_Gypsy","LTR_Pao","PLE"))

temp<-c("#156082","#EB7131",met.brewer("Derain",6))

tiff('TBL_num0.tiff', units="in", width=12, height=4, res=600, compression = 'lzw')

ggplot(combi_3.5, aes(x=Copy_Number,y=Median,size=num0,col=Family)) +
  geom_hline(yintercept=0.056, colour="red",linetype="dashed") + 
  scale_colour_manual(values = temp)  +
  geom_point(shape = 19) + 
  scale_x_continuous(trans="log10") + 
  scale_y_continuous() + 
  scale_size(range=c(3,20))+ 
  xlab("Copy Number")

dev.off()

##Figure 3A
setwd("/Users/tobybrann/Documents/HGT/polished/new_library/3A/")
schisto_tbl <- read.delim(file="all_schisto_v2.tbl",header=FALSE,sep="\t")

schisto_tbl$V2 <- gsub("japonicum", "schistosoma_japonicum",schisto_tbl$V2)
schisto_tbl$V2 <- gsub("turkestanicum", "schistosoma_turkestanicum",schisto_tbl$V2)
schisto_tbl$V2 <- gsub("mansoni", "schistosoma_mansoni",schisto_tbl$V2)
schisto_tbl$V2 <- gsub("rodhaini", "schistosoma_rodhaini",schisto_tbl$V2)
schisto_tbl$V2 <- gsub("spindale", "schistosoma_spindale",schisto_tbl$V2)
schisto_tbl$V2 <- gsub("haematobium", "schistosoma_haematobium",schisto_tbl$V2)

schisto_tbl$V2<-factor(schisto_tbl$V2,levels=c("schistosoma_japonicum","schistosoma_turkestanicum",
                                               "schistosoma_mansoni","schistosoma_rodhaini","schistosoma_spindale",
                                               "schistosoma_haematobium"))

tiff('schisto_tbls.tiff', units="in", width=8, height=3, res=600, compression = 'lzw')

ggplot(schisto_tbl, aes(x=V2,y=V1,fill=V2)) + 
  geom_violin() + 
  stat_summary(fun = "mean",
               geom = "point",
               color = "black")+
  ylab("Terminal Branch Length") + 
  xlab("Schistosoma species") + 
  theme(panel.background = element_blank()) +
  #ylim(0,1) + 
  geom_boxplot(width=0.2, outlier.shape = NA, coef = 0)  + 
  stat_summary(fun=mean,geom="point", shape=23, size=3) + 
  theme(legend.position='none')  + 
  ylab("Terminal Branch Length") + 
  xlab("Schistosoma sp.") + 
  ylim(0,0.4) + 
  theme(axis.text = element_text(size=12), axis.title = element_text(size=12, face="bold")) + 
  #scale_y_continuous(trans='log10') +
  scale_fill_manual(values = met.brewer("Derain",8)) 
#ylim(0,0.5) +

kruskal.test(V1~V2, data=schisto_tbl)

# Perform Dunn's test for post-hoc analysis
schisto_tbl$V2 <- gsub("schistosoma_", "",schisto_tbl$V2)

# Bonferroni correction for multiple comparisons
dunn.test(schisto_tbl$V1,schisto_tbl$V2, method = "bonferroni")  

dev.off()

##Figure 3B
setwd("/Users/tobybrann/Documents/HGT/polished/Figure_2/distmat")
snail_schisto<-read.delim(file="all_TEs.out",sep="\t",header=FALSE)

dist_mat <- snail_schisto

dist_mat<-dist_mat[,2:37]
rownames(dist_mat) <- dist_mat[,36]
dist_mat<-dist_mat[,c(1:34)]
dist_mat<-dist_mat[-1,]

library(pheatmap)

for(i in 1:nrow(dist_mat)){
  for(j in 1:ncol(dist_mat)){
    if(i>j){
      dist_mat[i,j]<-dist_mat[j,i]
    }
  }
}
colnames(dist_mat) <- rownames(dist_mat)
##Note to viewer, heatmap was chopped in half for visualisation and row and col heights of trees were edited for visualisation
tiff("all_wTEs_tree_schisto.tiff", units="in", width=8, height=8, res=200, compression = 'lzw')

pheatmap(dist_mat, cluster_cols=TRUE, cluster_rows=TRUE, border_color = "black",
         treeheight_row = 300, treeheight_col = 300, display_numbers = FALSE,
         show_rownames=FALSE,show_colnames=FALSE, legend=TRUE)

dev.off()

mean(dist_mat[,13])
dist_mat[1:13,1:13]










