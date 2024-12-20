library(ggplot2)
library(MetBrewer)
library(ggpubr)



setwd("/Users/tobybrann/Documents/Subtelomeres/thesis/rep_homo/")
TEs<-read.delim(file="TE_parsed.out",sep="\t",header=FALSE)
WEs<-read.delim(file="WE_parsed.out",sep="\t",header=FALSE)
TRF<-read.delim(file="TRF_parsed.out",sep="\t",header=FALSE)

combi<-rbind(TEs,WEs,TRF)

combi$V2 <-factor(combi$V2, levels =c("TE","WE","TRF"))

tiff('deeptools.tiff', units="in", width=10, height=12, res=600, compression = 'lzw')

ggplot(combi,aes(x=V3,y=V1,col=V2)) + 
  scale_color_manual(values=c(met.brewer("Egypt",4))) + 
  facet_wrap(.~V2, scales = "free_y",nrow=3) +
  theme_Publication(base_family = "Arial") + 
  theme(legend.position="none",
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_blank()) + 
  geom_hline(
    data = subset(combi, V2 == "TRF"), # Subset data for specific facet
    aes(yintercept = 0.00891),
    linetype = "dashed",
    color = "red"
  ) + 
  geom_hline(
    data = subset(combi, V2 == "TE"), # Subset data for specific facet
    aes(yintercept = 0.335),
    linetype = "dashed",
    color = "red"
  ) +
  geom_hline(
    data = subset(combi, V2 == "WE"), # Subset data for specific facet
    aes(yintercept = 0.091),
    linetype = "dashed",
    color = "red"
  ) +
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 700) +
  geom_line(linewidth=1.5) 

dev.off()




theme_Publication <- function(base_size=14, base_family="helvetica") {
  library(grid)
  library(ggthemes)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.margin = unit(0, "cm"),
            legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))
  
}

##MA Plots
##Gene defining


###Hypothetical
setwd("/Users/tobybrann/Documents/Subtelomeres/RNASeq/tags")
hypo<-read.delim(file="egg.txt", header=FALSE)   

hypo <- which(row.names(returned_data) %in% hypo$V1)
for(i in 1:length(hypo)){
  num<-(hypo[i])
  #writeLines(num)
  returned_data_2$Is_Sub[num] <- "Hypothetical Protein"
}

annex<-read.delim(file="annexin.txt", header=FALSE)   

annex <- which(row.names(returned_data) %in% annex$V1)
for(i in 1:length(annex)){
  num<-(annex[i])
  #writeLines(num)
  returned_data_2$Is_Sub[num] <- "Annexin"
}

EG<-read.delim(file="EG.txt", header=FALSE)   

EG <- which(row.names(returned_data) %in% EG$V1)
for(i in 1:length(EG)){
  num<-(EG[i])
  #writeLines(num)
  returned_data_2$Is_Sub[num] <- "Endoglycoceramidase"
}


###MA Plots
setwd("/Users/tobybrann/Documents/Subtelomeres/manuscript/Figure_3/expression/sarah_BB/MA_plots/")

bbcounts <- read.delim(file="BB_counts.txt",header=TRUE,sep="\t")
coldata <- read.delim(file="coldata.txt",header=TRUE,sep="\t")

###Merge Sporocyst and Intramammalian
coldata$Timepoint[coldata$Timepoint=="1d_Sporocysts"] <- "Sporocyst"
coldata$Timepoint[coldata$Timepoint=="5d_Sporocysts"] <- "Sporocyst"
coldata$Timepoint[coldata$Timepoint=="32d_Sporocysts"] <- "Sporocyst"

coldata$Timepoint[coldata$Timepoint=="2d_Somules"] <- "Intramammalian"
coldata$Timepoint[coldata$Timepoint=="26d_Juveniles"] <- "Intramammalian"

##Fix column names
test<-colnames(bbcounts)
counter=1
for(item in test){
  test2<-(strsplit(item, "[_]"[1]))
  #test3<-(strsplit(test2, "[_]"[1]))
  colnames(bbcounts)[counter]<-test2[[1]][1]
  counter = counter + 1
}

##Fix rownames
ctsdata2 <- bbcounts[,-1]
rownames(ctsdata2) <- bbcounts[,1]
bbcounts<-ctsdata2
ctsdata2<- bbcounts[-grep('sma', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
ctsdata2<-bbcounts[-grep('3f41e371-e1ce-441d-ad6e-22d58b6d60fd', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
ctsdata2<-bbcounts[-grep('57585205-152a-490f-a316', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
###Save length
lens <- bbcounts$Length
bbcounts <- subset(bbcounts, select=-c(Chr,Start,End,Strand,Length))

##Fix rwnames of coldata
coldata2 <- data.frame(coldata[,-1])
rownames(coldata2) <- coldata[,1]
coldata<-coldata2
###might need to fix the col names
colnames(coldata)[1] <- "condition"
colnames(coldata)[2] <- "sex"


library(DESeq2)

dds<- DESeqDataSetFromMatrix(countData = bbcounts,
                             colData=coldata,
                             design = ~ condition,
)


dds$condition <- relevel(dds$condition,ref="Sporocyst")


keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

dds <- DESeq(dds)
res <- results(dds,contrast=c("condition","Intramammalian","Sporocyst"))

plotMA(res, ylim=c(-5,5))

##implements the method of Independent Hypothesis Weighting (Ignatiadis et al. 2016)
##Removes low copy number signifcant values
resLFC <- lfcShrink(dds, coef="condition_Intramammalian_vs_Sporocyst", type="apeglm")
returned_data <- plotMA(resLFC, ylim=c(-5,5), returnData=TRUE)

colnames(returned_data)[3] <- "Significance"

ggplot(returned_data, aes(x=mean, y=lfc,col=Significance)) + geom_point(cex=1) + 
  scale_x_continuous(trans="log10") +
  ylab("Fold Change (Log2)") + xlab("Mean of Normalised Counts") + 
  theme(axis.text = element_text(size=14), axis.title = element_text(size=16, face="bold"), legend.title = element_text(size=16), 
        legend.text = element_text(size=14), strip.text.x = element_text(size = 16)) + 
  geom_hline(yintercept=0,col="grey",lty="dashed")


returned_data_2 <- cbind(returned_data,"Genes")
colnames(returned_data_2)[4] <- "Is_Sub"

returned_data_2$Is_Sub[returned_data_2$Is_Sub=="Genes"] <- "All Genes"
returned_data_2$Is_Sub[returned_data_2$Is_Sub=="Hypothetical Protein"] <- "All Genes"
returned_data_2$Is_Sub[returned_data_2$Is_Sub=="Annexin"] <- "All Genes"
returned_data_2$Is_Sub[returned_data_2$Is_Sub=="Endoglycoceramidase"] <- "All Genes"


#EG_comp<-returned_data_2
#annex_comp<-returned_data_2
#hypo_comp<-returned_data_2

#setwd("/Users/tobybrann/Documents/Subtelomeres/manuscript/Figure_3/expression/sarah_BB/MA_plots")

#tiff('hypo_mammal_vs_sporo.png', units="in", width=14, height=6.5, res=1200, compression = 'lzw')


#dev.off()
setwd("/Users/tobybrann/Documents/Subtelomeres/thesis/MA")

#annex_comp <- returned_data_2
#hypo_comp <- returned_data_2
#EG_comp <- returned_data_2

write.csv(annex_comp,"/Users/tobybrann/Documents/Subtelomeres/thesis/MA/annex.tbl", row.names = FALSE)
write.csv(hypo_comp,"/Users/tobybrann/Documents/Subtelomeres/thesis/MA/hypo.tbl", row.names = FALSE)
write.csv(EG_comp,"/Users/tobybrann/Documents/Subtelomeres/thesis/MA/EG.tbl", row.names = FALSE)

setwd("/Users/tobybrann/Documents/Subtelomeres/thesis/MA/")
annex_comp<-read.delim(file="annex.tbl",sep=",",header=TRUE)
hypo_comp<-read.delim(file="hypo.tbl",sep=",",header=TRUE)
EG_comp<-read.delim(file="EG.tbl",sep=",",header=TRUE)



###Hypothetical
setwd("/Users/tobybrann/Documents/Subtelomeres/RNASeq/tags")
hypo<-read.delim(file="egg.txt", header=FALSE)   

hypo <- which(row.names(hypo_comp) %in% hypo$V1)
for(i in 1:length(hypo)){
  num<-(hypo[i])
  #writeLines(num)
  hypo_comp$Is_Sub[num] <- "Hypothetical Protein"
}

annex<-read.delim(file="annexin.txt", header=FALSE)   

annex <- which(row.names(annex_comp) %in% annex$V1)
for(i in 1:length(annex)){
  num<-(annex[i])
  #writeLines(num)
  annex_comp$Is_Sub[num] <- "Annexin"
}

EG<-read.delim(file="EG.txt", header=FALSE)   

EG <- which(row.names(EG_comp) %in% EG$V1)
for(i in 1:length(EG)){
  num<-(EG[i])
  #writeLines(num)
  EG_comp$Is_Sub[num] <- "Endoglycoceramidase"
}


library(MetBrewer)

all<-annex_comp%>%filter(Is_Sub=="All Genes") 
supp<-annex_comp%>%filter(Is_Sub!="All Genes") 

library(cowplot)
ggplot(all, aes(x=mean, y=lfc,col=Significance)) + geom_point(cex=1.25) + 
  scale_x_continuous(trans="log10") +
  ylab("Fold Change (Log2)") + xlab("Mean of Normalised Counts") + 
  theme(axis.text = element_text(size=20), axis.title = element_text(size=24, face="bold"), legend.title = element_text(size=14), 
        legend.text = element_text(size=14), strip.text.x = element_text(size = 20), plot.title=element_text(size=24)) + 
  geom_hline(yintercept=0,col="black",lty="dashed",size=1.25) + 
  ggtitle("Intramammalian vs Sporocyst") + 
  scale_color_manual(values=c(met.brewer("Johnson",2))) + 
  theme_minimal_grid() +
  theme(axis.line.x = element_line(color="black", size = 0.5),
        axis.line.y = element_line(color="black", size = 0.5)) + 
  geom_point(data=supp,aes(x=mean,y=lfc),col="black",size=3)
  facet_wrap(.~factor(Is_Sub))


tiff('annex_mammal_vs_sporo.png', units="in", width=14, height=6.5, res=400, compression = 'lzw')

ggplot(annex_comp, aes(x=mean, y=lfc,col=Significance)) + geom_point(cex=1.25) + 
  scale_x_continuous(trans="log10") +
  ylab("Fold Change (Log2)") + xlab("Mean of Normalised Counts") + 
  theme(axis.text = element_text(size=20), axis.title = element_text(size=24, face="bold"), legend.title = element_text(size=14), 
        legend.text = element_text(size=14), strip.text.x = element_text(size = 20), plot.title=element_text(size=24)) + 
  geom_hline(yintercept=0,col="black",lty="dashed",size=1.25) + 
  ggtitle("Intramammalian vs Sporocyst") + 
  scale_color_manual(values=c(met.brewer("Johnson",2))) + 
  facet_wrap(.~factor(Is_Sub))

dev.off()
tiff('hypo_mammal_vs_sporo.png', units="in", width=14, height=6.5, res=400, compression = 'lzw')

ggplot(hypo_comp, aes(x=mean, y=lfc,col=Significance)) + geom_point(cex=1.25) + 
  scale_x_continuous(trans="log10") +
  ylab("Fold Change (Log2)") + xlab("Mean of Normalised Counts") + 
  theme(axis.text = element_text(size=20), axis.title = element_text(size=24, face="bold"), legend.title = element_text(size=14), 
        legend.text = element_text(size=14), strip.text.x = element_text(size = 20), plot.title=element_text(size=24)) + 
  geom_hline(yintercept=0,col="black",lty="dashed",size=1.25) + 
  ggtitle("Intramammalian vs Sporocyst") + 
  scale_color_manual(values=c(met.brewer("Johnson",2))) + 
  facet_wrap(.~factor(Is_Sub))

dev.off()
tiff('eg_mammal_vs_sporo.png', units="in", width=14, height=6.5, res=400, compression = 'lzw')

ggplot(EG_comp, aes(x=mean, y=lfc,col=Significance)) + geom_point(cex=1.25) + 
  scale_x_continuous(trans="log10") +
  ylab("Fold Change (Log2)") + xlab("Mean of Normalised Counts") + 
  theme(axis.text = element_text(size=20), axis.title = element_text(size=24, face="bold"), legend.title = element_text(size=14), 
        legend.text = element_text(size=14), strip.text.x = element_text(size = 20), plot.title=element_text(size=24)) + 
  geom_hline(yintercept=0,col="black",lty="dashed",size=1.25) + 
  ggtitle("Intramammalian vs Eggs") + 
  scale_color_manual(values=c(met.brewer("Johnson",2))) + 
  facet_wrap(.~factor(Is_Sub))
dev.off()




##Volcano plots
if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)

dds <- DESeqDataSet(airway, design = ~ cell + dex)
dds <- DESeq(dds, betaPrior=FALSE)
res <- results(dds,
               contrast = c('dex','trt','untrt'))
res <- lfcShrink(dds,
                 contrast = c('dex','trt','untrt'), res=res, type = 'normal')

####Eggs
setwd("/Users/tobybrann/Documents/Subtelomeres/manuscript/Figure_3/expression/sarah_BB/MA_plots/")

bbcounts <- read.delim(file="BB_counts.txt",header=TRUE,sep="\t")
coldata <- read.delim(file="coldata.txt",header=TRUE,sep="\t")

###Merge Sporocyst and Intramammalian
coldata$Timepoint[coldata$Timepoint=="1d_Sporocysts"] <- "Sporocyst"
coldata$Timepoint[coldata$Timepoint=="5d_Sporocysts"] <- "Sporocyst"
coldata$Timepoint[coldata$Timepoint=="32d_Sporocysts"] <- "Sporocyst"

coldata$Timepoint[coldata$Timepoint=="2d_Somules"] <- "Intramammalian"
coldata$Timepoint[coldata$Timepoint=="26d_Juveniles"] <- "Intramammalian"

##Fix column names
test<-colnames(bbcounts)
counter=1
for(item in test){
  test2<-(strsplit(item, "[_]"[1]))
  #test3<-(strsplit(test2, "[_]"[1]))
  colnames(bbcounts)[counter]<-test2[[1]][1]
  counter = counter + 1
}

##Fix rownames
ctsdata2 <- bbcounts[,-1]
rownames(ctsdata2) <- bbcounts[,1]
bbcounts<-ctsdata2
ctsdata2<- bbcounts[-grep('sma', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
ctsdata2<-bbcounts[-grep('3f41e371-e1ce-441d-ad6e-22d58b6d60fd', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
ctsdata2<-bbcounts[-grep('57585205-152a-490f-a316', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
###Save length
lens <- bbcounts$Length
bbcounts <- subset(bbcounts, select=-c(Chr,Start,End,Strand,Length))

##Fix rwnames of coldata
coldata2 <- data.frame(coldata[,-1])
rownames(coldata2) <- coldata[,1]
coldata<-coldata2
###might need to fix the col names
colnames(coldata)[1] <- "condition"
colnames(coldata)[2] <- "sex"


library(DESeq2)

dds<- DESeqDataSetFromMatrix(countData = bbcounts,
                             colData=coldata,
                             design = ~ condition,
)
#dds$condition <- relevel(dds$condition,ref="Intramammalian")

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

dds <- DESeq(dds)
#res <- results(dds)
res <- results(dds,contrast=c("condition","Sporocyst","Intramammalian"))

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements

keyvals.shape <- ifelse(
  rownames(res) %in% egg_list, 16,17)

names(keyvals.shape)[keyvals.shape == 17] <- 'Other'
names(keyvals.shape)[keyvals.shape == 16] <- 'Egg'

keyvals.colour <- ifelse(abs(res$log2FoldChange) < 1.0 & res$pvalue > pcutoff, 'grey',
  ifelse(abs(res$log2FoldChange) > 1.0 & res$pvalue > pcutoff, 'green',
         ifelse(abs(res$log2FoldChange) < 1.0 & res$pvalue < pcutoff, 'blue',
         ifelse(abs(res$log2FoldChange) > 1.0 & res$pvalue < pcutoff, 'red',
         'black'))))

keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'green'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'low'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'low2'

matching_indices <- which(rownames(res) %in% egg_list)

for (index in matching_indices) {
  keyvals.colour[[index]] <- "black"
}

names(keyvals.colour)[keyvals.colour == 'black'] <- 'egg'

##Remove labels in bad spots
to_remove<-c("Smp_350320.1","Smp_165590.1","Smp_192050.1","Smp_349650.1","Smp_350230.1","Smp_165600.1","Smp_350150.1","Smp_349210.1","Smp_349200.1","Smp_192040.1","Smp_349250.1",
             "Smp_350220.1","Smp_349580.1","Smp_093070.1","Smp_349720.1","Smp_349080.1","Smp_102020.1","Smp_349240.1","Smp_349710.1","Smp_349170.1",
             "Smp_349180.1","Smp_086630.1","Smp_349150.1","Smp_349730.1","Smp_179390.1","Smp_350510.1","Smp_349740.1","Smp_154590.1","Smp_350160.1",
             "Smp_349700.1","Smp_349570.1","Smp_349520.1","Smp_350240.1","Smp_349530.1","Smp_349100.1","Smp_349090.1","Smp_349160.1","Smp_179970.1")
egg_list2 <- egg_list[!egg_list %in% to_remove]


setwd("/Users/tobybrann/Documents/Subtelomeres/thesis/volcano")
tiff('egg_intra_vs_sporo.tiff', units="in", width=10, height=6, res=400, compression = 'lzw')

EnhancedVolcano(res,
                lab = rownames(res),
                selectLab = egg_list2,
                x = 'log2FoldChange',
                y = 'pvalue',
                shapeCustom = keyvals.shape,
                pointSize = 1.5,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                pCutoff = pcutoff,
                FCcutoff = 1.0,
                legendPosition = 'none',
                #legendLabSize = 14,
                #legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE,
                colConnectors = 'black',
                colCustom = keyvals.colour)
                #colCustom = egg_list)

dev.off()


####Annexin
setwd("/Users/tobybrann/Documents/Subtelomeres/manuscript/Figure_3/expression/sarah_BB/MA_plots/")

bbcounts <- read.delim(file="BB_counts.txt",header=TRUE,sep="\t")
coldata <- read.delim(file="coldata.txt",header=TRUE,sep="\t")

###Merge Sporocyst and Intramammalian
coldata$Timepoint[coldata$Timepoint=="1d_Sporocysts"] <- "Sporocyst"
coldata$Timepoint[coldata$Timepoint=="5d_Sporocysts"] <- "Sporocyst"
coldata$Timepoint[coldata$Timepoint=="32d_Sporocysts"] <- "Sporocyst"

coldata$Timepoint[coldata$Timepoint=="2d_Somules"] <- "Intramammalian"
coldata$Timepoint[coldata$Timepoint=="26d_Juveniles"] <- "Intramammalian"

##Fix column names
test<-colnames(bbcounts)
counter=1
for(item in test){
  test2<-(strsplit(item, "[_]"[1]))
  #test3<-(strsplit(test2, "[_]"[1]))
  colnames(bbcounts)[counter]<-test2[[1]][1]
  counter = counter + 1
}

##Fix rownames
ctsdata2 <- bbcounts[,-1]
rownames(ctsdata2) <- bbcounts[,1]
bbcounts<-ctsdata2
ctsdata2<- bbcounts[-grep('sma', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
ctsdata2<-bbcounts[-grep('3f41e371-e1ce-441d-ad6e-22d58b6d60fd', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
ctsdata2<-bbcounts[-grep('57585205-152a-490f-a316', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
###Save length
lens <- bbcounts$Length
bbcounts <- subset(bbcounts, select=-c(Chr,Start,End,Strand,Length))

##Fix rwnames of coldata
coldata2 <- data.frame(coldata[,-1])
rownames(coldata2) <- coldata[,1]
coldata<-coldata2
###might need to fix the col names
colnames(coldata)[1] <- "condition"
colnames(coldata)[2] <- "sex"


library(DESeq2)

dds<- DESeqDataSetFromMatrix(countData = bbcounts,
                             colData=coldata,
                             design = ~ condition,
)
#dds$condition <- relevel(dds$condition,ref="Intramammalian")

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

dds <- DESeq(dds)
#res <- results(dds)
res <- results(dds,contrast=c("condition","Sporocyst","Intramammalian"))

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements

annex_list<-annex[,1]

keyvals.shape <- ifelse(
  rownames(res) %in% annex_list, 16,17)

names(keyvals.shape)[keyvals.shape == 17] <- 'Other'
names(keyvals.shape)[keyvals.shape == 16] <- 'Annexin'

keyvals.colour <- ifelse(abs(res$log2FoldChange) < 1.0 & res$pvalue > pcutoff, 'grey',
                         ifelse(abs(res$log2FoldChange) > 1.0 & res$pvalue > pcutoff, 'green',
                                ifelse(abs(res$log2FoldChange) < 1.0 & res$pvalue < pcutoff, 'blue',
                                       ifelse(abs(res$log2FoldChange) > 1.0 & res$pvalue < pcutoff, 'red',
                                              'black'))))

keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'green'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'low'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'low2'

matching_indices <- which(rownames(res) %in% annex_list)

for (index in matching_indices) {
  keyvals.colour[[index]] <- "black"
}

names(keyvals.colour)[keyvals.colour == 'black'] <- 'annexin'

##Remove labels in bad spots
to_remove<-c("Smp_337100.2","Smp_334820.1","Smp_324450.1","Smp_329050.1","Smp_341290.1","Smp_334670.1","Smp_328930.1","Smp_311250.1")
annex_list2 <- annex_list[!annex_list %in% to_remove]


setwd("/Users/tobybrann/Documents/Subtelomeres/thesis/volcano")
tiff('annexin_intra_vs_sporo.tiff', units="in", width=10, height=6, res=400, compression = 'lzw')

EnhancedVolcano(res,
                lab = rownames(res),
                selectLab = annex_list2,
                x = 'log2FoldChange',
                y = 'pvalue',
                shapeCustom = keyvals.shape,
                pointSize = 1.5,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                boxedLabels = TRUE,
                colAlpha = 4/5,
                pCutoff = pcutoff,
                FCcutoff = 1.0,
                legendPosition = 'none',
                #legendLabSize = 14,
                #legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE,
                colConnectors = 'black',
                colCustom = keyvals.colour)
#colCustom = egg_list)

dev.off()


####Endoglycoceramidases
setwd("/Users/tobybrann/Documents/Subtelomeres/manuscript/Figure_3/expression/sarah_BB/MA_plots/")

bbcounts <- read.delim(file="BB_counts.txt",header=TRUE,sep="\t")
coldata <- read.delim(file="coldata.txt",header=TRUE,sep="\t")

###Merge Sporocyst and Intramammalian
coldata$Timepoint[coldata$Timepoint=="1d_Sporocysts"] <- "Sporocyst"
coldata$Timepoint[coldata$Timepoint=="5d_Sporocysts"] <- "Sporocyst"
coldata$Timepoint[coldata$Timepoint=="32d_Sporocysts"] <- "Sporocyst"

coldata$Timepoint[coldata$Timepoint=="2d_Somules"] <- "Intramammalian"
coldata$Timepoint[coldata$Timepoint=="26d_Juveniles"] <- "Intramammalian"

##Fix column names
test<-colnames(bbcounts)
counter=1
for(item in test){
  test2<-(strsplit(item, "[_]"[1]))
  #test3<-(strsplit(test2, "[_]"[1]))
  colnames(bbcounts)[counter]<-test2[[1]][1]
  counter = counter + 1
}

##Fix rownames
ctsdata2 <- bbcounts[,-1]
rownames(ctsdata2) <- bbcounts[,1]
bbcounts<-ctsdata2
ctsdata2<- bbcounts[-grep('sma', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
ctsdata2<-bbcounts[-grep('3f41e371-e1ce-441d-ad6e-22d58b6d60fd', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
ctsdata2<-bbcounts[-grep('57585205-152a-490f-a316', row.names(bbcounts)),, drop = FALSE]
bbcounts<-ctsdata2
###Save length
lens <- bbcounts$Length
bbcounts <- subset(bbcounts, select=-c(Chr,Start,End,Strand,Length))

##Fix rwnames of coldata
coldata2 <- data.frame(coldata[,-1])
rownames(coldata2) <- coldata[,1]
coldata<-coldata2
###might need to fix the col names
colnames(coldata)[1] <- "condition"
colnames(coldata)[2] <- "sex"


library(DESeq2)

dds<- DESeqDataSetFromMatrix(countData = bbcounts,
                             colData=coldata,
                             design = ~ condition,
)
#dds$condition <- relevel(dds$condition,ref="Intramammalian")

keep <- rowSums(counts(dds)) >=10
dds <- dds[keep,]

dds <- DESeq(dds)
#res <- results(dds)
res <- results(dds,contrast=c("condition","Eggs","Intramammalian"))

# create custom key-value pairs for different cell-types
# this can be achieved with nested ifelse statements

EG_list<-EG[,1]

keyvals.shape <- ifelse(
  rownames(res) %in% EG_list, 16,17)

names(keyvals.shape)[keyvals.shape == 17] <- 'Other'
names(keyvals.shape)[keyvals.shape == 16] <- 'Endoglycoceramidase'

keyvals.colour <- ifelse(abs(res$log2FoldChange) < 1.0 & res$pvalue > pcutoff, 'grey',
                         ifelse(abs(res$log2FoldChange) > 1.0 & res$pvalue > pcutoff, 'green',
                                ifelse(abs(res$log2FoldChange) < 1.0 & res$pvalue < pcutoff, 'blue',
                                       ifelse(abs(res$log2FoldChange) > 1.0 & res$pvalue < pcutoff, 'red',
                                              'black'))))

keyvals.colour[is.na(keyvals.colour)] <- 'black'
names(keyvals.colour)[keyvals.colour == 'grey'] <- 'high'
names(keyvals.colour)[keyvals.colour == 'green'] <- 'mid'
names(keyvals.colour)[keyvals.colour == 'blue'] <- 'low'
names(keyvals.colour)[keyvals.colour == 'red'] <- 'low2'

matching_indices <- which(rownames(res) %in% EG_list)

for (index in matching_indices) {
  keyvals.colour[[index]] <- "black"
}

names(keyvals.colour)[keyvals.colour == 'black'] <- 'Endoglycoceramidase'

##Remove labels in bad spots
to_remove<-c("Smp_344740.1","Smp_317230.1","Smp_317470.2","Smp_329030.1","Smp_317220.1")
#to_remove<-c("")
EG_list2 <- EG_list[!EG_list %in% to_remove]


setwd("/Users/tobybrann/Documents/Subtelomeres/thesis/volcano")
tiff('EG_intra_vs_sporo.tiff', units="in", width=10, height=6, res=400, compression = 'lzw')

pcutoff=1e-5

EnhancedVolcano(res,
                lab = rownames(res),
                selectLab = EG_list2,
                x = 'log2FoldChange',
                y = 'pvalue',
                shapeCustom = keyvals.shape,
                pointSize = 1.5,
                labSize = 3.0,
                labCol = 'black',
                labFace = 'bold',
                maxoverlapsConnectors = Inf,
                boxedLabels = TRUE,
                colAlpha = 4/5,
                pCutoff = pcutoff,
                FCcutoff = 1.0,
                legendPosition = 'none',
                #legendLabSize = 14,
                #legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.75,
                arrowheads = FALSE,
                colConnectors = 'black',
                colCustom = keyvals.colour)
#colCustom = egg_list)

dev.off()







###Schisto biser hits
####Figure 2, hits from biser
setwd("/Users/tobybrann/Documents/Subtelomeres/thesis/other_schisto/")

library(circlize)

all_pairs <- read.delim("jap_pairs.tbl",sep="\t",header=FALSE)

#Cytoband v2
cytoband.file = system.file(package = "circlize", "extdata", "cytoBand.txt")
cytoband.df = read.table("/Users/tobybrann/Documents/Subtelomeres/thesis/other_schisto/orig_jap.genome", colClasses = c("character", "numeric","numeric", "character", "character"), sep = "\t")

LHS_all <- all_pairs[,c(1:3,7)]
RHS_all <- all_pairs[,c(4:7)]

colnames(LHS_all)[4]<-"V4"
colnames(RHS_all)[1]<-"V1"
colnames(RHS_all)[2]<-"V2"
colnames(RHS_all)[3]<-"V3"
colnames(RHS_all)[4]<-"V4"

##Subset
LHS_sub <- LHS_all %>% filter(V4=="Subtelomere")
RHS_sub <- RHS_all %>% filter(V4=="Subtelomere")

LHS_no <- LHS_all %>% filter(V4=="No")
RHS_no <- RHS_all %>% filter(V4=="No")

colgen <- "#72bcd5"

#Function for genic colours 
colgen <- LHS_all$V4
Eles<-unique(LHS_all$V4)
X=length(Eles)
col = rep(c(values = met.brewer("Hiroshige",X)))
fors=length(LHS_all$V4)
for (i in 1:X) {
  #colgen <- gsub(Eles[val], col[val],colgen)
  colgen <- gsub(Eles[i], adjustcolor(col[i], alpha.f=0.8),colgen)
  
}

###Filter chromosomes
temp<-unique(c(LHS_all$V1,RHS_all$V1))

df<-cytoband.df %>% filter(V1 %in% temp)
segmented_df <- data.frame()

###Add tags to cytoband
# Assuming df is your data frame, with V3 as the total length column

# Set the 3M limit
segment_limit <- 3000000

# Create an empty data frame to store the segmented data
segmented_df <- data.frame()

# Loop through each row and create three segments if total length > 6M
for (i in 1:nrow(df)) {
  # Extract row data
  V1 <- df$V1[i]
  V2 <- df$V2[i]
  total_length <- df$V3[i]
  V4 <- df$V4[i]
  V5 <- df$V5[i]
  
  # Check if total length is greater than 6,000,000
  if (total_length > 6000000) {
    # Define the segments
    first_segment <- data.frame(V1 = V1, V2 = 0, V3 = segment_limit, V4 = V4, V5 = V5)
    middle_segment <- data.frame(V1 = V1, V2 = segment_limit, V3 = total_length - segment_limit, V4 = V4, V5 = V5)
    last_segment <- data.frame(V1 = V1, V2 = total_length - segment_limit, V3 = total_length, V4 = V4, V5 = V5)
    
    # Add the three segments to the segmented data frame
    segmented_df <- rbind(segmented_df, first_segment, middle_segment, last_segment)
  } else {
    # If total length <= 6,000,000, keep the original row
    segmented_df <- rbind(segmented_df, df[i, ])
  }
}

for (i in 1:nrow(segmented_df)) {
  if(segmented_df$V3[i]-segmented_df$V2[i]<6000000){
    segmented_df$V5[i] <- "gpos25"
  }
}

cytoband.df<-segmented_df



tiff('jap_hits_body.tiff', units="in", width=5, height=4, res=400, compression = 'lzw')

circos.clear()
circos.par("track.height" = 0.25, start.degree=90,gap.degree=c(1))
circos.initializeWithIdeogram(cytoband.df, ideogram.height = 0.04, plotType=c("labels", "ideogram"))
circos.genomicLink(LHS_sub, RHS_sub, col= "#EF8A47CC"
)
title("S. japonicum")
dev.off()



##Haem
all_pairs <- read.delim("haem_pairs.tbl",sep="\t",header=FALSE)


#Cytoband v2
cytoband.file = system.file(package = "circlize", "extdata", "cytoBand.txt")
cytoband.df = read.table("/Users/tobybrann/Documents/Subtelomeres/thesis/other_schisto/orig_haem.genome", colClasses = c("character", "numeric","numeric", "character", "character"), sep = "\t")

LHS_all <- all_pairs[,c(1:3,7)]
RHS_all <- all_pairs[,c(4:7)]

colnames(LHS_all)[4]<-"V4"
colnames(RHS_all)[1]<-"V1"
colnames(RHS_all)[2]<-"V2"
colnames(RHS_all)[3]<-"V3"
colnames(RHS_all)[4]<-"V4"

##Subset
LHS_sub <- LHS_all %>% filter(V4=="Subtelomere")
RHS_sub <- RHS_all %>% filter(V4=="Subtelomere")

LHS_no <- LHS_all %>% filter(V4=="No")
RHS_no <- RHS_all %>% filter(V4=="No")

colgen <- "#72bcd5"

#Function for genic colours 
colgen <- LHS_all$V4
Eles<-unique(LHS_all$V4)
X=length(Eles)
col = rep(c(values = met.brewer("Hiroshige",X)))
fors=length(LHS_all$V4)
for (i in 1:X) {
  #colgen <- gsub(Eles[val], col[val],colgen)
  colgen <- gsub(Eles[i], adjustcolor(col[i], alpha.f=0.8),colgen)
  
}

###Filter chromosomes
temp<-unique(c(LHS_all$V1,RHS_all$V1))

df<-cytoband.df %>% filter(V1 %in% temp)
segmented_df <- data.frame()

###Add tags to cytoband
# Assuming df is your data frame, with V3 as the total length column

# Set the 3M limit
segment_limit <- 3000000

# Create an empty data frame to store the segmented data
segmented_df <- data.frame()

# Loop through each row and create three segments if total length > 6M
for (i in 1:nrow(df)) {
  # Extract row data
  V1 <- df$V1[i]
  V2 <- df$V2[i]
  total_length <- df$V3[i]
  V4 <- df$V4[i]
  V5 <- df$V5[i]
  
  # Check if total length is greater than 6,000,000
  if (total_length > 6000000) {
    # Define the segments
    first_segment <- data.frame(V1 = V1, V2 = 0, V3 = segment_limit, V4 = V4, V5 = V5)
    middle_segment <- data.frame(V1 = V1, V2 = segment_limit, V3 = total_length - segment_limit, V4 = V4, V5 = V5)
    last_segment <- data.frame(V1 = V1, V2 = total_length - segment_limit, V3 = total_length, V4 = V4, V5 = V5)
    
    # Add the three segments to the segmented data frame
    segmented_df <- rbind(segmented_df, first_segment, middle_segment, last_segment)
  } else {
    # If total length <= 6,000,000, keep the original row
    segmented_df <- rbind(segmented_df, df[i, ])
  }
}

for (i in 1:nrow(segmented_df)) {
  if(segmented_df$V3[i]-segmented_df$V2[i]<6000000){
    segmented_df$V5[i] <- "gpos25"
  }
}

cytoband.df<-segmented_df



tiff('haem_hits_body.tiff', units="in", width=5, height=4, res=400, compression = 'lzw')

circos.clear()
circos.par("track.height" = 0.25, start.degree=90,gap.degree=c(1))
circos.initializeWithIdeogram(cytoband.df, ideogram.height = 0.04, plotType=c("labels", "ideogram"))
circos.genomicLink(LHS_sub, RHS_sub, col= "#EF8A47CC"
)
title("S. haematobium")
dev.off()



####Enrichment analysis
setwd("/Users/tobybrann/Documents/Subtelomeres/thesis/enrich/")
non_sub<-read.delim(file="non_subtelomeric_TEs.gff",sep="\t",header=FALSE)
sub<-read.delim(file="subtelomeric_TEs.gff",sep="\t",header=FALSE)

non_sub$V9 <- sapply(strsplit(as.character(non_sub$V9), " "), `[`, 2)
sub$V9 <- sapply(strsplit(as.character(sub$V9), " "), `[`, 2)

sums_non_sub <- non_sub %>%
  mutate(Difference = V5 - V4) %>%  # Calculate the difference for each row
  group_by(V9) %>%                 # Group by unique entries in V9
  summarise(SumDifference = sum(Difference, na.rm = TRUE))  # Sum the differences for each group

sums_sub <- sub %>%
  mutate(Difference = V5 - V4) %>%  # Calculate the difference for each row
  group_by(V9) %>%                 # Group by unique entries in V9
  summarise(SumDifference = sum(Difference, na.rm = TRUE))  # Sum the differences for each group

sums_non_sub<-cbind(sums_non_sub,"Body")
sums_sub<-cbind(sums_sub,"Subtelomere")

colnames(sums_sub)[3]<- "V3"
colnames(sums_non_sub)[3]<- "V3"

combi_TEs<-rbind(sums_non_sub,sums_sub)




data_bar <- data.frame(
  stringsAsFactors = F,
  Sample = rep(c("A", "B"), each = 10),
  Percentage = runif(20),
  Taxon = rep(1:10, by = 2)
)

comb_TEs_perc<-combi_TEs %>%
  group_by(V3) %>%
  mutate(perc = SumDifference* 100/sum(SumDifference)) 

library(ggalluvial)

comb_TEs_perc$V9 <- substr(comb_TEs_perc$V9, 7, nchar(comb_TEs_perc$V9))

TE_families <- read.delim(file="/Users/tobybrann/Documents/Subtelomeres/thesis/enrich/TE_families.txt",sep="\t",header=FALSE)

 

colnames(comb_TEs_perc)[1]<-"V1"
comb_TEs_perc_wfam<-comb_TEs_perc %>% inner_join(TE_families, by="V1")

comb_TEs_perc_filt <- comb_TEs_perc %>%
  group_by(V9) %>%                                      # Group by unique V9 values
  filter(all(SumDifference >= 100000)) %>%                           # Keep groups where all V6 >= 1000
  ungroup()   

ggplot(comb_TEs_perc_wfam, aes(y = perc, x = V3, fill = V1)) +
  geom_flow(aes(alluvium = V1), alpha= .5, color = "white",
            curve_type = "linear", 
            width = .5,) +
  geom_col(width = .5, color = "white") +
  #scale_fill_brewer(palette = "RdBu")+
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  facet_wrap(~ V2, scales = "fixed")

comb_TEs_perc_wfam_filt2<-comb_TEs_perc_wfam %>% filter(V2!="DNA_CACTA") %>%
  filter(V2!="DNA_Merlin") %>%
  filter(V2!="DNA_Mut") %>%
  filter(V2!="LTR_Unknown") %>%
  filter(V2!="Unknown") %>%
  filter(V2!="SINE") 

ggplot(comb_TEs_perc_wfam_filt2, aes(y = perc, x = V3, fill = V1)) +
  geom_flow(aes(alluvium = V1), alpha= .5, color = "white",
            curve_type = "linear", 
            width = .5,) +
  geom_col(width = .5, color = "white") +
  #scale_fill_brewer(palette = "RdBu")+
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) +
  facet_wrap(~ V2, scales = "fixed")

comb_TEs_perc_wfam_filt2$V2[comb_TEs_perc_wfam_filt2$V1 == "Perere-3"] <- "Perere-3"
comb_TEs_perc_wfam_filt2$V2[comb_TEs_perc_wfam_filt2$V1 == "Sr3"] <- "Sr3"

comb_TEs_perc_wfam_filt2_group <- comb_TEs_perc_wfam_filt2 %>%
  group_by(V2, V3) %>%                          # Group by V2 and V3
  summarise(TotalSum = sum(perc, na.rm = TRUE), .groups = "drop")  # Sum the SumDifference column



comb_TEs_perc_wfam_filt2_group$V2 <-factor(comb_TEs_perc_wfam_filt2_group$V2, levels =c("LINE_Jockey","LINE_R2","LTR_Gypsy","LTR_Pao","PLE","LINE_RTE","Perere-3","Sr3")) 

Pr3<-"#156082"
Sr3<-"#EB7131"

RTE<-met.brewer("Tam",8)[1]
col<-rev(met.brewer("Monet",5))
col_array<-c(col,RTE,Sr3,Pr3)

##Appendix

#"LINE_Jockey" "LINE_RTE"    "LTR_Gypsy"   "PLE"        
#"LTR_Pao"     "Perere-3"    "LINE_R2"     "Sr3"


##Subsets
subset<-comb_TEs_perc_wfam_filt2 %>% filter(V2=="LINE_R2")

col<-rev(met.brewer("Monet",length(unique(subset$V1))))

tiff('R2.tiff', units="in", width=8, height=6, res=600, compression = 'lzw')

ggplot(subset, aes(y = perc, x = V3, fill = V1)) +
  geom_flow(aes(alluvium = V1), alpha= .5, color = "white",
            curve_type = "linear", 
            width = .5,) +
  geom_col(width = .5, color = "white") +
  scale_fill_manual(values = col)+
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) 

dev.off()


setwd("/Users/tobybrann/Documents/Subtelomeres/thesis/enrich")

tiff('all_fam.tiff', units="in", width=8, height=6, res=600, compression = 'lzw')

ggplot(comb_TEs_perc_wfam_filt2_group, aes(y = TotalSum, x = V3, fill = V2)) +
  geom_flow(aes(alluvium = V2), alpha= .5, color = "white",
            curve_type = "linear", 
            width = .5,) +
  geom_col(width = .5, color = "white") +
  scale_fill_manual(values = col_array)+
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) 

dev.off()

Jockey_cols <- met.brewer("Redon",12)[1:10]

tiff('jockey.tiff', units="in", width=8, height=6, res=600, compression = 'lzw')

comb_TEs_perc_wfam_filt2 %>% filter(V2=="LINE_Jockey") %>%
ggplot(aes(y = perc, x = V3, fill = V1)) +
  geom_flow(aes(alluvium = V1), alpha= .5, color = "white",
            curve_type = "linear", 
            width = .5,) +
  geom_col(width = .5, color = "white") +
  scale_fill_manual(values = Jockey_cols)+
  scale_y_continuous(NULL, expand = c(0,0)) +
  cowplot::theme_minimal_hgrid() +
  theme(panel.grid.major = element_blank(), 
        axis.text.y = element_blank(), 
        axis.ticks.y = element_blank()) 

dev.off()

Jockeys_to_calc<-comb_TEs_perc_wfam_filt2 %>% filter(V2=="LINE_Jockey") 


###Stats
TE_subfam<-read.delim(file="/Users/tobybrann/Documents/Subtelomeres/thesis/enrich/TEs_for_calc.txt",sep="\t",header=TRUE)

t_test_result<-t.test(TE_subfam$Body,TE_subfam$Subtelomere,paired = FALSE, var.equal = TRUE)
print(t_test_result)

comb_TEs_perc_wfam_filt2_group %>%
  group_by(V2, V3) %>%                          # Group by V2 and V3
  summarise(TotalSum = sum(perc, na.rm = TRUE), .groups = "drop")  # Sum the SumDifference column



