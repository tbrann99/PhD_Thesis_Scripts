library(see)
library(ggplot2)
library(ggsignif)
library(MetBrewer)
library(dunn.test)
library(ggpubr)
library(circlize)
library(dplyr)


###Circos
setwd("/Users/tobybrann/Documents/FYR/circos/")
GC<-read.delim(file="GC.tbl",sep="\t",header=TRUE)
TE<-read.delim(file="TEs.tbl",sep="\t",header=FALSE)
WE<-read.delim(file="WEs.tbl",sep="\t",header=FALSE)

TE$x <- as.integer((TE$V2 + TE$V3))/2
names(TE)[names(TE) == "V1"] <- "sectors"
names(TE)[names(TE) == "V4"] <- "y"
TE$sectors <- gsub("SM_V10_", "",TE$sectors)

GC$x <- as.integer((GC$X2_usercol + GC$X3_usercol))/2
names(GC)[names(GC) == "X.1_usercol"] <- "sectors"
names(GC)[names(GC) == "X4_pct_at"] <- "y"
GC$sectors <- gsub("SM_V10_", "",GC$sectors)
GC$y <- 1-GC$y

WE$x <- as.integer((WE$V2 + WE$V3))/2
names(WE)[names(WE) == "V1"] <- "sectors"
names(WE)[names(WE) == "V4"] <- "y"
WE$sectors <- gsub("SM_V10_", "",WE$sectors)

TE_mean<-mean(TE$y)
GC_mean<-mean(GC$y)
WE_mean<-mean(WE$y)

circos.initialize(TE$sectors, x=TE$x)

circos.track(TE$sectors, x = TE$x, y = TE$y,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(7), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
               cell.xlim = get.cell.meta.data("cell.xlim")
               circos.lines(cell.xlim,  c(TE_mean, TE_mean),  col="black", lwd=1.5, lty="dashed")
             })
circos.trackLines(TE$sectors,TE$x,TE$y, col=col, cex=10, lwd=2)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = '1')

circos.track(GC$sectors, x = GC$x, y = GC$y,
             panel.fun = function(x, y) {
               cell.xlim = get.cell.meta.data("cell.xlim")
               circos.lines(cell.xlim,  c(GC_mean, GC_mean),  col="black", lwd=1.5, lty="dashed")
             })
circos.trackPoints(GC$sectors,GC$x,GC$y, col=col, cex=0.15)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = '1')

circos.track(WE$sectors, x = WE$x, y = WE$y,
             panel.fun = function(x, y) {
               cell.xlim = get.cell.meta.data("cell.xlim")
               circos.lines(cell.xlim,  c(WE_mean, WE_mean),  col="black", lwd=1.5, lty="dashed")
             })
circos.trackLines(WE$sectors,WE$x,WE$y, col=col,lwd=2)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = '1')
circos.clear


TE_auto<-TE %>% filter(sectors!="Z" , sectors!="WSR")
GC_auto<-GC %>% filter(sectors!="Z" , sectors!="WSR")
WE_auto<-WE %>% filter(sectors!="Z" , sectors!="WSR")

#col = rep(c('#8b4513', '#006400','#4682b4','#4b0082','#ff0000','#ffd700','#7cfc00','#00ffff','#0000ff','#ffe4b5','#ff69b4'), 1)
col<-met.brewer("Austria",7)

tiff('circos_autosomes.tiff', units="in", width=8, height=8, res=600, compression = 'lzw')

circos.clear

circos.par("track.height" = 0.22, start.degree=90, gap.degree=c(1,1,1,1,1,1,12))
circos.initialize(TE_auto$sectors, x=TE_auto$x)
circos.track(TE_auto$sectors, x = TE_auto$x, y = TE_auto$y,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(7), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
               cell.xlim = get.cell.meta.data("cell.xlim")
               circos.lines(cell.xlim,  c(TE_mean, TE_mean),  col="black", lwd=1.5, lty="dashed")
             })
circos.trackLines(TE_auto$sectors,TE_auto$x,TE_auto$y, col=col, cex=10, lwd=2)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = '1')

circos.track(GC_auto$sectors, x = GC_auto$x, y = GC_auto$y,
             panel.fun = function(x, y) {
               cell.xlim = get.cell.meta.data("cell.xlim")
               circos.lines(cell.xlim,  c(GC_mean, GC_mean),  col="black", lwd=1.5, lty="dashed")
             })
circos.trackPoints(GC_auto$sectors,GC_auto$x,GC_auto$y, col=col, cex=0.15)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = '1')

circos.track(WE_auto$sectors, x = WE_auto$x, y = WE_auto$y,
             panel.fun = function(x, y) {
               cell.xlim = get.cell.meta.data("cell.xlim")
               circos.lines(cell.xlim,  c(WE_mean, WE_mean),  col="black", lwd=1.5, lty="dashed")
             })
circos.trackLines(WE_auto$sectors,WE_auto$x,WE_auto$y, col=col,lwd=2)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = '1')
circos.clear

dev.off()

tiff('circos_sex.tiff', units="in", width=8, height=8, res=600, compression = 'lzw')

TE_sex<-TE %>% filter(sectors=="Z" | sectors=="WSR")
GC_sex<-GC %>% filter(sectors=="Z" | sectors=="WSR")
WE_sex<-WE %>% filter(sectors=="Z" | sectors=="WSR")

col<-met.brewer("Egypt",2)

circos.clear

circos.par("track.height" = 0.22, start.degree=90, gap.degree=c(3,12))
circos.initialize(TE_sex$sectors, x=TE_sex$x)
circos.track(TE_sex$sectors, x = TE_sex$x, y = TE_sex$y,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(7), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
               cell.xlim = get.cell.meta.data("cell.xlim")
               circos.lines(cell.xlim,  c(TE_mean, TE_mean),  col="black", lwd=1.5, lty="dashed")
             })
circos.trackLines(TE_sex$sectors,TE_sex$x,TE_sex$y, col=col, cex=10, lwd=2)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = 'WSR')

circos.track(GC_sex$sectors, x = GC_sex$x, y = GC_sex$y,
             panel.fun = function(x, y) {
               cell.xlim = get.cell.meta.data("cell.xlim")
               circos.lines(cell.xlim,  c(GC_mean, GC_mean),  col="black", lwd=1.5, lty="dashed")
             })
circos.trackPoints(GC_sex$sectors,GC_sex$x,GC_sex$y, col=col, cex=0.15)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = 'WSR')

circos.track(WE_sex$sectors, x = WE_sex$x, y = WE_sex$y,
             panel.fun = function(x, y) {
               cell.xlim = get.cell.meta.data("cell.xlim")
               circos.lines(cell.xlim,  c(WE_mean, WE_mean),  col="black", lwd=1.5, lty="dashed")
             })
circos.trackLines(WE_sex$sectors,WE_sex$x,WE_sex$y, col=col,lwd=2)
circos.yaxis(labels.cex=0.4, side = 'left', tick = T, sector.index = 'WSR')
circos.clear

dev.off()

##GC Cont
setwd("/Users/tobybrann/Documents/FYR/norm_plots/")

GC_mask<-read.delim(file="GCs_masked_norm.tbl",sep="\t",header=FALSE)
GC_unmask<-read.delim(file="GCs_unmasked_norm.tbl",sep="\t",header=FALSE)
TEs_norm<-read.delim(file="TEs_norm.tbl",sep="\t",header=FALSE)
WEs_norm<-read.delim(file="WEs_norm.tbl",sep="\t",header=FALSE)

TEs_norm<-cbind(TEs_norm,"TEs")
WEs_norm<-cbind(WEs_norm,"WEs")

colnames(TEs_norm)[4]<-"V4"
colnames(WEs_norm)[4]<-"V4"

REs_norm<-rbind(TEs_norm,WEs_norm)

  
GC_mask<-cbind(GC_mask,"GC_mask")
GC_unmask<-cbind(GC_unmask,"GC_unmask")

colnames(GC_mask)[4]<-"V4"
colnames(GC_unmask)[4]<-"V4"

GC_norm<-rbind(GC_mask,GC_unmask)

col<-met.brewer("Austria",7)
ggplot(GC_norm, aes(x=V2,y=V3,fill=V1)) + 
  geom_smooth()+
  facet_grid(.~V4) + ylim(0.3,0.5) 
  

GC_norm_2<-GC_norm %>% filter(V1!="SM_V10_WSR") %>%filter(V1!="SM_V10_Z")

ggplot(GC_norm_2, aes(x=V2,y=V3,col=V1)) + 
  geom_hline(yintercept=0.3544052,linetype="dashed",col="black") +
  geom_smooth(se=FALSE)+
  facet_grid(.~V4) + 
  ylab("GC Content") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=c(met.brewer("Austria",7)))

REs_norm_2<-REs_norm %>% filter(V1!="SM_V10_WSR") %>%filter(V1!="SM_V10_Z")

ggplot(REs_norm_2, aes(x=V2,y=V3,col=V1)) + 
  geom_smooth(se=FALSE)+
  facet_grid(.~V4)+ 
  ylab("Repeat Content") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=c(met.brewer("Austria",7)))
  
GC_norm_2 %>% group_by(V4) %>% summarise(mean=mean(V3))

GC_norm_3 <- GC_norm_2 %>%
  group_by(V4) %>%
  mutate(delta = V3 - mean(V3)) %>%
  ungroup()

ggplot(GC_norm_3, aes(x=V2,y=delta,col=V1)) + 
  #geom_hline(yintercept=0.3544052,linetype="dashed",col="black") +
  geom_smooth(se=FALSE)+
  #geom_line()+
  facet_grid(.~V4,scales="free_y") + 
  ylab("GC Content") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_color_manual(values=c(met.brewer("Austria",7)))

TEs_v2 <- REs_norm_2 %>% filter(V4=="TEs") %>% filter(V1!="SM_V10_WSR") %>%filter(V1!="SM_V10_Z")
WEs_v2 <- REs_norm_2 %>% filter(V4=="WEs") %>% filter(V1!="SM_V10_WSR") %>%filter(V1!="SM_V10_Z")
GC_un_v2 <- GC_norm_2 %>% filter(V4=="GC_unmask") %>% filter(V1!="SM_V10_WSR") %>%filter(V1!="SM_V10_Z")
GC_ma_v2 <- GC_norm_2 %>% filter(V4=="GC_mask") %>% filter(V1!="SM_V10_WSR") %>%filter(V1!="SM_V10_Z")

setwd("/Users/tobybrann/Documents/FYR/norm_plots/figs")

tiff('TE.tiff', units="in", width=6, height=4, res=600, compression = 'lzw')

ggplot(TEs_v2, aes(x=V2,y=V3,col=V1)) + 
  geom_smooth(se=FALSE)+
  facet_grid(.~V4)+ 
  ylab("Repeat Content") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept=0.3357814,linetype="dashed",col="black") +
  scale_color_manual(values=c(met.brewer("Austria",7)))

dev.off()
tiff('WE.tiff', units="in", width=6, height=4, res=600, compression = 'lzw')

ggplot(WEs_v2, aes(x=V2,y=V3,col=V1)) + 
  #geom_point() +
  geom_smooth(se=FALSE)+
  facet_grid(.~V4)+ 
  ylab("Repeat Content") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept=0.09858602,linetype="dashed",col="black") +
  scale_color_manual(values=c(met.brewer("Austria",7)))

dev.off()
tiff('GC_un.tiff', units="in", width=6, height=4, res=600, compression = 'lzw')

ggplot(GC_un_v2, aes(x=V2,y=V3,col=V1)) + 
  geom_smooth(se=FALSE)+
  facet_grid(.~V4)+ 
  ylab("Repeat Content") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept=0.3544052,linetype="dashed",col="black") +
  scale_color_manual(values=c(met.brewer("Austria",7)))

dev.off()
tiff('GC_ma.tiff', units="in", width=6, height=4, res=600, compression = 'lzw')

ggplot(GC_ma_v2, aes(x=V2,y=V3,col=V1)) + 
  geom_smooth(se=FALSE)+
  facet_grid(.~V4)+ 
  ylab("Repeat Content") + xlab("Chromosome") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) + 
  geom_hline(yintercept=0.3544052,linetype="dashed",col="black") +
  scale_color_manual(values=c(met.brewer("Austria",7)))

dev.off()

###GC V2
GC_V2<-as.data.frame(cbind(GC_mask$V1,GC_mask$V2,GC_unmask$V3-GC_mask$V3,"DeltaGC"))

GC_V2$V2<-as.numeric(GC_V2$V2)
GC_V2$V3<-as.numeric(GC_V2$V3)

ggplot(GC_V2, aes(x=V2,y=V3,col=V1)) + 
  geom_smooth()
  

setwd("/Users/tobybrann/Documents/FYR/GC_cont/length_filt/nucs_out/")
all_nucs<-read.delim(file="nucs_combi.out",sep="\t",header=FALSE)

#genomic "#454a74"
#intron "#6f9969"
#exon "#97c684"
#
all_nucs$V2 <- gsub(".bed", "",all_nucs$V2)
all_nucs$V2 <- gsub("genomic_shuffs", "Genomic",all_nucs$V2)
all_nucs$V2 <- gsub("introns", "Intron",all_nucs$V2)
all_nucs$V2 <- gsub("exon", "Exon",all_nucs$V2)

all_nucs$V2 <-factor(all_nucs$V2, levels =c("Genomic","Intron","Exon","TEs","WEs"))

ggplot(all_nucs,aes(x=V2,y=V1,fill=V2)) +
  geom_jitter(alpha=0.005) +
  geom_boxplot(outlier.shape = NA) + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ylab("GC Content") + 
  xlab("Genomic Feature") + 
  scale_fill_manual(values=c("#454a74","#6f9969","#97c684","#856C61","#FFD1DC")) + 
  ylim(0.15,0.75) + 
  theme(axis.text = element_text(size = 15),legend.text = element_text(size=12), legend.title = element_text(size=15)) + 
  theme(axis.title = element_text(size = 15)) + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) + 
  geom_signif(
    comparisons = list(c("Genomic","Intron"),c("Intron","Exon"),c("Exon","TEs"),c("TEs","WEs")),
    map_signif_level = TRUE,
    step_increase = 0.05,
    y_position = c(0.55, 0.565,0.58,0.595)
  ) +
  stat_summary(fun=mean,geom="point", shape=23, size=3,position=position_dodge(width=0.5)) 

anova <- aov(V1~V2, data=all_nucs)
tukey <- TukeyHSD(anova)

###New genomics
setwd("/Users/tobybrann/Documents/FYR/GC_cont/v2/")
gen_TE<-read.delim(file="genomic_TE.bed",sep="\t",header=TRUE)
gen_WE<-read.delim(file="genomic_WE.bed",sep="\t",header=TRUE)
gen_ex<-read.delim(file="genomic_exon.bed",sep="\t",header=TRUE)
gen_int<-read.delim(file="genomic_intron.bed",sep="\t",header=TRUE)

gen_WE<-data.frame(cbind(gen_WE[,5],"WEs","Genomic"))
gen_TE<-data.frame(cbind(gen_TE[,5],"TEs","Genomic"))
gen_ex<-data.frame(cbind(gen_ex[,5],"Exon","Genomic"))
gen_int<-data.frame(cbind(gen_int[,5],"Intron","Genomic"))

all_gen<-rbind(gen_WE,gen_TE,gen_ex,gen_int)


all_nucs<-cbind(all_nucs %>% filter(V2!="Genomic"),"Actual")

colnames(all_gen)<-c("V1","V2","V3")
colnames(all_nucs)<-c("V1","V2","V3")

all_GC<-rbind(all_gen,all_nucs)

all_GC$V1<-as.numeric(all_GC$V1)

all_GC$V2 <-factor(all_GC$V2, levels =c("Intron","Exon","TEs","WEs"))

tiff('GC_comp.tiff', units="in", width=9, height=4, res=600, compression = 'lzw')

ggplot(all_GC,aes(x=V2,y=V1,fill=interaction(V2,V3))) +
  #geom_jitter(alpha=0.005) +
  geom_boxplot(outlier.shape = NA) + 
  theme(
    legend.position="none",
    plot.title = element_text(size=11)
  ) +
  ylab("GC Content") + 
  xlab("Genomic Feature") + 
  scale_fill_manual(values=c(
      "#e04d4d", "#d04e00",
      "#f6c200", "#0086a8",
      "#fbd9d9", "#ffddcc",
      "#fff2cc", "#cceff7"
    )) + 
  coord_cartesian(ylim=c(0.15,0.65)) + 
  theme(axis.text = element_text(size = 15),legend.text = element_text(size=12), legend.title = element_text(size=15)) + 
  theme(axis.title = element_text(size = 15)) + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) + 
  stat_summary(fun=mean,geom="point", shape=23, size=3,position=position_dodge(width=0.75)) +
  geom_signif(
    xmin = 0.8125, 
    xmax = 1.1875, 
    annotations = "***",
    y_position = 0.48,
    tip_length = 0.005
  ) +
  geom_signif(
    xmin = 1.8125, 
    xmax = 2.1875, 
    annotations = "***",
    y_position = 0.525,
    tip_length = 0.005
  ) +
  geom_signif(
    xmin = 2.8125, 
    xmax = 3.1875, 
    annotations = "***",
    y_position = 0.585,
    tip_length = 0.005
  ) +
  geom_signif(
    xmin = 3.8125, 
    xmax = 4.1875, 
    annotations = "***",
    y_position = 0.575,
    tip_length = 0.005
  ) +
  geom_signif(
    xmin = 0.8125, 
    xmax = 1.8125, 
    annotations = "***",
    y_position = 0.575,
    tip_length = 0.005
  ) +
  geom_signif(
    xmin = 1.8125, 
    xmax = 2.8125, 
    annotations = "***",
    y_position = 0.6,
    tip_length = 0.005
  ) +
  geom_signif(
    xmin = 2.8125, 
    xmax = 3.8125, 
    annotations = "***",
    y_position = 0.625,
    tip_length = 0.005
  )

dev.off()
###Stats
  
 # Perform independent t-tests for each V1 group
anova_results <- all_GC %>%
  group_by(V2) %>%
  do(anova = summary(aov(V1 ~ V3, data = .))) %>%
  summarise(
    V2_group = V2,
    p_value = anova[[1]]$`Pr(>F)`[1]
  )
just_actual <- all_GC %>%
  filter(V3 == "Actual") 
  
anova_results_2<-aov(V1~V2, data=just_actual)
tukey_result <- TukeyHSD(anova_results_2)
summary(tukey_result)

###T.test 
categ<-"WEs"
gen_set<-all_GC %>% filter(V2==categ) %>% filter(V3=="Genomic")
act_set<-all_GC %>% filter(V2==categ) %>% filter(V3=="Actual")
gen_set<-gen_set[,1]
act_set<-act_set[,1]
t.test(act_set,gen_set,var.equal=TRUE)

###Sims work
setwd("/Users/tobybrann/Documents/FYR/overlaps/outputs/")
TEs_act<-read.delim(file="TEs_act.tbl",sep="\t",header=FALSE)
WEs_act<-read.delim(file="WEs_act.tbl",sep="\t",header=FALSE)

TEs_sim<-read.delim(file="TE_all.tbl",sep="\t",header=FALSE)
WEs_sim<-read.delim(file="WE_all.tbl",sep="\t",header=FALSE)

TEs_sim<-cbind(TEs_sim,"TEs","Sim")
WEs_sim<-cbind(WEs_sim,"WEs","Sim")

TEs_sim<-TEs_sim[,-c(3)]
WEs_sim<-WEs_sim[,-c(3)]
colnames(TEs_sim)<-colnames(TEs_act)
colnames(WEs_sim)<-colnames(WEs_act)

combi_overlap<-rbind(TEs_act,WEs_act,TEs_sim,WEs_sim)

combi_overlap<-combi_overlap %>% filter(V1!="Exon")
combi_overlap$V1 <- gsub("5'UTR", "5' UTR",combi_overlap$V1)
combi_overlap$V1 <- gsub("3'UTR", "3' UTR",combi_overlap$V1)

sims <- combi_overlap %>% filter(V4=="Sim")
actual <- combi_overlap %>% filter(V4!="Sim")

sims$V1 <-factor(sims$V1, levels =c("5' UTR","CDS","Intron","3' UTR","Gene"))
actual$V1 <-factor(actual$V1, levels =c("5' UTR","CDS","Intron","3' UTR","Gene"))

ggplot(sims, aes(x=V1,y=V2, fill=V1)) +
  geom_violin() + 
  geom_point(data=actual,aes(x=V1,y=V2,col=V1),size=3,pch=23) + 
  scale_fill_manual(values=met.brewer("Derain",6)) +
  scale_color_manual(values=met.brewer("Derain",6)) +
  ylab("Overlap Length") + xlab("Genomic Feature") + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) + 
  theme(legend.position='none') +
  facet_wrap(.~factor(V3)) 

sims %>% filter(V1=="CDS") %>% filter(V3=="WEs") %>% summarise(mean_value = mean(V2))

###Proportions
setwd("/Users/tobybrann/Documents/FYR/overlaps/props/")
#covs_both <- read.delim(file="catted_covs_both.tbl",sep=" ",header=FALSE)

TE_cov<-read.delim(file="TEs_cov.tbl",sep="\t",header=FALSE)
WE_cov<-read.delim(file="WEs_cov.tbl",sep="\t",header=FALSE)

covs_both<-rbind(TE_cov,WE_cov)

covs_both<-covs_both %>% filter(V1!=0)
covs_both<-covs_both %>% filter(V2!="Exon")

covs_both$V3 <- gsub("Opp", "Opposite",covs_both$V3)
covs_both$V2 <- gsub("5UTR", "5' UTR",covs_both$V2)
covs_both$V2 <- gsub("ThreeUTR", "3' UTR",covs_both$V2)

TE_cov<-covs_both %>% filter(V4=="TE")
WE_cov<-covs_both %>% filter(V4=="WE")

TE_cov$V2 <-factor(TE_cov$V2, levels =c("5' UTR","CDS","Intron","3' UTR","Gene","Intergene"))
WE_cov$V2 <-factor(WE_cov$V2, levels =c("5' UTR","CDS","Intron","3' UTR","Gene","Intergene"))

TE_cov_out<- ggplot(TE_cov, aes(x=V2,y=V1,fill=V2)) + 
  geom_boxplot(width=0.5,outlier.shape = NA) +
  stat_summary(fun=mean,geom="point", shape=23, size=3§,position=position_dodge(width=0.5)) + 
  ylab("Overlap Length") + xlab("Genomic Feature") + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) + 
  scale_fill_manual(values=met.brewer("Derain",6)) +
  theme(legend.position='none') + 
  geom_signif(
  comparisons = list(c("3' UTR","Gene")),
  map_signif_level = FALSE,
  step_increase = 0.05,
  y_position = 0.85,
  tip_length = 0.01,
  annotation = "NS") +
  facet_wrap(V4 ~ .) 

WE_cov_out<-ggplot(WE_cov, aes(x=V2,y=V1,fill=V2)) + 
  geom_boxplot(width=0.5,outlier.shape = NA) +
  #ylim(0,1) +
  stat_summary(fun=mean,geom="point", shape=23, size=3,position=position_dodge(width=0.5)) + 
  ylab("Overlap Length") + xlab("Genomic Feature") + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) + 
  scale_fill_manual(values=met.brewer("Derain",6)) +
  theme(legend.position='none') + 
  geom_signif(
    comparisons = list(c("CDS","5' UTR"),c("Intron","3' UTR"),c("3' UTR","5' UTR"),c("3' UTR","CDS")),
    map_signif_level = FALSE,
    tip_length = 0.01,
    #step_increase = 0.05,
    y_position = c(0.43,0.61,0.49,0.55),
    annotation = c("NS","NS","NS","*")) +
  facet_wrap(V4 ~ .) 

ggarrange(TE_cov_out, WE_cov_out,  ncol = 2)

TE_cov %>% filter(V2=="CDS") %>% summarise(means=median(V1))

###Stats
temp <- covs_both %>% filter(V4=="WE")
anova <- aov(V1~V2, data=temp)
tukey <- TukeyHSD(anova)
tukey


setwd("/Users/tobybrann/Documents/FYR/overlaps/props/")
covs <- read.delim(file="catted_covs.tbl",sep=" ",header=FALSE)

covs_no0<-covs %>% filter(V1!=0)
covs_no0<-covs_no0 %>% filter(V2!="Exon")

covs_no0$V3 <- gsub("Opp", "Opposite",covs_no0$V3)
covs_no0$V2 <- gsub("5UTR", "5' UTR",covs_no0$V2)
covs_no0$V2 <- gsub("3UTR", "3' UTR",covs_no0$V2)

covs_no0$V2 <-factor(covs_no0$V2, levels =c("5' UTR","CDS","Intron","3' UTR","Gene","Intergene"))

ggplot(covs_no0, aes(x=V2,y=V1,fill=V2,col=V3)) + 
  geom_boxplot(width=0.5,outlier.shape=NA) +
  stat_summary(fun=mean,geom="point", shape=23, size=3,position=position_dodge(width=0.5)) + 
  ylab("Overlap Length") + xlab("Genomic Feature") + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) + 
  scale_color_manual(values = c("#000000","#000000")) +  
  scale_fill_manual(values = met.brewer("Derain", 6)) +   
  labs(col = "Strand")

###Intron Occupancy
setwd("/Users/tobybrann/Documents/FYR/introns/intron_occ/")
dros <- read.delim(file="dros_introns_TEs.tbl",sep="\t",header=FALSE)
eleg <- read.delim(file="eleg_introns_TEs.tbl",sep="\t",header=FALSE)
sman <- read.delim(file="sman_introns_TEs.tbl",sep="\t",header=FALSE)
thal <- read.delim(file="thal_introns_TEs.tbl",sep="\t",header=FALSE)

combi_introcc <- rbind(dros,eleg,sman,thal)

combi_introcc$V3 <- gsub("dros", "D. melanogaster",combi_introcc$V3)
combi_introcc$V3 <- gsub("eleg", "C. elegans",combi_introcc$V3)
combi_introcc$V3 <- gsub("thal", "A. thaliana",combi_introcc$V3)
combi_introcc$V3 <- gsub("sman", "S. mansoni",combi_introcc$V3)

cols<-c(met.brewer("Egypt",4)[1],"#bca89f",met.brewer("Egypt",4)[2],"#bca89f",met.brewer("Egypt",4)[3],"#bca89f",met.brewer("Egypt",4)[4],"#bca89f")

tiff('introns.tiff', units="in", width=12, height=6, res=600, compression = 'lzw')

ggplot(combi_introcc,aes(x=V1, fill = interaction(V3,V2, sep="-", lex.order=TRUE))) +
  geom_histogram(color ="#000000", center=-150,binwidth = 300) + 
  xlim(0,30000) +
  scale_fill_manual(values=cols)+
  xlab("Feature Length (bp)") + ylab("Frequency") +theme(legend.position='none') +
  theme(axis.text = element_text(size=14)) + theme(axis.title = element_text(size=20)) + 
  theme(legend.text = element_text(face="italic"))  +scale_y_continuous(trans="log10") +
  theme(strip.text = element_text(size =18, face="italic")) + labs(fill="Species") + facet_grid(V3 ~ .,scales = "free_x") + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) 

dev.off()

LHS<-ggplot(combi_introcc,aes(x=V1, fill = interaction(V3,V2, sep="-", lex.order=TRUE))) +
  geom_histogram(color ="#000000", center=-150,binwidth = 100) + 
  xlim(0,5000) +
  scale_fill_manual(values=cols)+
  xlab("Feature Length (bp)") + ylab("Frequency") +theme(legend.position='none') +
  theme(axis.text = element_text(size=14)) + theme(axis.title = element_text(size=20)) + 
  theme(legend.text = element_text(face="italic"))  +scale_y_continuous(trans="log10") +
  theme(strip.text = element_text(size =18, face="italic")) + labs(fill="Species") + facet_grid(V3 ~ .,scales = "free_x") + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

RHS<-ggplot(combi_introcc,aes(x=V1, fill = interaction(V3,V2, sep="-", lex.order=TRUE))) +
  geom_histogram(color ="#000000", center=-150,binwidth = 300) + 
  xlim(0,30000) +
  scale_fill_manual(values=cols)+
  xlab("Feature Length (bp)") + ylab("Frequency") +theme(legend.position='none') +
  theme(axis.text = element_text(size=14)) + theme(axis.title = element_text(size=20)) + 
  theme(legend.text = element_text(face="italic"))  +scale_y_continuous(trans="log10") +
  theme(strip.text = element_text(size =18, face="italic")) + labs(fill="Species") + facet_grid(V3 ~ .,scales = "free_x") + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) +
  theme(axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

tiff('introns_v2.tiff', units="in", width=12, height=6, res=600, compression = 'lzw')
ggarrange(LHS, RHS,  ncol = 2,widths = c(0.75, 1))
dev.off()
##Std.freq
setwd("/Users/tobybrann/Documents/FYR/introns/intron_freq/calc_std/")
all_freq<-read.delim(file="all_freq.tbl",header=FALSE,sep="\t")

library(dplyr)
all_freq <- all_freq %>%
  group_by(V2, V3) %>%
  mutate(std_dev_V1 = sd(V1))

colnames(all_freq)[4] <- "V4"

all_freq$V3 <- gsub("dros", "D. melanogaster",all_freq$V3)
all_freq$V3 <- gsub("eleg", "C. elegans",all_freq$V3)
all_freq$V3 <- gsub("thal", "A. thaliana",all_freq$V3)
all_freq$V3 <- gsub("sman", "S. mansoni",all_freq$V3)

all_freq_summary <- all_freq %>%
  group_by(V2, V3) %>%
  summarize(mean_V1 = mean(V1), std_dev_V1 = sd(V1))
colnames(all_freq_summary)[3] <- "V1"
colnames(all_freq_summary)[4] <- "V4"

all_freq_summary$V3 <-factor(all_freq_summary$V3, levels =c("A. thaliana","C. elegans","D. melanogaster","S. mansoni"))

ggplot(all_freq_summary,aes(x=V2,y=V1,col=V3)) + 
  scale_color_manual(values=met.brewer("Egypt",4)) +
  geom_hline(yintercept=0, linetype="dashed", linewidth=0.5) + 
  geom_line(size=0.75) +
  geom_errorbar(aes(ymax = V1 + V4, ymin = V1 - V4), 
                width = 10, size=1) +
  ylab("Standardised TE Frequency") +
  xlab("Distance to Exon (bp)") + 
  theme(axis.text = element_text(size = 15),legend.text = element_text(size=12,face="italic"), legend.title = element_text(size=15)) + 
  theme(axis.title = element_text(size = 15)) + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) 

setwd("/Users/tobybrann/Documents/FYR/introns/intron_pos/tbl_out/")
intron_pos <- read.delim(file="tbl_proc.tbl",sep="\t",header=FALSE)

intron_pos$V2 <- gsub(".gff.bed", "",intron_pos$V2)
intron_pos$V2 <- gsub("first", "First",intron_pos$V2)
intron_pos$V2 <- gsub("two", "Second",intron_pos$V2)
intron_pos$V2 <- gsub("three", "Third",intron_pos$V2)
intron_pos$V2 <- gsub("four", "Fourth",intron_pos$V2)
intron_pos$V2 <- gsub("five", "Fifth",intron_pos$V2)
intron_pos$V2 <- gsub("six", "Sixth",intron_pos$V2)
intron_pos$V2 <- gsub("seven", "Seventh",intron_pos$V2)
intron_pos$V2 <- gsub("eight", "Eighth Plus",intron_pos$V2)

intron_pos$V2 <-factor(intron_pos$V2, levels =c("First","Second","Third","Fourth","Fifth","Sixth","Seventh","Eighth Plus"))

ggplot(intron_pos, aes(x=V3,y=V1,col=V2)) + 
  geom_vline(xintercept = 20, linetype="dashed", col="red") +
  geom_vline(xintercept = 220,linetype="dashed", col="red") +
  geom_line(size=2) + 
  scale_color_manual(values=met.brewer("Demuth",8)) + 
  ylab("TE Density") +
  xlab("") + 
  theme(axis.text = element_text(size = 15),legend.text = element_text(size=12), legend.title = element_text(size=15)) + 
  theme(axis.title = element_text(size = 15)) + 
  theme(axis.text = element_text(size=18), axis.title = element_text(size=22, face="bold")) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())





