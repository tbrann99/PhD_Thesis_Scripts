setwd("/Users/tobybrann/Documents/HGT/thesis/host_coverage/outputs")

##Smansoni

sman_Sr3 <- read.delim(file="sman_Sr3.tbl",sep=" ",header=FALSE)
sman_colo <- read.delim(file="sman_colo.tbl",sep=" ",header=FALSE)
sman_columb <- read.delim(file="sman_columbina.tbl",sep=" ",header=FALSE)
sman_saci2 <- read.delim(file="sman_saci2.tbl",sep=" ",header=FALSE)


sman_Sr3<-cbind(sman_Sr3,"mansoni")
colnames(sman_Sr3)[3]<-"V3"

##Bulinus

trun_Sr3 <- read.delim(file="trunc_Sr3.tbl",sep=" ",header=FALSE)

trun_P3<-cbind(trun_P3,"mansoni")
colnames(trun_P3)[3]<-"V3"
trun_Sr3<-cbind(trun_Sr3,"mansoni")
colnames(trun_Sr3)[3]<-"V3"

##P3s
sman_P3 <- read.delim(file="sman_P3.tbl",sep=" ",header=FALSE)
trun_P3 <- read.delim(file="trunc_P3.tbl",sep=" ",header=FALSE)
lymn_P3 <- read.delim(file="lymnaea_P3.tbl",sep=" ",header=FALSE)
stram_P3 <- read.delim(file="stram_P3.tbl",sep=" ",header=FALSE)

sman_P3<-cbind(sman_P3,"mansoni")
colnames(sman_P3)[3]<-"V3"
trun_P3<-cbind(trun_P3,"truncatus")
colnames(trun_P3)[3]<-"V3"
lymn_P3<-cbind(lymn_P3,"lymnaea")
colnames(lymn_P3)[3]<-"V3"
stram_P3<-cbind(stram_P3,"straminaea")
colnames(stram_P3)[3]<-"V3"

all_P3<-rbind(sman_P3,trun_P3,lymn_P3,stram_P3)

ggplot(data=all_P3, aes(x=V1, y=V2, col=V3)) + geom_line(linewidth=1.1) +
  ylab("Coverage (copies)") + xlab("Consensus Sequence Location (bp)") +
  scale_y_continuous(trans="log10")
  #geom_line(data=s_acon, aes(x=V1, y=V2, col='Homology Hit'), linewidth=1.1)


setwd("/Users/tobybrann/Documents/HGT/coverage")


library(ggplot2)
##bioms
b_acon <- read.delim(file="biom_aconcagua_cov.tbl", sep=" ", header=FALSE)
b_pere <- read.delim(file="biom_perere-3_cov.tbl", sep=" ", header=FALSE)
b_sr3 <- read.delim(file="biom_sr3_cov.tbl", sep=" ", header=FALSE)

##smans
s_acon <- read.delim(file="sman_aconcagua_cov.tbl", sep=" ", header=FALSE)
s_pere <- read.delim(file="sman_perere-3_cov.tbl", sep=" ", header=FALSE)
s_sr3 <- read.delim(file="sman_sr3_cov.tbl", sep=" ", header=FALSE)

##combine
c_acon<-cbind(b_acon,s_acon$V2)
colnames(c_acon)[3] <- "V3"

c_pere<-cbind(b_pere,s_pere$V2)
colnames(c_pere)[3] <- "V3"

c_sr3<-cbind(b_sr3,s_sr3$V2)
colnames(c_sr3)[3] <- "V3"

b_acon<-cbind(b_acon,"biom_a")
b_pere<-cbind(b_pere,"biom_p")
b_sr3<-cbind(b_sr3,"biom_s")

s_acon<-cbind(s_acon,"sman_a")
s_pere<-cbind(s_pere,"sman_p")
s_sr3<-cbind(s_sr3,"sman_s")

colnames(b_acon)[3] <- "V3"
colnames(b_pere)[3] <- "V3"
colnames(b_sr3)[3] <- "V3"
colnames(s_acon)[3] <- "V3"
colnames(s_pere)[3] <- "V3"
colnames(s_sr3)[3] <- "V3"

c_acon<-rbind(b_acon,s_acon)
c_pere<-rbind(b_pere,s_pere)
c_sr3<-rbind(b_sr3,s_sr3)

c_acon<-cbind(c_acon,"aconcagua")
c_pere<-cbind(c_pere,"perere-3")
c_sr3<-cbind(c_sr3,"sr3")

colnames(c_acon)[4] <- "V4"
colnames(c_pere)[4] <- "V4"
colnames(c_sr3)[4] <- "V4"

combi_to_assess<-rbind(c_acon,c_pere,c_sr3)

combi_to_assess$V4 <- gsub("aconcagua", "Aconcagua",combi_to_assess$V4)
combi_to_assess$V4 <- gsub("perere-3", "Perere-3",combi_to_assess$V4)
combi_to_assess$V4 <- gsub("sr3", "Sr3",combi_to_assess$V4)

col_array2<-c(Pr3,Sr3,rep(green_col,5))

col_array2<-c("#94b78e","#3587a2","#f49567","#4b7045","#0d425c","#b15324")
##Greenish
#Lighter: #94b78e
#  Darker: #4b7045
#  2. #156082 (Teal-Blue):
#Lighter: #3587a2
#  Darker: #0d425c
#  3. #EB7131 (Orange):
#Lighter: #f49567
#  Darker: #b15324

setwd("/Users/tobybrann/Documents/HGT/thesis/pre_figs")

tiff('B_stram_cov.tiff', units="in", width=12, height=4, res=600, compression = 'lzw')

ggplot(data=combi_to_assess, aes(x=V1, y=V2, col=V3)) + geom_line(linewidth=1.1) +
  ylab("Coverage (copies)") + xlab("Consensus Sequence Location (bp)") +
  scale_color_manual(values=col_array2) + 
  scale_y_continuous(trans="log10") + 
  theme_minimal() +
  facet_wrap(.~factor(V4, levels=c("Aconcagua","Perere-3","Sr3")), scales="free",ncol=3) +
  theme(legend.position='none')

dev.off()  

###By Line Type
combi_to_assess<-cbind(combi_to_assess,combi_to_assess$V3)
colnames(combi_to_assess)[5]<-"V5"
combi_to_assess$V3 <- gsub("biom_a", "Biom",combi_to_assess$V3)
combi_to_assess$V3 <- gsub("biom_p", "Biom",combi_to_assess$V3)
combi_to_assess$V3 <- gsub("biom_s", "Biom",combi_to_assess$V3)
combi_to_assess$V3 <- gsub("sman_a", "Sman",combi_to_assess$V3)
combi_to_assess$V3 <- gsub("sman_p", "Sman",combi_to_assess$V3)
combi_to_assess$V3 <- gsub("sman_s", "Sman",combi_to_assess$V3)

setwd("/Users/tobybrann/Documents/HGT/thesis/pre_figs")

tiff('B_stram_cov_v2.tiff', units="in", width=12, height=4, res=600, compression = 'lzw')

ggplot(data=combi_to_assess, aes(x=V1, y=V2, col=V5, linetype=V3)) + 
  geom_line(linewidth=0.8) +
  ylab("Coverage (copies)") + xlab("Consensus Sequence Location (bp)") +
  scale_color_manual(values=col_array2) + 
  scale_y_continuous(trans="log10") + 
  theme_minimal() +
  facet_wrap(.~factor(V4, levels=c("Aconcagua","Perere-3","Sr3")), scales="free",ncol=3) +
  theme(legend.position='none')

dev.off()

##Plot
ggplot(data=b_acon, aes(x=V1, y=V2, col='Genome')) + geom_line(linewidth=1.1) +
  ylab("Coverage (copies)") + xlab("Consensus Sequence Location (bp)") +
  geom_line(data=s_acon, aes(x=V1, y=V2, col='Homology Hit'), linewidth=1.1) +
  labs(col='')

ggplot(data=b_pere, aes(x=V1, y=V2, col='Genome')) + geom_line(linewidth=1.1) +
  ylab("Coverage (copies)") + xlab("Consensus Sequence Location (bp)") +
  geom_line(data=s_pere, aes(x=V1, y=V2, col='Homology Hit'), linewidth=1.1) +
  labs(col='') + scale_y_continuous(trans="log10")

ggplot(data=b_sr3, aes(x=V1, y=V2, col='Genome')) + geom_line(linewidth=1.1) +
  ylab("Coverage (copies)") + xlab("Consensus Sequence Location (bp)") +
  geom_line(data=s_sr3, aes(x=V1, y=V2, col='Homology Hit'), linewidth=1.1) +
  labs(col='') + scale_y_continuous(trans="log10")


##
setwd("/Users/tobybrann/Documents/HGT")
TE <- read.delim(file="TE_counts.txt", header=FALSE, sep="\t")

#TE$V2 <- factor(TE$V2 , levels=TE[order(TE$V1),], ordered=TRUE)

TE2 <- TE %>% filter(V1>=100)


Pr3<-"#156082"
Sr3<-"#EB7131"
green_col<-met.brewer("Derain",6)[3]

col_array<-c(Pr3,Sr3,rep(green_col,7))

TE2$V2 <- factor(TE2$V2, levels = TE2$V2[order(-TE2$V1)])

setwd("/Users/tobybrann/Documents/HGT/thesis/pre_figs")

tiff('B_stram_counts.tiff', units="in", width=8, height=4, res=600, compression = 'lzw')

ggplot(TE2, aes(x=reorder(V2, -V1), y=V1,fill=V2)) + geom_bar(stat="identity") +
  theme(axis.text.x = element_text(angle = 75, vjust = 1, hjust=1)) +
  xlab("Transposable Element") + ylab("Number of Copies in B. straminaea") + 
  theme_minimal() + 
  scale_fill_manual(values=col_array) + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  theme(legend.position='none') 

dev.off()



##TBL
#setwd("/Users/tobybrann/Documents/HGT/polished/new_library/2E")
setwd("/Users/tobybrann/Documents/HGT/polished/Figure_1/term_branch_length")
TEs.tbl <- read.delim(file="temp.txt",sep=" ",header=FALSE)

TEs.tbl$V2 <- gsub("\\.tbl$", "", TEs.tbl$V2)

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

combi_3.5[36,4]<-"Perere-3"
combi_3.5[61,4]<-"Sr3"


temp<-c("#156082","#EB7131",met.brewer("Derain",6)[3])

combi_3.5$Family <- gsub("LINE_RTE", "Other",combi_3.5$Family)
combi_3.5$Family <- gsub("LINE_Jockey", "Other",combi_3.5$Family)
combi_3.5$Family <- gsub("LINE_R2", "Other",combi_3.5$Family)
combi_3.5$Family <- gsub("LTR_Gypsy", "Other",combi_3.5$Family)
combi_3.5$Family <- gsub("LTR_Pao", "Other",combi_3.5$Family)
combi_3.5$Family <- gsub("PLE", "Other",combi_3.5$Family)

combi_3.5$Family<-factor(combi_3.5$Family,levels=c("Perere-3","Sr3","Other"))
setwd("/Users/tobybrann/Documents/HGT/thesis/tbl")
tiff('TBL_num0.tiff', units="in", width=12, height=4, res=600, compression = 'lzw')

ggplot(combi_3.5, aes(x=Copy_Number,y=Median,size=num0,col=Family)) +
  geom_hline(yintercept=0.056, colour="red",linetype="dashed") + 
  scale_colour_manual(values = temp)  +
  geom_point(shape = 19) + 
  scale_x_continuous(trans="log10") + 
  ylim(-0.01,0.27)+
  scale_size(range=c(3,20))+ 
  xlab("Copy Number") + 
  theme_minimal()

dev.off()



library(treemapify)
library(ggplot2)
library(viridis)



setwd("/Users/tobybrann/Documents/HGT/polished/Figure_1/genome_coverage")
treemap <- read.delim(file="treemap_gen_v2.out",sep="\t",header=FALSE)


col_arr<-rep("#6f9969",80)
setwd("/Users/tobybrann/Documents/HGT/thesis/treemap")
tiff('gen_cov.tiff', units="in", width=10, height=8, res=600, compression = 'lzw')

ggplot(treemap, aes(area = (V3/391395382)*100, fill = V1, label=V1)) +
  geom_treemap(color = "black", size=2) +
  geom_treemap_text(colour = "white",
                    place = "centre",
                    size = 25,
                    grow=TRUE,
                    fontface="bold")  +
  labs(colour="Percentage of Genome") + 
  scale_fill_manual(values=col_arr) + 
  theme(legend.position='none') 

dev.off()




###Expression
setwd("/Users/tobybrann/Documents/HGT/basics/expression/")
TE_lens <- read.delim(file="/Users/tobybrann/Documents/HGT/basics/expression/TE_lengths.txt", sep="\t", header=FALSE)

raw2 <- read.delim(file="TE_expression.tbl", sep="\t",header=FALSE)

colnames(raw2) <- raw2[2,]
raw2<-raw2[-c(1,2),]
raw2$Length<-as.numeric(raw2$Length)

colnames(TE_lens)[1] <- "Geneid"
raw_wlen<-merge(raw2, TE_lens, by = "Geneid")

raw_wlen<-cbind(raw_wlen$Geneid,raw_wlen$V2,raw_wlen)

raw_wlen<-raw_wlen[,c(1,2,5:49)]

colnames(raw_wlen)[1]<-"Geneid"
colnames(raw_wlen)[2]<-"TE_Length"

FPKM_output<-raw_wlen

##FPKM per cell
for(i in seq(from=3, to=ncol(raw_wlen))){
  for(j in seq(from=1,to=nrow(raw_wlen))){
    yield<-sum(as.numeric(raw_wlen[,i]))
    total_length<-raw_wlen[j,2]
    FPKM<-(as.numeric(raw_wlen[j,i])*1e9)/(yield*total_length)
    FPKM_output[j,i]<-FPKM
  }
}
TPM_output<-FPKM_output
#Calculate TPMs
for (i in seq(from = 3, to = ncol(FPKM_output))) {
  FPKM_output[, i]<-as.numeric(FPKM_output[, i])
  sum_fpkm <- sum(FPKM_output[, i])
  TPM_output[, i] <- (FPKM_output[, i] / sum_fpkm) * 1e6
}


colnames(TPM_output)<-c("TE","Total_Length","F_26d_Juveniles","F_26d_Juveniles", "F_26d_Juveniles", "F_26d_Juveniles", "F_26d_Juveniles",
                      "M_26d_Juveniles","M_26d_Juveniles","M_26d_Juveniles","M_26d_Juveniles","M_26d_Juveniles",
                      "X_1d_Sporocysts","X_1d_Sporocysts","X_1d_Sporocysts","X_1d_Sporocysts","X_1d_Sporocysts",
                      "X_2d_Somules","X_2d_Somules","X_2d_Somules","X_2d_Somules","X_2d_Somules",
                      "X_32d_Sporocysts","X_32d_Sporocysts","X_32d_Sporocysts","X_32d_Sporocysts","X_32d_Sporocysts",
                      "X_5d_Sporocysts", "X_5d_Sporocysts","X_5d_Sporocysts","X_5d_Sporocysts","X_5d_Sporocysts",
                      "X_Cercariae","X_Cercariae","X_Cercariae","X_Cercariae","X_Cercariae",
                      "X_Eggs","X_Eggs","X_Eggs","X_Eggs","X_Eggs",
                      "X_Miracidia","X_Miracidia","X_Miracidia","X_Miracidia","X_Miracidia")

TPM_output<-TPM_output[,-2]
row.names(TPM_output) <- TPM_output[,1]
TPM_output<-TPM_output[,-1]

TPM_output<-cbind(rownames(TPM_output),TPM_output)
colnames(TPM_output)[1]<-"TE"
t<-melt(TPM_output, id.vars = "TE")

t$variable<-as.character(t$variable)

t$variable<-gsub(".1","",t$variable)
t$variable<-gsub(".2","",t$variable)
t$variable<-gsub(".3","",t$variable)
t$variable<-gsub(".4","",t$variable)

t$value<-as.numeric(t$value)

TE_family <- read.delim(file="/Users/tobybrann/Documents/HGT/basics/expression/TE_family.nums", sep="\t", header=FALSE)

colnames(TE_family)[1]<-"TE"

t<-merge(t, TE_family, by = "TE")

mean_values <- t %>%
  group_by(TE,variable) %>%
  summarise(Mean_Value = mean(value))

t2<-merge(mean_values, TE_family, by = "TE")

t3<-t2 %>% filter(V2!="SINE")
t3<-t3 %>% filter(V2!="Unknown")

for(i in seq(from=1, to=nrow(t3))){
  if(t3[i,1]=="Perere-3"){
    t3[i,4]<-"Perere-3"
  }
}

for(i in seq(from=1, to=nrow(t3))){
  if(t3[i,1]=="Sr3"){
    t3[i,4]<-"Sr3"
  }
}

t3$V2<-gsub("_CACTA","",t3$V2)
t3$V2<-gsub("_Merlin","",t3$V2)
t3$V2<-gsub("_Mut","",t3$V2)
t3$V2<-gsub("_Gypsy","",t3$V2)
t3$V2<-gsub("_Pao","",t3$V2)
t3$V2<-gsub("_Unknown","",t3$V2)
t3$V2<-gsub("_Jockey","",t3$V2)
t3$V2<-gsub("_R2","",t3$V2)
t3$V2<-gsub("DNA","Other",t3$V2)
t3$V2<-gsub("LTR","Other",t3$V2)
t3$V2<-gsub("PLE","Other",t3$V2)
t3$V2<-gsub("LINE","Other",t3$V2)
t3$V2<-gsub("Other_RTE","LINE_RTE",t3$V2)
#t3$V2<-gsub("LINE_RTE","Other",t3$V2)

t3$variable<-gsub("X_","",t3$variable)

library(MetBrewer)

t3$variable<-gsub("Xd_Somules","2d_Somules",t3$variable)
t3$variable<-gsub("5d_Sporocyst","temp",t3$variable)
t3$variable<-gsub("Xd_Sporocyst","tran",t3$variable)
t3$variable<-gsub("d_Sporocyst","2d_Sporocyst",t3$variable)
t3$variable<-gsub("temp","5d_Sporocyst",t3$variable)
t3$variable<-gsub("tran","32d_Sporocyst",t3$variable)

t3$variable <- factor(t3$variable, levels=c("Eggs", "Miracidia","2d_Sporocysts","5d_Sporocysts","32d_Sporocysts","Cercariae","2d_Somules","F6d_Juveniles","M6d_Juveniles"))
t3$V2 <- factor(t3$V2, levels=c("Perere-3","Sr3","LINE_RTE","Other"))

Pr3<-"#156082"
Sr3<-"#EB7131"
green_col<-met.brewer("Derain",6)[3]
other_col<-met.brewer("Derain",6)[1]

col_array<-c(Pr3,Sr3,other_col,rep(green_col,1))

setwd("/Users/tobybrann/Documents/HGT/thesis/express")

tiff('FPKM_TEs.tiff', units="in", width=12, height=6, res=600, compression = 'lzw')

ggplot(t3, aes(x=variable,y=Mean_Value,col=V2)) + 
  geom_jitter(size=3) +
  #scale_y_continuous(trans="log10") + 
  ylab("FPKM") + xlab("Time Stage") + 
  scale_color_manual(values = col_array) +
  theme_minimal() + 
  theme(legend.position='none') 

dev.off()

t3 %>%
  filter(V2 %in% c("Perere-3", "Sr3")) %>%
  group_by(variable) %>%             
  summarize(sum = sum((Mean_Value), na.rm = TRUE), .groups = "drop")

t3 %>%
  filter(V2 %in% c("Perere-3", "Sr3","LINE_RTE")) %>%
  group_by(variable) %>%             
  summarize(sum = sum((Mean_Value), na.rm = TRUE), .groups = "drop")

setwd("/Users/tobybrann/Documents/HGT/thesis/comp_P3_Sr3/")

blastout<-read.delim(file="blastn_summary_namefix_dupefix_remove.tbl",sep="\t",header = FALSE)

library(tidyr)

data_split <- blastout %>%
  separate(V5, into = c("v5_1", "v5_2", "v5_3","v5_4"), sep = ",", fill = "right") %>%
  separate(V6, into = c("v6_1", "v6_2", "v6_3","v6_4"), sep = ",", fill = "right")

data_split <- data_split %>%
  mutate(v5_4 = ifelse(is.na(v5_4), 0, v5_4))

data_split <- data_split %>%
  mutate(v6_4 = ifelse(is.na(v6_4), 0, v6_4))

data_split <- data_split %>%
  mutate(v5_2 = ifelse(v5_2 == 0, NA, v5_2))

data_split <- data_split %>%
  mutate(v6_2 = ifelse(v6_2 == 0, NA, v6_2))


max_len_1<-data_split[,6]
max_len_2<-data_split[,10]

max_len<-as.data.frame(cbind(max_len_1,max_len_2))

max_len$max_len_1<-as.numeric(max_len$max_len_1)
max_len$max_len_2<-as.numeric(max_len$max_len_2)

ggplot(max_len,aes(x=max_len_1,y=max_len_2)) + 
  geom_point() 
  

##Number over 1,000bp
over_1000_1<-data_split[,8]
over_1000_2<-data_split[,12]

over_1000<-as.data.frame(cbind(over_1000_1,over_1000_2))

over_1000$over_1000_1<-as.numeric(over_1000$over_1000_1)
over_1000$over_1000_2<-as.numeric(over_1000$over_1000_2)

over_1000_non0<-over_1000 %>% filter(over_1000_1>0 & over_1000_2>0)

model <- lm(over_1000_non0$over_1000_1 ~ over_1000_non0$over_1000_2)
summary(model)$r.squared
cor_xy2 <- cor(over_1000_non0$over_1000_1, over_1000_non0$over_1000_2, method = "pearson")

tiff('Hits_1000.tiff', units="in", width=5, height=4, res=600, compression = 'lzw')

ggplot(over_1000_non0,aes(x=over_1000_1,y=over_1000_2)) + 
  scale_y_continuous(trans="log10")+
  scale_x_continuous(trans="log10")+ 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_point(size=2, col="#856C61") + 
  ylab("Number of Sr3 Hits Over 1kbp") + 
  xlab("Number of Perere-3 Hits Over 1kbp") + 
  theme_minimal() + 
  ggtitle("0.9983706")
  
dev.off()

##Evalue
eval_1<-data_split[,7]
eval_2<-data_split[,11]

eval<-as.data.frame(cbind(eval_1,eval_2))

eval_non0<-eval %>% filter(!is.na(eval_1) & !is.na(eval_2))

eval_non0$eval_1<-as.numeric(eval_non0$eval_1)
eval_non0$eval_2<-as.numeric(eval_non0$eval_2)

eval_non0<-eval_non0 %>% filter(eval_1!=0 & eval_2!=0)

eval_non0 <- eval_non0 %>%
  mutate(across(everything(), ~ sub(".*e-", "", .)))

eval_non0$eval_1<-as.numeric(eval_non0$eval_1)
eval_non0$eval_2<-as.numeric(eval_non0$eval_2)

model <- lm(eval_non0$eval_1 ~ eval_non0$eval_2)
summary(model)$r.squared
cor_xy <- cor(eval_non0$eval_1, eval_non0$eval_2, method = "pearson")

tiff('eval_hits.tiff', units="in", width=5, height=4, res=600, compression = 'lzw')

ggplot(eval_non0,aes(x=eval_1,y=eval_2)) + 
  geom_smooth(method = "lm", se = FALSE, color = "black") +
  geom_point(size=2, col="#856C61") + 
  ylab("E Value (e-) of Perere-3 Hit") + 
  xlab("E Value (e-) of Sr3 Hit") +
  theme_minimal() + 
  ggtitle("0.8932749")

dev.off()


###Expression v2
setwd("/Users/tobybrann/Documents/HGT/basics/expression/")
TE_lens <- read.delim(file="/Users/tobybrann/Documents/HGT/basics/expression/TE_lengths.txt", sep="\t", header=FALSE)

raw2 <- read.delim(file="TE_expression.tbl", sep="\t",header=FALSE)

colnames(raw2) <- raw2[2,]
raw2<-raw2[-c(1,2),]
raw2$Length<-as.numeric(raw2$Length)

colnames(TE_lens)[1] <- "Geneid"
raw_wlen<-merge(raw2, TE_lens, by = "Geneid")

raw_wlen<-cbind(raw_wlen$Geneid,raw_wlen$V2,raw_wlen)

raw_wlen<-raw_wlen[,c(1,2,5:49)]

colnames(raw_wlen)[1]<-"Geneid"
colnames(raw_wlen)[2]<-"TE_Length"

FPKM_output<-raw_wlen
fractions <- 

##FPKM per cell
for(i in seq(from=3, to=ncol(raw_wlen))){
    yield<-sum(as.numeric(raw_wlen[,i]))
    fraction<-as.numeric(raw_wlen[,i])/yield
    FPKM_output[,i]<-fraction*1e6
}


colnames(FPKM_output)<-c("TE","Total_Length","F_26d_Juveniles","F_26d_Juveniles", "F_26d_Juveniles", "F_26d_Juveniles", "F_26d_Juveniles",
                        "M_26d_Juveniles","M_26d_Juveniles","M_26d_Juveniles","M_26d_Juveniles","M_26d_Juveniles",
                        "X_1d_Sporocysts","X_1d_Sporocysts","X_1d_Sporocysts","X_1d_Sporocysts","X_1d_Sporocysts",
                        "X_2d_Somules","X_2d_Somules","X_2d_Somules","X_2d_Somules","X_2d_Somules",
                        "X_32d_Sporocysts","X_32d_Sporocysts","X_32d_Sporocysts","X_32d_Sporocysts","X_32d_Sporocysts",
                        "X_5d_Sporocysts", "X_5d_Sporocysts","X_5d_Sporocysts","X_5d_Sporocysts","X_5d_Sporocysts",
                        "X_Cercariae","X_Cercariae","X_Cercariae","X_Cercariae","X_Cercariae",
                        "X_Eggs","X_Eggs","X_Eggs","X_Eggs","X_Eggs",
                        "X_Miracidia","X_Miracidia","X_Miracidia","X_Miracidia","X_Miracidia")

TPM_output <- FPKM_output

TPM_output<-TPM_output[,-2]
row.names(TPM_output) <- TPM_output[,1]
TPM_output<-TPM_output[,-1]

TPM_output<-cbind(rownames(TPM_output),TPM_output)
colnames(TPM_output)[1]<-"TE"
t<-melt(TPM_output, id.vars = "TE")

t$variable<-as.character(t$variable)

t$variable<-gsub(".1","",t$variable)
t$variable<-gsub(".2","",t$variable)
t$variable<-gsub(".3","",t$variable)
t$variable<-gsub(".4","",t$variable)

t$value<-as.numeric(t$value)

TE_family <- read.delim(file="/Users/tobybrann/Documents/HGT/basics/expression/TE_family.nums", sep="\t", header=FALSE)

colnames(TE_family)[1]<-"TE"

t<-merge(t, TE_family, by = "TE")

mean_values <- t %>%
  group_by(TE,variable) %>%
  summarise(Mean_Value = mean(value))

t2<-merge(mean_values, TE_family, by = "TE")

t3<-t2 %>% filter(V2!="SINE")
t3<-t3 %>% filter(V2!="Unknown")

for(i in seq(from=1, to=nrow(t3))){
  if(t3[i,1]=="Perere-3"){
    t3[i,4]<-"Perere-3"
  }
}

for(i in seq(from=1, to=nrow(t3))){
  if(t3[i,1]=="Sr3"){
    t3[i,4]<-"Sr3"
  }
}

t3$V2<-gsub("_CACTA","",t3$V2)
t3$V2<-gsub("_Merlin","",t3$V2)
t3$V2<-gsub("_Mut","",t3$V2)
t3$V2<-gsub("_Gypsy","",t3$V2)
t3$V2<-gsub("_Pao","",t3$V2)
t3$V2<-gsub("_Unknown","",t3$V2)
t3$V2<-gsub("_Jockey","",t3$V2)
t3$V2<-gsub("_R2","",t3$V2)
t3$V2<-gsub("DNA","Other",t3$V2)
t3$V2<-gsub("LTR","Other",t3$V2)
t3$V2<-gsub("PLE","Other",t3$V2)
t3$V2<-gsub("LINE","Other",t3$V2)
t3$V2<-gsub("Other_RTE","LINE_RTE",t3$V2)
#t3$V2<-gsub("LINE_RTE","Other",t3$V2)

t3$variable<-gsub("X_","",t3$variable)

library(MetBrewer)

t3$variable<-gsub("Xd_Somules","2d_Somules",t3$variable)
t3$variable<-gsub("5d_Sporocyst","temp",t3$variable)
t3$variable<-gsub("Xd_Sporocyst","tran",t3$variable)
t3$variable<-gsub("d_Sporocyst","2d_Sporocyst",t3$variable)
t3$variable<-gsub("temp","5d_Sporocyst",t3$variable)
t3$variable<-gsub("tran","32d_Sporocyst",t3$variable)

t3$variable <- factor(t3$variable, levels=c("Eggs", "Miracidia","2d_Sporocysts","5d_Sporocysts","32d_Sporocysts","Cercariae","2d_Somules","F6d_Juveniles","M6d_Juveniles"))
t3$V2 <- factor(t3$V2, levels=c("Perere-3","Sr3","LINE_RTE","Other"))

Pr3<-"#156082"
Sr3<-"#EB7131"
green_col<-met.brewer("Derain",6)[3]
other_col<-met.brewer("Derain",6)[1]

col_array<-c(Pr3,Sr3,other_col,rep(green_col,1))

setwd("/Users/tobybrann/Documents/HGT/thesis/express")

options(scipen=999)

tiff('FPKM_TEs_v2.tiff', units="in", width=12, height=6, res=600, compression = 'lzw')

ggplot(t3, aes(x=variable,y=Mean_Value,col=V2)) + 
  geom_jitter(size=3,width = 0.2) +
  #scale_y_continuous(trans="log10") + 
  ylab("FPKM") + xlab("Time Stage") + 
  scale_color_manual(values = col_array) +
  theme_minimal() + 
  theme(legend.position='none') 

dev.off()

t3 %>%
  filter(V2 %in% c("Perere-3", "Sr3")) %>%
  group_by(variable) %>%             
  summarize(sum = sum((Mean_Value), na.rm = TRUE), .groups = "drop")

t3 %>%
  filter(V2 %in% c("Perere-3", "Sr3","LINE_RTE")) %>%
  group_by(variable) %>%             
  summarize(sum = sum((Mean_Value), na.rm = TRUE), .groups = "drop")














