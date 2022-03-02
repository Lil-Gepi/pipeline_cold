#!/usr/bin/env Rscript
# install.packages("~/Downloads/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")
#?xtable
setwd("~/Dropbox (PopGen)/Yiwen/RS/data/")
library(poolSeq)
library(ggplot2)
library(factoextra)
library(reshape2)
library(data.table)
#Functions----
PCA <- function(df,c,s,coloring_factor)
{
  pcadata<-na.omit(df)#remove missing data
  pca.res <- prcomp(pcadata, retx=TRUE, center = c, scale. =s)
  
  print(fviz_pca_ind(X = pca.res, 
                     geom = "point",
                     col.ind=coloring_factor,
                     repel = TRUE,   # Avoid text overlapping
                     mean.point = FALSE,
                     addEllipses=FALSE)) # Default parameters for ellipses
}



replicate <- c(10,10,10,10,10,10,10,10,20,20,20,20,20,20,20,20,30,30,30,30,30,30,
                30,30,30,30,30,30,30,30,30,30,30,30,30,30,40,40,40,40,40,40,40,40,
                50,50,50,50,50,50,50,50,60,60,60,60,80,80,80,80,80,80,80,80,0,0,0,0)
generation <- c(1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,1,2,3,4,5,6,7,8,1,2,3,4,5,6,7,8,1,2,3,4,1,2,3,4,5,6,7,8,1,2,3,4)
#Data handling----
my.sync <- read.sync(file="./grand_PCA_all_generations_all_chr_mq20_bq20_Neda_filtered.sync", 
                 repl=replicate, gen=generation)
my.af <- as.data.frame(af(my.sync, repl=replicate, gen=generation))

my.af <- my.af[which(apply(my.af,1,var)!=0),]

colnames(my.af) = c("F0-R11.freq","F0-R12.freq","F0-R13.freq","F0-R14.freq",
                    "F10-R11S1.freq","F10-R11S2.freq","F10-R12S1.freq","F10-R12S2.freq",
                    "F10-R13S1.freq","F10-R13S2.freq","F10-R14S1.freq","F10-R14S2.freq",
                    "F20-R11S1.freq","F20-R11S2.freq","F20-R12S1.freq","F20-R12S2.freq",
                    "F20-R13S1.freq","F20-R13S2.freq","F20-R14S1.freq","F20-R14S2.freq",
                    "F30-R11S1.freq","F30-R11S2.freq","F30-R11S3.freq","F30-R11S4.freq","F30-R11S5.freq",
                    "F30-R12S1.freq","F30-R12S2.freq","F30-R12S3.freq","F30-R12S4.freq","F30-R12S5.freq",
                    "F30-R13S1.freq","F30-R13S2.freq","F30-R13S3.freq","F30-R13S4.freq","F30-R13S5.freq",
                    "F30-R14S1.freq","F30-R14S2.freq","F30-R14S3.freq","F30-R14S4.freq","F30-R14S5.freq",
                    "F40-R11S1.freq","F40-R11S2.freq","F40-R12S1.freq","F40-R12S2.freq",
                    "F40-R13S1.freq","F40-R13S2.freq","F40-R14S1.freq","F40-R14S2.freq",
                    "F50-R11S1.freq","F50-R11S2.freq","F50-R12S1.freq","F50-R12S2.freq",
                    "F50-R13S1.freq","F50-R13S2.freq","F50-R14S1.freq","F50-R14S2.freq",
                    "F60-R11S1.freq","F60-R11S2.freq","F60-R12S1.freq","F60-R12S2.freq",
                    "F80-R11S1.freq","F80-R11S2.freq","F80-R12S1.freq","F80-R12S2.freq",
                    "F80-R13S1.freq","F80-R13S2.freq","F80-R14S1.freq","F80-R14S2.freq")
#frequencies in dat must be transformed
my.af = 2*asin(sqrt(my.af))
#meta information
clmnames <- colnames(my.af)
sampleIDs <- sub(".freq","",clmnames)
splits <- unlist(strsplit(sampleIDs, "[-]"))
gen <- as.factor(splits[grep("F",splits)])
rep <- as.factor(splits[grep("R",splits)])
#transpose data frame
t_my.af <- as.data.frame(t(my.af))
#set sampleIDs as rownames for labeling in the PCA
row.names(t_my.af) <- sampleIDs
t_my.af <- na.omit(t_my.af)#remove missing data


pdf("grand_PCA_arcsin_test2.pdf")
PCA(t_my.af,TRUE,TRUE,rep)
PCA(t_my.af,TRUE,TRUE,gen)
dev.off()
rm(clmnames, gen, rep, sampleIDs, splits)

####r11----
r11_af <- my.af[,c(1,5,6,13,14,21,22,23,24,25,41,42,49,50,57,58,61,62)]
r11_af <- r11_af[which(apply(r11_af,1,var)!=0),]
#meta information
clmnames <- colnames(r11_af)
sampleIDs <- sub(".freq","",clmnames)
splits <- unlist(strsplit(sampleIDs, "[-]"))
gen <- as.factor(splits[grep("F",splits)])
rep <- as.factor(splits[grep("R",splits)])
#transpose data frame
t_r11_af <- as.data.frame(t(r11_af))
#set sampleIDs as rownames for labeling in the PCA
row.names(t_r11_af) <- sampleIDs
pdf("PCA_r11_all_generation_arcsin.pdf")
PCA(t_r11_af,TRUE,TRUE,rep)
PCA(t_r11_af,TRUE,TRUE,gen)
dev.off()
rm(clmnames, gen, rep, sampleIDs, splits)

####r12----
r12_af <- my.af[,c(2,7,8,15,16,26,27,28,29,30,43,44,51,52,59,60,63,64)]
r12_af <- r12_af[which(apply(r12_af,1,var)!=0),]
#meta information
clmnames <- colnames(r12_af)
sampleIDs <- sub(".freq","",clmnames)
splits <- unlist(strsplit(sampleIDs, "[-]"))
gen <- as.factor(splits[grep("F",splits)])
rep <- as.factor(splits[grep("R",splits)])
#transpose data frame
t_r12_af <- as.data.frame(t(r12_af))
#set sampleIDs as rownames for labeling in the PCA
row.names(t_r12_af) <- sampleIDs
pdf("PCA_r12_all_generation_arcsin.pdf")
PCA(t_r12_af,TRUE,TRUE,rep)
PCA(t_r12_af,TRUE,TRUE,gen)
dev.off()
rm(clmnames, gen, rep, sampleIDs, splits)


####r13----
r13_af <- my.af[,c(3,9,10,17,18,31,32,33,34,35,45,46,53,54,65,66)]
r13_af <- r13_af[which(apply(r13_af,1,var)!=0),]
#meta information
clmnames <- colnames(r13_af)
sampleIDs <- sub(".freq","",clmnames)
splits <- unlist(strsplit(sampleIDs, "[-]"))
gen <- as.factor(splits[grep("F",splits)])
rep <- as.factor(splits[grep("R",splits)])
#transpose data frame
t_r13_af <- as.data.frame(t(r13_af))
#set sampleIDs as rownames for labeling in the PCA
row.names(t_r13_af) <- sampleIDs
pdf("PCA_r13_all_generation_arcsin.pdf")
PCA(t_r13_af,TRUE,TRUE,rep)
PCA(t_r13_af,TRUE,TRUE,gen)
dev.off()
rm(clmnames, gen, rep, sampleIDs, splits)


####r14----
r14_af <- my.af[,c(4,11,12,19,20,36,37,38,39,40,47,48,55,56,67,68)]
r14_af <- r14_af[which(apply(r14_af,1,var)!=0),]
#meta information
clmnames <- colnames(r14_af)
sampleIDs <- sub(".freq","",clmnames)
splits <- unlist(strsplit(sampleIDs, "[-]"))
gen <- as.factor(splits[grep("F",splits)])
rep <- as.factor(splits[grep("R",splits)])
#transpose data frame
t_r14_af <- as.data.frame(t(r14_af))
#set sampleIDs as rownames for labeling in the PCA
row.names(t_r14_af) <- sampleIDs
pdf("PCA_r14_all_generation_arcsin.pdf")
PCA(t_r14_af,TRUE,TRUE,rep)
PCA(t_r14_af,TRUE,TRUE,gen)
dev.off()
rm(clmnames, gen, rep, sampleIDs, splits)



