#!/usr/bin/env Rscript
# install.packages("~/Downloads/poolSeq-0.3.5.tar.gz", repos=NULL, type="source")
#?xtable
setwd("~//RS/F30_all/result")
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
  #Visualize eigenvalues (scree plot).
  #print(fviz_eig(pca.res))

  #Show the percentage of variances explained by each principal component.
  #print(fviz_contrib(pca.res, choice = "var", axes = 1:2, top = 10))

  #if mean.point = TRUE group means are added for each group, specified by coloring_factor
  print(fviz_pca_ind(pca.res, geom = "point",
                     col.ind=coloring_factor,
                     repel = TRUE,   # Avoid text overlapping
                     mean.point = FALSE,
                     addEllipses=FALSE)) # Default parameters for ellipses
}

PCAcos <- function(df,c,s,coloring_factor)
{
  pcadata <- na.omit(df)#remove missing data
  pca.res <- prcomp(pcadata, retx=TRUE, center = c, scale. =s)
  #Visualize eigenvalues (scree plot).
  #print(fviz_eig(pca.res))

  #Show the percentage of variances explained by each principal component.
  #print(fviz_contrib(pca.res, choice = "var", axes = 1:2, top = 10))

  #if mean.point = TRUE group means are added for each group, specified by coloring_factor
  #Graph of individuals. Individuals with a similar profile are grouped together.
  print(fviz_pca_ind(pca.res, col.ind = "cos2", geom = "point",
                     gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                     repel = TRUE # Avoid text overlapping (slow if many points)
  ))
}




#Data handling----
F30 <- read.sync(file="./F30_all_chr_mq20_bq20_Neda_filtered.sync", 
                 repl=rep(c(1:5), each=4), gen=rep(c(1:4),times=5))
F30_af <- as.data.frame(af(F30, repl=rep(c(1:5), each=4), gen=rep(c(1:4),times=5)))
F30_af <- F30_af[which(apply(F30_af,1,var)!=0),]
colnames(F30_af) = c("R1-S1.freq","R1-S2.freq","R1-S3.freq","R1-S4.freq","R1-S5.freq",
                     "R2-S1.freq","R2-S2.freq","R2-S3.freq","R2-S4.freq","R2-S5.freq",
                     "R3-S1.freq","R3-S2.freq","R3-S3.freq","R3-S4.freq","R3-S5.freq",
                     "R4-S1.freq","R4-S2.freq","R4-S3.freq","R4-S4.freq","R4-S5.freq")
#frequencies in dat must be transformed
F30_af = 2*asin(sqrt(F30_af))
#meta information
clmnames <- colnames(F30_af)
sampleIDs <- sub(".freq","",clmnames)
splits <- unlist(strsplit(sampleIDs, "[-]"))
repl <- as.factor(splits[grep("R",splits)])
subr <- as.factor(splits[grep("S",splits)])
#transpose data frame
t_F30_af <- as.data.frame(t(F30_af))
#set sampleIDs as rownames for labeling in the PCA
row.names(t_F30_af) <- sampleIDs
pdf("./PCA_rising_cons.pdf")
PCA(t_F30_af,TRUE,TRUE,repl)
#PCA(t_F30_af,TRUE,TRUE,subr)
dev.off()





# #Coverage----
# F30_cov <- coverage(F30, repl=rep(c(1:5), each=4), gen=rep(c(1:4),times=5))
# (meani <- apply(F30_cov,2,mean))
# hist(meani, breaks = 30, xlim = c(20,100))
#
# (meanl <- apply(F30_cov,1,mean))
# hist(meanl, breaks = 1000, xlim = c(0,100))
# meanl.good <- meanl[which(meanl > 20)]
# ## needs to be done here. Bad loci needs to be excluded from the analyses
# ## the name of meanl.good can be retrieved by names(meanl.good)
# ## in we can get all teh loci's name in F30 as F30@alleles$posID
# ## but I don;t know how to do this filtering
