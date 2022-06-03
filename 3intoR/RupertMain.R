
setwd("~/Dropbox (PopGen)/Yiwen/RS/pipeline/3intoR")
source("RupertCalcPACER.R")
source("NeEstimate.R")
source("RupertPolarizing.R")
source("RupertReadFile.R")
source("RupertRemoveMulti.R")
library(matrixStats)
require(data.table)

#the polarizing function wants it in this format: "gen_rep" format, e.g. "10_r05"
gens <- rep(c(0,10,20,30,40,50,60), each=10)
reps <- rep(c(1:10), times=7)

mycols <- paste(paste(gens, "r", sep = "_"), reps, sep = "")
# read allele frequencies (what=AF) or coverage depths (what=DP)
# from a given variant file using bcftools-query
# input VCF file
vcf <- "/Users/ychen/COLD/result/vcf/cold_piled_MAX_flt.010.vcf.gz"
# full SNP set
dt.af <- read_bcf(vcf,what="AF")
head(dt.af)
dt.dp <- read_bcf(vcf,what="DP")
head(dt.dp)
colnames(dt.af) <- c(colnames(dt.af[,1:5]), mycols)
colnames(dt.dp) <- c(colnames(dt.dp[,1:2]), mycols)
# prune non-dominant ALTs by freq
dt.af.multi <- prune_multi(dt.af,dt.dp,how="most_freq")
# polarize
dt.pol <- polarize(dt.af.multi,direct=TRUE)
# 
# colnames(dt.pol) <- c(colnames(dt.pol[,1:5]), mycols)
# ne estimate
# ne.est <- NeEstimate(dt.pol = dt.pol, dt.dp = dt.dp, whichChrom = "autosomes", nrSNPs = 1000, nrIteration = 100, timepoints = c(0,60), replicates = 1:10, poolSize = c(600,600), census = 1250)
# ne.est$median_ne <- as.integer(ne.est$median_ne)

nb_SNPs <- 1000 #number of SNPs
nb_rounds <- 100 #number of sampling trials
times <- c(0, 60) #2 time points
census <- 600 #census size
poolSize <- rep(census, times=2)
replicates <- 1:10
ne_estimates <- NULL
nb_replicates <- 10
generations <- c(0, 60) #generations sequenced
tp <- length(generations) #number of time points sequenced
library(dplyr)
####autosomes----
dt.dp.auto <- dt.dp[which(dt.dp$CHROM != "X"),]
dt.af.auto <- dt.pol[which(dt.pol$chr != "X"),]
for (r in replicates){ 
  print(r)
  pref_i <- paste( times[1], "_r", r, sep = "") #time point i
  pref_j <- paste(times[2], "_r", r, sep = "") #time point j with j>i
  
  for(j in seq_len(nb_rounds)){
    cov_trial <- sample_n(dt.dp.auto, size = nb_SNPs)
    af_trial <- sample_n(dt.af.auto, size = nb_SNPs)
    pi <- unlist(subset(af_trial, select = pref_i))
    pj <- unlist(subset(af_trial, select = pref_j))
    covi <- unlist(subset(cov_trial, select = pref_i))
    covj <- unlist(subset(cov_trial, select = pref_j))
    ne <- estimateNe(p0 = pi, pt = pj, cov0 = covi, covt = covj, t = 60, 
                     ploidy = 2, truncAF = 0.05, method = "P.planI", poolSize = poolSize, Ncensus = census)
    ne_estimates <- rbind(ne_estimates, data.frame(replicate = r, start = times[1], end = times[2], trial = j,  ne = ne))
  } 
}

ne_estimates <- na.omit(ne_estimates)
ne <- ne_estimates$ne; if(length(which(ne<0))>0){ne_estimates <- ne_estimates[-which(ne<0)]}
ne_estimates$replicate <- as.factor(ne_estimates$replicate)
median_ne <- t(sapply(replicates, function(x) {idx <- which(ne_estimates$replicate == x);
val <- sort(ne_estimates$ne[idx]);
med <- median(val); # median
up <- val[qbinom(1-0.025, length(idx), 0.5)+1]; # CI 95% of a median
down <- val[qbinom(0.025, length(idx), 0.5)];  # CI 95% of a median
round(c(down, med, up))}))
median_auto <- data.frame(replicate = replicates, type = rep("autosomes", nb_replicates), 
                          start = rep(times[1], nb_replicates), end = rep(times[2], nb_replicates), 
                          CI_95_down_median_ne = median_ne[, 1], median_ne = median_ne[, 2], CI_95_up_median_ne = median_ne[, 3])
median_auto$replicate <- as.factor(median_auto$replicate)
rm(median_ne, ne_estimates, covi, covj, j)

####X chromosome----
ne_estimates <- NULL
dt.dp.X <- dt.dp[which(dt.dp$CHROM == "X"),]
dt.af.X <- dt.pol[which(dt.pol$chr == "X"),]
for (r in replicates){ 
  print(r)
  pref_i <- paste( times[1], "_r", r, sep = "") #time point i
  pref_j <- paste(times[2], "_r", r, sep = "") #time point j with j>i
  
  for(j in seq_len(nb_rounds)){
    cov_trial <- sample_n(dt.dp.X, size = nb_SNPs)
    af_trial <- sample_n(dt.af.X, size = nb_SNPs)
    pi <- unlist(subset(af_trial, select = pref_i))
    pj <- unlist(subset(af_trial, select = pref_j))
    covi <- unlist(subset(cov_trial, select = pref_i))
    covj <- unlist(subset(cov_trial, select = pref_j))
    ne <- estimateNe(p0 = pi, pt = pj, cov0 = covi, covt = covj, t = 60, 
                     ploidy = 2, truncAF = 0.05, method = "P.planI", poolSize = poolSize, Ncensus = census)
    ne_estimates <- rbind(ne_estimates, data.frame(replicate = r, start = times[1], end = times[2], trial = j,  ne = ne))
  } 
}

ne_estimates <- na.omit(ne_estimates)
ne <- ne_estimates$ne; if(length(which(ne<0))>0){ne_estimates <- ne_estimates[-which(ne<0)]}
ne_estimates$replicate <- as.factor(ne_estimates$replicate)
median_ne <- t(sapply(replicates, function(x) {idx <- which(ne_estimates$replicate == x);
val <- sort(ne_estimates$ne[idx]);
med <- median(val); # median
up <- val[qbinom(1-0.025, length(idx), 0.5)+1]; # CI 95% of a median
down <- val[qbinom(0.025, length(idx), 0.5)];  # CI 95% of a median
round(c(down, med, up))}))
median_X <- data.frame(replicate = replicates, type = rep("X", nb_replicates), 
                          start = rep(times[1], nb_replicates), end = rep(times[2], nb_replicates), 
                          CI_95_down_median_ne = median_ne[, 1], median_ne = median_ne[, 2], CI_95_up_median_ne = median_ne[, 3])
median_X$replicate <- as.factor(median_X$replicate)
rm(census, covi, covj, generations, j, nb_replicates, nb_rounds, nb_SNPs, ne, pi, pj, poolSize,median_ne, ne_estimates, pref_i, pref_j,cov_freq_trial,
   r, replicates, times, tp)















# calculate p-values
dt.pval <- calc_pvals(dt.pol,dt.dp,
                      par.rep = c(1:10), 
                      par.gen = c(0,10,20,30,40,50,60), 
                      par.Ne = median_auto$median_ne,
                      par.poolSize = rep(625,70))
colnames(dt.pval) <- c("CHR", "BP", "P")


dt.pval[which(dt.pval[,"CHR"] == "2L"), "CHR"]<- 2
dt.pval[which(dt.pval[,"CHR"] == "2R"), "CHR"]<- 3
dt.pval[which(dt.pval[,"CHR"] == "4"), "CHR"]<- 6
dt.pval[which(dt.pval[,"CHR"] == "3L"), "CHR"]<- 4
dt.pval[which(dt.pval[,"CHR"] == "3R"), "CHR"]<- 5
dt.pval[which(dt.pval[,"CHR"] == "X"), "CHR"]<- 1
dt.pval$CHR <- as.numeric(dt.pval$CHR)

saveRDS(dt.pval, "cold_cmh_result")
saveRDS(dt.af, "cold_af")
saveRDS(dt.dp, "cold_dp")

require(ggplot2)
thres_cmh_genome <- 5e-8
thres_cmh_suggest <- 5e-6
result <- dt.pval
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) 

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
tiff("./cold_CMH_adapted_0-60.tiff", width = 24, height = 16, res = 900, units = "cm")
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("genomic signatures in cold adaptation contrasting F0 and F60")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_cmh_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_cmh_suggest), linetype="dashed", color = "blue") +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

# pca
library(poolSeq)
library(ggplot2)
library(factoextra)
library(reshape2)
library(data.table)
PCA <- function(df,c,s,coloring_factor)
{
  pcadata<-na.omit(df)#remove missing data
  pca.res <- prcomp(pcadata, retx=TRUE, center = c, scale. =s)
  print(fviz_pca_ind(X = pca.res, 
                     geom = "point",
                     col.ind=coloring_factor,
                     palette="aaas",
                     repel = TRUE,   # Avoid text overlapping
                     mean.point = FALSE,
                     addEllipses=FALSE)) # Default parameters for ellipses
}
SNP <- paste(dt.af.multi$CHROM, dt.af.multi$POS, sep = ".")
my.af <- dt.af.multi[,6:ncol(dt.af.multi)]
rownames(my.af) <- SNP
my.af <- 2*asin(sqrt(my.af))

t_my.af <- as.data.frame(t(my.af))
#set sampleIDs as rownames for labeling in the PCA
row.names(t_my.af) <- mycols
t_my.af <- na.omit(t_my.af)#remove missing data

tiff("cold_PCA_arcsined.tiff", width = 24, height = 24, res = 450, units = "cm")
PCA(t_my.af,TRUE,TRUE,coloring_factor = reps)
PCA(t_my.af,TRUE,TRUE,gens)
dev.off()

# manhattan
dt.pval$SNP <- c(1:nrow(dt.pval)) #go to another script


#enrichment of "population specific SNP's"
#threshold of -log10(0.0001) = 4

# prune non-dominant ALTs by freq
dt.af.bases.2 <- prune_multi(dt.af.base,dt.dp.base,how="most_freq")
dt.af.base[1:50, ]
dt.pval00001 <- dt.pval[which(dt.pval$PVAL < 0.00005), ]
nrow(dt.pval00001)

backgroundCHROMPOS <- paste0(dt.pval$CHROM, dt.pval$POS)

FLSpesific <- dt.af.base[which(dt.af.base$BasePT < 0.05 & dt.af.base$BaseFL > 0.05) , ]
FLSpesificCHROMPOS <- paste0(FLSpesific$CHROM, FLSpesific$POS)
length(FLSpesificCHROMPOS)

PTSpesific <- dt.af.base[which(dt.af.base$BaseFL < 0.05 & dt.af.base$BasePT > 0.05), ]
PTSpesificCHROMPOS <- paste0(PTSpesific$CHROM, PTSpesific$POS)
length(PTSpesificCHROMPOS)

dt.pval00001CHROMPOS <- paste0(dt.pval00001$CHROM, dt.pval00001$POS)

cont_table=function(query,background,classifyer){
  p1=length(Reduce(intersect,list(query,background,classifyer)))
  q1=length(intersect(query,background))-p1
  p0=length(setdiff(intersect(background,classifyer),intersect(query,classifyer)))
  q0=length(setdiff(background,query))-p0
  return(matrix(c(p1,p0,q1,q0),2,2))
}

fisher.test(cont_table(dt.pval00001CHROMPOS,backgroundCHROMPOS,FLSpesificCHROMPOS),alternative = "greater") 
fisher.test(cont_table(dt.pval00001CHROMPOS,backgroundCHROMPOS,PTSpesificCHROMPOS),alternative = "greater") 




# dt.af <- fread(cmd="bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%XF]\n' filt.reheader.FLPT23.flt.vcf",header=FALSE,col.names=mycols,key=c("CHROM","POS"))
# dt.dp <- unique(fread(cmd="bcftools query -f '%CHROM\t%POS[\t%SAD]\n' filt.reheader.FLPT23.flt.vcf",header=FALSE,col.names=mycols2,key=c("CHROM","POS")))
#
# dt.af.base <- dt.af[,c(1:15)] #extracting bases 
# head(dt.af.base)
# colnames(dt.af.base) <- c(colnames(dt.af.base[,1:5]), c("BasePT", "BaseFL", "BasePTFL"))
# dt.af <- dt.af[,-c(6:8)]
# head(dt.af)
# setcolorder(x = dt.af, neworder = c(colnames(dt.af)[1:5], sort(colnames(dt.af)[6:47])))
# dt.dp.base <- dt.dp[,c(1:12)] #extracting bases
# head(dt.dp.base)
# colnames(dt.dp.base) <- c(colnames(dt.dp.base[,1:2]), c("BasePT", "BaseFL", "BasePTFL"))
# head(dt.dp.base)
# dt.dp <- dt.dp[,-c(3:5)]
# head(dt.dp)
# setcolorder(x = (dt.dp), neworder = c(colnames((dt.dp))[1:2], sort(colnames((dt.dp))[3:42])))
