setwd("~/Dropbox (PopGen)/Yiwen/RS/data/")
require(poolSeq)
require(parallel)
require(ACER)
library(viridis)
require(ggplot2)
library(factoextra)
library(reshape2)
library(data.table)
require(gridExtra)
require(ggrepel)
require(dplyr)
require(permute)

#Data handling----
####read r14----
Repl <- rep(c(1:5), times=2)
Gen <- rep(c(0,30),each=5)
Sync.data <- read.sync(file="./cold_F30_r14_all_chr_mq20_bq20_Neda_filtered_base_duplicated_nozerocov.sync", 
                       repl=Repl, gen=Gen, polarization = "rising",keepOnlyBiallelic = F)
Allele_f <- as.data.frame(af(Sync.data, repl=Repl, gen=Gen))
Coverage <- as.data.frame(coverage(Sync.data, repl=Repl, gen=Gen))
Alleles <- as.data.frame(alleles(Sync.data))
row.names(Alleles) <- row.names(Coverage)
Allele_f <- na.omit(Allele_f)
Coverage <- Coverage[rownames(Allele_f),]

#### filter af and cov----
quan_divide <- function(x) quantile(x, probs = c(0.01, 0.02, 0.25, 0.50, 0.75, 0.98, 0.99))
largest_99per_cov <- max(apply(Coverage, 2, quan_divide))

lowest_cov <- apply(Coverage,2,function(x) quantile(x, probs = c(0.01)))
for (i in 1:ncol(Coverage)) {
  Coverage <- Coverage[Coverage[,i]> lowest_cov[i],]
}
Coverage <- Coverage[apply(Coverage, 1, function(x) !any(x>=largest_99per_cov)),]

Allele_f <- Allele_f[rownames(Coverage),]
Alleles <- Alleles[rownames(Coverage),]
cov_freq <- cbind(Allele_f, Coverage, Alleles)

rm(Sync.data, largest_99per_cov, lowest_cov, Allele_f, Alleles, Coverage, Gen, Repl, quan_divide, quanti)
# NE estimate----

nb_SNPs <- 1000 #number of SNPs
nb_rounds <- 100 #number of sampling trials
times <- c(0, 30) #2 time points
census <- 600 #census size
poolSize <- rep(census, times=2)
replicates <- 1:5
ne_estimates <- NULL
nb_replicates <- 5
generations <- c(0, 30) #generations sequenced
tp <- length(generations) #number of time points sequenced

####autosomes----
cov_freq_auto <- cov_freq[which(cov_freq$chr != "X"),]
for (r in replicates){ 
  print(r)
  pref_i <- paste("F", times[1], ".R", r, sep = "") #time point i
  pref_j <- paste("F", times[2], ".R", r, sep = "") #time point j with j>i
  
  for(j in seq_len(nb_rounds)){
    cov_freq_trial <- sample_n(cov_freq_auto, size = nb_SNPs)
    pi <- unlist(subset(cov_freq_trial, select = paste(pref_i, ".freq", sep = "")))
    pj <- unlist(subset(cov_freq_trial, select = paste(pref_j, ".freq", sep = "")))
    covi <- unlist(subset(cov_freq_trial, select = paste(pref_i, ".cov", sep = "")))
    covj <- unlist(subset(cov_freq_trial, select = paste(pref_j, ".cov", sep = "")))
    ne <- estimateNe(p0 = pi, pt = pj, cov0 = covi, covt = covj, t = 30, 
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
rm(median_ne, ne_estimates, covi, covj, i, j)

####X chromosome----
ne_estimates <- NULL
cov_freq_x <- cov_freq[which(cov_freq$chr == "X"),]
for (r in replicates){ 
  print(r)
  pref_i <- paste("F", times[1], ".R", r, sep = "") #time point i
  pref_j <- paste("F", times[2], ".R", r, sep = "") #time point j with j>i
  
  for(j in seq_len(nb_rounds)){
    cov_freq_trial <- sample_n(cov_freq_x, size = nb_SNPs)
    pi <- unlist(subset(cov_freq_trial, select = paste(pref_i, ".freq", sep = "")))
    pj <- unlist(subset(cov_freq_trial, select = paste(pref_j, ".freq", sep = "")))
    covi <- unlist(subset(cov_freq_trial, select = paste(pref_i, ".cov", sep = "")))
    covj <- unlist(subset(cov_freq_trial, select = paste(pref_j, ".cov", sep = "")))
    ne <- estimateNe(p0 = pi, pt = pj, cov0 = covi, covt = covj, t = 30, 
                     ploidy = 2, truncAF = 0.05, method = "P.planI", poolSize = poolSize, Ncensus = census)
    ne_estimates <- rbind(ne_estimates, data.frame(replicate = r, start = times[1], end = times[2], trial = j,  ne = ne))
  } 
}

# Check the presence of negative / NA Ne values and take the median over the trials 
ne_estimates <- na.omit(ne_estimates)
ne <- ne_estimates$ne; if(length(which(ne<0))>0){ne_estimates <- ne_estimates[-which(ne<0)]}
ne_estimates$replicate <- as.factor(ne_estimates$replicate)

# Compute median and CI 95% of the median
median_ne <- t(sapply(replicates, function(x) {idx <- which(ne_estimates$replicate == x);
val <- sort(ne_estimates$ne[idx]);
med <- median(val); # median
up <- val[qbinom(1-0.025, length(idx), 0.5)+1]; # CI 95% of a median
down <- val[qbinom(0.025, length(idx), 0.5)];  # CI 95% of a median
round(c(down, med, up))}))

median_x <- data.frame(replicate = replicates, type = rep("X chromosome", nb_replicates), 
                     start = rep(times[1], nb_replicates), end = rep(times[2], nb_replicates), 
                     CI_95_down_median_ne = median_ne[, 1], median_ne = median_ne[, 2], CI_95_up_median_ne = median_ne[, 3])
median_x$replicate <- as.factor(median_x$replicate)
rm(census, covi, covj, generations, j, nb_replicates, nb_rounds, nb_SNPs, ne, pi, pj, poolSize,median_ne, ne_estimates, pref_i, pref_j,cov_freq_trial,
   r, replicates, times, tp)
#FET----
####r14s1----
#####autosome----
freq_auto_r14_s1 <- cov_freq_auto[,c(1,2)]
cov_auto_r14_s1 <-  cov_freq_auto[,c(11,12)]
ne_auto_r14_s1 <- median_auto[1,6]
p_auto_r14_s1<- ACER::adapted.chisq.test(freq=freq_auto_r14_s1, coverage=cov_auto_r14_s1, Ne=ne_auto_r14_s1, gen=c(0,30),poolSize = rep(600, 2), 
                              RetVal = 0, IntGen = FALSE)
#####x chromosome----
freq_x_r14_s1 <- cov_freq_x[,c(1,2)]
cov_x_r14_s1 <-  cov_freq_x[,c(11,12)]
ne_x_r14_s1 <- median_x[1,6]
p_x_r14_s1<- ACER::adapted.chisq.test(freq=freq_x_r14_s1, coverage=cov_x_r14_s1, Ne=ne_x_r14_s1, gen=c(0,30),poolSize = rep(600, 2), 
                                         RetVal = 0, IntGen = FALSE)


res_r14_s1_auto <- as.data.frame(cbind(CHR = cov_freq_auto[,21], BP = cov_freq_auto[,22], P = p.adjust(p_auto_r14_s1, method = "fdr"), SNP = row.names(cov_freq_auto)))
res_r14_s1_x <- as.data.frame(cbind(CHR = cov_freq_x[,21], BP = cov_freq_x[,22], P = p.adjust(p_x_r14_s1, method = "fdr"), SNP = row.names(cov_freq_x)))
res_r14_s1 <- rbind(res_r14_s1_x, res_r14_s1_auto)
res_r14_s1$BP <- as.numeric(res_r14_s1$BP)
res_r14_s1$P <- as.numeric(res_r14_s1$P)
####r14s2----
#####autosome----
freq_auto_r14_s2 <- cov_freq_auto[,c(3,4)]
cov_auto_r14_s2 <-  cov_freq_auto[,c(13,14)]
ne_auto_r14_s2 <- median_auto[2,6]
p_auto_r14_s2<- ACER::adapted.chisq.test(freq=freq_auto_r14_s2, coverage=cov_auto_r14_s2, Ne=ne_auto_r14_s2, gen=c(0,30),poolSize = rep(600, 2), 
                                         RetVal = 0, IntGen = FALSE)
#####x chromosome----
freq_x_r14_s2 <- cov_freq_x[,c(3,4)]
cov_x_r14_s2 <-  cov_freq_x[,c(13,14)]
ne_x_r14_s2 <- median_x[2,6]
p_x_r14_s2<- ACER::adapted.chisq.test(freq=freq_x_r14_s2, coverage=cov_x_r14_s2, Ne=ne_x_r14_s2, gen=c(0,30),poolSize = rep(600, 2), 
                                      RetVal = 0, IntGen = FALSE)


res_r14_s2_auto <- as.data.frame(cbind(CHR = cov_freq_auto[,21], BP = cov_freq_auto[,22], P = p.adjust(p_auto_r14_s2, method = "fdr"), SNP = row.names(cov_freq_auto)))
res_r14_s2_x <- as.data.frame(cbind(CHR = cov_freq_x[,21], BP = cov_freq_x[,22], P = p.adjust(p_x_r14_s2, method = "fdr"), SNP = row.names(cov_freq_x)))
res_r14_s2 <- rbind(res_r14_s2_x, res_r14_s2_auto)
res_r14_s2$BP <- as.numeric(res_r14_s2$BP)
res_r14_s2$P <- as.numeric(res_r14_s2$P)


####r14s3----
#####autosome----
freq_auto_r14_s3 <- cov_freq_auto[,c(5,6)]
cov_auto_r14_s3 <-  cov_freq_auto[,c(15,16)]
ne_auto_r14_s3 <- median_auto[3,6]
p_auto_r14_s3<- ACER::adapted.chisq.test(freq=freq_auto_r14_s3, coverage=cov_auto_r14_s3, Ne=ne_auto_r14_s3, gen=c(0,30),poolSize = rep(600, 2), 
                                         RetVal = 0, IntGen = FALSE)
#####x chromosome----
freq_x_r14_s3 <- cov_freq_x[,c(5,6)]
cov_x_r14_s3 <-  cov_freq_x[,c(15,16)]
ne_x_r14_s3 <- median_x[3,6]
p_x_r14_s3<- ACER::adapted.chisq.test(freq=freq_x_r14_s3, coverage=cov_x_r14_s3, Ne=ne_x_r14_s3, gen=c(0,30),poolSize = rep(600, 2), 
                                      RetVal = 0, IntGen = FALSE)


res_r14_s3_auto <- as.data.frame(cbind(CHR = cov_freq_auto[,21], BP = cov_freq_auto[,22], P = p.adjust(p_auto_r14_s3, method = "fdr"), SNP = row.names(cov_freq_auto)))
res_r14_s3_x <- as.data.frame(cbind(CHR = cov_freq_x[,21], BP = cov_freq_x[,22], P = p.adjust(p_x_r14_s3, method = "fdr"), SNP = row.names(cov_freq_x)))
res_r14_s3 <- rbind(res_r14_s3_x, res_r14_s3_auto)
res_r14_s3$BP <- as.numeric(res_r14_s3$BP)
res_r14_s3$P <- as.numeric(res_r14_s3$P)


####r14s4----
#####autosome----
freq_auto_r14_s4 <- cov_freq_auto[,c(7,8)]
cov_auto_r14_s4 <-  cov_freq_auto[,c(17,18)]
ne_auto_r14_s4 <- median_auto[4,6]
p_auto_r14_s4<- ACER::adapted.chisq.test(freq=freq_auto_r14_s4, coverage=cov_auto_r14_s4, Ne=ne_auto_r14_s4, gen=c(0,30),poolSize = rep(600, 2), 
                                         RetVal = 0, IntGen = FALSE)
#####x chromosome----
freq_x_r14_s4 <- cov_freq_x[,c(7,8)]
cov_x_r14_s4 <-  cov_freq_x[,c(17,18)]
ne_x_r14_s4 <- median_x[4,6]
p_x_r14_s4<- ACER::adapted.chisq.test(freq=freq_x_r14_s4, coverage=cov_x_r14_s4, Ne=ne_x_r14_s4, gen=c(0,30),poolSize = rep(600, 2), 
                                      RetVal = 0, IntGen = FALSE)


res_r14_s4_auto <- as.data.frame(cbind(CHR = cov_freq_auto[,21], BP = cov_freq_auto[,22], P = p.adjust(p_auto_r14_s4, method = "fdr"), SNP = row.names(cov_freq_auto)))
res_r14_s4_x <- as.data.frame(cbind(CHR = cov_freq_x[,21], BP = cov_freq_x[,22], P = p.adjust(p_x_r14_s4, method = "fdr"), SNP = row.names(cov_freq_x)))
res_r14_s4 <- rbind(res_r14_s4_x, res_r14_s4_auto)
res_r14_s4$BP <- as.numeric(res_r14_s4$BP)
res_r14_s4$P <- as.numeric(res_r14_s4$P)


####r14s5----
#####autosome----
freq_auto_r14_s5 <- cov_freq_auto[,c(9,10)]
cov_auto_r14_s5 <-  cov_freq_auto[,c(19,20)]
ne_auto_r14_s5 <- median_auto[5,6]
p_auto_r14_s5<- ACER::adapted.chisq.test(freq=freq_auto_r14_s5, coverage=cov_auto_r14_s5, Ne=ne_auto_r14_s5, gen=c(0,30),poolSize = rep(600, 2), 
                                         RetVal = 0, IntGen = FALSE)
#####x chromosome----
freq_x_r14_s5 <- cov_freq_x[,c(9,10)]
cov_x_r14_s5 <-  cov_freq_x[,c(19,20)]
ne_x_r14_s5 <- median_x[5,6]
p_x_r14_s5<- ACER::adapted.chisq.test(freq=freq_x_r14_s5, coverage=cov_x_r14_s5, Ne=ne_x_r14_s5, gen=c(0,30),poolSize = rep(600, 2), 
                                      RetVal = 0, IntGen = FALSE)


res_r14_s5_auto <- as.data.frame(cbind(CHR = cov_freq_auto[,21], BP = cov_freq_auto[,22], P = p.adjust(p_auto_r14_s5, method = "fdr"), SNP = row.names(cov_freq_auto)))
res_r14_s5_x <- as.data.frame(cbind(CHR = cov_freq_x[,21], BP = cov_freq_x[,22], P = p.adjust(p_x_r14_s5, method = "fdr"), SNP = row.names(cov_freq_x)))
res_r14_s5 <- rbind(res_r14_s5_x, res_r14_s5_auto)
res_r14_s5$BP <- as.numeric(res_r14_s5$BP)
res_r14_s5$P <- as.numeric(res_r14_s5$P)

####save results----
saveRDS(res_r14_s1, file = "./res_r14_s1")
saveRDS(res_r14_s2, file = "./res_r14_s2")
saveRDS(res_r14_s3, file = "./res_r14_s3")
saveRDS(res_r14_s4, file = "./res_r14_s4")
saveRDS(res_r14_s5, file = "./res_r14_s5")