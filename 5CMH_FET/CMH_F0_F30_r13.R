#!/usr/bin/env Rscript
#### loading packages####
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
# for r13
Repl <- rep(c(1:5), times=2)
Gen <- rep(c(0,30),each=5)
Sync.data <- read.sync(file="./cold_F30_r13_all_chr_mq20_bq20_Neda_filtered_base_duplicated_nozerocov.sync", 
                       repl=Repl, gen=Gen, polarization = "rising")
Allele_f <- as.data.frame(af(Sync.data, repl=Repl, gen=Gen))
Coverage <- as.data.frame(coverage(Sync.data, repl=Repl, gen=Gen))
Alleles <- as.data.frame(alleles(Sync.data))
row.names(Alleles) <- row.names(Coverage)
Allele_f <- na.omit(Allele_f)
Coverage <- Coverage[rownames(Allele_f),]

#### filter af and cov####
quan_divide <- function(x) quantile(x, probs = c(0.01, 0.02, 0.25, 0.50, 0.75, 0.98, 0.99))
largest_99per_cov <- max(apply(Coverage, 2, quan_divide))

lowest_cov <- apply(Coverage,2,function(x) quantile(x, probs = c(0.01)))
for (i in 1:ncol(Coverage)) {
  Coverage <- Coverage[Coverage[,i]> lowest_cov[i],]
}
Coverage <- Coverage[apply(Coverage, 1, function(x) !any(x>=largest_99per_cov)),]

Allele_f <- Allele_f[rownames(Coverage),]
Alleles <- Alleles[rownames(Coverage),]
cov_freq_r13 <- cbind(Allele_f, Coverage, Alleles)
rm(Sync.data, largest_99per_cov, lowest_cov,i, Allele_f, Alleles, Coverage, Gen, Repl, quan_divide)

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
cov_freq_auto_r13 <- cov_freq_r13[which(cov_freq_r13$chr != "X"),]
for (r in replicates){ 
  print(r)
  pref_i <- paste("F", times[1], ".R", r, sep = "") #time point i
  pref_j <- paste("F", times[2], ".R", r, sep = "") #time point j with j>i
  
  for(j in seq_len(nb_rounds)){
    cov_freq_trial <- sample_n(cov_freq_auto_r13, size = nb_SNPs)
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
rm(median_ne, ne_estimates, covi, covj, j)

####X chromosome----
ne_estimates <- NULL
cov_freq_x_r13 <- cov_freq_r13[which(cov_freq_r13$chr == "X"),]
for (r in replicates){ 
  print(r)
  pref_i <- paste("F", times[1], ".R", r, sep = "") #time point i
  pref_j <- paste("F", times[2], ".R", r, sep = "") #time point j with j>i
  
  for(j in seq_len(nb_rounds)){
    cov_freq_trial <- sample_n(cov_freq_x_r13, size = nb_SNPs)
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

# CMH----
#####autosome----
freq_auto_r13 <- as.matrix(cov_freq_auto_r13[,1:10])
cov_auto_r13 <-  as.matrix(cov_freq_auto_r13[,11:20])
ne_auto_r13 <- as.integer(median_auto[,6])
p_auto_r13_cmh <- ACER::adapted.cmh.test(freq=freq_auto_r13, coverage=cov_auto_r13, 
                                         Ne=ne_auto_r13, gen=c(0,30), repl = 1:5,
                                         IntGen = F, order = 0)
#####x chromosome----
freq_x_r13 <- as.matrix(cov_freq_x_r13[,1:10])
cov_x_r13 <-  as.matrix(cov_freq_x_r13[,11:20])
ne_x_r13 <- as.integer(median_x[,6])
p_x_r13_cmh <- ACER::adapted.cmh.test(freq=freq_x_r13, coverage=cov_x_r13, 
                                         Ne=ne_x_r13, gen=c(0,30), repl = 1:5,
                                         IntGen = F, order = 0)

res_r13_auto <- as.data.frame(cbind(CHR = cov_freq_auto_r13["chr"], BP = cov_freq_auto_r13["pos"], P = p.adjust(p_auto_r13_cmh, method = "fdr"), SNP = row.names(cov_freq_auto_r13)))
res_r13_x <- as.data.frame(cbind(CHR = cov_freq_x_r13["chr"], BP = cov_freq_x_r13["pos"], P = p.adjust(p_x_r13_cmh, method = "fdr"), SNP = row.names(cov_freq_x_r13)))
res_r13 <- rbind(res_r13_x, res_r13_auto)

res_r13 <- res_r13[which(res_r13$P > 1e-50),]
colnames(res_r13) <- c("CHR","BP","P","SNP")
res_r13$BP <- as.numeric(res_r13$BP); res_r13$P <- as.numeric(res_r13$P)
res_r13["CHR"][res_r13["CHR"] == "2L"] <- 2
res_r13["CHR"][res_r13["CHR"] == "2R"] <- 3
res_r13["CHR"][res_r13["CHR"] == "4"] <- 6
res_r13["CHR"][res_r13["CHR"] == "3L"] <- 4
res_r13["CHR"][res_r13["CHR"] == "3R"] <- 5
res_r13["CHR"][res_r13["CHR"] == "X"] <- 1
res_r13$CHR <- as.numeric(res_r13$CHR)

saveRDS(object = res_r13, file = "./res_r13")
# res_r13 <- readRDS("./res_r13")

thres_genome <- 1e-7
thres_suggest <- 5e-7

sig_SNP_cmh <- res_r13[which(res_r13$P < thres_genome), 4]
saveRDS(sig_SNP_cmh, file = "./sig_SNP_cmh_r13")
top_SNP_cmh <- res_r13[order(res_r13$P, decreasing = F),][1:100, 4]
saveRDS(top_SNP_cmh, file = "./top_SNP_cmh_r13")

# Plotting----
result <- res_r13
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) #%>%
  # Add highlight and annotation information
  #mutate( is_highlight=ifelse(SNP %in% top_SNP_cmh, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("./Plots/r13_plots/CMH_adapted_r13_F0-F30_fdr_adjusted_rising.pdf", width = 10, height = 6)
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted CMH test among 5 subreplicates in r13 contrasting F0 and F30")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_suggest), linetype="dashed", color = "blue") +
  #geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=1) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()











