setwd("~/Dropbox (PopGen)/Yiwen/RS/data/")
require(ggplot2)

#manhattan plot----
sig_SNP_cmh <- readRDS(file = "./sig_SNP_cmh_r12")
thres_genome <- 0.05
thres_suggest <- 0.1
####r12s1----
res_r12_s1 <- readRDS(file = "./res_r12_s1")
res_r12_s1 <- res_r12_s1[which(res_r12_s1$P > 1e-6),]

result <- res_r12_s1
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_SNP_cmh & P < thres_genome, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("./Plots/r12_plots/FET_adapted_r12_s1_F0-F30_fdr_adjusted_rising_cmh_highlight.pdf", width = 10, height = 6, compress = F )
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted chi-square test in r12_s1 contrasting F0 and F30")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=1) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()
####r12s2----
res_r12_s2 <- readRDS(file = "./res_r12_s2")
res_r12_s2 <- res_r12_s2[which(res_r12_s2$P > 1e-6),]

result <- res_r12_s2
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_SNP_cmh & P < thres_genome, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("FET_adapted_r12_s2_F0-F30_fdr_adjusted_rising_cmh_highlight.pdf", width = 10, height = 6, compress = F )
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted chi-square test in r12_s2 contrasting F0 and F30")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1 ) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=1 ) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()
####r12s3----
res_r12_s3 <- readRDS(file = "./res_r12_s3")
res_r12_s3 <- res_r12_s3[which(res_r12_s3$P > 1e-6),]

result <- res_r12_s3
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_SNP_cmh & P < thres_genome, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("FET_adapted_r12_s3_F0-F30_fdr_adjusted_rising_cmh_highlight.pdf", width = 10, height = 6, compress = F )
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted chi-square test in r12_s3 contrasting F0 and F30")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1 ) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=1 ) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()

####r12s4----
res_r12_s4 <- readRDS(file = "./res_r12_s4")
res_r12_s4 <- res_r12_s4[which(res_r12_s4$P > 1e-6),]

result <- res_r12_s4
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_SNP_cmh & P < thres_genome, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("FET_adapted_r12_s4_F0-F30_fdr_adjusted_rising_cmh_highlight.pdf", width = 10, height = 6, compress = F )
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted chi-square test in r12_s4 contrasting F0 and F30")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1 ) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=1 ) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()
####r12s5----
res_r12_s5 <- readRDS(file = "./res_r12_s5")
res_r12_s5 <- res_r12_s5[which(res_r12_s5$P > 1e-6),]

result <- res_r12_s5
don <- result %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  mutate(tot=cumsum(chr_len)-chr_len) %>%
  select(-chr_len) %>%
  left_join(result, ., by=c("CHR"="CHR")) %>%
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  # Add highlight and annotation information
  mutate( is_highlight=ifelse(SNP %in% sig_SNP_cmh & P < thres_genome, "yes", "no"))

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
pdf("FET_adapted_r12_s5_F0-F30_fdr_adjusted_rising_cmh_highlight.pdf", width = 10, height = 6, compress = F )
ggplot(don, aes(x=BPcum, y=-log10(P))) +
  ggtitle("Adapted chi-square test in r12_s5 contrasting F0 and F30")+
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1 ) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  scale_x_continuous(name = "Chromosome", label = c("X", "2L","2R","3L", "3R", "4"), breaks= axisdf$center ) +
  scale_y_continuous(name = "-log10(q)", expand = c(0, 0) ) +     
  geom_hline(yintercept=-log10(thres_genome), linetype="dashed", color = "red") +
  geom_hline(yintercept=-log10(thres_suggest), linetype="dashed", color = "blue") +
  geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=1 ) +
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )
dev.off()
