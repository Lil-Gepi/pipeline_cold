res_r11_s1 <- readRDS(file = "./res_r11_s1")
res_r11_s1 <- res_r11_s1[which(res_r11_s1$P > 1e-6),]
top_s1 <- res_r11_s1[which(res_r11_s1$P < quantile(x = res_r11_s1$P, probs = 0.001)),]
row.names(top_s1) <- top_s1$SNP

res_r11_s2 <- readRDS(file = "./res_r11_s2")
res_r11_s2 <- res_r11_s2[which(res_r11_s2$P > 1e-6),]
top_s2 <- res_r11_s2[which(res_r11_s2$P < quantile(x = res_r11_s2$P, probs = 0.001)),]
row.names(top_s2) <- top_s2$SNP

res_r11_s3 <- readRDS(file = "./res_r11_s3")
res_r11_s3 <- res_r11_s3[which(res_r11_s3$P > 1e-6),]
top_s3 <- res_r11_s3[which(res_r11_s3$P < quantile(x = res_r11_s3$P, probs = 0.001)),]
row.names(top_s3) <- top_s3$SNP

res_r11_s4 <- readRDS(file = "./res_r11_s4")
res_r11_s4 <- res_r11_s4[which(res_r11_s4$P > 1e-6),]
top_s4 <- res_r11_s4[which(res_r11_s4$P < quantile(x = res_r11_s4$P, probs = 0.001)),]
row.names(top_s4) <- top_s4$SNP

res_r11_s5 <- readRDS(file = "./res_r11_s5")
res_r11_s5 <- res_r11_s5[which(res_r11_s5$P > 1e-6),]
top_s5 <- res_r11_s5[which(res_r11_s5$P < quantile(x = res_r11_s5$P, probs = 0.001)),]
row.names(top_s5) <- top_s5$SNP

top_across <- top_s1[rownames(top_s2),];top_across <- na.omit(top_across)
top_across <- top_across[rownames(top_s3),];top_across <- na.omit(top_across)
top_across <- top_across[rownames(top_s4),];top_across <- na.omit(top_across)
top_across <- top_across[rownames(top_s5),];top_across <- na.omit(top_across)
top_across <- top_across[,3]
saveRDS(top_across, file = "./top_SNPs_acorss_subreplicate")
