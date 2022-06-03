#to estimate Ne usting PoolSeq

NeEstimate <- function(dt.pol, 
                       dt.dp, 
                       whichChrom = "autosomes", 
                       nrSNPs = 100, 
                       nrIteration = 100, 
                       timepoints = c(1,60), 
                       poolSize = c(600,600), 
                       replicates = 1:6,
                       census = 1250) 
  {
  require(poolSeq)
  ne_estimates <- NULL
  if(whichChrom == "autosomes") {
    auto.dt.pol <- dt.pol[which(dt.pol$chr != "X"), ]
    auto.dt.pol$randomSamp  <- 1:nrow(auto.dt.pol)
    auto.dt.dp <- dt.dp[which(dt.dp$CHROM != "X"), ]
    auto.dt.dp$randomSamp  <- 1:nrow(auto.dt.dp)
    for (r in 1:length(replicates)){ 
      print(r)
      pref_i <- paste(timepoints[1], "_r", r, sep = "") #time point i
      pref_j <- paste(timepoints[2], "_r", r, sep = "") #time point j with j>i
      for(j in 1:nrIteration){
        randomSNPs <- sample(1:nrow(auto.dt.pol), size = nrSNPs)
        randomSNPs <- sort(randomSNPs, decreasing = F)
        freqi <- unlist(subset(auto.dt.pol, randomSamp %in% randomSNPs, select = paste(pref_i)), use.names = F)
        freqj <- unlist(subset(auto.dt.pol, randomSamp %in% randomSNPs, select = paste(pref_j)), use.names = F)
        covi <- unlist(subset(auto.dt.dp, randomSamp %in% randomSNPs, select = paste(pref_i)), use.names = F)
        covj <- unlist(subset(auto.dt.dp, randomSamp %in% randomSNPs, select = paste(pref_j)), use.names = F)
        ne <- estimateNe(p0 = freqi, pt = freqj, cov0 = covi, covt = covj, t = timepoints[2]-timepoints[1], 
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
    median_auto
  }
  else if(whichChrom == "X") {
    auto.dt.pol <- dt.pol[which(dt.pol$chr == "X"), ]
    auto.dt.pol$randomSamp  <- 1:nrow(auto.dt.pol)
    auto.dt.dp <- dt.dp[which(dt.dp$CHROM == "X"), ]
    auto.dt.dp$randomSamp  <- 1:nrow(auto.dt.dp)
    for (r in 1:length(replicates)){ 
      print(r)
      pref_i <- paste(timepoints[1], "_r", r, sep = "") #time point i
      pref_j <- paste(timepoints[2], "_r", r, sep = "") #time point j with j>i
      for(j in 1:length(nrIteration)){
        randomSNPs <- sample(1:nrow(auto.dt.pol), size = nrSNPs)
        randomSNPs <- sort(randomSNPs, decreasing = F)
        freqi <- unlist(subset(auto.dt.pol, randomSamp %in% randomSNPs, select = paste(pref_i)), use.names = F)
        freqj <- unlist(subset(auto.dt.pol, randomSamp %in% randomSNPs, select = paste(pref_j)), use.names = F)
        covi <- unlist(subset(auto.dt.dp, randomSamp %in% randomSNPs, select = paste(pref_i)), use.names = F)
        covj <- unlist(subset(auto.dt.dp, randomSamp %in% randomSNPs, select = paste(pref_j)), use.names = F)
        ne <- estimateNe(p0 = freqi, pt = freqj, cov0 = covi, covt = covj, t = timepoints[2]-timepoints[1], 
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
    median_auto <- data.frame(replicate = replicates, type = rep("autosomes", ), 
                              start = rep(times[1], nb_replicates), end = rep(times[2], nb_replicates), 
                              CI_95_down_median_ne = median_ne[, 1], median_ne = median_ne[, 2], CI_95_up_median_ne = median_ne[, 3])
    median_auto$replicate <- as.factor(median_auto$replicate)
    median_auto
  }
}
