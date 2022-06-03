
# calculate p-values using ACER
calc_pvals <- function(dt.af,dt.dp,
                       drop.af.cols=5L,
                       drop.dp.cols=2L,
                       par.rep=1L:10L,
                       par.gen=10L*0L:6L,
                       par.Ne=rep(500L,60L),
                       par.poolSize=rep(852L,60L),
                       correct.multi=TRUE) {
  require(ACER)
  stopifnot(all(key(dt.dp)==c("CHROM","POS")))
  af.mat <- as.matrix(dt.af[,-(1:drop.af.cols)])
  dp.mat <- as.matrix(dt.dp[dt.af[,1:2]][,-(1:drop.dp.cols)])
  # ACER does not like 0 coverage
  dp.mat[dp.mat==0L] <- 1L
  dt.pval <- dt.af[,1:2][
    ,PVAL:=ACER::adapted.cmh.test(freq=af.mat,
                                  coverage=dp.mat,
                                  Ne=par.Ne,
                                  gen=par.gen,
                                  repl=par.rep,
                                  poolSize=par.poolSize,
                                  IntGen=TRUE,
                                  order=1)]
  names(dt.pval)[1:2] <- c("CHROM","POS")
  # for multi-allelic sites
  if (correct.multi) {
    dt.pval[,c("N","I"):=.(.N,seq_len(.N)),by=.(CHROM,POS)]
    dt.pval[N>1,PVAL:=p.adjust(PVAL,method="fdr"),by=.(CHROM,POS)]
    dt.pval[,N:=NULL]
  }
  # multiple testing across chromosome
  aux <- dt.pval[,.(N=length(unique(POS))),keyby=CHROM]
  dt.pval[,PVAL:=p.adjust(PVAL,method="fdr",n=aux[CHROM,N]),by=.(CHROM,I)][
    ,I:=NULL]
}
