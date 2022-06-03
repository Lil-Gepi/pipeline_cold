

# polarize for rising alleles
polarize <- function(dt.af,drop.af.cols=5L,direct=TRUE) {
  if (direct) {
    af.mat <- t(as.matrix(dt.af[,-(1:drop.af.cols)]))
    gens <- as.integer(sub("_.*$","",rownames(af.mat)))
    mlm.obj <- lm(af.mat~gens)
    # keep SNPs with a rising trend
    keep <- coef(mlm.obj)["gens",]>0
    setkey(rbind(
      cbind(dt.af[keep,.(chr=CHROM,pos=POS,ref=REF,riseallele=ALT,fallallele=REF)],
            dt.af[keep,-(1:drop.af.cols)]),
      cbind(dt.af[!keep,.(chr=CHROM,pos=POS,ref=REF,riseallele=REF,fallallele=ALT)],
            as.data.table(t(1-af.mat[,!keep])))),chr,pos,riseallele)
  } else {
    require(foreach)
    stopifnot(all(names(dt.af)[1:5]==c("CHROM","POS","REF","ALT","N_ALT")))
    stopifnot(drop.af.cols==5L)
    # REF frequencies
    dt.ref <- cbind(setkey(unique(dt.af[,.(CHROM,POS,REF,ALT=REF,N_ALT)]),N_ALT),
                    foreach (i=1:3,.combine=rbind) %do% {
                      which <- dt.af[,N_ALT==i]
                      as.data.table(1-colSums(array(unname(unlist(dt.af[which,-(1:drop.af.cols)])),
                                                    dim=c(i,sum(which)/i,ncol(dt.af)-5),
                                                    dimnames=list(NULL, # paste0("ALT",seq_len(i))
                                                                  NULL, # dt.af[which,paste0(CHROM,":",POS)]
                                                                  names(dt.af)[-(1:drop.af.cols)])),
                                              dims=1))
                    })
    # combine with dt.af
    dt.pol <- rbind(dt.af,dt.ref)
    # linear model
    af.mat <- t(as.matrix(dt.pol[,-(1:drop.af.cols)]))
    gens <- as.integer(sub("_.*$","",rownames(af.mat)))
    mlm.obj <- lm(af.mat~gens)
    # keep SNPs with a rising trend
    keep <- coef(mlm.obj)["gens",]>0
    setkey(cbind(dt.pol[keep,.(chr=CHROM,pos=POS,ref=REF,riseallele=ALT,fallallele="N")],
                 dt.pol[keep,-(1:drop.af.cols)]),
           chr,pos,riseallele)
  }
}
