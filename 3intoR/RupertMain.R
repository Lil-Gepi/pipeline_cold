
setwd("/Volumes/Data/Dagny/Dropbox (PopGen)/DagnyPhD/Projects/PopulationSpecificAdaptation/scripts/")

samples <- read.table("../output/bcfListNames.txt")
head(samples)
sampSplit <- strsplit(samples$V2, split = "_")
gens <- c(0, 0, 0)
reps <- c("BasePT", "BaseFL", "BasePTFL")
for(i in 4:length(sampSplit)) {
  gens <- c(gens, substring(sampSplit[[i]][6], first = 2))
  reps <- c(reps, substring(sampSplit[[i]][7], first = 2))
  
}

#the polarizing function wants it in this format: "gen_rep" format, e.g. "10_r05"

mycols <- paste(paste(gens, "r", sep = "_"), reps, sep = "")
setwd("/Volumes/Temp2/preBcfHybridData/FLPT/")

require(data.table)
dt.af <- fread(cmd="bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%XF]\n' filt.reheader.FLPT23.flt.vcf",header=FALSE,col.names=mycols,key=c("CHROM","POS"))
dt.dp <- unique(fread(cmd="bcftools query -f '%CHROM\t%POS[\t%SAD]\n' filt.reheader.FLPT23.flt.vcf",header=FALSE,col.names=mycols2,key=c("CHROM","POS")))


# read allele frequencies (what=AF) or coverage depths (what=DP)
# from a given variant file using bcftools-query


setwd("/Volumes/Temp2/preBcfHybridData/FLPT/")

# input VCF file
vcf <- "filt.reheader.FLPT23.flt.vcf"

# full SNP set
dt.af <- read_bcf(vcf,what="AF")
head(dt.af)
dt.dp <- read_bcf(vcf,what="DP")
head(dt.dp)

dt.af.multi <- prune_multi(dt.af,dt.dp,how="most_freq")

dt.af.base <- dt.af[,c(1:8)] #extracting bases 
head(dt.af.base)
colnames(dt.af.base) <- c(colnames(dt.af.base[,1:5]), c("BasePT", "BaseFL", "BasePTFL"))

dt.af <- dt.af[,-c(6:8)]
head(dt.af)
setcolorder(x = dt.af, neworder = c(colnames(dt.af)[1:5], sort(colnames(dt.af)[6:47])))


dt.dp.base <- dt.dp[,c(1:5)] #extracting bases
head(dt.dp.base)
colnames(dt.dp.base) <- c(colnames(dt.dp.base[,1:2]), c("BasePT", "BaseFL", "BasePTFL"))
head(dt.dp.base)

dt.dp <- dt.dp[,-c(3:5)]
head(dt.dp)
setcolorder(x = (dt.dp), neworder = c(colnames((dt.dp))[1:2], sort(colnames((dt.dp))[3:42])))


setwd("/Volumes/Data/Dagny/Dropbox (PopGen)/DagnyPhD/Projects/PopulationSpecificAdaptation/scripts/")


source("RupertCalcPACER.R")
source("RupertPolarizing.R")
source("RupertReadFile.R")
source("RupertRemoveMulti.R")




# prune non-dominant ALTs by freq

# polarize
dt.pol <- polarize(dt.af,direct=TRUE)

# ne estimate
ne.est <- NeEstimate(dt.pol, dt.dp)
ne.est$median_ne <- as.integer(ne.est$median_ne)

# calculate p-values
dt.pval <- calc_pvals(dt.pol,dt.dp,
                      par.rep = c(1:6), 
                      par.gen = c(1,10,20,30,40,50,60), 
                      par.Ne = ne.est$median_ne,
                      par.poolSize = rep(625,42))

# pca


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




