#! /usr/bin/env Rscript
# Ting, 2012.7.6
# get statistic significance from outputs, the same with the function "post" in the main script, used for integrating parallel null_GSCoL results in the end

re <- as.matrix(read.table('null_GSCoL.txt',sep='\t'));
real.tmp <- as.matrix(read.table("observe_GSCoL.txt"));
obs <- as.matrix(as.numeric(real.tmp[,2]));
rownames(obs)<-real.tmp[,1]
pairname.tmp <- rownames(obs);
real <- obs;
perm <- ncol(re);
n <- nrow(re)/2;

#normalize GSCoL
cat('normalize GSCoL\n');
data.mean <- apply(re, 1, mean);
data.sd <- apply(re, 1, sd);
re.data.norm <- (re - data.mean)/data.sd;
real.data.norm <- (real - data.mean)/data.sd;

i.seq <- seq(nrow(re));
#real.mean <- (real[which(i.seq%%2==0),] + real[which(i.seq%%2==1),])/2;
re.score <- (re.data.norm[which(i.seq%%2==0),] + re.data.norm[which(i.seq%%2==1),])/2;
real.score <- (real.data.norm[which(i.seq%%2==0)] + real.data.norm[which(i.seq%%2==1)])/2;
pairname <- pairname.tmp[which(i.seq%%2==0)];

# pvalue
pvalue <- vector(length=n, mode="numeric");
cat("calculate pvalue.\n");
for(i in 1:n){
  pvalue[i] <- sum( re.score[i,] >= real.score[i] ) / perm;
}

cat("calculate qvalue with fdrtool.\n");
library("fdrtool");
qvalue.fdrtool <- fdrtool(pvalue, statistic="pvalue",plot=FALSE)$qval;

# p.adjust
cat("calculate FDR with p.adjust.\n");
fdr.adjust <- p.adjust(pvalue, method="fdr");

pvalue_fdr <- paste(pairname, real.score, pvalue, qvalue.fdrtool, fdr.adjust, sep="\t");
pvalue_fdr <- c(as.character("pairid\tGSCoL\tpvalue\tqvalue.fdrtool\tfdr.padjust"), pvalue_fdr);
writeLines(pvalue_fdr, "output_pvalue_fdr.txt");
