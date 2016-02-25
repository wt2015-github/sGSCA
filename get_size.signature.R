#! /usr/bin/env Rscript
# Ting, 2013.5.7
# extract signature genes from outputs and get size and overlap information

setpair <- as.matrix(read.table('observe_GSCoL.txt',sep='\t'));
setpair <- setpair[,1];
iseq <- seq(length(setpair));
setpair <- setpair[which(iseq%%2==1)];
rm(iseq);

uv <- strsplit(readLines('obs_uv.txt'),'\t');
uv.new <- list();
for(i in 1:length(setpair)){
  u1 <- as.numeric(uv[[8*i-6]]);
  names(u1) <- uv[[(8*i-7)]];
  v1 <- as.numeric(uv[[8*i-4]]);
  names(v1) <- uv[[(8*i-5)]];
  u2 <- as.numeric(uv[[8*i-2]]);
  names(u2) <- uv[[(8*i-3)]];
  v2 <- as.numeric(uv[[8*i]]);
  names(v2) <- uv[[(8*i-1)]];
  uv.new[[setpair[i]]] <- list(u1=u1,v1=v1,u2=u2,v2=v2);
}

output <- c('pairid\toverlap.set\tsignature.set1\tsignature.set2');
for(i in 1:length(uv.new)){
  eg <- uv.new[[i]];
  u1 <- eg$u1;
  v1 <- eg$v1;
  u2 <- eg$u2;
  v2 <- eg$v2;
  
  u1<-u1[order(abs(u1),decreasing=T)];
  v1<-v1[order(abs(v1),decreasing=T)];
  u2<-u2[order(abs(u2),decreasing=T)];
  v2<-v2[order(abs(v2),decreasing=T)];
  
  u <- union(names(u1),names(u2));
  v <- union(names(v1),names(v2));
  
  gene.u <- union(names(u1)[u1!=0], names(u2)[u2!=0]);
  gene.v <- union(names(v1)[v1!=0], names(v2)[v2!=0]);
  
  overlap <- paste(length(u),length(v),length(intersect(u,v)),sep=',');
  u.sig <- paste(length(gene.u),'(',paste(gene.u,collapse=','),')',sep='');
  v.sig <- paste(length(gene.v),'(',paste(gene.v,collapse=','),')',sep='');
  tmp <- c(names(uv.new)[i], overlap, u.sig, v.sig);
  output <- c(output, paste(tmp, sep='\t', collapse='\t'));
}
writeLines(output, "output_size.signatures.txt");
