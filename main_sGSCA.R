#! /usr/bin/env Rscript
#Ting, 2013.4.5
#use sCCA to calculate pathway-pathway correlation with the fold-change matrix, use gene permutation to get null distributions, output results and weight vectors u,v

# inputs:
#	perm: number of permutation, set # according to parallel computation
#	fc.file: fold-change matrix (or other expression profile), more samples are better
#	geneset.file: gene set matrix, row is gene, column is geneset, 0/1 represent containing relationship
#	LM.filter: TRUE or FALSE
#	LMgene.file: gene filtering file, use gene names, so at least 2 columns (gene name colume and gene score column) with column names
#	ppi.file: ppi file with 3 columns, the third column is source
# outputs:
#	geneset_filtered.txt: geneset after filtering with expression data 
#	geneset_pair.txt: candidate geneset pairs
#	observe_GSCoL.txt: observed/foreground GSCoL
#	obs_uv.txt: u,v vectors from sCCA
#	null_GSCoL.txt: null/background GSCoL
#	output_pvalue_fdr.txt: results


main <- function(perm=10, fc.file="./data_exp.gse22058/fc_matrix.txt", geneset.file="./data_exp.gse22058/geneset.new.txt", LM.filter=TRUE, LMgene.file="./data_exp.gse22058/ttest_genes.txt", ppi.file="./data_exp.gse22058/incorporated_gene_relations.txt"){
	library(plyr);
	library(PMA);
  
	# read fold-change matrix file
	cat('load expression data\n');
	fc <- as.matrix(read.table(fc.file, header=TRUE, sep="\t", row.names=1));
  
	# read geneset file, construct RPN
	GS <- init.GeneSetPair(LM.filter=LM.filter, fc.gene=row.names(fc), LMgene.file=LMgene.file, geneset.file=geneset.file, ppi.file=ppi.file);
	geneset <- GS[[1]];
	setpair <- GS[[2]];
  
	# observed GSCoL
	obs <- observe.GSCoL(fc, setpair, geneset);
  
	# permutation for null GSCoL distribution
	GSCoL.permd <- perm.GSCoL(fc, setpair, geneset, perm);
  
	# result
	post(obs, GSCoL.permd);
}


init.GeneSetPair <- function(LM.filter, fc.gene, LMgene.file, geneset.file, ppi.file){
	ppi <- matrix(unlist (strsplit (readLines(ppi.file), "\t")), ncol=3, byrow=TRUE);
	ppi.tmp <- ppi[,c(2,1,3)];
	ppi <- rbind(ppi, ppi.tmp);
	ppi <- paste(ppi[,1], ppi[,2], sep="\t");
  
	genes <- fc.gene;
	if(LM.filter){
		cat('filter genes with predefined list\n');
		LMgenes <- read.table(LMgene.file, header=TRUE, sep="\t", row.names=1);
		genes <- intersect(genes, row.names(LMgenes));
	}
	
	geneset <- as.matrix(read.table(geneset.file, header=TRUE, sep='\t',row.names=1));
	set <- intersect(genes, rownames(geneset));
	geneset.filted <- geneset[set,];
	write.table(geneset.filted, 'geneset_filtered.txt',sep='\t');
	
	cat('get candidate geneset pairs\n');
	set.pair <- c();
	for(i in 1:(ncol(geneset.filted)-1)){
		set1 <- rownames(geneset.filted)[geneset.filted[,i] == 1];
		for(j in (i+1):ncol(geneset.filted)){
			set2 <- rownames(geneset.filted)[geneset.filted[,j] ==1 ];
			if(contained(set1, set2)){
				next;
			}
			relation.set12 <- paste(set1, rep(set2, each=length(set1)), sep="\t");
			if(length(intersect(set1, set2)) >0 || length(intersect(relation.set12, ppi)) > 0){ # filter geneset pairs with ppi
				set.pair <- rbind(set.pair, c(colnames(geneset.filted)[i], colnames(geneset.filted)[j]));
			}
		}
	}
	write.table(set.pair, 'geneset_pair.txt',sep='\t', row.names=FALSE, col.names=FALSE);
	
	return ( list(geneset.filted, set.pair) );
}

# inclusion relation of two sets
contained<-function(a,b){
	if(length(setdiff(a,b)) < 2  || length(setdiff(b,a)) < 2) return(TRUE) else return(FALSE);
}

observe.GSCoL<-function(fc, setpair, geneset){
	cat("get observation score.\n");
	obs <- c();
	#setpair.notNA <- c();
	uv <- c();
	for (i in 1:nrow(setpair)){
		cat('setpair',i,'\n');
		set1 <- rownames(geneset)[geneset[,setpair[i,1]] == 1];
		set2 <- rownames(geneset)[geneset[,setpair[i,2]] == 1];
		tmp <- GSCoL(fc, set1, set2); # every two continuous output from this function are for one setpair because of the potential overlap
		
		obs <- c(obs, tmp[[1]]$GSCoL, tmp[[2]]$GSCoL);
		#setpair.notNA <- rbind(setpair.notNA,setpair[j,]);
		
		#output the weight vectors u,v, note every 8 lines are for one setpair because of the overlap
		uv <- c(uv, paste(tmp[[1]]$set1, sep="", collapse="\t"));
		uv <- c(uv, paste(tmp[[1]]$u, sep="", collapse="\t"));
		uv <- c(uv, paste(tmp[[1]]$set2, sep="", collapse="\t"));
		uv <- c(uv, paste(tmp[[1]]$v, sep="", collapse="\t"));
		uv <- c(uv, paste(tmp[[2]]$set1, sep="", collapse="\t"));
		uv <- c(uv, paste(tmp[[2]]$u, sep="", collapse="\t"));
		uv <- c(uv, paste(tmp[[2]]$set2, sep="", collapse="\t"));
		uv <- c(uv, paste(tmp[[2]]$v, sep="", collapse="\t"));
	}
	pairname <- rep (paste(setpair[,1], setpair[,2], sep="---"), each=2);
	obs <- matrix(obs, ncol=1, dimnames=list(pairname, NULL));
	write.table(obs, "observe_GSCoL.txt",row.names = TRUE, col.names = FALSE, sep="\t");
	writeLines(uv, "obs_uv.txt");
	return(obs);
}

perm.GSCoL <- function(fc, setpair, geneset, perm){
	N <- nrow(fc);
	tmp.n <- nrow(geneset);
	re <- c();
	
	for (i in 1:perm){
		cat("permutation",i,"\n");
		geneset.perm <- geneset;
		rownames(geneset.perm) <- rownames(fc)[sample(N, tmp.n)];
		
		tmp <- c();
		for (j in 1:nrow(setpair)){
			set1 <- rownames(geneset.perm)[geneset.perm[,setpair[j,1]] == 1];
			set2 <- rownames(geneset.perm)[geneset.perm[,setpair[j,2]] == 1];
			
			tmp2 <- GSCoL(fc, set1, set2);
			tmp <- c(tmp, tmp2[[1]]$GSCoL, tmp2[[2]]$GSCoL);
		}
		re <- cbind(re, tmp);
	}
	
	write.table(re, 'null_GSCoL.txt', row.names=FALSE, col.names=FALSE, sep="\t");
	return (re);
}

GSCoL <- function(fc, set1.geneID, set2.geneID){
	if(contained(set1.geneID, set2.geneID)){
		return( list(NA, NA) );
	} else{
		return(GSCoL.inter(fc, set1.geneID, set2.geneID));
	}
}

GSCoL.inter <- function(fc, set1.geneID, set2.geneID){
	overlap.id<-intersect(set1.geneID,set2.geneID);
	set1.geneID.diff<-setdiff(set1.geneID, overlap.id);
	set2.geneID.diff<-setdiff(set2.geneID, overlap.id);
	
	Ex.set1<-t(fc[set1.geneID.diff,]);
	Ex.set2<-t(fc[set2.geneID,]);
	GSCoL1 <- GSCoL.compu(Ex.set1, Ex.set2);
	
	Ex.set1<-t(fc[set1.geneID,]);
	Ex.set2<-t(fc[set2.geneID.diff,]);
	GSCoL2 <- GSCoL.compu(Ex.set1, Ex.set2);
	
	list(GSCoL1, GSCoL2)
}

# differential correlation of pathway, output the weight vectors u,v
GSCoL.compu <- function(Ex.set1, Ex.set2){
	CCA.output <- CCA.compu(Ex.set1, Ex.set2);
	DC <- cor(scale(Ex.set1,T,T)%*%CCA.output$u, scale(Ex.set2,T,T)%*%CCA.output$v);
	list(set1=colnames(Ex.set1), set2=colnames(Ex.set2), u=CCA.output$u, v=CCA.output$v, GSCoL=DC)
}

CCA.compu<-function(x,z){
	#tmp.num<-nrow(x)-2;
	#penaltyx<-sqrt(min(ncol(x),tmp.num))/sqrt(ncol(x));
	#penaltyz<-sqrt(min(ncol(z),tmp.num))/sqrt(ncol(z));
	return(CCA(x,z, typex="standard", typez="standard", penaltyx=NULL, penaltyz=NULL, trace=FALSE));
}

post <- function(obs, re){
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
}

main()
