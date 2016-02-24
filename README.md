# sGSCA
Signature-based gene set co-expression analysis

R codes (Version 2013.4.5)

[Homepage @ Github](http://wt2015-github.github.io/sGSCA/) | [Homepage @ Jin Gu's Lab](http://bioinfo.au.tsinghua.edu.cn/member/jgu/sGSCA/) | [Source Code](https://github.com/wt2015-github/sGSCA)

## Introduction
sGSCA uses sparse Canonical Correlation Analysis (sCCA) to calculate correlation between gene sets with the fold-change matrix and identify significantly correlated gene set pairs as well as inclusive signature genes.

## Usage
* Main script **main_sGSCA.R** uses sCCA to calculate correlation between gene sets with the fold-change matrix, and uses gene permutation to get null distributions, outputs results and weight vectors u,v
```
source('main_sGSCA.R')
main(perm=10, fc.file="./data_exp.gse22058/fc_matrix.txt", geneset.file="./data_exp.gse22058/geneset.new.txt", LM.filter=TRUE, LMgene.file="./data_exp.gse22058/ttest_genes.txt", ppi.file="./data_exp.gse22058/incorporated_gene_relations.txt")
```
* Script **get_size.signature.R** extracts signature genes from outputs and gets the size and overlap information
* Script **get_pvalue.R** gets statistic significance from outputs, the same with the function "post" in the main script, used for integrating parallel null_GSCoL results in the end

## Arguments in the main script
* **perm**: number of permutation, set according to parallel computation
*	**fc.file**: fold-change matrix (or other expression profile), more samples are better
*	**geneset.file**: gene set matrix, row is gene, column is geneset, 0/1 represent containing relationship
*	**LM.filter**: TRUE or FALSE, filter genes according to Leterature Mining results or not
*	**LMgene.file**: gene filtering file, use gene names, so at least 2 columns (gene name colume and gene score column) with column names
*	**ppi.file**: ppi file with 3 columns, the third column is source

## Output files
* **geneset_filtered.txt**: geneset after filtering with expression data 
*	**geneset_pair.txt**: candidate geneset pairs
*	**observe_GSCoL.txt**: observed/foreground GSCoL
*	**obs_uv.txt**: u,v vectors from sCCA
*	**null_GSCoL.txt**: null/background GSCoL
*	**output_pvalue_fdr.txt**: results

## Publication
* Ting Wang, Jin Gu, Jun Yuan, Ran Tao, Yanda Li and Shao Li. Inferring pathway crosstalk networks using gene set co-expression signatures. *Molecular BioSystems* 2013, 9(7):1822-1888. PMID: 23591523.

## Contact
[Ting Wang](http://wt2015-github.github.io/) ([email](wang9ting@gmail.com)), [Jin Gu](http://bioinfo.au.tsinghua.edu.cn/member/jgu/) ([email](wellgoo@gmail.com))
