## Table of Contents
- [Compare with CBA and PB](#compare-with-cba-and-pb)
- [compare with FHS](#compare-with-fhs)
- [compare with BEST](#compare-with-best)


#### Compare with CBA and PB
```R
#----------
# compare with CBA and PB
#----------

AA_thr = 2.267195e-4
BMC2014 = read.csv("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/compare/BMCgenomics2014/12864_2013_5906_MOESM1_ESM.csv")
BMC_PB = BMC2014[which(BMC2014$PB_p_holms==TRUE),]
BMC_CBA = BMC2014[which(BMC2014$CBA_p_holms==TRUE),]

lifted_coor = read.table("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/compare/BMCgenomics2014/hglft_genome_280db_1b3a50.bed")
colnames(lifted_coor) = c("chr","newpos","newpos_")
BMC2014 = cbind(BMC2014,lifted_coor)

BMC_PB = BMC2014[which(BMC2014$PB_p_holms==TRUE),]
BMC_CBA = BMC2014[which(BMC2014$CBA_p_holms==TRUE),]

# PB
common_num = c()
significant_num = c()
common_pairs = list()
for(chr in 1:22){
  print(chr)
  load(paste0(path_analysis,"/dfcpg_chr_",chr,".RData"))
  BMC_PB_chr = BMC_PB[which(BMC_PB$Chr==chr),]
  index_rs = which(as.character(dfcpg_chr$ps) %in% as.character(BMC_PB_chr$newpos ))

  if(length(index_rs)>0){
    filter_rs = dfcpg_chr[index_rs , ]
    index_cpg = which(as.character(filter_rs$cpg) %in% as.character(BMC_PB$CpG))
    if(length(index_cpg)>0){
		filter_rs$cpg_snppos = paste0(filter_rs$cpg,filter_rs$ps)
		BMC_PB$cpg_snppos = paste0(BMC_PB$CpG,BMC_PB$newpos)
		common_pairs[[chr]] = merge(filter_rs, BMC_PB,by="cpg_snppos")
		common_num[chr] = dim(common_pairs[[chr]])[1]
		significant_num[chr] = sum(common_pairs[[chr]]$significant==1)
    }
  }
}



# CBA
common_num = c()
significant_num = c()
common_pairs = list()
for(chr in 1:22){
  print(chr)
  load(paste0(path_analysis,"/dfcpg_chr_",chr,".RData"))
  BMC_CBA_chr = BMC_CBA[which(BMC_CBA$Chr==chr),]
  index_rs = which(as.character(dfcpg_chr$ps) %in% as.character(BMC_CBA_chr$newpos ))

  if(length(index_rs)>0){
    filter_rs = dfcpg_chr[index_rs , ]
    index_cpg = which(as.character(filter_rs$cpg) %in% as.character(BMC_CBA$CpG))
    if(length(index_cpg)>0){
		filter_rs$cpg_snppos = paste0(filter_rs$cpg,filter_rs$ps)
		BMC_CBA$cpg_snppos = paste0(BMC_CBA$CpG,BMC_CBA$newpos)
		common_pairs[[chr]] = merge(filter_rs, BMC_CBA,by="cpg_snppos")
		common_num[chr] = dim(common_pairs[[chr]])[1]
		significant_num[chr] = sum(common_pairs[[chr]]$significant==1)
    }
  }
}
```

#### compare with FHS
```R

#---------------------
# compare with FHS
#---------------------

AA_thr = 2.267195e-4
FHS_thr = 2e-11

# all FHS significant pairs
27235697

common_pairs_num = c()
also_significant_num = c()
FHS_pairs_num_50kb = c()
cpg_diff = c()
# common_pairs_num_50kb = c()
# also_significant_num_50kb = c()
for(chr in 1:22){
	print(chr)
	load(paste0(path_analysis,"/dfcpg_chr_",chr,".RData"))
	load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/compare/FHS/FHS_chr",chr,".RDta"))
	FHS$cpg=FHS$CpG
	FHS$ps = FHS$SNP_Pos_hg19
	common_data = merge(dfcpg_chr, FHS, by=c("cpg","ps"))
	common_pairs_num[chr] = dim(common_data)[1]
	also_significant_num[chr] = sum(common_data$significant==1)

	FHS_50kb = FHS[which(abs(FHS$DNAm_Pos_hg19-FHS$ps)<=50001),]
	FHS_pairs_num_50kb[chr] = dim(FHS_50kb)[1]
	cpg_diff = c(cpg_diff, setdiff(FHS$cpg, dfcpg_chr$cpg))
	# common_data_50kb = merge(dfcpg_chr, FHS_50kb, by=c("cpg","ps"))
	# common_pairs_num_50kb[chr] = dim(common_data_50kb)[1]
	# also_significant_num_50kb[chr] = sum(common_data_50kb$significant==1)
}

sum(common_pairs_num)
sum(also_significant_num)
sum(also_significant_num)/sum(common_pairs_num)
sum(FHS_pairs_num_50kb)
		
```

#### compare with BEST
```R

#---------------------
# compare with BEST
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.hq68q
#---------------------
library(purrr)
library(dplyr)
# need to map with positions
# (SNPs with P<0.05)
cooccurring_cismqtl = read.table("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/compare/cooccuring/genome_wide_cis_mQTL_p_ltp05.txt")
cooccurring = cooccurring_cismqtl
cooccurring_cismqtl=cooccurring[-1,]
colnames(cooccurring_cismqtl) = c( "SNP"   ,    "cpg" ,  "beta" , "tstat"  ,  "pvalue")
rsnum = unlist(lapply(as.character(cooccurring_cismqtl$SNP), function(x) as.numeric(strsplit(x, '[._]')[[1]][2])))
cooccurring_cismqtl$pos=as.character(rsnum)
chr = unlist(lapply(as.character(cooccurring_cismqtl$SNP), function(x) as.numeric(strsplit(x, '[X._]')[[1]][2])))
cooccurring_cismqtl$chr=as.character(chr)
cooccurring_cismqtl$pvalue = as.numeric(as.character(cooccurring_cismqtl$pvalue))
BEST_cutoff= cooccurring_cismqtl %>% group_by(cpg) %>% slice(which.min(pvalue))
cooccurring_cismqtl_significant = cooccurring_cismqtl[which(cooccurring_cismqtl$pvalue <=4.138e-06),]
cooccurring_cismqtl_significant$cpg_rspos = unlist(pmap(cooccurring_cismqtl_significant[c(2,6)], paste, sep = '_'))
#save(cooccurring_cismqtl_significant, file = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v5/cooccurring_cismqtl_significant.RData")


load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v5/cooccurring_cismqtl_significant.RData")
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/Illumina_anno450k.RData")
cooccurring_cismqtl_significant$DNAm_Pos_hg19 = ann450k$pos[match(cooccurring_cismqtl_significant$cpg, ann450k$cpg)]
cooccurring_cismqtl_significant$pos = as.numeric(cooccurring_cismqtl_significant$pos)
common_pairs_num_BEST = c()
also_significant_num_BEST = c()
BEST_pairs_num_50kb = c()
cpg_diff_BEST = c()
# common_pairs_num_50kb = c()
# also_significant_num_50kb = c()
for(chr in 1:22){
	print(chr)
	load(paste0(path_analysis,"/dfcpg_chr_",chr,".RData"))
	cooccurring_cismqtl_significant_chr = cooccurring_cismqtl_significant[which(cooccurring_cismqtl_significant$chr == chr),]
	cooccurring_cismqtl_significant_chr$ps = cooccurring_cismqtl_significant_chr$pos
	common_data = merge(dfcpg_chr, cooccurring_cismqtl_significant_chr, by=c("cpg","ps"))
	common_pairs_num_BEST[chr] = dim(common_data)[1]
	also_significant_num_BEST[chr] = sum(common_data$significant==1)

# need cpg position in 450k annotation
	BEST_50kb = cooccurring_cismqtl_significant_chr[which(abs(cooccurring_cismqtl_significant_chr$DNAm_Pos_hg19-cooccurring_cismqtl_significant_chr$pos)<=50001),]
	BEST_pairs_num_50kb[chr] = dim(BEST_50kb)[1]
	cpg_diff_BEST = c(cpg_diff_BEST, setdiff(cooccurring_cismqtl_significant_chr$cpg, dfcpg_chr$cpg))
}





```
