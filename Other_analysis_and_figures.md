


```R

# histogram of PVE
pdf(paste0(path_figure,"/Histogram_pve.pdf"),width=8, height=8)
ggplot(topcpg_pve, aes(x=pve)) +
  geom_histogram(bins = 150,color="purple", fill="lightpurple")+
  geom_vline(aes(xintercept=median(pve)),linetype="dashed")+ 
  theme_minimal(base_size = 22)+
  geom_density(alpha=0.6)
dev.off()

# histogram of PVE in significant cpgs
pdf(paste0(path_figure,"/Histogram_pve_mCPG.pdf"),width=8, height=8)
ggplot(topcpg_pve_mQTL, aes(x=pve)) +
  geom_histogram(bins = 150,color="purple", fill="lightpurple")+
  geom_vline(aes(xintercept=median(pve)),linetype="dashed")+ 
  theme_minimal(base_size = 22)+
  geom_density(alpha=0.6)
dev.off()


# cis-PVE vs trans-PVE

library(ggplot2)

topcpg_pve$por_cis_combined = topcpg_pve$pve * topcpg_pve$pge
topcpg_pve$por_trans_combined = topcpg_pve$pve * (1 - topcpg_pve$pge)
cis = c(topcpg_pve$por_cis_combined, topcpg_pve$por_cis_combined[topcpg_pve$significant==0],topcpg_pve$por_cis_combined[topcpg_pve$significant==1])
trans = c(topcpg_pve$por_trans_combined, topcpg_pve$por_trans_combined[topcpg_pve$significant==0],topcpg_pve$por_trans_combined[topcpg_pve$significant==1])
additive = c(topcpg_pve$pve, topcpg_pve$pve[topcpg_pve$significant==0],topcpg_pve$pve[topcpg_pve$significant==1])
effect = c(cis, trans, additive)
effect_type = c(rep("cis-PVE",length(cis)), rep("trans-PVE",length(trans)), rep("total PVE",length(additive)))
effect_type <- factor(effect_type,levels = c('cis-PVE','trans-PVE','total PVE'),ordered = TRUE)
tmp = c(rep("All cpg sites", length(topcpg_pve$por_cis_combined)), rep("non cis-mQTL cpgs", sum(topcpg_pve$significant==0)),rep("cis-mQTL cpgs", sum(topcpg_pve$significant==1)))
effect_class = factor(rep(tmp, 3),levels = c('All cpg sites','non cis-mQTL cpgs','cis-mQTL cpgs'),ordered = TRUE)
dat = data.frame(effect,effect_type,effect_class)



pdf(paste0("PVE.pdf"),width=18, height=12)
ggplot(dat, aes(x = effect_class, y = effect,fill = effect_type)) +
scale_y_continuous(name = "Proportion of variance explained",breaks = seq(0, 1, 0.1),limits=c(0, 1)) +
scale_x_discrete(name = "") + 
geom_violin(trim = FALSE, alpha=0.5, scale = "width",show.legend = FALSE) + 
geom_boxplot( aes(x = effect_class, y = effect,fill = effect_type),width = 0.2,position=position_dodge(0.9)) +
ggtitle("African American") +
theme_bw(base_size = 30) +
theme(plot.title = element_text(size = 40,face = "bold"),text = element_text(size = 40),axis.title = element_text(face="bold"),axis.text.x=element_text(size = 30)) +
scale_fill_brewer(palette = "Accent")+
theme(legend.position="bottom")
dev.off()


# enrichment of meQTL signals around CpG site

snppos_0 = top_cpg_significant$ps 
snp_cpg_start = top_cpg_significant$start
snp_cpg_end = top_cpg_significant$end
snp_cpg_len = top_cpg_significant$end - top_cpg_significant$start

snpdat = data.frame(snppos_0, snp_cpg_start, snp_cpg_end,snp_cpg_len)
dist_to_start = (snpdat$snp_cpg_start - snpdat$snppos_0)/1000
snpdat$dist_to_start = dist_to_start

pdf(paste0(path_figure,"/enrichment_mqtlcpg.pdf"),width=10, height=10)
a <- ggplot(snpdat, aes(x = dist_to_start))
a + geom_histogram(aes(y = ..density..), 
                   colour="black", fill="#0073C2FF",bins = 100) +
  geom_density(alpha = 0.2, fill = "#FF6666") +
  theme_minimal(base_size = 22)+
   labs(title="+/- 50kb",x="Distance to CpG starting site (kb)", y = "Density")
dev.off()

snppos_0 = topcpg$ps 
snp_cpg_start = topcpg$start
snp_cpg_end = topcpg$end
snp_cpg_len = topcpg$end - topcpg$start

snpdat = data.frame(snppos_0, snp_cpg_start, snp_cpg_end,snp_cpg_len)
dist_to_start = (snpdat$snp_cpg_start - snpdat$snppos_0)/1000
snpdat$dist_to_start = dist_to_start

pdf(paste0(path_figure,"/enrichment_allcpg.pdf"),width=10, height=10)
a <- ggplot(snpdat, aes(x = dist_to_start))
a + geom_histogram(aes(y = ..density..), 
                   colour="black", fill="#0073C2FF",bins = 250) +
  geom_density(alpha = 0.2, fill = "#FF6666") +
  theme_minimal(base_size = 22)+
   labs(title="+/- 50kb",x="Distance to CpG starting site (kb)", y = "Density")
dev.off()

#---------------------------
# functional enrichment
#---------------------------



testEnrichment = function(Goi, gsets, background,  mindiffexp=2) {
  scores = do.call(rbind.data.frame, lapply(names(gsets), function(gsetid){
    gset = gsets[[gsetid]]

    # much faster than using table
    tp = length(intersect(Goi, gset))
    fn = length(Goi) - tp
    fp = length(gset) - tp
    tn = length(background) - tp - fp - fn

    contingency_table = matrix(c(tp, fp, fn, tn), 2, 2)

    if (contingency_table[1,1] < mindiffexp){
      return(NULL)
    }

    fisher_result = stats::fisher.test(contingency_table, alternative="two.sided",conf.int = TRUE)
    list(pval=fisher_result$p.value, odds=fisher_result$estimate, odds_conf_low = fisher_result$conf.int[1],odds_conf_high = fisher_result$conf.int[2],found=contingency_table[1,1], gsetid=gsetid)
  }))
  if (nrow(scores) > 0) {
    scores$qval = stats::p.adjust(scores$pval, method="fdr")
  }
  scores
}

### for hypomethylated cpg sites:

annodata$refgene_group[which(is.na(annodata$refgene_group))] = "Intergenic"
Goi = as.character(annodata$cpg)[which(annodata$meanBeta<0.5)]
gsets = list("Island"= as.character(annodata$cpg)[which(annodata$Relation_to_Island=="Island")], 
       "N_Shelf" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="N_Shelf")], 
       "N_Shore" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="N_Shore")], 
       "OpenSea" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="OpenSea")], 
       "S_Shelf" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="S_Shelf")], 
       "S_Shore" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="S_Shore")],
       "OpenSea" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="OpenSea")])
background = as.character(annodata$cpg)
res_fisher1 = testEnrichment(Goi, gsets, background)

Goi = as.character(annodata$cpg)[which(annodata$meanBeta<0.5)]
gsets = list("1stExon"= as.character(annodata$cpg)[which(annodata$refgene_group=="1stExon")], 
       "3'UTR" = as.character(annodata$cpg)[which(annodata$refgene_group=="3'UTR")], 
       "5'UTR" = as.character(annodata$cpg)[which(annodata$refgene_group=="5'UTR")], 
       "Body" = as.character(annodata$cpg)[which(annodata$refgene_group=="Body")], 
       "ExonBnd" = as.character(annodata$cpg)[which(annodata$refgene_group=="ExonBnd")], 
       "TSS1500" = as.character(annodata$cpg)[which(annodata$refgene_group=="TSS1500")], 
       "TSS200" = as.character(annodata$cpg)[which(annodata$refgene_group=="TSS200")],
       "Intergenic" = as.character(annodata$cpg)[which(annodata$refgene_group=="Intergenic")])
background = as.character(annodata$cpg)
res_fisher2 = testEnrichment(Goi, gsets, background)

res_fisher = rbind(res_fisher1, res_fisher2)

cbp1=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
res_fisher$gsetid = factor(res_fisher$gsetid, levels=c("Island","S_Shelf","N_Shelf","N_Shore","S_Shore","OpenSea",
  "TSS200","1stExon","5'UTR","3'UTR","TSS1500","Body","Intergenic"),order=T)
res_fisher$l = res_fisher$gsetid 
pdf("v3_hypomethylated_enrichment.pdf", width=8,height=5)
ggplot(res_fisher, aes(y= odds, x = l,color=l)) +
geom_point(size=3) +
geom_errorbar(aes(ymin=odds_conf_low, ymax=odds_conf_high), width=0.5,size=1) +
#scale_y_log10(breaks=ticks, labels = ticks) +
geom_hline(yintercept = 1, linetype=2) +
coord_flip() +
# ylim(0.4,1.6)+
labs(x="", y = "Odds Ratio") +
theme_linedraw()+
scale_color_manual(values = cbp1)+
theme(plot.title = element_text(size = 15,  face = "bold"),
              text = element_text(size = 15),
              axis.title = element_text(),
              axis.text.x=element_text(size = 15) ,
              legend.position = "none",
              plot.margin = margin(1, 1, 1, 1, "cm")) 
dev.off()


### for hypermethylated cpg sites:

annodata$refgene_group[which(is.na(annodata$refgene_group))] = "Intergenic"
Goi = as.character(annodata$cpg)[which(annodata$meanBeta>0.5)]
gsets = list("Island"= as.character(annodata$cpg)[which(annodata$Relation_to_Island=="Island")], 
       "N_Shelf" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="N_Shelf")], 
       "N_Shore" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="N_Shore")], 
       "OpenSea" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="OpenSea")], 
       "S_Shelf" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="S_Shelf")], 
       "S_Shore" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="S_Shore")],
       "OpenSea" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="OpenSea")])
background = as.character(annodata$cpg)
res_fisher1 = testEnrichment(Goi, gsets, background)

Goi = as.character(annodata$cpg)[which(annodata$meanBeta>0.5)]
gsets = list("1stExon"= as.character(annodata$cpg)[which(annodata$refgene_group=="1stExon")], 
       "3'UTR" = as.character(annodata$cpg)[which(annodata$refgene_group=="3'UTR")], 
       "5'UTR" = as.character(annodata$cpg)[which(annodata$refgene_group=="5'UTR")], 
       "Body" = as.character(annodata$cpg)[which(annodata$refgene_group=="Body")], 
       "ExonBnd" = as.character(annodata$cpg)[which(annodata$refgene_group=="ExonBnd")], 
       "TSS1500" = as.character(annodata$cpg)[which(annodata$refgene_group=="TSS1500")], 
       "TSS200" = as.character(annodata$cpg)[which(annodata$refgene_group=="TSS200")],
       "Intergenic" = as.character(annodata$cpg)[which(annodata$refgene_group=="Intergenic")])
background = as.character(annodata$cpg)
res_fisher2 = testEnrichment(Goi, gsets, background)

res_fisher = rbind(res_fisher1, res_fisher2)

cbp1=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
res_fisher$gsetid = factor(res_fisher$gsetid, levels=c("Island","S_Shelf","N_Shelf","N_Shore","S_Shore","OpenSea",
  "TSS200","1stExon","5'UTR","3'UTR","TSS1500","Body","Intergenic"),order=T)
res_fisher$l = res_fisher$gsetid 
pdf("v3_hypermethylated_enrichment.pdf", width=8,height=5)
ggplot(res_fisher, aes(y= odds, x = l,color=l)) +
geom_point(size=3) +
geom_errorbar(aes(ymin=odds_conf_low, ymax=odds_conf_high), width=0.5,size=1) +
#scale_y_log10(breaks=ticks, labels = ticks) +
geom_hline(yintercept = 1, linetype=2) +
coord_flip() +
# ylim(0.4,1.6)+
labs(x="", y = "Odds Ratio") +
theme_linedraw()+
scale_color_manual(values = cbp1)+
theme(plot.title = element_text(size = 15,  face = "bold"),
              text = element_text(size = 15),
              axis.title = element_text(),
              axis.text.x=element_text(size = 15) ,
              legend.position = "none",
              plot.margin = margin(1, 1, 1, 1, "cm")) 
dev.off()


#---------------------------
# for colocalized eQTL and meQTL pairs:
# direction of association of eQTL and meQTL
# direction of correlation of methylation and expression
#---------------------------
library(plyr)
library(dplyr)

load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/CpG_gene_table.RData")
path_mQTL = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu"
path_eQTL = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc"
path_analysis = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/analysis"
path_mediation = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/mediation"

load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/gene_AA_anno_order_correct.RData")
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/CpG_gene_table.RData")

percent_p12 = 0.75
load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/mQTL_leadSNPgeno/Version1_p12_percent",percent_p12*100,".RData"))
p12=0.0008403957*percent_p12
p1=0.0008403957-p12
p2=0.03115616-p12
out_p12_percent = data.frame(do.call(rbind, out))
out_p12_percent$prob = out_p12_percent$PP.H4.abf / (out_p12_percent$PP.H3.abf +out_p12_percent$PP.H4.abf )
tmpp=na.omit(out_p12_percent)

genenames = c()
for(i in 1:dim(tmpp)[1]){
  print(i)
  load(paste0(path_mediation,"/result/Gene_",i,".RData"))
  genenames[i] = colnames(res_expr)
}


gene_names = c()
cpg_names = c()
SNP_names = c()
correlations = c()
associations = c()
correlations_pval = c()
for(i in 1:dim(tmpp)[1]){
  print(i)
  gene = as.character(tmpp$gene_names[i])
  if(length(which(genenames_4994 %in% gene))>0){
  load(paste0(path_mediation,"/result/Gene_",which(genenames_4994 %in% gene),".RData"))
  expr = res_expr[,1]
  if(length(which(colnames(res_CpG_mat) %in% as.character(tmpp$cpg_names[i]) ))>0){
  methy = res_CpG_mat[,which(colnames(res_CpG_mat) %in% as.character(tmpp$cpg_names[i]) )]
  correlations[i] = cor(expr, methy,method="spearman")
  correlations_pval[i] = cor.test(expr, methy,method="spearman")$p.value
  associations[i] = tmpp$associations_eqtl[i] * tmpp$associations_mqtl[i]
}
}
}


coloc_result = tmpp
coloc_result$correlation = correlations
coloc_result$associations = associations
coloc_result$correlations_pval = correlations_pval
coloc_result = na.omit(coloc_result)

# focus on promoter region association and correlation
coloc_result$genesymbol = c()
for(i in 1:dim(coloc_result)[1]){
	coloc_result$genesymbol[i] = as.character(gene_AA_anno_order$gene)[which(gene_AA_anno_order$GENE == coloc_result$gene_names[i])]
}

# annotate gene and cpg relation

coloc_result$region = NA
for(i in 1:dim(coloc_result)[1]){
	print(i)
	ind = intersect(which(CpG_gene_table$allCpG_vec == coloc_result$cpg_names[i]) ,which(CpG_gene_table$allgenes_vec == coloc_result$genesymbol[i]))
	if(length(ind)!=0){
	coloc_result$region[i] = CpG_gene_table$allregion_vec[ind]
	print(coloc_result$region[i])
}
}


ind1 = which(coloc_result$region == "5'UTR")
ind2 = which(coloc_result$region == "TSS200")
ind3 = which(coloc_result$region == "TSS1500")
ind4 = which(coloc_result$region == "1stExon")


#------- promoter region
coloc_result_promoter = coloc_result[unique(c(ind1,ind2,ind3,ind4)),]
asso = c(sum(coloc_result_promoter$associations<0)/dim(coloc_result_promoter)[1],
sum(coloc_result_promoter$associations<0 & coloc_result_promoter$PP.H4.abf>0.8)/sum(coloc_result_promoter$PP.H4.abf>0.8),
sum(coloc_result_promoter$associations<0 & coloc_result_promoter$PP.H4.abf>0.9)/sum(coloc_result_promoter$PP.H4.abf>0.9),
sum(coloc_result_promoter$associations<0 & coloc_result_promoter$PP.H4.abf>0.95)/sum(coloc_result_promoter$PP.H4.abf>0.95),
sum(coloc_result_promoter$associations<0 & coloc_result_promoter$PP.H4.abf>0.99)/sum(coloc_result_promoter$PP.H4.abf>0.99))

PP4 = rep(c("All","0.8","0.9","0.95","0.99"),2)
PP4 = factor(PP4, levels = c("All","0.8","0.9","0.95","0.99"), order=T)
neg = round(asso,2)*100
pos = 100-neg
Percent = c(pos,neg)
class = c(rep("Same direction",length(pos)),rep("Different direction",length(neg)))
dat = data.frame(PP4,Percent,class)
pdf("coloc_Direction_eQTL_mQTL_promoter_region.pdf",width=8,height=6)
ggplot(data=dat, aes(x=PP4, y=Percent, fill=class)) +
  geom_bar(stat="identity", color="black",position=position_dodge(),width=.8)+
  scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values=c( "#FFDB6D", "#4E84C4"))+
  labs(title="Promoter region",
  		subtitle="Direction of association: eQTL and meQTL",
  				x="PP4 in coloc", y = "Percentage (%)")+
  theme_bw()+
  ylim(0,100)+
  geom_text(aes(label=Percent), vjust=-1, color="black",
            position = position_dodge(0.9), size=6)+
  theme(plot.title = element_text(size = 20,  face = "bold"),
              text = element_text(size = 20),
              axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 20) ,
              legend.position = "bottom") 

dev.off()

coloc_result_nonpromoter = coloc_result[-unique(c(ind1,ind2,ind3,ind4)),]
dim(coloc_result_nonpromoter)
asso = c(sum(coloc_result_nonpromoter$associations<0)/dim(coloc_result_nonpromoter)[1],
sum(coloc_result_nonpromoter$associations<0 & coloc_result_nonpromoter$PP.H4.abf>0.8)/sum(coloc_result_nonpromoter$PP.H4.abf>0.8),
sum(coloc_result_nonpromoter$associations<0 & coloc_result_nonpromoter$PP.H4.abf>0.9)/sum(coloc_result_nonpromoter$PP.H4.abf>0.9),
sum(coloc_result_nonpromoter$associations<0 & coloc_result_nonpromoter$PP.H4.abf>0.95)/sum(coloc_result_nonpromoter$PP.H4.abf>0.95),
sum(coloc_result_nonpromoter$associations<0 & coloc_result_nonpromoter$PP.H4.abf>0.99)/sum(coloc_result_nonpromoter$PP.H4.abf>0.99))

dim(coloc_result_nonpromoter)[1]
sum(coloc_result_nonpromoter$PP.H4.abf>0.8)
sum(coloc_result_nonpromoter$PP.H4.abf>0.9)
sum(coloc_result_nonpromoter$PP.H4.abf>0.95)
sum(coloc_result_nonpromoter$PP.H4.abf>0.99)


PP4 = rep(c("All","0.8","0.9","0.95","0.99"),2)
PP4 = factor(PP4, levels = c("All","0.8","0.9","0.95","0.99"), order=T)
neg = round(asso,2)*100
pos = 100-neg
Percent = c(pos,neg)
class = c(rep("Same direction",length(pos)),rep("Different direction",length(neg)))
dat = data.frame(PP4,Percent,class)
pdf("coloc_Direction_eQTL_mQTL_non_promoter_region.pdf",width=8,height=6)
ggplot(data=dat, aes(x=PP4, y=Percent, fill=class)) +
  geom_bar(stat="identity", color="black",position=position_dodge(),width=.8)+
  scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values=c( "#FFDB6D", "#4E84C4"))+
  labs(title="Non-Promoter region",
  		subtitle="Direction of association: eQTL and meQTL",
  				x="PP4 in coloc", y = "Percentage (%)")+
  theme_bw()+
  ylim(0,100)+
  geom_text(aes(label=Percent), vjust=-1, color="black",
            position = position_dodge(0.9), size=6)+
  theme(plot.title = element_text(size = 20,  face = "bold"),
              text = element_text(size = 20),
              axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 20) ,
              legend.position = "bottom") 

dev.off()

//////////////////////
///// correlaion 
//////////////////////

# ---promoter

sum(coloc_result_promoter$correlations_pval<0.05)
# > sum(coloc_result_promoter$correlations_pval<0.05)
# [1] 510
coloc_result_promoter_sig = coloc_result_promoter[which(coloc_result_promoter$correlations_pval<0.05),]


library(ggplot2)
corr = c(sum(coloc_result_promoter_sig$correlation<0)/dim(coloc_result_promoter_sig)[1],
sum(coloc_result_promoter_sig$correlation<0 & coloc_result_promoter_sig$PP.H4.abf>0.8)/sum(coloc_result_promoter_sig$PP.H4.abf>0.8),
sum(coloc_result_promoter_sig$correlation<0 & coloc_result_promoter_sig$PP.H4.abf>0.9)/sum(coloc_result_promoter_sig$PP.H4.abf>0.9),
sum(coloc_result_promoter_sig$correlation<0 & coloc_result_promoter_sig$PP.H4.abf>0.95)/sum(coloc_result_promoter_sig$PP.H4.abf>0.95),
sum(coloc_result_promoter_sig$correlation<0 & coloc_result_promoter_sig$PP.H4.abf>0.99)/sum(coloc_result_promoter_sig$PP.H4.abf>0.99))

dim(coloc_result_promoter_sig)[1]
sum(coloc_result_promoter_sig$PP.H4.abf>0.8)
sum(coloc_result_promoter_sig$PP.H4.abf>0.9)
sum(coloc_result_promoter_sig$PP.H4.abf>0.95)
sum(coloc_result_promoter_sig$PP.H4.abf>0.99)


PP4 = rep(c("All","0.8","0.9","0.95","0.99"),2)
PP4 = factor(PP4, levels = c("All","0.8","0.9","0.95","0.99"), order=T)
neg = round(corr,2)*100
pos = 100-neg
Percent = c(pos,neg)
class = c(rep("Same direction",length(pos)),rep("Different direction",length(neg)))
dat = data.frame(PP4,Percent,class)
pdf("coloc_Direction_expression_methylation_promoter_region.pdf",width=8,height=6)
ggplot(data=dat, aes(x=PP4, y=Percent, fill=class)) +
  geom_bar(stat="identity", color="black",position=position_dodge(),width=.8)+
  scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values=c( "#CC79A7", "#C3D7A4"))+
  labs(title="Promoter region",
  	subtitle="Direction of correlation: expression and methylation", x="", y = "Percentage (%)")+
  theme_bw()+
  ylim(0,100)+
  geom_text(aes(label=Percent), vjust=-1, color="black",
            position = position_dodge(0.9), size=6)+
  theme(plot.title = element_text(size = 20,  face = "bold"),
              text = element_text(size = 20),
              axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 20) ,
              legend.position = "bottom") 

dev.off()

# non-promoter region


coloc_result_nonpromoter = coloc_result[-unique(c(ind1,ind2,ind3,ind4)),]

coloc_result_nonpromoter_sig = coloc_result_nonpromoter[which(coloc_result_nonpromoter$correlations_pval<0.05),]


library(ggplot2)
corr = c(sum(coloc_result_nonpromoter_sig$correlation<0)/dim(coloc_result_nonpromoter_sig)[1],
sum(coloc_result_nonpromoter_sig$correlation<0 & coloc_result_nonpromoter_sig$PP.H4.abf>0.8)/sum(coloc_result_nonpromoter_sig$PP.H4.abf>0.8),
sum(coloc_result_nonpromoter_sig$correlation<0 & coloc_result_nonpromoter_sig$PP.H4.abf>0.9)/sum(coloc_result_nonpromoter_sig$PP.H4.abf>0.9),
sum(coloc_result_nonpromoter_sig$correlation<0 & coloc_result_nonpromoter_sig$PP.H4.abf>0.95)/sum(coloc_result_nonpromoter_sig$PP.H4.abf>0.95),
sum(coloc_result_nonpromoter_sig$correlation<0 & coloc_result_nonpromoter_sig$PP.H4.abf>0.99)/sum(coloc_result_nonpromoter_sig$PP.H4.abf>0.99))

dim(coloc_result_nonpromoter_sig)[1]
sum(coloc_result_nonpromoter_sig$PP.H4.abf>0.8)
sum(coloc_result_nonpromoter_sig$PP.H4.abf>0.9)
sum(coloc_result_nonpromoter_sig$PP.H4.abf>0.95)
sum(coloc_result_nonpromoter_sig$PP.H4.abf>0.99)


PP4 = rep(c("All","0.8","0.9","0.95","0.99"),2)
PP4 = factor(PP4, levels = c("All","0.8","0.9","0.95","0.99"), order=T)
neg = round(corr,2)*100
pos = 100-neg
Percent = c(pos,neg)
class = c(rep("Same direction",length(pos)),rep("Different direction",length(neg)))
dat = data.frame(PP4,Percent,class)
pdf("Direction_correlation_expression_methylation_nonpromoter_region.pdf",width=8,height=6)
ggplot(data=dat, aes(x=PP4, y=Percent, fill=class)) +
  geom_bar(stat="identity", color="black",position=position_dodge(),width=.8)+
  scale_fill_brewer(palette="Paired")+
  scale_fill_manual(values=c( "#CC79A7", "#C3D7A4"))+
  labs(title="Non-promoter region",
  		subtitle="Direction of correlation: expression and methylation",
  				x="", y = "Percentage (%)")+
  theme_bw()+
  ylim(0,100)+
  geom_text(aes(label=Percent), vjust=-1, color="black",
            position = position_dodge(0.9), size=6)+
  theme(plot.title = element_text(size = 20,  face = "bold"),
              text = element_text(size = 20),
              axis.title = element_text(face="bold"),
              axis.text.x=element_text(size = 20) ,
              legend.position = "bottom") 

dev.off()

#------------------------------------
# SME and SEM p value correlations 
#------------------------------------

coloc_result_PP4$log10_SEM = -log10(coloc_result_PP4$sobel_p_SEM)
coloc_result_PP4$log10_SME = -log10(coloc_result_PP4$sobel_p_SME)
pdf("mediation_SEM_SME_pval.pdf")
p <- ggplot(coloc_result_PP4, aes(x=log10_SEM, y=log10_SME))+  
labs(title="",
         x="SEM: -log10(Sobel p)",
         y="SME: -log10(Sobel p)")+
  geom_point(size=2,color = "cornflowerblue") + 
  geom_smooth(method="lm",size=0.5,colour="red")+
  theme_bw(base_size = 22)
p
dev.off() 


cor.test(coloc_result_PP4$log10_SEM,coloc_result_PP4$log10_SME,method="spearman")$p.value
cor.test(coloc_result_PP4$log10_SEM,coloc_result_PP4$log10_SME,method="spearman")$estimate

#------------------------------------
# SEM and SME examples
#------------------------------------

load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v1/Figure4/coloc_result_PP4.RData)
coloc_result_PP4[which(coloc_result_PP4$sobel_q_SME<0.05 & coloc_result_PP4$sobel_q_SEM>=0.05 ),c(1,2,3,20)]

which(coloc_result_PP4$gene_names=="ENSG00000117280" & coloc_result_PP4$cpg_names=="cg06442372")
library(rtracklayer)
# install.packages("coMET")
# library("coMET")


load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v1/Figure5/topeQTL_mQTL_anno.RData")
geneENSG = "ENSG00000117280"


geneENSGchr = gene_AA_anno_order[which(gene_AA_anno_order$GENE==geneENSG),]$chr
cpgCPG = "cg06442372"
snpSNP = "rs1775151"
snppos=gene_AA_anno_order[which(gene_AA_anno_order$GENE==geneENSG),]$topSNP_ps

path_res = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/data/chunk/result"
res = list()
topres = list()
bslmm = list()
for(chunk in 1:100){
  print(chunk)
   load(paste0(path_res,"/chr_",geneENSGchr,"_chunk_",chunk,"_pheno_1_result.RData"))
   if(length(which(as.character(df$cpg)==cpgCPG)>0)){
    break
   }
}

load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/AA_full_result.RData")


# all QTL version
cpg_dat = df[which(as.character(df$cpg)==cpgCPG),]
gene_dat = AA_full_result[which(AA_full_result$GENE == geneENSG),]

# > gene_dat[302,]
#         chr        rs        ps n_miss allele1 allele0    af      beta
# 1192743   1 rs1775151 205745684      1       T       C 0.377 0.2270112
#                 se   l_remle       p_wald            GENE
# 1192743 0.04497418 0.1842419 5.287617e-07 ENSG00000117280


p_vals = c(gene_dat$p_wald, cpg_dat$p_wald)
positions = c(gene_dat$ps, cpg_dat$ps)
types = c(rep("Gene-SNP",length(gene_dat$ps)),rep("CpG-SNP",length(cpg_dat$ps)))
types = as.factor(types)
library(ggplot2)
data = data.frame(p_vals,positions,  types)
data$log10p = -log10(p_vals)
data$positionsMB = data$positions/1000000
data$l = data$types

data = data[which(data$positions>=(snppos-50000) &  data$positions<=(snppos+50000)),]

min(data$positions)
max(data$positions)


# ENSG00000117280:
> min(data$positions)
[1] 205695695
> max(data$positions)
[1] 205795281


# pval_thr =  0.001241127

# pval_thr = coloc_result_PP4$sobel_p_SME[which(coloc_result_PP4$sobel_q_SME==max(coloc_result_PP4$sobel_q_SME[which(coloc_result_PP4$sobel_q_SME<=0.05)]))]
# > pval_thr
# [1] 0.002763829
genename = geneENSG
ind = which(topeQTL_mQTL_anno$GENE == genename)
TargetID = topeQTL_mQTL_anno$cpg[ind]
CHR = topeQTL_mQTL_anno$chr.x[ind]
MAPINFO = topeQTL_mQTL_anno$cpgstart[ind]
Pval = topeQTL_mQTL_anno$sobel_p_SME[ind]

infotable = data.frame(TargetID, CHR, MAPINFO, Pval)
write.table(infotable, paste0("infotable_",genename,".txt"), row.names=F,quote=F,sep = "\t")

path_mediation="/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/mediation"
methylation_table = matrix(0,800,length(TargetID))
for(i in 1:length(TargetID)){
methylation = read.table(paste0(path_mediation,"/mediator_cpg/covariate_cpg/",TargetID[i],".txt"))$V2
methylation_table[,i] = methylation
}

colnames(methylation_table) = TargetID
write.table(methylation_table, paste0("methylation_table_",genename,".txt"), row.names=F,quote=F,sep = "\t")

correlation = cor(methylation_table)
colnames(correlation) = TargetID
write.table(correlation, paste0("correlation_",genename,".txt"), row.names=F,quote=F,sep = "\t")


---mulan

library(rtracklayer)
library(coMET)

genename = "ENSG00000117280"


setwd("/net/mulan/home/shanglu/GENOA_mQTL/manuscript/v1/Figure5")

extdata <- system.file("extdata", package="coMET",mustWork=TRUE) 
configfile <- file.path( paste0("/net/mulan/home/shanglu/GENOA_mQTL/manuscript/v1/Figure5/format_",genename,".txt"))
myinfofile <- file.path( paste0("/net/mulan/home/shanglu/GENOA_mQTL/manuscript/v1/Figure5/infotable_",genename,".txt") )
mycorrelation <- file.path( paste0("/net/mulan/home/shanglu/GENOA_mQTL/manuscript/v1/Figure5/methylation_table_",genename,".txt"))

# ENSG00000117280:
chrom <- "chr1" 
start <-  205695695
end <- 205795281 


gen <- "hg19" 
strand <- "*"
BROWSER.SESSION="UCSC"
mySession <- browserSession(BROWSER.SESSION) 
genome(mySession) <- gen
genetrack <-genes_ENSEMBL(gen,chrom,start,end,showId=TRUE) 
snptrack <- snpBiomart_ENSEMBL(gen,chrom, start, end,dataset="hsapiens_snp_som",showId=FALSE)
cpgIstrack <- cpgIslands_UCSC(gen,chrom,start,end)
# gwastrack <-GWAScatalog_UCSC(gen,chrom,start,end)

listgviz <- list(genetrack,snptrack,cpgIstrack)

pdf(paste0(genename,"_cpgIstrack.pdf"))
comet(config.file=configfile, 
	mydata.file=myinfofile,
	mydata.type="file",
	cormatrix.file=mycorrelation,
	cormatrix.type="listfile",
	#mydata.large.file=myexpressfile, 
	#mydata.large.type="listfile",
	fontsize.gviz =12,
tracks.gviz=listgviz,verbose=FALSE, print.image=FALSE)
dev.off()


regulatorytrack = regulatoryFeaturesBiomart_ENSEMBL(gen,chrom,start,end)
listgviz <- list(genetrack,regulatorytrack,cpgIstrack)
pdf(paste0(genename,"_regulatorytrack.pdf"))
comet(config.file=configfile, 
	mydata.file=myinfofile,
	mydata.type="file",
	cormatrix.file=mycorrelation,
	cormatrix.type="listfile",
	#mydata.large.file=myexpressfile, 
	#mydata.large.type="listfile",
	fontsize.gviz =12,
tracks.gviz=listgviz,verbose=FALSE, print.image=TRUE)
dev.off()




geneENSG = "ENSG00000117280"
cpgCPG = "cg06442372"
snpSNP = "rs1775151"
path_phenotype = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/data/phenotype"
load(paste0(path_phenotype,"/M_value_EPIC_962.RData"))
load(paste0(path_phenotype,"/sample_meta.RData"))
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_GE/preprocessing_zxu_March2017/Lulu_Mar2018/AJHG_revision/data/expr_AA_proteingene.RData")
methy_col = M_value_EPIC[which(rownames(M_value_EPIC)=="cg06442372"),na.omit(match( sample_meta$Sample_name,colnames(M_value_EPIC)))]
names = as.character(sample_meta$PIN)[match(sample_meta$Sample_name, names(methy_col))]
dat1 = data.frame(names,methy_col)
expr_col = expr_AA_proteingene[which(rownames(expr_AA_proteingene)=="ENSG00000117280"),]
dat2 = data.frame("names"=colnames(expr_AA_proteingene),"expr_col"=unlist(expr_col))
datt = merge(dat1,dat2,by="names")

cor.test(datt[,2],datt[,3])


#------------------------------------
# Disruptive SNPs
#------------------------------------

# map to meSNP results
disrupted_cpg = c()
disrupted_significant = c()
disrupted_is = c()
disrupted_pval = c()
for(chr in 1:22){

    print(chr)

    load(paste0(path_analysis,"/dfcpg_chr_",chr,".RData"))
    dfcpg_chr$disrupted = rep(0,dim(dfcpg_chr)[1])
    ind = which(dfcpg_chr$ps >= dfcpg_chr$cpgstart & dfcpg_chr$ps <= dfcpg_chr$cpgend )
    dfcpg_chr$disrupted[ind]=1

    disrupted_cpg = c(disrupted_cpg, dfcpg_chr$cpg)
    disrupted_significant = c(disrupted_significant, dfcpg_chr$significant)
    disrupted_is = c(disrupted_is, dfcpg_chr$disrupted)
    disrupted_pval = c(disrupted_pval, dfcpg_chr$p_wald)

}


disrupted_pval[which(disrupted_pval==0)]=2.23e-324
sum(disrupted_is==1)
sum(disrupted_is==0)




# proportion of disrupted CpGs to be meCPGs
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v3/disrupted_cpg/anno_disrupt_CpG_meSNP_version.RData")
length(unique(disrupted_cpg[disrupted_is==1]))
[1] 21399
length(unique(disrupted_cpg[disrupted_is==1 & disrupted_significant==1]))
[1] 18303
> 18303/21399
[1] 0.8553203

# proportion of disrupting SNPs to be meSNPs
length(unique(SNPs_all[disrupted_is==1]))
[1] 21545
length(unique(SNPs_all[disrupted_is==1 & disrupted_significant==1]))
[1] 18383
> 18383/21545
[1] 0.8532374

# disrupted but not significant:
summary(disrupted_pval[which(disrupted_is==1 & disrupted_significant==0)])
summary(disrupted_pval[which(disrupted_is==1 & disrupted_significant==1)])
> summary(disrupted_pval[which(disrupted_is==1 & disrupted_significant==0)])
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.0002275 0.0079570 0.0954118 0.2455413 0.4272044 0.9994601 
> summary(disrupted_pval[which(disrupted_is==1 & disrupted_significant==1)])
     Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
0.000e+00 0.000e+00 0.000e+00 2.076e-06 0.000e+00 2.263e-04 



Goi = as.character(disrupted_cpg[which(disrupted_is==1)])
gsets = list("meSNP"= as.character(disrupted_cpg)[which(disrupted_significant==1)], 
       "non meSNP" = as.character(disrupted_cpg)[which(disrupted_significant==0)])
background = as.character(disrupted_cpg)
res_fisher = testEnrichment(Goi, gsets, background)

# > res_fisher
#    pval       odds odds_conf_low odds_conf_high found    gsetid qval
# 2     0 157.692393    148.549273     166.573196 20093     meSNP    0
# 21    0   7.767561      6.819204       8.886893 21399 non meSNP    0

cbp1=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
res_fisher$gsetid = factor(res_fisher$gsetid, levels=c("meSNP","non meSNP"),order=T)
res_fisher$l = res_fisher$gsetid 

pdf("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v3/disrupted_cpg/v3_disrupted_CpG_enrichment_meSNP.pdf", width=8,height=5)
ggplot(res_fisher, aes(y= odds, x = l,color=l)) +
geom_point(size=3) +
geom_errorbar(aes(ymin=odds_conf_low, ymax=odds_conf_high), width=0.5,size=1) +
#scale_y_log10(breaks=ticks, labels = ticks) +
geom_hline(yintercept = 1, linetype=2) +
coord_flip() +
#ylim(0.5,2)+
labs(x="", y = "Odds Ratio",title="Disrupted CpG sites") +
theme_linedraw()+
scale_color_manual(values = cbp1)+
theme(plot.title = element_text(size = 22,  face = "bold"),
              text = element_text(size = 22),
              axis.title = element_text(),
              axis.text.x=element_text(size = 22) ,
              legend.position = "none",
              plot.margin = margin(1, 1, 1, 1, "cm")) 
dev.off()



# functional enrichment
annodata = annodata_meta
annodata$refgene_group[which(is.na(annodata$refgene_group))] = "Intergenic"
Goi = as.character(disrupted_cpg[which(disrupted_is==1)])
gsets = list("Island"= as.character(annodata$cpg)[which(annodata$Relation_to_Island=="Island")], 
       "N_Shelf" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="N_Shelf")], 
       "N_Shore" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="N_Shore")], 
       "OpenSea" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="OpenSea")], 
       "S_Shelf" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="S_Shelf")], 
       "S_Shore" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="S_Shore")],
       "OpenSea" = as.character(annodata$cpg)[which(annodata$Relation_to_Island=="OpenSea")])
background = as.character(annodata$cpg)
res_fisher1 = testEnrichment(Goi, gsets, background)
> res_fisher1
            pval      odds odds_conf_low odds_conf_high found  gsetid
2   0.000000e+00 0.2522336     0.2377358      0.2674414  1202  Island
21  2.281485e-07 1.1969968     1.1187954      1.2794797   938 N_Shelf
3  1.653184e-120 0.5226879     0.4919026      0.5549838  1150 N_Shore
4   0.000000e+00 2.3038436     2.2334909      2.3767376 16147 OpenSea
5   4.256065e-07 1.1998075     1.1187682      1.2853850   874 S_Shelf
6   1.307641e-97 0.5350299     0.5012934      0.5705031   995 S_Shore
7   0.000000e+00 2.3038436     2.2334909      2.3767376 16147 OpenSea
            qval
2   0.000000e+00
21  2.661732e-07
3  2.893071e-120
4   0.000000e+00
5   4.256065e-07
6   1.830698e-97
7   0.000000e+00

gsets = list("1stExon"= as.character(annodata$cpg)[which(annodata$refgene_group=="1stExon")], 
       "3'UTR" = as.character(annodata$cpg)[which(annodata$refgene_group=="3'UTR")], 
       "5'UTR" = as.character(annodata$cpg)[which(annodata$refgene_group=="5'UTR")], 
       "Body" = as.character(annodata$cpg)[which(annodata$refgene_group=="Body")], 
       "ExonBnd" = as.character(annodata$cpg)[which(annodata$refgene_group=="ExonBnd")], 
       "TSS1500" = as.character(annodata$cpg)[which(annodata$refgene_group=="TSS1500")], 
       "TSS200" = as.character(annodata$cpg)[which(annodata$refgene_group=="TSS200")],
       "Intergenic" = as.character(annodata$cpg)[which(annodata$refgene_group=="Intergenic")])
background = as.character(annodata$cpg)
res_fisher2 = testEnrichment(Goi, gsets, background)
> res_fisher2
            pval      odds odds_conf_low odds_conf_high found     gsetid
2   1.662065e-23 0.4659832     0.3903177      0.5522071   136    1stExon
21  1.134537e-02 1.1244292     1.0259976      1.2299974   505      3'UTR
3   6.039950e-14 0.8326701     0.7924466      0.8744402  1795      5'UTR
4   7.402809e-23 1.1504669     1.1188101      1.1829578  8376       Body
5  6.015077e-103 0.6105757     0.5817229      0.6405757  1869    TSS1500
6  7.802848e-194 0.4058361     0.3785927      0.4345837   867     TSS200
7  1.209008e-114 1.3954280     1.3564232      1.4354827  7758 Intergenic
            qval
2   2.908614e-23
21  1.134537e-02
3   7.046608e-14
4   1.036393e-22
5  1.403518e-102
6  5.461994e-193
7  4.231527e-114

res_fisher = rbind(res_fisher1, res_fisher2)


cbp1=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
res_fisher$gsetid = factor(res_fisher$gsetid, levels=c("Island","S_Shelf","N_Shelf","N_Shore","S_Shore","OpenSea",
  "TSS200","1stExon","5'UTR","3'UTR","TSS1500","Body","Intergenic"),order=T)

res_fisher$l = res_fisher$gsetid 

pdf("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v3/disrupted_cpg/v3_disrupted_CpG_allpairs_functional_enrichment.pdf", width=8,height=5)
ggplot(res_fisher, aes(y= odds, x = l,color=l)) +
geom_point(size=3) +
geom_errorbar(aes(ymin=odds_conf_low, ymax=odds_conf_high), width=0.5,size=1) +
#scale_y_log10(breaks=ticks, labels = ticks) +
geom_hline(yintercept = 1, linetype=2) +
coord_flip() +
#ylim(0.5,2)+
labs(x="", y = "Odds Ratio",title="Disrupted CpG sites") +
theme_linedraw()+
scale_color_manual(values = cbp1)+
theme(plot.title = element_text(size = 22,  face = "bold"),
              text = element_text(size = 22),
              axis.title = element_text(),
              axis.text.x=element_text(size = 22) ,
              legend.position = "none",
              plot.margin = margin(1, 1, 1, 1, "cm")) 
dev.off()

#---------------
# disrupted TFBS
#---------------

motifmap = read.table("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/TFBS/HUMAN_hg19_BBLS_1_00_FDR_0_10.bed",header = FALSE, sep="\t",stringsAsFactors=FALSE, quote="")
> dim(motifmap)
[1] 4474877       6
> head(motifmap)
     V1        V2        V3            V4       V5 V6
1 chr17  58239126  58239144 LM1_RFX1=RFX1 1.763390  +
2 chr17  58239126  58239144 LM1_RFX1=RFX1 1.763377  -
3  chr7 128800229 128800247 LM1_RFX1=RFX1 1.815382  +
4  chr7 128800229 128800247 LM1_RFX1=RFX1 1.815382  -
5  chr1   6133012   6133030 LM1_RFX1=RFX1 1.650912  -
6 chr12 104996389 104996407 LM1_RFX1=RFX1 1.247502  -
> dim(unique(motifmap[,1:3]))
[1] 3961042       3
> length(unique(motifmap[,4]))
[1] 607


library(parallel)
library(ggplot2)

fx_motif = function(i,motifmap_chr_in){
  ind_j = which(dfcpg_chr$ps>= motifmap_chr_in$V2[i]  & dfcpg_chr$ps<=motifmap_chr_in$V3[i]  )
  if(length(ind_j!=0)){
    return(list(ind_j,i))
  }
}

# map to meSNP results
# too large, need to get each chromosome seperately
for(chr in 1:22){

    print(chr)

    load(paste0(path_analysis,"/dfcpg_chr_",chr,".RData"))
    motifmap_chr = motifmap[which(as.character(motifmap$V1)==paste0("chr",chr)),]
    dfcpg_chr$TF = rep(0,dim(dfcpg_chr)[1])

    fx_motif = function(i,motifmap_chr_in){
      ind_j = which(dfcpg_chr$ps>= motifmap_chr_in$V2[i]  & dfcpg_chr$ps<=motifmap_chr_in$V3[i]  )
      if(length(ind_j!=0)){
        return(list(ind_j,i))
      }
    }

    motifmap_chr = unique(motifmap_chr[,1:4])


    site_name = c()
    for(k in 1:round(dim(motifmap_chr)[1]/1000+1)){
        print(k/round(dim(motifmap_chr)[1]/1000+1))
        if(k<dim(motifmap_chr)[1]/1000){
          start = 1000*(k-1)+1
          end = 1000*k
        }else{
          start = 1000*(k-1)+1
          end = dim(motifmap_chr)[1]
        }
        results = mclapply(start:end, fx_motif, motifmap_chr_in=motifmap_chr, mc.cores = 20)
        dfcpg_chr$TF[unlist(lapply(results, `[[`, 1))]=1
        site_name = c(site_name,rownames(motifmap_chr)[unlist(lapply(results, `[[`, 2))])
    }

    motif_is = dfcpg_chr$TF
    motif_cpg = dfcpg_chr$cpg
    motif_significant = dfcpg_chr$significant
    motif_pval = dfcpg_chr$p_wald
    save(motif_cpg, motif_significant, motif_is, motif_pval,file = paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v3/Motif/anno_disrupt_Motif_meSNP_chr_",chr,".RData"))

}


##########  too slow on local, submit jobs

cd /home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v3/Motif
for ((myline=1;myline<=22;myline++)); do
qsub -cwd -N motif  -v FOO=$myline -l vf=5G  motif_run.sh
done

##########
# collect results
motif_is_all = c()
motif_cpg_all = c()
motif_significant_all = c()
motif_pval_all = c()
site_name_all = c()
for(chr in 1:22){
  print(chr)
  load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v3/Motif/anno_disrupt_Motif_meSNP_chr_",chr,".RData"))
  motif_is_all = c(motif_is_all, motif_is)
  motif_cpg_all = c(motif_cpg_all, motif_cpg)
  motif_significant_all = c(motif_significant_all, motif_significant)
  motif_pval_all = c(motif_pval_all, motif_pval)
  site_name_all = c(site_name_all, site_name)
}


SNPs_all = c()
for(chr in 1:22){
    print(chr)
    load(paste0(path_analysis,"/dfcpg_chr_",chr,".RData"))
    SNPs_all = c(SNPs_all, dfcpg_chr$rs)
}


motif_pval_all[which(motif_pval_all==0)]=2.23e-308

# mapped to  motif sites
length(site_name_all)
[1] 87077
# SNPs disrupting motifs
sum(motif_is_all==1)
[1] 2386480
sum(motif_is_all==0)
[1] 239325283

# for motif disrupting SNPs, see if p values are smaller than non disrupting ones
dat = data.frame("minus_log10_pval"=-log10(motif_pval_all), "Disrupted" = as.character(motif_is_all))

summary(dat$minus_log10_pval[which(dat$Disrupted==1)])
# median 0.4600548
summary(dat$minus_log10_pval[which(dat$Disrupted==0)])
# median 0.4694

# These SNPs could directly disrupt motif binding and are indeed enriched with meSNPs 
# (OR= 1.802, p-value <2.23e-308) than in non-meSNPs (OR= 0.029, p-value<2.23e-308) 
Goi = as.character(motif_cpg_all[which(motif_is_all==1)]) # CpG names
gsets = list("meSNP"= as.character(motif_cpg_all)[which(motif_significant_all==1)], # CpG with SNP significant
       "non meSNP" = as.character(motif_cpg_all)[which(motif_significant_all==0)])  # CpG with SNP non-significant
background = as.character(motif_cpg_all)
res_fisher = testEnrichment(Goi, gsets, background)

> res_fisher
   pval       odds odds_conf_low odds_conf_high  found    gsetid qval
2     0 1.80210158     1.7952042     1.80896543 308887     meSNP    0
21    0 0.02946881     0.0293082     0.02956757 668575 non meSNP    0


cbp1=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
res_fisher$gsetid = factor(res_fisher$gsetid, levels=c("meSNP","non meSNP"),order=T)
res_fisher$l = res_fisher$gsetid 

pdf("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v3/Motif/v3_disrupted_Motif_enrichment_meSNP.pdf", width=8,height=5)
ggplot(res_fisher, aes(y= odds, x = l,color=l)) +
geom_point(size=3) +
geom_errorbar(aes(ymin=odds_conf_low, ymax=odds_conf_high), width=0.5,size=1) +
#scale_y_log10(breaks=ticks, labels = ticks) +
geom_hline(yintercept = 1, linetype=2) +
coord_flip() +
#ylim(0.5,2)+
labs(x="", y = "Odds Ratio",title="Disrupted Motif") +
theme_linedraw()+
scale_color_manual(values = cbp1)+
theme(plot.title = element_text(size = 22,  face = "bold"),
              text = element_text(size = 22),
              axis.title = element_text(),
              axis.text.x=element_text(size = 22) ,
              legend.position = "none",
              plot.margin = margin(1, 1, 1, 1, "cm")) 
dev.off()

# number of pairs, disrupted 
disrupted_cpg = unique(motif_cpg_all[which(motif_is_all==1)])
length(disrupted_cpg)
[1] 668576
# number of pairs, disrupted & significant:
disrupted_meCpG = unique(motif_cpg_all[which(motif_significant_all==1 & motif_is_all==1)])
length(disrupted_meCpG)
[1] 103825


# number of SNPs disrupting 
disrupting_SNP = unique(SNPs_all[which(motif_is_all==1)])
length(disrupting_SNP)
[1] 78940
# number of SNPs disrupting & significant
disrupting_meSNP = unique(SNPs_all[which(motif_significant_all==1 & motif_is_all==1)])
length(disrupting_meSNP)
[1] 46384


SNPs_all_unique = unique(SNPs_all)


# test if meSNPs are enriched in SNPs that disrupting motifs
Goi = unique(SNPs_all[which(motif_significant_all==1)])
gsets = list("SNPs Disrupting Motifs"= unique(as.character(SNPs_all)[which(motif_is_all==1)]), 
       "SNPs Not Disrupting Motifs" = unique(as.character(SNPs_all)[which(motif_is_all==0)])  )
background = SNPs_all_unique
res_fisher = testEnrichment(Goi, gsets, background)
res_fisher
> res_fisher
        pval     odds odds_conf_low odds_conf_high   found         gsetid
2  0.0000000 1.887737     1.8579628       1.918085   57861     disrupting
21 0.3177117 1.009369     0.9909982       1.028025 3490015 non disrupting
        qval
2  0.0000000
21 0.3177117

cbp1=c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
  "#FFDB6D", "#C4961A", "#F4EDCA", "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")
res_fisher$gsetid = factor(res_fisher$gsetid, levels=c("SNPs Disrupting Motifs","SNPs Not Disrupting Motifs"),order=T)
res_fisher$l = res_fisher$gsetid 

pdf("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v3/Motif/v3_meSNPs_enrich_in_SNPs_disrupting_Motifs.pdf", width=10,height=4)
ggplot(res_fisher, aes(y= odds, x = l,color=l)) +
geom_point(size=3) +
geom_errorbar(aes(ymin=odds_conf_low, ymax=odds_conf_high), width=0.5,size=1) +
#scale_y_log10(breaks=ticks, labels = ticks) +
geom_hline(yintercept = 1, linetype=2) +
coord_flip() +
#ylim(0.5,2)+
labs(x="", y = "Odds Ratio",title="meSNP enrichment") +
theme_linedraw()+
scale_color_manual(values = cbp1)+
theme(plot.title = element_text(size = 22,  face = "bold"),
              text = element_text(size = 22),
              axis.title = element_text(),
              axis.text.x=element_text(size = 22) ,
              legend.position = "none",
              plot.margin = margin(1, 1, 1, 1, "cm")) 
dev.off()



#-----------------------------------
# make coMet figures for SME and SEM
#-----------------------------------

library(rtracklayer)
library(coMET)
load("/net/mulan/home/shanglu/GENOA_mQTL/manuscript/v1/Figure5/coMet/min_max_loc_only_SME_use.RData")
load("/net/mulan/home/shanglu/GENOA_mQTL/manuscript/v1/Figure5/coMet/only_SME_use.RData")

k = 0

k = k + 1
genename = as.character(only_SME_use$gene_names[k])
extdata <- system.file("extdata", package="coMET",mustWork=TRUE) 
configfile <- file.path( paste0("/net/mulan/home/shanglu/GENOA_mQTL/manuscript/v1/Figure5/coMet/format_",genename,".txt"))
myinfofile <- file.path( paste0("/net/mulan/home/shanglu/GENOA_mQTL/manuscript/v1/Figure5/coMet/infotable_",genename,".txt") )
mycorrelation <- file.path( paste0("/net/mulan/home/shanglu/GENOA_mQTL/manuscript/v1/Figure5/coMet/methylation_table_",genename,".txt"))
chrom <- only_SME_use$chr[k]
start <-  min_loc[k]
end <- max_loc[k] 
gen <- "hg19" 
strand <- "*"
BROWSER.SESSION="UCSC"
mySession <- browserSession(BROWSER.SESSION) 
genome(mySession) <- gen
genetrack <-genes_ENSEMBL(gen,chrom,start,end,showId=TRUE) 
snptrack <- snpBiomart_ENSEMBL(gen,chrom, start, end,dataset="hsapiens_snp_som",showId=FALSE)
cpgIstrack <- cpgIslands_UCSC(gen,chrom,start,end)
listgviz <- list(genetrack,snptrack,cpgIstrack)
pdf(paste0(genename,"_cpgIstrack_k",k,".pdf"))
comet(config.file=configfile, 
  mydata.file=myinfofile,
  mydata.type="file",
  cormatrix.file=mycorrelation,
  cormatrix.type="listfile",
  #mydata.large.file=myexpressfile, 
  #mydata.large.type="listfile",
  fontsize.gviz =12,
tracks.gviz=listgviz,verbose=FALSE, print.image=FALSE)
dev.off()



regulatorytrack = regulatoryFeaturesBiomart_ENSEMBL(gen,chrom,start,end)
chromosome(regulatorytrack)=chrom
listgviz <- list(genetrack,regulatorytrack,cpgIstrack)
pdf(paste0(genename,"_regulatorytrack_k",k,".pdf"))
comet(config.file=configfile, 
  mydata.file=myinfofile,
  mydata.type="file",
  cormatrix.file=mycorrelation,
  cormatrix.type="listfile",
  #mydata.large.file=myexpressfile, 
  #mydata.large.type="listfile",
  fontsize.gviz =12,
tracks.gviz=listgviz,verbose=FALSE, print.image=TRUE)
dev.off()

```








