


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


# beta vs MAF

top_cpg_significant = topcpg[which(topcpg$significant==1),]
tiff(paste0(path_figure,"/mQTL_beta_vs_MAF.tiff"), units="in", width=6, height=6, res=150)
ggplot(top_cpg_significant, aes(x=af, y=abs(beta))) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=1)+
  geom_smooth(method=lm, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  labs(title=paste0("mQTL"),
       x="MAF", y = "abs(effect size)")
dev.off()

pdf(paste0(path_figure,"/mQTL_beta_vs_MAF_loess.pdf"))
ggplot(top_cpg_significant, aes(x=af, y=abs(beta))) +
  geom_point(shape=19, fill="#56B4E9", color="#56B4E9", size=1)+
  geom_smooth(method=loess, se=FALSE, fullrange=TRUE,color="tomato1")+
  theme_bw(base_size=25)+
  labs(title=paste0("mQTL"),
       x="MAF", y = "abs(effect size)")
dev.off()


# enrichment

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


```








