
#### molocalization and colocalization

```R

args <- as.numeric(commandArgs(TRUE))
traitnum = args[1]
count = args[2]
print(traitnum)
print(count)


#########################################

# load data

#########################################


load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/mcoloc_prepare_dat_",count,".RData"))
coloc_dat=datt
chrnum = coloc_dat$coloc_input[1,][1,1]

if(traitnum==1){
	load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/SBP_trait_chr_",chrnum,".RData"))
	myGWAS_summary=SBP_trait_chr
	trait_type="quant"
	GWASsamplesize = 29378

}else if(traitnum==2){
	load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/BMI_trait_chr_",chrnum,".RData"))
	myGWAS_summary=BMI_trait_chr
	trait_type="quant"
	GWASsamplesize = 23095

}else if(traitnum==3){
	load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/DBP_trait_chr_",chrnum,".RData"))
	myGWAS_summary=DBP_trait_chr
	trait_type="quant"
	GWASsamplesize = 29378

}else if(traitnum==4){
	load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/HTN_trait_chr_",chrnum,".RData"))
	myGWAS_summary=HTN_trait_chr
	trait_type="quant"
	GWASsamplesize = 29378

}else if(traitnum==5){
	load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/PP_trait_chr_",chrnum,".RData"))
	myGWAS_summary=PP_trait_chr
	trait_type="quant"
	GWASsamplesize = 29378

}else if(traitnum==6){
	load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/T2D_trait_chr_",chrnum,".RData"))
	myGWAS_summary=T2D_trait_chr
	trait_type="cc"
	GWASsamplesize = 23827
	s=0.3476728
}


path_mQTL = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu"
path_eQTL = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc"
path_analysis = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/analysis"
path_mediation = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/mediation"
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/annodata.RData")
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/CpG_gene_table.RData")
library(plyr)
library(dplyr)
#load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/gene_AA_anno_order_correct.RData")
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/CpG_gene_table.RData")
thr_mQTL = 0.0002267195
thr_eQTL = 6.245907e-05


#########################################

# load codes

#########################################


library(data.table)
library(coloc)
library(moloc)
library(susieR)

source("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/coloc_susie_package_code.R")


#########################################

# run

#########################################



if(coloc_dat$SNPnum>0){

	coloc_dat$coloc_input$p_wald_eQTL[which(coloc_dat$coloc_input$p_wald_eQTL==0)] = 1e-300
	coloc_dat$coloc_input$p_wald_mQTL[which(coloc_dat$coloc_input$p_wald_mQTL==0)] = 1e-300
	coloc_dat$coloc_input$z_eqtl = coloc_dat$coloc_input$beta_eQTL/coloc_dat$coloc_input$se_eQTL
	coloc_dat$coloc_input$N_eqtl = 1032
	eQTL_prep = coloc_dat$coloc_input[,c(2,8,9,26,11,27,7)]
	coloc_dat$coloc_input$z_mqtl = coloc_dat$coloc_input$beta_mQTL/coloc_dat$coloc_input$se_mQTL
	coloc_dat$coloc_input$N_mqtl = 961
	mQTL_prep = coloc_dat$coloc_input[,c(2,19,20,28,21,29,18)]

	colnames(eQTL_prep) = c( "SNP", "BETA" , "SE" ,"z" , "PVAL" ,"N" , "MAF")
	colnames(mQTL_prep) = c( "SNP", "BETA" , "SE" ,"z" , "PVAL" ,"N" , "MAF")

	ind_match_gwas = na.omit(match(as.character(eQTL_prep$SNP), myGWAS_summary$ID))

	if(length(ind_match_gwas)>0){

		print(length(ind_match_gwas))
		myGWAS_temp = myGWAS_summary[ind_match_gwas,]
		coloc_temp = coloc_dat$coloc_input[match(myGWAS_temp$ID,as.character(coloc_dat$coloc_input$rs_eQTL),),]

		ind_ALLELE = which(myGWAS_temp$Allele2 == coloc_temp$allele1_eQTL)
		myGWAS_temp$Freq1[ind_ALLELE] = 1-myGWAS_temp$Freq1[ind_ALLELE] 
		myGWAS_temp$Effect[ind_ALLELE] = -myGWAS_temp$Effect[ind_ALLELE]
		myGWAS_temp$z[ind_ALLELE] = -myGWAS_temp$z[ind_ALLELE]
		GWAS_summary = myGWAS_temp[,c(1,7,8,11,9,10,6)]
		colnames(GWAS_summary) = c( "SNP", "BETA" , "SE" ,"z" , "PVAL" ,"N" , "MAF")
		eQTL_summary = eQTL_prep[match(GWAS_summary$SNP, eQTL_prep$SNP),]
		mQTL_summary = mQTL_prep[match(GWAS_summary$SNP, mQTL_prep$SNP),]

		data_single = list(GWAS_summary,eQTL_summary,mQTL_summary )
		moloc <- moloc_test(data_single, prior_var=c(0.01, 0.1, 0.5), priors=c(1e-04, 1e-06, 1e-07))
		print(moloc$priors_lkl_ppa[14,])
		PPA=moloc$priors_lkl_ppa[14,4]

		#------------ coloc, eqtl/mqtl vs gwas

		coloc_dat$coloc_input$p_wald_eQTL[which(coloc_dat$coloc_input$p_wald_eQTL==0)] = 1e-300
		coloc_dat$coloc_input$p_wald_mQTL[which(coloc_dat$coloc_input$p_wald_mQTL==0)] = 1e-300
		dataset1 = list("beta"=unlist(eQTL_summary$BETA),"varbeta" = unlist(eQTL_summary$SE)^2,"snp"= as.character(unlist(eQTL_summary$SNP)),"pvalues"=unlist(eQTL_summary$PVAL) ,"sdY"=1,type="quant","N"=1032)
		dataset2 = list("beta"=unlist(mQTL_summary$BETA),"varbeta" = unlist(mQTL_summary$SE)^2,"snp"= as.character(unlist(mQTL_summary$SNP)),"pvalues"=unlist(mQTL_summary$PVAL) ,"sdY"=1,type="quant","N"=961)
		
		if(traitnum !=6){
			dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type, "N"=GWASsamplesize) 
		}else if(traitnum==6){
			dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type, "N"=GWASsamplesize,"s"=s) 
		}
		
		
		# coloc.abf:
		res_eqtl_gwas <- coloc.abf(dataset1, dataset3,MAF=eQTL_summary$MAF) # use default setting for p1, p2, p12
		res_mqtl_gwas <- coloc.abf(dataset2, dataset3,MAF=mQTL_summary$MAF) # use default setting for p1, p2, p12
		res_eqtl_mqtl_default <- coloc.abf(dataset1, dataset2,MAF=eQTL_summary$MAF) # use default setting for p1, p2, p12


		percent_p12 = 0.75
		p12=0.0008403957*percent_p12
		p1=0.0008403957-p12
		p2=0.03115616-p12
		res_eqtl_mqtl_guo <- coloc.abf(dataset1, dataset2,MAF=eQTL_summary$MAF,p1,p2,p12) # use guo et al. select p1,p2,p12
		PP4_res_eqtl_gwas = res_eqtl_gwas$summary
		PP4_res_mqtl_gwas = res_mqtl_gwas$summary
		PP4_res_eqtl_mqtl_default = res_eqtl_mqtl_default$summary
		PP4_res_eqtl_mqtl_guo = res_eqtl_mqtl_guo$summary

		# coloc.susie
		dataset1 = list("beta"=unlist(eQTL_summary$BETA),"varbeta" = unlist(eQTL_summary$SE)^2,"snp"= as.character(unlist(eQTL_summary$SNP)),"pvalues"=unlist(eQTL_summary$PVAL) ,"sdY"=1,type="quant")
		dataset2 = list("beta"=unlist(mQTL_summary$BETA),"varbeta" = unlist(mQTL_summary$SE)^2,"snp"= as.character(unlist(mQTL_summary$SNP)),"pvalues"=unlist(mQTL_summary$PVAL) ,"sdY"=1,type="quant")
		# dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type)
		if(traitnum !=6){
			dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type, "N"=GWASsamplesize) 
		}else if(traitnum==6){
			dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type, "N"=GWASsamplesize,"s"=s) 
		}
		LD=cor(datt$geno[,match(unlist(as.character(eQTL_summary$SNP)),colnames(datt$geno))])
		dataset1$LD=LD
		dataset2$LD=LD
		dataset3$LD=LD

		S1=runsusie(dataset1)
		S2=runsusie(dataset2)
		S3=runsusie(dataset3)

		if(sum(dim(summary(S1)$cs)[1]>0) & sum(dim(summary(S3)$cs)[1]>0)){
			res_eqtl_gwas=coloc.susie(S1,S3)
			PP4_res_eqtl_gwas_susie = res_eqtl_gwas$summary
			print(PP4_res_eqtl_gwas_susie)
		}else{
			PP4_res_eqtl_gwas_susie=NA
		}

		if(sum(dim(summary(S2)$cs)[1]>0) & sum(dim(summary(S3)$cs)[1]>0)){
			res_mqtl_gwas=coloc.susie(S2,S3)
			PP4_res_mqtl_gwas_susie = res_mqtl_gwas$summary
			print(PP4_res_mqtl_gwas_susie)
		}else{
			PP4_res_mqtl_gwas_susie=NA
		}

		if(sum(dim(summary(S1)$cs)[1]>0) & sum(dim(summary(S2)$cs)[1]>0)){
			res_eqtl_mqtl=coloc.susie(S1,S2)
			PP4_res_eqtk_mqtl_susie = res_eqtl_mqtl$summary
			print(PP4_res_eqtk_mqtl_susie)
		}else{
			PP4_res_eqtk_mqtl_susie=NA
		}

		gene_names = coloc_dat$eGene[1]
		cpg_names = coloc_dat$mCPG
		SNP_names = coloc_dat$eSNP
		ind_tmp = which(coloc_dat[[1]]$rs_eQTL==SNP_names & coloc_dat[[1]]$GENE_eQTL==gene_names & coloc_dat[[1]]$cpg_mQTL==cpg_names)
		if(length(ind_tmp)>=1){
		print(count)
		associations_eqtl = coloc_dat[[1]]$beta_eQTL[ind_tmp]
		associations_mqtl = coloc_dat[[1]]$beta_mQTL[ind_tmp]

		out = list("gene_names"=gene_names, "cpg_names"=cpg_names, "SNP_names"=SNP_names, "associations_eqtl"=associations_eqtl,"associations_mqtl"=associations_mqtl, 
									"PPA"=PPA,"PP4_res_eqtl_gwas"=PP4_res_eqtl_gwas,"PP4_res_mqtl_gwas"=PP4_res_mqtl_gwas,
									"PP4_res_eqtl_mqtl_default"=PP4_res_eqtl_mqtl_default,"PP4_res_eqtl_mqtl_guo"=PP4_res_eqtl_mqtl_guo,
									"PP4_res_eqtl_gwas_susie"=PP4_res_eqtl_gwas_susie,"PP4_res_mqtl_gwas_susie"=PP4_res_mqtl_gwas_susie,
									"PP4_res_eqtk_mqtl_susie"=PP4_res_eqtk_mqtl_susie)
		print(PP4_res_eqtl_gwas_susie)
		print(PP4_res_mqtl_gwas_susie)
		print(PP4_res_eqtk_mqtl_susie)
		print(PP4_res_eqtl_gwas)
		print(PP4_res_mqtl_gwas)
		print(PP4_res_eqtl_mqtl_default)
		print(PP4_res_eqtl_mqtl_guo)	
	}
}

save(out, file = paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/result/Trait",traitnum,"_gene_",count,".RData"))

}

###########################################################################################


path_mQTL = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu"
path_eQTL = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc"
path_analysis = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/analysis"
path_mediation = "/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/mediation"
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/annodata.RData")
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/CpG_gene_table.RData")
library(plyr)
library(dplyr)
#load("/net/mulan/home/shanglu/GENOA/data/AA/gene_AA_anno_order_correct.RData")
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/gene_AA_anno_order_correct.RData")
load("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/coloc/CpG_gene_table.RData")
thr_mQTL = 0.0002267195
thr_eQTL = 6.245907e-05


traitnum=6
		print(traitnum)
	#for(count in c( 198 , 571, 1028 ,2664 ,2679 ,2738 ,2739, 2801 ,3126 ,3478 ,3512, 3529 ,3572 ,3573 ,3646 ,3674, 3751, 3832 ,3836 ,3856 ,4253 ,4530 ,4723, 4779, 4909)){
	for(count in missing_file){
		print(count)
		load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/mcoloc_prepare_dat_",count,".RData"))
		coloc_dat=datt
		chrnum = coloc_dat$coloc_input[1,][1,1]

		if(traitnum==1){
			load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/SBP_trait_chr_",chrnum,".RData"))
			myGWAS_summary=SBP_trait_chr
			trait_type="quant"
			GWASsamplesize = 29378

		}else if(traitnum==2){
			load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/BMI_trait_chr_",chrnum,".RData"))
			myGWAS_summary=BMI_trait_chr
			trait_type="quant"
			GWASsamplesize = 23095

		}else if(traitnum==3){
			load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/DBP_trait_chr_",chrnum,".RData"))
			myGWAS_summary=DBP_trait_chr
			trait_type="quant"
			GWASsamplesize = 29378

		}else if(traitnum==4){
			load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/HTN_trait_chr_",chrnum,".RData"))
			myGWAS_summary=HTN_trait_chr
			trait_type="quant"
			GWASsamplesize = 29378

		}else if(traitnum==5){
			load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/PP_trait_chr_",chrnum,".RData"))
			myGWAS_summary=PP_trait_chr
			trait_type="quant"
			GWASsamplesize = 29378

		}else if(traitnum==6){
			load(paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/data/T2D_trait_chr_",chrnum,".RData"))
			myGWAS_summary=T2D_trait_chr
			trait_type="cc"
			GWASsamplesize = 23827
			s=0.3476728
		}
		if(coloc_dat$SNPnum>0){

			coloc_dat$coloc_input$p_wald_eQTL[which(coloc_dat$coloc_input$p_wald_eQTL==0)] = 1e-300
			coloc_dat$coloc_input$p_wald_mQTL[which(coloc_dat$coloc_input$p_wald_mQTL==0)] = 1e-300
			coloc_dat$coloc_input$z_eqtl = coloc_dat$coloc_input$beta_eQTL/coloc_dat$coloc_input$se_eQTL
			coloc_dat$coloc_input$N_eqtl = 1032
			eQTL_prep = coloc_dat$coloc_input[,c(2,8,9,26,11,27,7)]
			coloc_dat$coloc_input$z_mqtl = coloc_dat$coloc_input$beta_mQTL/coloc_dat$coloc_input$se_mQTL
			coloc_dat$coloc_input$N_mqtl = 961
			mQTL_prep = coloc_dat$coloc_input[,c(2,19,20,28,21,29,18)]

			colnames(eQTL_prep) = c( "SNP", "BETA" , "SE" ,"z" , "PVAL" ,"N" , "MAF")
			colnames(mQTL_prep) = c( "SNP", "BETA" , "SE" ,"z" , "PVAL" ,"N" , "MAF")

			ind_match_gwas = na.omit(match(as.character(eQTL_prep$SNP), myGWAS_summary$ID))

			if(length(ind_match_gwas)>0){

				print(length(ind_match_gwas))
				myGWAS_temp = myGWAS_summary[ind_match_gwas,]
				coloc_temp = coloc_dat$coloc_input[match(myGWAS_temp$ID,as.character(coloc_dat$coloc_input$rs_eQTL),),]

				ind_ALLELE = which(myGWAS_temp$Allele2 == coloc_temp$allele1_eQTL)
				myGWAS_temp$Freq1[ind_ALLELE] = 1-myGWAS_temp$Freq1[ind_ALLELE] 
				myGWAS_temp$Effect[ind_ALLELE] = -myGWAS_temp$Effect[ind_ALLELE]
				myGWAS_temp$z[ind_ALLELE] = -myGWAS_temp$z[ind_ALLELE]
				GWAS_summary = myGWAS_temp[,c(1,7,8,11,9,10,6)]
				colnames(GWAS_summary) = c( "SNP", "BETA" , "SE" ,"z" , "PVAL" ,"N" , "MAF")
				eQTL_summary = eQTL_prep[match(GWAS_summary$SNP, eQTL_prep$SNP),]
				mQTL_summary = mQTL_prep[match(GWAS_summary$SNP, mQTL_prep$SNP),]

				data_single = list(GWAS_summary,eQTL_summary,mQTL_summary )
				if(dim(GWAS_summary)[1]>1){
				moloc <- moloc_test(data_single, prior_var=c(0.01, 0.1, 0.5), priors=c(1e-04, 1e-06, 1e-07))
				print(moloc$priors_lkl_ppa[14,])
				PPA=moloc$priors_lkl_ppa[14,4]
				}else{
				PPA=NA
				}
				#------------ coloc, eqtl/mqtl vs gwas

				coloc_dat$coloc_input$p_wald_eQTL[which(coloc_dat$coloc_input$p_wald_eQTL==0)] = 1e-300
				coloc_dat$coloc_input$p_wald_mQTL[which(coloc_dat$coloc_input$p_wald_mQTL==0)] = 1e-300
				dataset1 = list("beta"=unlist(eQTL_summary$BETA),"varbeta" = unlist(eQTL_summary$SE)^2,"snp"= as.character(unlist(eQTL_summary$SNP)),"pvalues"=unlist(eQTL_summary$PVAL) ,"sdY"=1,type="quant","N"=1032)
				dataset2 = list("beta"=unlist(mQTL_summary$BETA),"varbeta" = unlist(mQTL_summary$SE)^2,"snp"= as.character(unlist(mQTL_summary$SNP)),"pvalues"=unlist(mQTL_summary$PVAL) ,"sdY"=1,type="quant","N"=961)
				
				if(traitnum !=6){
					dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type, "N"=GWASsamplesize) 
				}else if(traitnum==6){
					dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type, "N"=GWASsamplesize,"s"=s) 
				}
				
				# coloc.abf:
				res_eqtl_gwas <- coloc.abf(dataset1, dataset3,MAF=eQTL_summary$MAF) # use default setting for p1, p2, p12
				res_mqtl_gwas <- coloc.abf(dataset2, dataset3,MAF=mQTL_summary$MAF) # use default setting for p1, p2, p12
				res_eqtl_mqtl_default <- coloc.abf(dataset1, dataset2,MAF=eQTL_summary$MAF) # use default setting for p1, p2, p12

				percent_p12 = 0.75
				p12=0.0008403957*percent_p12
				p1=0.0008403957-p12
				p2=0.03115616-p12
				res_eqtl_mqtl_guo <- coloc.abf(dataset1, dataset2,MAF=eQTL_summary$MAF,p1,p2,p12) # use guo et al. select p1,p2,p12
				PP4_res_eqtl_gwas = res_eqtl_gwas$summary
				PP4_res_mqtl_gwas = res_mqtl_gwas$summary
				PP4_res_eqtl_mqtl_default = res_eqtl_mqtl_default$summary
				PP4_res_eqtl_mqtl_guo = res_eqtl_mqtl_guo$summary

				# # coloc.susie
				# dataset1 = list("beta"=unlist(eQTL_summary$BETA),"varbeta" = unlist(eQTL_summary$SE)^2,"snp"= as.character(unlist(eQTL_summary$SNP)),"pvalues"=unlist(eQTL_summary$PVAL) ,"sdY"=1,type="quant")
				# dataset2 = list("beta"=unlist(mQTL_summary$BETA),"varbeta" = unlist(mQTL_summary$SE)^2,"snp"= as.character(unlist(mQTL_summary$SNP)),"pvalues"=unlist(mQTL_summary$PVAL) ,"sdY"=1,type="quant")
				# # dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type)
				# if(traitnum !=6){
				# 	dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type, "N"=GWASsamplesize) 
				# }else if(traitnum==6){
				# 	dataset3 = list("beta"=unlist(GWAS_summary$BETA),"varbeta" = unlist(GWAS_summary$SE)^2,"snp"= as.character(unlist(GWAS_summary$SNP)),"pvalues"=unlist(GWAS_summary$PVAL) ,"sdY"=1,type=trait_type, "N"=GWASsamplesize,"s"=s) 
				# }
				# LD=cor(datt$geno[,match(unlist(as.character(eQTL_summary$SNP)),colnames(datt$geno))])
				# dataset1$LD=LD
				# dataset2$LD=LD
				# dataset3$LD=LD

				# S1=runsusie(dataset1)
				# S2=runsusie(dataset2)
				# S3=runsusie(dataset3)

				# if(sum(dim(summary(S1)$cs)[1]>0) & sum(dim(summary(S3)$cs)[1]>0)){
				# 	res_eqtl_gwas=coloc.susie(S1,S3)
				# 	PP4_res_eqtl_gwas_susie = res_eqtl_gwas$summary
				# 	print(PP4_res_eqtl_gwas_susie)
				# }else{
					PP4_res_eqtl_gwas_susie=NA
				# }

				# if(sum(dim(summary(S2)$cs)[1]>0) & sum(dim(summary(S3)$cs)[1]>0)){
				# 	res_mqtl_gwas=coloc.susie(S2,S3)
				# 	PP4_res_mqtl_gwas_susie = res_mqtl_gwas$summary
				# 	print(PP4_res_mqtl_gwas_susie)
				# }else{
					PP4_res_mqtl_gwas_susie=NA
				# }

				# if(sum(dim(summary(S1)$cs)[1]>0) & sum(dim(summary(S2)$cs)[1]>0)){
				# 	res_eqtl_mqtl=coloc.susie(S1,S2)
				# 	PP4_res_eqtk_mqtl_susie = res_eqtl_mqtl$summary
				# 	print(PP4_res_eqtk_mqtl_susie)
				# }else{
					PP4_res_eqtk_mqtl_susie=NA
				#}

				gene_names = coloc_dat$eGene[1]
				cpg_names = coloc_dat$mCPG
				SNP_names = coloc_dat$eSNP
				ind_tmp = which(coloc_dat[[1]]$rs_eQTL==SNP_names & coloc_dat[[1]]$GENE_eQTL==gene_names & coloc_dat[[1]]$cpg_mQTL==cpg_names)
				if(length(ind_tmp)>=1){
					print(count)
					associations_eqtl = coloc_dat[[1]]$beta_eQTL[ind_tmp]
					associations_mqtl = coloc_dat[[1]]$beta_mQTL[ind_tmp]

					out = list("gene_names"=gene_names, "cpg_names"=cpg_names, "SNP_names"=SNP_names, "associations_eqtl"=associations_eqtl,"associations_mqtl"=associations_mqtl, 
												"PPA"=PPA,"PP4_res_eqtl_gwas"=PP4_res_eqtl_gwas,"PP4_res_mqtl_gwas"=PP4_res_mqtl_gwas,
												"PP4_res_eqtl_mqtl_default"=PP4_res_eqtl_mqtl_default,"PP4_res_eqtl_mqtl_guo"=PP4_res_eqtl_mqtl_guo,
												"PP4_res_eqtl_gwas_susie"=PP4_res_eqtl_gwas_susie,"PP4_res_mqtl_gwas_susie"=PP4_res_mqtl_gwas_susie,
												"PP4_res_eqtk_mqtl_susie"=PP4_res_eqtk_mqtl_susie)
					print(PP4_res_eqtl_gwas_susie)
					print(PP4_res_mqtl_gwas_susie)
					print(PP4_res_eqtk_mqtl_susie)
					print(PP4_res_eqtl_gwas)
					print(PP4_res_mqtl_gwas)
					print(PP4_res_eqtl_mqtl_default)
					print(PP4_res_eqtl_mqtl_guo)	}
			}

		save(out, file = paste0("/home/skardia_lab/clubhouse/research/projects/GENOA_AA_mQTL/lulu/manuscript/v2/coloc_moloc/result/Trait",traitnum,"_gene_",count,".RData"))

	}

}




```
