# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Analysis of NMD-determinant exons ----


# ******************** ----
# ******************** Proportions of poison and essential exons ----
# ******************** ----


library(DOSE)
library(igraph)
library(org.Hs.eg.db)
library(tidyverse)
library(magrittr)
library(plyr)
library(clusterProfiler)
library(msigdbr)

dexVAR=readRDS('Rdata/dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds')


pll = load('Rdata/KEGG_list.148.SPLICING_RELATED.Rdata')


fromType="SYMBOL"
toType=c("UNIPROT", "ENSEMBL", 'ENTREZID')
OrgDb="org.Hs.eg.db"
typeFormat='ENTREZID'
ids <- bitr(unique(pl$SPLICING_RELATED), fromType=fromType, toType=toType, OrgDb=OrgDb)

sum(pl$SPLICING_RELATED%in%dexVAR$symbol) # 50
sum(ids$ENSEMBL%in%dexVAR$ensemble)

dexVAR$SRG='other'
dexVAR$SRG[which(dexVAR$ensemble%in%ids$ENSEMBL)]='SRG'


tab=table(dexVAR$pois_ess,dexVAR$FOXA1_sign)

tab[1,2]/colSums(tab)[2]
tab[1,1]/colSums(tab)[1]


tab[2,2]/colSums(tab)[2]
tab[2,1]/colSums(tab)[1]


to_test_poi=t(matrix(c(tab[2,],colSums(tab[c(1,3),])),nrow = 2))
colnames(to_test_poi)=c('vASE','FOXA1_regulated')
rownames(to_test_poi)=c('POISON','OTHER')
fisher.test(to_test_poi)

to_test_ess=t(matrix(c(tab[1,],colSums(tab[c(2,3),])),nrow = 2))
colnames(to_test_ess)=c('vASE','FOXA1_regulated')
rownames(to_test_ess)=c('ESSENTIAL','OTHER')
fisher.test(to_test_ess)

pdf(file='pois_ess_enrichment_TCGA.pdf', height = unit(5,'cm'), width = unit(6, 'cm'), useDingbats = F)
barplot(c(tab[2,2]/colSums(tab)[2],
          tab[2,1]/colSums(tab)[1],
          tab[1,1]/colSums(tab)[1],
          tab[1,2]/colSums(tab)[2]),col = c('darkgreen','white','white','darkgreen'))
dev.off()


# SRGs only ----

tab=table(subset(dexVAR,SRG=='SRG')$pois_ess,subset(dexVAR,SRG=='SRG')$FOXA1_sign)

tab[1,2]/colSums(tab)[2]
tab[1,1]/colSums(tab)[1]


tab[2,2]/colSums(tab)[2]
tab[2,1]/colSums(tab)[1]


to_test_poi=t(matrix(c(tab[2,],colSums(tab[c(1,3),])),nrow = 2))
colnames(to_test_poi)=c('vASE','FOXA1_regulated')
rownames(to_test_poi)=c('POISON','OTHER')
fisher.test(to_test_poi)

to_test_ess=t(matrix(c(tab[1,],colSums(tab[c(2,3),])),nrow = 2))
colnames(to_test_ess)=c('vASE','FOXA1_regulated')
rownames(to_test_ess)=c('ESSENTIAL','OTHER')
fisher.test(to_test_ess)

pdf(file='pois_ess_enrichment_in_SRG_only_TCGA.pdf', height = unit(5,'cm'), width = unit(6, 'cm'), useDingbats = F)
barplot(c(tab[2,2]/colSums(tab)[2],
          tab[2,1]/colSums(tab)[1],
          tab[1,1]/colSums(tab)[1],
          tab[1,2]/colSums(tab)[2]),col = c('darkgreen','white','white','darkgreen'))
dev.off()


# SRGs only - ENHANCED ----

tab=table(subset(dexVAR,SRG=='SRG' & delta.mean>0)$pois_ess,subset(dexVAR,SRG=='SRG' & delta.mean>0)$FOXA1_sign)

tab[1,2]/colSums(tab)[2]
tab[1,1]/colSums(tab)[1]


tab[2,2]/colSums(tab)[2]
tab[2,1]/colSums(tab)[1]


to_test_poi=t(matrix(c(tab[2,],colSums(tab[c(1,3),])),nrow = 2))
colnames(to_test_poi)=c('vASE','FOXA1_regulated')
rownames(to_test_poi)=c('POISON','OTHER')
fisher.test(to_test_poi)

to_test_ess=t(matrix(c(tab[1,],colSums(tab[c(2,3),])),nrow = 2))
colnames(to_test_ess)=c('vASE','FOXA1_regulated')
rownames(to_test_ess)=c('ESSENTIAL','OTHER')
fisher.test(to_test_ess)

pdf(file='pois_ess_ENHANCED_enrichment_in_SRG_only_TCGA.pdf', height = unit(5,'cm'), width = unit(6, 'cm'), useDingbats = F)
barplot(c(tab[2,2]/colSums(tab)[2],
          tab[2,1]/colSums(tab)[1],
          tab[1,1]/colSums(tab)[1],
          tab[1,2]/colSums(tab)[2]),col = c('darkgreen','white','white','darkgreen'))
dev.off()


# SRGs only - SILENCED ----

tab=table(subset(dexVAR,SRG=='SRG' & delta.mean<0)$pois_ess,subset(dexVAR,SRG=='SRG' & delta.mean<0)$FOXA1_sign)

tab[1,2]/colSums(tab)[2]
tab[1,1]/colSums(tab)[1]


tab[2,2]/colSums(tab)[2]
tab[2,1]/colSums(tab)[1]


to_test_poi=t(matrix(c(tab[2,],colSums(tab[c(1,3),])),nrow = 2))
colnames(to_test_poi)=c('vASE','FOXA1_regulated')
rownames(to_test_poi)=c('POISON','OTHER')
fisher.test(to_test_poi)

to_test_ess=t(matrix(c(tab[1,],colSums(tab[c(2,3),])),nrow = 2))
colnames(to_test_ess)=c('vASE','FOXA1_regulated')
rownames(to_test_ess)=c('ESSENTIAL','OTHER')
fisher.test(to_test_ess)

pdf(file='pois_ess_SILENCED_enrichment_in_SRG_only_TCGA.pdf', height = unit(5,'cm'), width = unit(6, 'cm'), useDingbats = F)
barplot(c(tab[2,2]/colSums(tab)[2],
          tab[2,1]/colSums(tab)[1],
          tab[1,1]/colSums(tab)[1],
          tab[1,2]/colSums(tab)[2]),col = c('darkgreen','white','white','darkgreen'))
dev.off()




# ******************** ----
# ******************** Delta mean PSI distributions ----
# ******************** ----


# *** TCGA ----

library(ggplot2)
library(ggpubr)
library(reshape2)

# NMD knockdown annotation on TCGA ----

dexVAR=readRDS('Rdata/dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds')
dexVAR$key=paste0(dexVAR$pois_ess,'_',dexVAR$FOXA1_sign)
tmp=subset(dexVAR,key%in%c('ESSENTIAL_FALSE','ESSENTIAL_TRUE','POISON_FALSE','POISON_TRUE'))
tmp$key=factor(tmp$key,levels = c('ESSENTIAL_TRUE','ESSENTIAL_FALSE','POISON_FALSE','POISON_TRUE'))
my_comparisons <- list( c("ESSENTIAL_TRUE", "ESSENTIAL_FALSE"), c("ESSENTIAL_TRUE", "POISON_FALSE"), c("ESSENTIAL_TRUE", "POISON_TRUE"),
                        c("ESSENTIAL_FALSE","POISON_FALSE"), c("ESSENTIAL_FALSE","POISON_TRUE"),
                        c("POISON_FALSE","POISON_TRUE"))
labels = as.data.frame(as.character(table(tmp$key)))
labels$key = c('ESSENTIAL_TRUE','ESSENTIAL_FALSE','POISON_FALSE','POISON_TRUE')
colnames(labels)[1]='counts'

pdf(file='boxplot_poison_essential_FOXA1_svESs_and_varESs.pdf', width = 6, height = 7, useDingbats = F)
ggplot(tmp,aes(y=delta.mean,x=key,fill=key))+geom_boxplot(notch=T)+stat_compare_means(comparisons = my_comparisons)+theme_bw()+
  geom_text(data = labels, aes(x = key, y = (-10), label = counts), size = 6)
dev.off()

pdf(file='boxplot_poison_essential_FOXA1_svESs_and_varESs_bis.pdf', width = 6, height = 7, useDingbats = F)
ggplot(tmp,aes(y=delta.mean,x=key,fill=key))+geom_boxplot(notch=T)+theme_bw()
dev.off()



# *** Cell lines ----


# GeneStructureTools annotation on PC3 NOVASEQ ----

delta = readRDS('Rdata/Whippet_PC3_nextseq_novaseq_filter_0_2.rds')

risultati = load('Rdata/NMD_analysis_GeneStructureTools_PC3.Rdata')

sensiCE = subset(delta,Type=='CE' & !Complexity=='K0')
sensiCE$nmd_prob_normalTrans = orfChange$nmd_prob_bygroup_x[match(sensiCE$Coord,orfChange$id)]
sensiCE$nmd_prob_skippedExonTrans = orfChange$nmd_prob_bygroup_y[match(sensiCE$Coord,orfChange$id)]

sensiCE$delta_nmd_prob = sensiCE$nmd_prob_normalTrans-sensiCE$nmd_prob_skippedExonTrans
sensiCE$new_pois_ess = 'REST'
qt=quantile(sensiCE$delta_nmd_prob,na.rm=T,seq(0,1,0.05))
sensiCE$new_pois_ess[which(sensiCE$delta_nmd_prob<qt['15%'])] = 'ESSENTIAL'
sensiCE$new_pois_ess[which(sensiCE$delta_nmd_prob>qt['85%'])] = 'POISON'

sensiCE$significant = F
sensiCE$significant[which(sensiCE$key%in%subset(sensiCE,abs(DeltaPsi)>0.05 & Probability>0.9)$key)] = T

sensiCE$code=paste0(sensiCE$new_pois_ess,'_',sensiCE$significant)
tmp=subset(sensiCE,code%in%c('ESSENTIAL_FALSE','ESSENTIAL_TRUE','POISON_FALSE','POISON_TRUE'))
tmp$code=factor(tmp$code,levels = c('ESSENTIAL_TRUE','ESSENTIAL_FALSE','POISON_FALSE','POISON_TRUE'))
my_comparisons <- list( c("ESSENTIAL_TRUE", "ESSENTIAL_FALSE"), c("ESSENTIAL_TRUE", "POISON_FALSE"), c("ESSENTIAL_TRUE", "POISON_TRUE"),
                        c("ESSENTIAL_FALSE","POISON_FALSE"), c("ESSENTIAL_FALSE","POISON_TRUE"),
                        c("POISON_FALSE","POISON_TRUE"))
labels = as.data.frame(as.character(table(tmp$code)))
labels$code = c('ESSENTIAL_TRUE','ESSENTIAL_FALSE','POISON_FALSE','POISON_TRUE')
colnames(labels)[1]='counts'

pdf(file='Boxplot_distribution_deltaPSI_GeneStructureTools_POISON_ESSENTIAL__PC3.pdf', width = 6, height = 7, useDingbats = F)
ggplot(tmp,aes(y=DeltaPsi,x=code,fill=code))+geom_boxplot(notch=T)+stat_compare_means(comparisons = my_comparisons)+theme_bw()+
  geom_text(data = labels, aes(x = code, y = (-0.30), label = counts), size = 6)
dev.off()

pdf(file='Boxplot_distribution_deltaPSI_GeneStructureTools_POISON_ESSENTIAL__PC3_bis.pdf', width = 6, height = 7, useDingbats = F)
ggplot(tmp,aes(y=DeltaPsi,x=code,fill=code))+geom_boxplot(notch=T)+theme_bw()
dev.off()


# GeneStructureTools annotation on VCaP ----

delta = readRDS('Rdata/Whippet_VCAP_filter_0_2.rds')

risultati = load('Rdata/NMD_analysis_GeneStructureTools_SENSITIVE_VCaP.Rdata')
tmp = orfChange
risultati = load('Rdata/NMD_analysis_GeneStructureTools_COMPLEMENTARY_TO_SENSITIVE_VCaP.Rdata')
orfChange = rbind(orfChange,tmp)

sensiCE = subset(delta,Type=='CE' & !Complexity=='K0')
sensiCE$nmd_prob_normalTrans = orfChange$nmd_prob_bygroup_x[match(sensiCE$Coord,orfChange$id)]
sensiCE$nmd_prob_skippedExonTrans = orfChange$nmd_prob_bygroup_y[match(sensiCE$Coord,orfChange$id)]

sensiCE$delta_nmd_prob = sensiCE$nmd_prob_normalTrans-sensiCE$nmd_prob_skippedExonTrans
sensiCE$new_pois_ess = 'REST'
qt=quantile(sensiCE$delta_nmd_prob,na.rm=T,seq(0,1,0.05))
sensiCE$new_pois_ess[which(sensiCE$delta_nmd_prob<qt['15%'])] = 'ESSENTIAL'
sensiCE$new_pois_ess[which(sensiCE$delta_nmd_prob>qt['85%'])] = 'POISON'

sensiCE$significant = F
sensiCE$significant[which(sensiCE$key%in%subset(sensiCE,abs(DeltaPsi)>0.05 & Probability>0.9)$key)] = T

sensiCE$code=paste0(sensiCE$new_pois_ess,'_',sensiCE$significant)
tmp=subset(sensiCE,code%in%c('ESSENTIAL_FALSE','ESSENTIAL_TRUE','POISON_FALSE','POISON_TRUE'))
tmp$code=factor(tmp$code,levels = c('ESSENTIAL_TRUE','ESSENTIAL_FALSE','POISON_FALSE','POISON_TRUE'))
my_comparisons <- list( c("ESSENTIAL_TRUE", "ESSENTIAL_FALSE"), c("ESSENTIAL_TRUE", "POISON_FALSE"), c("ESSENTIAL_TRUE", "POISON_TRUE"),
                        c("ESSENTIAL_FALSE","POISON_FALSE"), c("ESSENTIAL_FALSE","POISON_TRUE"),
                        c("POISON_FALSE","POISON_TRUE"))
labels = as.data.frame(as.character(table(tmp$code)))
labels$code = c('ESSENTIAL_TRUE','ESSENTIAL_FALSE','POISON_FALSE','POISON_TRUE')
colnames(labels)[1]='counts'

pdf(file='Boxplot_distribution_deltaPSI_GeneStructureTools_POISON_ESSENTIAL__VCaP.pdf', width = 6, height = 7, useDingbats = F)
ggplot(tmp,aes(y=DeltaPsi,x=code,fill=code))+geom_boxplot(notch=T)+stat_compare_means(comparisons = my_comparisons)+theme_bw()+
  geom_text(data = labels, aes(x = code, y = (-0.30), label = counts), size = 6)
dev.off()

pdf(file='Boxplot_distribution_deltaPSI_GeneStructureTools_POISON_ESSENTIAL__VCaP_bis.pdf', width = 6, height = 7, useDingbats = F)
ggplot(tmp,aes(y=DeltaPsi,x=code,fill=code))+geom_boxplot(notch=T)+theme_bw()
dev.off()






