#!/usr/bin/env Rscript --vanilla


message('[****] SCRIPT: Splicing_analysis_NMD_determinant_exons_annotations.R')

print('[****] SCRIPT: Splicing_analysis_NMD_determinant_exons_annotations.R')


message('[*] Setting working path and creating directories ...')

repo = '/Volumes/Prule/repo/FOXA1_public/'
setwd(repo)

gene = 'FOXA1'

print(paste0('STRATA: ',gene))

work_dir     ='~/Downloads/FOXA1/Automatic_splicing_analysis_results/'
STRATA_code  = paste0(gene,'_HE') # 'PTEN_loss'
work_dir     = paste0(work_dir,STRATA_code)
Fig_out_dir  = paste0(work_dir,'/Figures/splicing_analysis_',STRATA_code)
Data_out_dir = paste0(work_dir,'/Rdata/splicing_analysis_',STRATA_code)


message('[*] Loading sources and libraries ...')

source('sources/config.R')
source('sources/config2.R')

library(reshape2)
library(plyr)
library(gtools)
library(ggpubr)


# POISON EXONS ----

message('[*] Annotating dexSIGN_ES with NMD...')

whole_coord_EXONS_unique <- readRDS("Rdata/whole_coord_EXONS_unique.rds")

dexSIGN_ES <- readRDS(paste0(Data_out_dir,'/dexSIGN_EXON_SKIPPING_only_',STRATA_code,'.rds'))
# dexSIGN_ES <- readRDS(paste0(repo, 'Rdata/dexSIGN_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds'))


dexSIGN_ES$strand = whole_coord_EXONS_unique$strand[match(dexSIGN_ES$gene_id,whole_coord_EXONS_unique$GeneName)]


dexSIGN_ES$start = sapply(strsplit(as.character(dexSIGN_ES$alt_region_coordinates),':'), `[`, 1)
dexSIGN_ES$end = sapply(strsplit(as.character(dexSIGN_ES$alt_region_coordinates),':'), `[`, 2)

dexSIGN_bed=as.data.frame(cbind(as.character(dexSIGN_ES$event_chr),as.character(dexSIGN_ES$start),as.character(dexSIGN_ES$end),as.character(dexSIGN_ES$event_id),as.character(dexSIGN_ES$overall.mean),as.character(dexSIGN_ES$strand)))

write.table(dexSIGN_bed,paste0(Data_out_dir,'/',STRATA_code,'_sign_regulated_ES_events.bed'),quote = F, sep = "\t", row.names = F, col.names = F)

comando = paste0('bedtools intersect -wa -a ',Data_out_dir,'/',STRATA_code,'_sign_regulated_ES_events.bed',' -b sources/bed_files/nmd_poison.bed > ',Data_out_dir,'/',STRATA_code,'_ES_nmd_POISON_intersection.bed')
system(comando)
comando = paste0('bedtools intersect -wa -a ',Data_out_dir,'/',STRATA_code,'_sign_regulated_ES_events.bed',' -b sources/bed_files/nmd_essential.bed > ',Data_out_dir,'/',STRATA_code,'_ES_nmd_ESSENTIAL_intersection.bed')
system(comando)

signES_poison=read.table(paste0(Data_out_dir,'/',STRATA_code,'_ES_nmd_POISON_intersection.bed'))
signES_essential=read.table(paste0(Data_out_dir,'/',STRATA_code,'_ES_nmd_ESSENTIAL_intersection.bed'))


length(unique(signES_poison$V4)) # 266
length(unique(signES_essential$V4)) # 277

dexSIGN_ES$POISON=F
dexSIGN_ES$POISON[which(dexSIGN_ES$event_id%in%signES_poison$V4)]=T
dexSIGN_ES$ESSENTIAL=F
dexSIGN_ES$ESSENTIAL[which(dexSIGN_ES$event_id%in%signES_essential$V4)]=T


message('[*] Annotating dexVAR_ES with NMD...')

dexVAR_ES <- readRDS(paste0(Data_out_dir,'/dexVAR_EXON_SKIPPING_only_',STRATA_code,'.rds'))
# dexVAR_ES <- readRDS(paste0(repo, 'Rdata/dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds'))

whole_coord_EXONS_unique <- readRDS("Rdata/whole_coord_EXONS_unique.rds")

dexVAR_ES$start = sapply(strsplit(as.character(dexVAR_ES$alt_region_coordinates),':'), `[`, 1)
dexVAR_ES$end = sapply(strsplit(as.character(dexVAR_ES$alt_region_coordinates),':'), `[`, 2)

dexVAR_bed=as.data.frame(cbind(as.character(dexVAR_ES$event_chr),as.character(dexVAR_ES$start),as.character(dexVAR_ES$end),as.character(dexVAR_ES$event_id),as.character(dexVAR_ES$overall.mean),as.character(dexVAR_ES$strand)))


write.table(dexVAR_bed,paste0(Data_out_dir,'/',STRATA_code,'_variable_ES_events.bed'),quote = F, sep = "\t", row.names = F, col.names = F)


comando = paste0('bedtools intersect -wa -a ',Data_out_dir,'/',STRATA_code,'_variable_ES_events.bed',' -b sources/bed_files/nmd_poison.bed > ',Data_out_dir,'/',STRATA_code,'_variable_ES_nmd_POISON_intersection.bed')
system(comando)
comando = paste0('bedtools intersect -wa -a ',Data_out_dir,'/',STRATA_code,'_variable_ES_events.bed',' -b sources/bed_files/nmd_essential.bed > ',Data_out_dir,'/',STRATA_code,'_variable_ES_nmd_ESSENTIAL_intersection.bed')
system(comando)


varES_poison=read.table(paste0(Data_out_dir,'/',STRATA_code,'_variable_ES_nmd_POISON_intersection.bed'))
varES_essential=read.table(paste0(Data_out_dir,'/',STRATA_code,'_variable_ES_nmd_ESSENTIAL_intersection.bed'))


length(unique(varES_poison$V4)) # 266
length(unique(varES_essential$V4)) # 277

dexVAR_ES$POISON=F
dexVAR_ES$POISON[which(dexVAR_ES$event_id%in%varES_poison$V4)]=T
dexVAR_ES$ESSENTIAL=F
dexVAR_ES$ESSENTIAL[which(dexVAR_ES$event_id%in%varES_essential$V4)]=T
dexVAR_ES$FOXA1_sign=F
dexVAR_ES$FOXA1_sign[which(dexVAR_ES$event_id%in%dexSIGN_ES$event_id)]=T


# FISHER TEST TO VERIFY THE ENRICHMENT OF POISON EXONS AMONG FOXA1-REGULATED EXONS WITH RESPECT TO VARIABLE ONES ---- 

message('[*] Computing enrichments...')

totest=matrix(c(nrow(subset(dexVAR_ES,event_type=='ES' & FOXA1_sign & POISON)),
                nrow(subset(dexVAR_ES,event_type=='ES' & FOXA1_sign & !POISON)),
                nrow(subset(dexVAR_ES,event_type=='ES' & !FOXA1_sign & POISON)),
                nrow(subset(dexVAR_ES,event_type=='ES' & !FOXA1_sign & !POISON))),nrow=2)
rownames(totest)=c('POISON','NO_POISON')
colnames(totest)=c('SIGN','VAR')
totest
totest[1,1]/sum(totest[,1])
totest[1,2]/sum(totest[,2])
pval_towrite=fisher.test(totest)$p.val
ypval=max(totest[1,1]/sum(totest[,1]),totest[1,2]/sum(totest[,2]))

pdf(file=paste0(Fig_out_dir,'/barplot_POISON_SIGN_and_VAR_',STRATA_code,'.pdf'), width = 4, height = 5, useDingbats = F)
barplot(c('svES events'=totest[1,1]/sum(totest[,1]),'varES events'=totest[1,2]/sum(totest[,2])),ylim = c(0,(ypval+0.05)))
text(1.2,ypval,as.character(round(pval_towrite,4)),pos=3,col = 'red')
dev.off()
fisher.test(totest)


# CONSIDER THE DELTA MEAN SIGN ---

totest=matrix(c(nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean<0 & FOXA1_sign & POISON)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean<0 & FOXA1_sign & !POISON)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean<0 & !FOXA1_sign & POISON)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean<0 & !FOXA1_sign & !POISON))),nrow=2)
rownames(totest)=c('POISON','NO_POISON')
colnames(totest)=c('SIGN','VAR')
totest
totest[1,1]/sum(totest[,1])
totest[1,2]/sum(totest[,2])
pval_towrite=fisher.test(totest)$p.val
ypval=max(totest[1,1]/sum(totest[,1]),totest[1,2]/sum(totest[,2]))
pdf(file=paste0(Fig_out_dir,'/barplot_POISON_SIGN_and_VAR_Negative_Delta_',STRATA_code,'.pdf'), width = 4, height = 5, useDingbats = F)
barplot(c('svES events'=totest[1,1]/sum(totest[,1]),'varES events'=totest[1,2]/sum(totest[,2])),ylim = c(0,(ypval+0.05)))
text(1.2,ypval,as.character(round(pval_towrite,4)),pos=3,col = 'red')
dev.off()
fisher.test(totest)


totest=matrix(c(nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean>0 & FOXA1_sign & POISON)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean>0 & FOXA1_sign & !POISON)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean>0 & !FOXA1_sign & POISON)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean>0 & !FOXA1_sign & !POISON))),nrow=2)
rownames(totest)=c('POISON','NO_POISON')
colnames(totest)=c('SIGN','VAR')
totest
totest[1,1]/sum(totest[,1])
totest[1,2]/sum(totest[,2])
pval_towrite=fisher.test(totest)$p.val
ypval=max(totest[1,1]/sum(totest[,1]),totest[1,2]/sum(totest[,2]))
pdf(file=paste0(Fig_out_dir,'/barplot_POISON_SIGN_and_VAR_Positive_Delta_',STRATA_code,'.pdf'), width = 4, height = 5, useDingbats = F)
barplot(c('svES events'=totest[1,1]/sum(totest[,1]),'varES events'=totest[1,2]/sum(totest[,2])),ylim = c(0,(ypval+0.05)))
text(1.2,ypval,as.character(round(pval_towrite,4)),pos=3,col = 'red')
dev.off()
fisher.test(totest)



# FISHER TEST TO VERIFY THE ENRICHMENT OF ESSENTIAL AND ESSENTIAL EXONS AMONG FOXA1-REGULATED EXONS WITH RESPECT TO VARIABLE ONES ---- 

totest=matrix(c(nrow(subset(dexVAR_ES,event_type=='ES' & FOXA1_sign & ESSENTIAL)),
                nrow(subset(dexVAR_ES,event_type=='ES' & FOXA1_sign & !ESSENTIAL)),
                nrow(subset(dexVAR_ES,event_type=='ES' & !FOXA1_sign & ESSENTIAL)),
                nrow(subset(dexVAR_ES,event_type=='ES' & !FOXA1_sign & !ESSENTIAL))),nrow=2)
rownames(totest)=c('ESSENTIAL','NO_ESSENTIAL')
colnames(totest)=c('SIGN','VAR')
totest
totest[1,1]/sum(totest[,1])
totest[1,2]/sum(totest[,2])
pval_towrite=fisher.test(totest)$p.val
ypval=max(totest[1,1]/sum(totest[,1]),totest[1,2]/sum(totest[,2]))
pdf(file=paste0(Fig_out_dir,'/barplot_ESSENTIAL_SIGN_and_VAR_',STRATA_code,'.pdf'), width = 4, height = 5, useDingbats = F)
barplot(c('svES events'=totest[1,1]/sum(totest[,1]),'varES events'=totest[1,2]/sum(totest[,2])),ylim = c(0,(ypval+0.05)))
text(1.2,ypval,as.character(round(pval_towrite,4)),pos=3,col = 'red')
dev.off()
fisher.test(totest)


# CONSIDER DELTA MEAN SIGN ---

totest=matrix(c(nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean<0 & FOXA1_sign & ESSENTIAL)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean<0 & FOXA1_sign & !ESSENTIAL)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean<0 & !FOXA1_sign & ESSENTIAL)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean<0 & !FOXA1_sign & !ESSENTIAL))),nrow=2)
rownames(totest)=c('ESSENTIAL','NO_ESSENTIAL')
colnames(totest)=c('SIGN','VAR')
totest
totest[1,1]/sum(totest[,1])
totest[1,2]/sum(totest[,2])
pval_towrite=fisher.test(totest)$p.val
ypval=max(totest[1,1]/sum(totest[,1]),totest[1,2]/sum(totest[,2]))
pdf(file=paste0(Fig_out_dir,'/barplot_ESSENTIAL_SIGN_and_VAR_Negative_Delta_',STRATA_code,'.pdf'), width = 4, height = 5, useDingbats = F)
barplot(c('svES events'=totest[1,1]/sum(totest[,1]),'varES events'=totest[1,2]/sum(totest[,2])),ylim = c(0,(ypval+0.05)))
text(1.2,ypval,as.character(round(pval_towrite,4)),pos=3,col = 'red')
dev.off()
fisher.test(totest)


totest=matrix(c(nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean>0 & FOXA1_sign & ESSENTIAL)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean>0 & FOXA1_sign & !ESSENTIAL)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean>0 & !FOXA1_sign & ESSENTIAL)),
                nrow(subset(dexVAR_ES,event_type=='ES' & delta.mean>0 & !FOXA1_sign & !ESSENTIAL))),nrow=2)
rownames(totest)=c('ESSENTIAL','NO_ESSENTIAL')
colnames(totest)=c('SIGN','VAR')
totest
totest[1,1]/sum(totest[,1])
totest[1,2]/sum(totest[,2])
pval_towrite=fisher.test(totest)$p.val
ypval=max(totest[1,1]/sum(totest[,1]),totest[1,2]/sum(totest[,2]))
pdf(file=paste0(Fig_out_dir,'/barplot_ESSENTIAL_SIGN_and_VAR_Positive_Delta_',STRATA_code,'.pdf'), width = 4, height = 5, useDingbats = F)
barplot(c('svES events'=totest[1,1]/sum(totest[,1]),'varES events'=totest[1,2]/sum(totest[,2])),ylim = c(0,(ypval+0.05)))
text(1.2,ypval,as.character(round(pval_towrite,4)),pos=3,col = 'red')
dev.off()
fisher.test(totest)



# Save dexSIGN_ES and dexVAR_ES ----

message('[*] Saving datasets...')

saveRDS(dexSIGN_ES, file=paste0(Data_out_dir,'/dexSIGN_EXON_SKIPPING_only_with_pfam_and_NMD_',STRATA_code,'.rds'))
saveRDS(dexVAR_ES, file=paste0(Data_out_dir,'/dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_',STRATA_code,'.rds'))


message('[*] Done')
