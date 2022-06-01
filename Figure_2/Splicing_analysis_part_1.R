#!/usr/bin/env Rscript --vanilla

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%  Analysis of splicing variability in Prostate Cancer TCGA patients %%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 1 # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message('[****] SCRIPT: Splicing_analysis_part_1.R')

print('[****] SCRIPT: Splicing_analysis_part_1.R')


# Set working path : ----

library(R.utils)

message('[*] Setting working path and creating directories ...')

gene = 'FOXA1'

print(paste0('STRATA: ',gene))

work_dir     = 'Automatic_splicing_analysis_results/'
STRATA_code  = paste0(gene,'_HE') 
work_dir     = paste0(work_dir,STRATA_code)
Fig_out_dir  = paste0(work_dir,'/Figures/splicing_analysis_',STRATA_code)
Data_out_dir = paste0(work_dir,'/Rdata/splicing_analysis_',STRATA_code)

mkdirs(work_dir)
mkdirs(Fig_out_dir)
mkdirs(Data_out_dir)


# Sources and libraries: ----

message('[*] Loading sources and libraries ...')

source('sources/config.R')
source('sources/config2.R')

library(reshape2)
library(plyr)
library(gtools)


# Load CASE samples IDs: ----

message('[*] Loading stratification ...')

rna <- readRDS("Rdata/RNA_TPM_tumor.rds")
qt=quantile(subset(rna, symbol==gene)$TPM, seq(0,1,.05))
cases = substr(unique(subset(rna, symbol==gene & TPM>=qt["75%"])$sample),1,12)



# XXX. Adapt "ex" dataset with the information about which samples are CASES ----

message('[*] Loading ex dataset, integrating with strata info and saving ...')

ex = readRDS('Rdata/PRAD_selected_exons.rds')

ex$patient = substr(ex$variable,1,12)
ex$cases = F
ex$cases[which(ex$patient%in%cases)] = T

# Save "ex" dataset with CASES info:
saveRDS(ex,paste0(Data_out_dir,'/PRAD_selected_exons_',STRATA_code,'.rds'))


# XXX. Create "dex" dataset with statistics related to CASES ----

ex_CASE = ex
rm(ex)

library(gtools)
library(plyr)

ex_CASE$type=ex_CASE$cases


message('[*] Building dex dataset ...')

dex_CASE        = ddply(ex_CASE, .(as_id), get_stats_ex_boolean_v2, .progress = 'text')
dex_CASE$BF     = p.adjust(dex_CASE$w, 'bonferroni', n=sum(!is.na(dex_CASE$w)))
dex_CASE$BF.ks  = p.adjust(dex_CASE$k, 'bonferroni', n=sum(!is.na(dex_CASE$k)))
dex_CASE$FDR    = p.adjust(dex_CASE$w, 'fdr', n=sum(!is.na(dex_CASE$w)))
dex_CASE$FDR.ks = p.adjust(dex_CASE$k, 'fdr', n=sum(!is.na(dex_CASE$k)))
dex_CASE_flig   = ddply(ex_CASE, .(as_id), summarise, fligner=fligner.test(list(value[type],value[!type]))$p.value, .progress = 'text')

dex_CASE        = cbind(dex_CASE,dex_CASE_flig$fligner)
colnames(dex_CASE)[length(colnames(dex_CASE))] = 'fligner'

saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))

for(i in c(3:5,7:9)) dex_CASE[,i]=dex_CASE[,i]*100

dex_CASE$delta.median           = with(dex_CASE, case.median - cntr.median)
dex_CASE$delta.mean             = with(dex_CASE, case.mean - cntr.mean)
dex_CASE$fc.median              = foldchange(dex_CASE$case.median,dex_CASE$cntr.median)
dex_CASE$fc.mean                = foldchange(dex_CASE$case.mean,dex_CASE$cntr.mean)

message('[*] Saving dex dataset ...')

saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))


# XXX. Add annotations to "dex" dataset ----

message('[*] Adding cancer gene annotations to dex dataset and saving ...')

dex_CASE = readRDS(paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))
events_info = readRDS("Rdata/events_info.rds")

dex_CASE = cbind.data.frame( events_info[match(dex_CASE$as_id,events_info$event_id),1:6], dex_CASE[,2:ncol(dex_CASE)] )
rm(events_info)

# prova = readRDS("/Volumes/LaCie/marco/FOXA1/ANALISI_POST_UV2_DISASTER/Analisi_splicing_automatica/sources/NCG6/NCG6_cancergenes.rds")
# prova2 = read.delim("/Volumes/LaCie/marco/FOXA1/ANALISI_POST_UV2_DISASTER/Analisi_splicing_automatica/sources/NCG6/NCG6_tsgoncogene.tsv", stringsAsFactors = F) # non ci sono però i candidate cancer genes
# prova2$symbol[578:580]=c('SEPT5', 'SEPT6', 'SEPT9')
# 
# prova$vogelstein=NA
# prova$NCG6_oncogene=NA
# prova$NCG6_tsg=NA
# prova$vogelstein=prova2$vogelstein_annotation[match(prova$symbol,prova2$symbol)]
# prova$NCG6_oncogene=prova2$NCG6_oncogene[match(prova$symbol,prova2$symbol)]
# prova$NCG6_tsg=prova2$NCG6_tsg[match(prova$symbol,prova2$symbol)]
# 
# dex_CASE$gene_id = sapply(strsplit(as.character(dex_CASE$gene_name),'\\.'),"[", 1 )
# 
# dex_CASE$gtype = 'rst'
# dex_CASE$gtype[which(dex_CASE$gene_id%in%prova$gene_id)] = 'cancer'
# 
# load("/Volumes/LaCie/marco/FOXA1/ANALISI_POST_UV2_DISASTER/Analisi_splicing_automatica/sources/falseCancerGenes.Rdata")
# dex_CASE$gtype[which(dex_CASE$gene_id%in%prova$gene_id[which(FP%in%prova$symbol)])] = 'rst'
# 
# dex_CASE$vogelstein=prova$vogelstein[match(dex_CASE$gene_id,prova$gene_id)]

saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))


# XXX. Add "overall.mean", "overall.sd" to "dex" dataset ----

message('[*] Adding overall.mean and overall.sd to dex dataset and saving ...')

nas = ddply(ex_CASE, .(as_id), summarise
            , good = sum(!is.na(value)) 
            , n    = sum(is.na(value))
            , l    = length(value)
            , md   = median(value, na.rm = T)
            , mn   = mean(value, na.rm = T)
            , sd   = sd(value, na.rm=T)
            , qt3  = quantile(value,seq(0,1,0.25), na.rm = T)["75%"]
            , .progress = 'text')

nas$splice_type = sapply(strsplit(as.character(nas$as_id),'_'),"[", 2 )

dex_CASE$overall.mean=nas$mn
dex_CASE$overall.sd=nas$sd

saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))


# XXX. Run script "simulations_splicing_pemp_succRate.R" for the simulations for pemp and success rate ----

message('[*] RUN BOOTSTRAPPING ON SEPARATE SCRIPT ...')