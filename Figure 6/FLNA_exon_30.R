# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Analyses on FLNA exon 30 ----



# ******************** ----
# ******************** Box plots of FLNA exon 30 inclusion with respect to the double stratification FOXA1-SRSF1 ----
# ******************** ----


library(ggplot2)
library(ggpubr)
library(plyr)
library(reshape2)
library(gridExtra)

rna <- readRDS("Rdata/RNA_TPM_tumor.rds")

ex <- readRDS("Rdata/PRAD_selected_exons.rds")
ex$cases=NULL
ex = subset(ex, splice_type=='ES')

qt = quantile(subset(rna,symbol=='FOXA1')$TPM)
rna$FOXA1_25_75 = NA
rna$FOXA1_25_75[which(rna$patient%in%unique(subset(rna,symbol=='FOXA1' & TPM>=qt['75%'])$patient))]='>=75'
rna$FOXA1_25_75[which(rna$patient%in%unique(subset(rna,symbol=='FOXA1' & TPM<=qt['25%'])$patient))]='<=25'

qt = quantile(subset(rna,symbol=='SFRS1')$TPM)
rna$SRSF1_25_75 = NA
rna$SRSF1_25_75[which(rna$patient%in%unique(subset(rna,symbol=='SFRS1' & TPM>=qt['75%'])$patient))]='>=75'
rna$SRSF1_25_75[which(rna$patient%in%unique(subset(rna,symbol=='SFRS1' & TPM<=qt['25%'])$patient))]='<=25'

ex$FOXA1_25_75='REST'
ex$FOXA1_25_75[which(ex$patient%in%subset(rna,FOXA1_25_75=='>=75')$patient)] = '>=75'
ex$FOXA1_25_75[which(ex$patient%in%subset(rna,FOXA1_25_75=='<=25')$patient)] = '<=25'
ex$SRSF1_25_75='REST'
ex$SRSF1_25_75[which(ex$patient%in%subset(rna,SRSF1_25_75=='>=75')$patient)] = '>=75'
ex$SRSF1_25_75[which(ex$patient%in%subset(rna,SRSF1_25_75=='<=25')$patient)] = '<=25'

ex$FOXA1_SRSF1 = paste0('FOXA1_',ex$FOXA1_25_75,'_SRSF1_',ex$SRSF1_25_75)

my.comparisons=list(c('FOXA1_<=25_SRSF1_<=25','FOXA1_<=25_SRSF1_>=75'),c('FOXA1_<=25_SRSF1_<=25','FOXA1_>=75_SRSF1_>=75'),c('FOXA1_<=25_SRSF1_>=75','FOXA1_>=75_SRSF1_>=75'))
pdf(file="FLNA_exon_PSI_distribution_double_segr_FOXA1_SRSF1_STRATA_25_75_JITTER.pdf", useDingbats = F, height = unit(7,'cm'), width = unit(6,'cm')) # paper = 'a4r', 
ggplot(subset(ex,as_id=='exon_skip_517627' & !FOXA1_25_75=='REST' & !SRSF1_25_75=='REST'),aes(x=FOXA1_SRSF1,y=value))+geom_boxplot(notch = T,outlier.shape = NA)+theme_bw()+
  stat_compare_means(comparisons = my.comparisons)+geom_jitter()
dev.off()







# ******************** ----
# ******************** GLM ----
# ******************** ----


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(viridis)
library(seriation)
library(caret)
library(plyr)

set.seed(30580)

load('Rdata/genes_expr_TPM.Rdata')

TFs  = genes_expr_TPM[, c('PTEN','FOXA1','AR','MYC','MYCN',"HOXB13",'ERG')]
SRPs = genes_expr_TPM[, !colnames(genes_expr_TPM) %in% colnames(TFs)]

rna <- readRDS("Rdata/RNA_TPM_tumor.rds")
rna = subset(rna,symbol=='FOXA1')

ex = readRDS(file="Rdata/ex_pois_ess_exons.rds")
ex$as_id = gsub('exon_skip',"ES", ex$as_id)

tmp <- reshape2::dcast( subset(ex, as_id %in% 'ES_517627')
                        , as_id ~ variable, value.var="value")
rownames(tmp) = tmp$as_id
tmp = tmp[,-1]
ex = as.matrix(t(tmp))
rownames(ex) = substr(rownames(ex), 1,12)


qtf = quantile(FOXA1_log2TPM)
ids = names(FOXA1_log2TPM[which(FOXA1_log2TPM>=qtf['75%'] | FOXA1_log2TPM<=qtf['25%'])])


es = predict(preProcess(ex, method = c("center", "scale","nzv"))
             , newdata = ex)


s = read.csv("Tables/deseq2_FC_PC3_nextseq_novaseq_VCAP.csv")
s = subset(s, cancer=='P' & same_direction & same_sign ) 
ss = ddply(s, .(symbol), summarise, n=length(cell_line) )

# SRPs Validated in both RNAseq 

SRP_selected = as.character(sort(subset(ss, n>1)$symbol))
# "DHX15"    "ESRP2"    "HNRNPK"   "HNRPLL"   "NCBP1"    "SF3B3"    "SNRNP200" "SNRNP40"  "SRSF1"    "WBP11" 

df = list()

df$srp_scaled  = predict( preProcess(SRPs[, SRP_selected], method = c("center", "scale", "YeoJohnson", "nzv"))
                          , newdata = SRPs[, SRP_selected] )

df = lapply(df, function(a,b) {a$FLNA_scaled    = b[rownames(a),'ES_517627']; return(a)}, b=es )


# GLM Generalized Linear Modelling ==========

library(relaimpo)

get_coef = function(f){
  coe = coef(f); coe  = coe[-1]; coe = sort(coe)
  pvalue = coef(summary(f))[,4]
  coe = data.frame(g = names(coe), coe = coe, pval = pvalue[names(coe)])
  coe$g = factor(coe$g, levels = coe$g)
  coe$pval[coe$pval>0.05] = 1
  print(ggplot(coe, aes(x=g, y=coe, fill = -log10(pval)) )+geom_bar(stat='identity')+theme_bw()+coord_flip())
  coe
}

FORMULA = as.formula("FLNA_scaled ~ SRSF1 + ESRP2 + HNRNPK + NCBP1 + SNRNP40 + SNRNP200 + SF3B3 + DHX15 + WBP11 + HNRPLL") # Only DE SRGs in TCGA with at least one active FOXA1 bs and validated in both VCaP and PC3 cells

fit = list()
fit$GLM_foxa1HL_patient_srp_scaled      = glm(data=df$srp_scaled[ids,], FORMULA )

lapply(fit, get_coef)

get_coef = function(f){
  coe = coef(f); coe  = coe[-1]; coe = sort(coe)
  pvalue = coef(summary(f))[,4]
  coe = data.frame(g = names(coe), coe = coe, pval = pvalue[names(coe)])
  coe$g = factor(coe$g, levels = coe$g)
  coe$pval[coe$pval>0.05] = 1
  print(ggplot(coe, aes(x=g, y=coe, fill = (coe+1), size=pval<0.05) )+
          geom_bar(stat='identity',color='black')+theme_bw()+coord_flip()+scale_fill_gradientn(values=scales::rescale(c(1.4,1.2,1,0.9,0.8)), colours=c("red", "red", "white", "blue", "blue"))+
          scale_size_manual(values = c(0,1)))
  coe
}

pdf(file=paste0("FLNA_SRP_scaled_glm_coeff.pdf"), height=unit(3,'cm'), width=unit(3, 'cm'), useDingbats = F )
coe = get_coef(fit$GLM_foxa1HL_patient_srp_scaled)
dev.off()

# Bootstrap Measures of Relative Importance (1000 samples)
set.seed(30580)
boot <- boot.relimp(fit$GLM_foxa1HL_patient_srp_scaled
                    , b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
b = booteval.relimp(boot,sort=TRUE) # print result

pdf(file=paste0("RelImpo_glm_FLNA.pdf"), height=unit(5,'cm'), width=unit(5,'cm'), useDingbats = F)
plot(b)
dev.off()



