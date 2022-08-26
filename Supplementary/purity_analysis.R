# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Analysis of the impact of tumor purity on FOXA1 regulation of SRGs and splicing ----


# Options and libraries ========
options(stringsAsFactors=F)
 
library(readxl)
library(ggplot2)
library(plyr)
library(ggpubr)
library(reshape2)
library(gridExtra)
library(caret)
library(relaimpo)
library(TCGAbiolinks)
library(SummarizedExperiment)
library(clusterProfiler)
library(org.Hs.eg.db)
library(UpSetR)
library(gtools)
library(R.utils)
library(gridExtra)
library(DOSE)
library(igraph)
library(org.Hs.eg.db)
library(tidyverse)
library(magrittr)
library(clusterProfiler)
library(msigdbr)
library(ggplotify)


setwd("")

FIG_DIR=""   # define folder

# Function =========

get_coef = function(f){
  coe = coef(f); coe  = coe[-1]; coe = sort(coe)
  pvalue = coef(summary(f))[,4]
  coe = data.frame(g = names(coe), coe = coe, pval = pvalue[names(coe)])
  coe$g = factor(coe$g, levels = coe$g)
  coe$pval[coe$pval>0.05] = 1
  print(ggplot(coe, aes(x=g, y=coe, fill = coe, size=pval<0.05) )+
          geom_bar(stat='identity',color='black')+theme_bw()+coord_flip()+scale_fill_gradientn(limits=c(0,0.5),colours = c('white','red'))+#scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red")+
          scale_size_manual(values = c(0,1)))
  coe
}

get_DESEQ = function(counts, samples, cond, ctrl){ 
  
  # Perform DESeq analysis
  require(DESeq2)
  require(BiocParallel)
  
  tmp=counts[,samples]
  colData <- data.frame(condition=factor(cond))
  colData$condition <- relevel(colData$condition, ref = ctrl) 
  dds <- DESeqDataSetFromMatrix(tmp, colData, formula(~condition))
  # mcols(dds)$basepairs=tlen
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds, parallel = T, BPPARAM = MulticoreParam(workers = 2))
  # dds <- DESeq(dds, parallel = F)
  res <- results(dds)
  res
  
}



simulazione=function(samples, ri, pl){
  ri$type="WT"
  ri$type[which(ri$sample%in%samples)]="HE"
  x = subset(ri, symbol %in% pl)
  t1 = ddply(x, .(symbol),
             summarise,
             ks=ks.test(TPM[type=="HE"], TPM[type=="WT"], alternative = 'two.sided' )$p.value # %%%%%%%%%%%%%%%%%%%%%% CAMBIATO
  )
  t1$BF.ks = p.adjust(t1$ks, 'bonferroni', n=sum(!is.na(t1$ks)))
  
  return(t1[,c('symbol','ks','BF.ks')])
}


get_wk_test <- function(x){ 
  w = with(x, wilcox.test(value[type], value[!type])$p.value)
  k = with(x, ks.test(    value[type], value[!type])$p.value)
  f = with(x, if(length(table(value))>1) { fligner.test(list(value[type],value[!type]))$p.value } else  {NA})
  return( c('w'=w,'k'=k,'f'=f) ) # c('w'=w,'k'=k,'f'=f)
}


execute_bootstrap <- function(samples, data, type='sample.size'){
  require(plyr)
  if(type=='p.empirical'){
    data$type <- samples$type[match(as.character(data$variable), as.character(samples$variable)) ]
  }else if(type=='sample.size'){
    data <- subset(data, as.character(variable)%in%samples)
  }
  rx  = ddply(data, .(as_id), get_wk_test)
  return(rx)
}



format_enrichr = function(x,  org=org.Hs.eg.db, key="ENTREZID"){
  y <- setReadable(x, OrgDb = org, keyType=key)
  y <- clusterProfiler::mutate(y, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  y <- clusterProfiler::mutate(y, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  y
}

get_enrich=function(x
                    , TERM2GENE,  pAdjustMethod = "BH", pvalueCutoff=1, qvalueCutoff=1, minGSSize = 10
                    ,fromType="ENSEMBL", toType=c("UNIPROT", "SYMBOL", 'ENTREZID','ALIAS'), OrgDb="org.Hs.eg.db",typeFormat='ENTREZID'){
  ids <- bitr(unique(x$ensembl), fromType=fromType, toType=toType, OrgDb=OrgDb)
  ego <- enricher(unique(ids$ENTREZID)  #, universe = unique(bg$ENTREZID)
                  , TERM2GENE = TERM2GENE
                  , pAdjustMethod = pAdjustMethod
                  , pvalueCutoff=pvalueCutoff
                  , qvalueCutoff=qvalueCutoff
                  , minGSSize = minGSSize )
  ego = format_enrichr(ego, org=OrgDb, key=typeFormat)
  list('id'=ids, 'go'=ego)
}



# 1. ASSIGN PATIENTS TO A TUMOR PURITY CATEGORY ----



purity=read_xlsx('Tables/Aran_Nat_Comm_2015/41467_2015_BFncomms9971_MOESM1236_ESM.xlsx',skip = 3)
purity=subset(purity,`Cancer type`=='PRAD')
purity$...8=NULL

rna=readRDS("Rdata/RNA_TPM_tumor.rds")
rna$Sample_ID=substr(rna$Tumor_Sample_Barcode,start = 1,stop = 16)

purity=subset(purity,`Sample ID`%in%unique(rna$Sample_ID))
rna$purity=purity$CPE[match(rna$Sample_ID,purity$`Sample ID`)]
rna$purity=as.numeric(rna$purity)


rna$purity_strata='low'
rna$purity_strata[which(rna$purity>=0.9)]='high'

tmp=subset(rna,symbol=='FOXA1')
tmp=subset(rna,symbol=='FOXA1' & purity_strata=='high')
qt=quantile(tmp$TPM)
tmp$FOXA1_overexpression=F
tmp$FOXA1_overexpression[which(tmp$TPM>=qt['75%'])]=T

tmp2=subset(rna,symbol=='FOXA1' & purity_strata=='low')
qt=quantile(tmp2$TPM)
tmp2$FOXA1_overexpression=F
tmp2$FOXA1_overexpression[which(tmp2$TPM>=qt['75%'])]=T

tmp=rbind(tmp,tmp2)

pdf(paste0(FIG_DIR, file ='TCGA_sample_purity_vs_FOXA1_boxplots.pdf'),width = unit(3,'cm'),height = unit(4,'cm'))
ggplot(tmp,aes(x=FOXA1_overexpression,y=purity))+theme_bw()+geom_boxplot(notch = T)+stat_compare_means()
dev.off()

rna$FOXA1_overexpression=F
rna$FOXA1_overexpression=tmp$FOXA1_overexpression[match(rna$Tumor_Sample_Barcode,tmp$Tumor_Sample_Barcode)]


p1=ggplot(tmp,aes(x=FOXA1_overexpression,y=purity))+theme_bw()+geom_boxplot(notch = T)+stat_compare_means()
p2=ggplot(tmp,aes(x=TPM,y=purity))+theme_bw()+geom_point()+geom_smooth(method = 'lm')+stat_cor(label.x = 300, label.y = 1.05)+xlab('FOXA1 (TPM)')

p3=ggplot(tmp,aes(x=purity_strata,y=TPM))+theme_bw()+geom_boxplot(notch = T)+stat_compare_means()+ylab('FOXA1 (TPM)')
p4=ggplot(tmp,aes(x=purity_strata,fill=FOXA1_overexpression))+theme_bw()+geom_bar()

pdf(paste0(FIG_DIR, file = 'TCGA_sample_purity_characterization__NEW_STRATA.pdf'),width = unit(12,'cm'),height = unit(12,'cm'))
grid.arrange(p1,p2,p3,p4)
dev.off()



pdf(paste0(FIG_DIR, file = 'TCGA_FOXA1_expr_wrt_sample_purity.pdf'),width = unit(4,'cm'),height = unit(5,'cm'))
ggplot(tmp,aes(x=purity_strata,fill=FOXA1_overexpression,y=TPM))+geom_boxplot(notch = T)+theme_bw()
dev.off()

saveRDS(rna, "Rdata/RNA_TPM_tumor.rds")

# 2. GLM ----

pl = readRDS("Rdata/KEGG_Genetic_information_process.GSECA.200105.rds")

rna=readRDS( "Rdata/RNA_TPM_tumor.rds")


lista=c('low', 'high')


for(purity_stratification in lista){
  
  print(purity_stratification)
  
  rna_sub=subset(rna,purity_strata==purity_stratification)
  
  srp.p = ddply(subset(rna_sub, symbol%in%pl$SPLICING_RELATED), .(sample),summarise, tpm = sum(TPM))
  
  p = dcast(subset(rna_sub, symbol %in%c("FOXA1",'AR',"ERG","MYC")), sample~symbol, value.var = 'TPM')
  
  
  p$SRP = srp.p$tpm[match(p$sample,srp.p$sample)]
  rownames(p) = p$sample; p = p[,-1]
  
  x   = predict(preProcess(p, method = c("center", "scale", "YeoJohnson", "nzv"))
                , newdata = p)
  
  
  set.seed(30580)
  
  pdf(file=paste0(FIG_DIR,"SRP_TF_primary_glm_coeff__purity_strata_",purity_stratification,".pdf"), height=unit(3,'cm'), width=unit(3, 'cm'), useDingbats = F )
  fit.p = glm(data=x, SRP~FOXA1+MYC+AR+ERG ); coe1 = get_coef(fit.p)
  dev.off()
  
  set.seed(30580)
  boot <- boot.relimp(fit.p, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
  b = booteval.relimp(boot,sort=TRUE) # print result
  
  pdf(file=paste0(FIG_DIR, "RelImpo_glm_primary__purity_strata_",purity_stratification,".pdf"), height=unit(5,'cm'), width=unit(5,'cm'), useDingbats = F)
  plot(b)
  dev.off()
}


# 3. DIFFERENTIAL EXPRESSION ----

source("sources/config.R")

rna=readRDS("Rdata/RNA_TPM_tumor.rds")

lista=c('low', 'high')


for(purity_stratification in lista){
  
print(purity_stratification)

load("Rdata/PRAD.Rdata")

load("Rdata/KEGG_list.148.SPLICING_RELATED.Rdata") 

# ∞∞∞ EDGER ======================

r2 = subset(rna,purity_strata==purity_stratification)


rc <- assay(data, "raw_count")
rc.norm <- TCGAanalyze_Normalization(tabDF = rc, geneInfo =  TCGAbiolinks::geneInfo)

# quantile filter of genes
rc.norm.filt <- TCGAanalyze_Filtering(tabDF = rc.norm,
                                      method = "quantile",
                                      qnt.cut =  0.1)

rc.norm.filt = rc.norm.filt[,colnames(rc.norm.filt)%in%unique(r2$Tumor_Sample_Barcode)] # %%%%%%%%%%%%


# selection of normal samples "NT"
cntr <- unique( subset(r2, !FOXA1_overexpression)$Tumor_Sample_Barcode)

# selection of tumor samples "TP"
case <- unique( subset(r2, FOXA1_overexpression)$Tumor_Sample_Barcode)

# Diff.expr.analysis (DEA)
dataDEGs <- TCGAanalyze_DEA(mat1 = rc.norm.filt[,cntr],
                            mat2 = rc.norm.filt[,case],
                            Cond1type = "CNTR",
                            Cond2type = "CASE",
                            method = "glmLRT")

kegg   =c("DDX39B","PABPC2","GM10110","QK", "RBFOX1","SAMD4", "SRSF1","SRSF10", "SRSF2","SRSF3","SRSF4","SRSF5","SRSF6","SRSF7",'SRSF8', "SRSF9","U2SURP","ZFP638",'CELF4','CELF5','CELF6','ELAVL1')
symbol =c('BAT1',  "PABPC1","PABPC4", "QKI","A2BP1", 'SAMD4A','SFRS1','SFRS13A','SFRS2','SFRS3','SFRS4','SFRS5','SFRS6','SFRS7','SFRS2B','SFRS9','SR140','ZNF638', "BRUNOL4","BRUNOL5","BRUNOL6","HUR")

names(kegg)=symbol
names(symbol)=kegg

id = which(!rownames(dataDEGs)%in%pl[['SPLICING_RELATED']] & rownames(dataDEGs)%in%symbol);
if(length(id)>0){
  tmp= kegg[rownames(dataDEGs)[id]]
  rownames(dataDEGs)[id]=tmp
}

id = which(!rownames(dataDEGs)%in%pl[['RBP_not_splicing']] & rownames(dataDEGs)%in%symbol);
if(length(id)>0){
  tmp= kegg[rownames(dataDEGs)[id]]
  rownames(dataDEGs)[id]=tmp
}

save(rc.norm, rc.norm.filt, dataDEGs, file=paste0("Rdata/PRAD_EDGR_FOXA1_75th__purity_strata",purity_stratification,".Rdata"))



# ∞∞∞ DESEQ ======================

load('Rdata/PRAD.Rdata')
prad=assays(data)$raw_count

prad <- reshape2::melt(prad)
colnames(prad) <- c("symbol","Tumor_Sample_Barcode","value")
prad$symbol=as.character(prad$symbol)
prad$Tumor_Sample_Barcode=as.character(prad$Tumor_Sample_Barcode)

r2 = subset(rna,purity_strata==purity_stratification)
prad = subset(prad,Tumor_Sample_Barcode%in%r2$Tumor_Sample_Barcode)

xxx=subset(r2,symbol=='FOXA1')

prad$type = "CNTR"
prad$type[which(prad$Tumor_Sample_Barcode%in%subset(xxx,FOXA1_overexpression)$Tumor_Sample_Barcode)] = "CASE"

prad$value=round(prad$value)

prad.rwc <- reshape2::dcast(prad, symbol ~ Tumor_Sample_Barcode, value.var = "value", fun.aggregate = mean)
rownames(prad.rwc) <- as.character(prad.rwc$symbol)
prad.rwc$symbol <- NULL

case = unique(subset(prad,type=='CASE')$Tumor_Sample_Barcode)
cntr = unique(subset(prad,type=='CNTR')$Tumor_Sample_Barcode)

pheno = c("CASE","CNTR")
pheno_nn = c(length(case),length(cntr))
prad.rwc = prad.rwc[,c(case,cntr)]


myc.deout2 <- get_DESEQ(counts = prad.rwc, samples = c(case,cntr), 
                        cond = c(rep(pheno[1],pheno_nn[1]),rep(pheno[2],pheno_nn[2])),
                        ctrl = pheno[2]) 

save(myc.deout2, file=paste0("Rdata/PRAD_DESEQ_FOXA1_75th__purity_strata",purity_stratification,".Rdata"))




# ∞∞∞ SUBSET ON SRGs ========


pl = readRDS('Rdata/KEGG_Genetic_information_process.GSECA.200105.rds')

edger=load(paste0("Rdata/PRAD_EDGR_FOXA1_75th__purity_strata",purity_stratification,".Rdata"))


pl$SPLICING_RELATED[!pl$SPLICING_RELATED%in%rownames(dataDEGs)] #"PABPC1"  "ZNF638"  "QKI"     "HNRNPLL" "SNU13"   "SRSF8" 
rownames(dataDEGs)[which(rownames(dataDEGs)=='SFRS8')]='SRSF8'
rownames(dataDEGs)[which(rownames(dataDEGs)=='QK')]='QKI'
rownames(dataDEGs)[which(rownames(dataDEGs)=='PABPC2')]='PABPC1'
rownames(dataDEGs)[which(rownames(dataDEGs)=='ZFP638')]='ZNF638'
rownames(dataDEGs)[which(rownames(dataDEGs)=='HNRPLL')]='HNRNPLL'
rownames(dataDEGs)[which(rownames(dataDEGs)=='NHP2L1')]='SNU13'
pl$SPLICING_RELATED[!pl$SPLICING_RELATED%in%rownames(dataDEGs)]

dataDEGs = subset(dataDEGs,rownames(dataDEGs)%in%pl$SPLICING_RELATED)


saveRDS(dataDEGs,paste0('Rdata/PRAD_EDGR_FOXA1_75th_148_SRPs__purity_strata',purity_stratification,'.rds'))



deseq=load(paste0('Rdata/PRAD_DESEQ_FOXA1_75th__purity_strata',purity_stratification,'.Rdata'))

dataDEGs = as.data.frame(myc.deout2)
dataDEGs$symbol = sapply(strsplit(rownames(dataDEGs),'[|]'),`[`,1)

pl$SPLICING_RELATED[!pl$SPLICING_RELATED%in%dataDEGs$symbol] # "RBFOX1"  "HNRNPLL" "SNU13"   "SRSF5"   "SRSF9"   "SRSF3"   "SRSF7"   "SRSF4"   "SRSF6"   "SRSF1"   "SRSF2"   "U2SURP"  "SRSF10"  "DDX39B"  "SRSF8" 
dataDEGs$symbol[which(dataDEGs$symbol=='A2BP1')]='RBFOX1'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS5')]='SRSF5'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS9')]='SRSF9'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS3')]='SRSF3'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS7')]='SRSF7'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS4')]='SRSF4'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS6')]='SRSF6'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS1')]='SRSF1'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS2')]='SRSF2'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS5')]='SRSF1'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS13A')]='SRSF10'
dataDEGs$symbol[which(dataDEGs$symbol=='SR140')]='U2SURP'
dataDEGs$symbol[which(dataDEGs$symbol=='BAT1')]='DDX39B'
dataDEGs$symbol[which(dataDEGs$symbol=='SFRS8')]='SRSF8'
dataDEGs$symbol[which(dataDEGs$symbol=='NHP2L1')]='SNU13'
dataDEGs$symbol[which(dataDEGs$symbol=='HNRPLL')]='HNRNPLL'


dataDEGs = subset(dataDEGs,symbol%in%pl$SPLICING_RELATED)

saveRDS(dataDEGs,paste0('Rdata/PRAD_DESEQ_FOXA1_75th_148_SRPs__purity_strata',purity_stratification,'.rds'))



# ∞∞∞ TCGA - simulations----

source('sources/config.R')

rna = rna 
pl = readRDS('Rdata/KEGG_Genetic_information_process.GSECA.200105.rds')
deseq=readRDS(paste0('Rdata/PRAD_DESEQ_FOXA1_75th_148_SRPs__purity_strata',purity_stratification,'.rds'))
dataDEGs=readRDS(paste0('Rdata/PRAD_EDGR_FOXA1_75th_148_SRPs__purity_strata',purity_stratification,'.rds'))

r2 = subset(rna,purity_strata==purity_stratification)



# ∞∞∞ TCGA - empirical pvalue ----

SRPs = pl[['SPLICING_RELATED']]
SRPs%in%unique(r2$symbol)
sum(SRPs%in%unique(r2$symbol))
SRPs[!SRPs%in%unique(r2$symbol)] #  "RBFOX1"  "HNRNPLL" "SNU13"   "SRSF5"   "SRSF9"   "SRSF3"   "SRSF7"   "SRSF4"   "SRSF6"   "SRSF1"   "SRSF2"   "U2SURP"  "SRSF10"  "DDX39B"  "SRSF8"

r2$alias = r2$symbol
r2$symbol[r2$symbol=='A2BP1'] = 'RBFOX1'
r2$symbol[r2$symbol=='SFRS5'] = 'SRSF5'
r2$symbol[r2$symbol=='SFRS9'] = 'SRSF9'
r2$symbol[r2$symbol=='SFRS3'] = 'SRSF3'
r2$symbol[r2$symbol=='SFRS7'] = 'SRSF7'
r2$symbol[r2$symbol=='SFRS4'] = 'SRSF4'
r2$symbol[r2$symbol=='SFRS6'] = 'SRSF6'
r2$symbol[r2$symbol=='SFRS1'] = 'SRSF1'
r2$symbol[r2$symbol=='SFRS2'] = 'SRSF2'
r2$symbol[r2$symbol=='SR140'] = 'U2SURP'
r2$symbol[r2$symbol=='SFRS13A'] = 'SRSF10'
r2$symbol[r2$symbol=='BAT1'] = 'DDX39B'
r2$symbol[r2$symbol=='SFRS2B'] = 'SRSF8'
r2$symbol[r2$symbol=='NHP2L1'] = 'SNU13'
r2$symbol[r2$symbol=='HNRPLL'] = 'HNRNPLL'




library(snow, verbose=F)
library(plyr)

r2 = r2[, c("sample", "symbol", "TPM", "FOXA1_overexpression")]


r2$type = r2$FOXA1_overexpression
r2$FOXA1_overexpression = NULL

pat  = unique(r2[, c("sample", "type")])

n_foxa1_he = sum(pat$type) # = 103

samples_boot = list()
set.seed(30580) 
for(i in 1:10000) samples_boot[[i]]=sample(pat$sample, n_foxa1_he) 

save(samples_boot, r2, SRPs, file=paste0("Rdata/RBP_simulations.montecarlo_TCGA__purity_strata",purity_stratification,".Rdata")) 


library(snow, verbose=F)
library(plyr)

load(paste0('Rdata/RBP_simulations.montecarlo_TCGA__purity_strata',purity_stratification,'.Rdata'))

nclust = 2
print(nclust)
clus <- makeCluster(nclust)

clusterEvalQ(clus, library(plyr) )

print("EXPORTING")
clusterExport(cl=clus, c("r2", "SRPs", "simulazione", "samples_boot"))
print("RUNNING")
sim = parLapply(clus, samples_boot, simulazione,   ri=r2, pl=SRPs)

stopCluster(clus)

print("SAVING RESULTS")

save(sim,file=paste0('Rdata/20220308_RBP_sim_1000.montecarlo.ks.test_TCGA_TWO_SIDED__purity_strata',purity_stratification,'.Rdata'))

print("DONE")


load(paste0('Rdata/RBP_simulations.montecarlo_TCGA__purity_strata',purity_stratification,'.Rdata'))

r2$FOXA1_overexpression = r2$type

x = subset(r2, symbol %in% SRPs)

summarySE <-
  function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
           conf.interval=.95, .drop=TRUE) {
    library(plyr)
    
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )
    
    # Rename the "mean" column
    colnames(datac)[which(colnames(datac)=='mean')]=measurevar
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }

t = summarySE(x, measurevar = "TPM", groupvars = c('symbol','FOXA1_overexpression'), na.rm = T)

t1 = ddply(x, .(symbol),
           summarise,
           median.TPM    = median(TPM[FOXA1_overexpression], na.rm=T),
           median.TPM.wt = median(TPM[!FOXA1_overexpression], na.rm=T),
           ks=ks.test(    TPM[FOXA1_overexpression], TPM[!FOXA1_overexpression], alternative = 'two.sided' )$p.value # %%%%%%%%%%%%%%%%%%%%%% CAMBIATO
)

t1$BF.ks  = p.adjust(t1$ks, 'bonferroni', n=sum(!is.na(t1$ks)))
t1$FDR.ks = p.adjust(t1$ks, 'fdr', n=sum(!is.na(t1$ks)))

load(paste0('Rdata/20220308_RBP_sim_1000.montecarlo.ks.test_TCGA_TWO_SIDED__purity_strata',purity_stratification,'.Rdata'))
ks = do.call(rbind,lapply(sim, "[[", 'ks'))
lks = list(); for (i in 1:ncol(ks)) lks[[i]]=ks[,i]
names(lks)=sim[[1]][,1]
pv = as.list(t1$ks); names(pv)=t1$symbol
pv=pv[names(lks)]
emp = mapply( function(x,y) sum(x<y)/10000, x=lks,y=pv, SIMPLIFY = F)
emp = unlist(emp)
t1$emp=emp[t1$symbol]

t = cbind(t, t1[match(t$symbol, t1$symbol),])
tc = subset(t,FOXA1_overexpression)
tc$TPM.wt = subset(t, !FOXA1_overexpression)$TPM
tc$se.wt  = subset(t, !FOXA1_overexpression)$se
tc$delta = with(tc, TPM-TPM.wt)
tc$L2R = log2(tc$median.TPM)-log2(tc$median.TPM.wt)
tc$FC  = tc$median.TPM/tc$median.TPM.wt

colnames(dataDEGs) = paste0("edger.", colnames(dataDEGs))

tc=cbind(tc, dataDEGs[tc$symbol,])

colnames(deseq) = paste0("deseq.", colnames(deseq))

rownames(deseq)=deseq$deseq.symbol
tc=cbind(tc, deseq[tc$symbol,])


tc$DEG.sign_FDR = with(tc, ( (
  abs(edger.logFC) >= 0.2 &
    edger.FDR<=0.01) |
    (
      abs(deseq.log2FoldChange)>=0.2 &
        deseq.padj<=0.01))
  & (FDR.ks<=0.01) &
    emp<=0.01
)


tc$direction = "DOWN"
tc$direction[which(tc$delta>0)] = "UP"

tc$DEG.sign_FDR[which(is.na(tc$DEG.sign_FDR))]=F

tc$concordant_sign=F
tc$concordant_sign[which(sign(tc$L2R)==sign(tc$edger.logFC) & sign(tc$edger.logFC)==sign(tc$deseq.log2FoldChange))]=T
tc$concordant_sign[which(!tc$DEG.sign_FDR)]='not_DE'

tc$DEG.sign_FDR[which(tc$concordant_sign=='FALSE')]=F

saveRDS(tc,paste0('Rdata/final_dataset_148_SRPs_TCGA__purity_strata',purity_stratification,'.rds'))


two_tailed_FDR=subset(tc,DEG.sign_FDR)

}

# 4. COMPARE DIFFERENTIAL GENE EXPRESSION RESULTS STRATIFICATION PURITY ----

tc_low=readRDS('Rdata/final_dataset_148_SRPs_TCGA__purity_stratalow.rds')
tc_low$dataset='purity_low'
tc_high=readRDS('Rdata/final_dataset_148_SRPs_TCGA__purity_stratahigh.rds')
tc_high$dataset='purity_high'

tc=readRDS('Rdata/final_dataset_148_SRPs_TCGA.rds')
tc$dataset='original'

tc_total=rbind(tc[,c('symbol','FDR.ks','L2R','emp','edger.logFC','edger.FDR','deseq.log2FoldChange','deseq.padj','DEG.sign_FDR','direction','dataset')],
               tc_low[,c('symbol','FDR.ks','L2R','emp','edger.logFC','edger.FDR','deseq.log2FoldChange','deseq.padj','DEG.sign_FDR','direction','dataset')],
               tc_high[,c('symbol','FDR.ks','L2R','emp','edger.logFC','edger.FDR','deseq.log2FoldChange','deseq.padj','DEG.sign_FDR','direction','dataset')])




listInput = list(original=subset(tc_total,DEG.sign_FDR & dataset=='original')$symbol,
                 purity_low=subset(tc_total,DEG.sign_FDR & dataset=='purity_low')$symbol,
                 purity_high=subset(tc_total,DEG.sign_FDR & dataset=='purity_high')$symbol)
pdf(paste0(FIG_DIR,file="upset_plot_DE_SRGs_in_purity_NEW_STRATA.pdf"), height=unit(3.5,'cm'), width=unit(5,'cm'), useDingbats = F)
upset(fromList(listInput), order.by = "freq")
dev.off()



# 5. DIFFERENTIAL SPLICING ----


# Set working path : ----

message('[*] Setting working path and creating directories ...')

args = commandArgs(trailingOnly=TRUE)

gene = 'FOXA1'

rna=readRDS("Rdata/RNA_TPM_tumor.rds")


# both to be run
purity_stratification='low'
purity_stratification='high'


  
print(purity_stratification)

print(paste0('STRATA: ',gene))

work_dir     = paste0(purity_stratification,'_purity/')
STRATA_code  = paste0(gene,'_HE') # 'PTEN_loss'
Data_out_dir    = paste0(work_dir,STRATA_code)
FIG_DIR2 = paste0(Data_out_dir,STRATA_code)

mkdirs(work_dir)
mkdirs(Data_out_dir)
mkdirs(FIG_DIR2)

# Sources and libraries: ----

message('[*] Loading sources and libraries ...')

source('sources/config.R')
source('sources/config2.R')


# Load CASE samples IDs: ----

message('[*] Loading stratification ...')

rna_sub <- subset(rna,purity_strata==purity_stratification)
cases = substr(unique(subset(rna_sub, symbol==gene & FOXA1_overexpression)$sample),1,12)

# Adapt "ex" dataset with the information about which samples are CASES ----

message('[*] Loading ex dataset, integrating with strata info and saving ...')

ex = readRDS('Rdata/PRAD_selected_exons.rds')

ex$patient = substr(ex$variable,1,12)
ex = subset(ex, variable%in%rna_sub$Tumor_Sample_Barcode)
ex$cases = F
ex$cases[which(ex$patient%in%cases)] = T

# Save "ex" dataset with CASES info:
saveRDS(ex,paste0(Data_out_dir,'/PRAD_selected_exons_',STRATA_code,'.rds'))

# Create "dex" dataset with statistics related to CASES ----

ex_CASE = ex
rm(ex)


ex_CASE$type=ex_CASE$cases

message('[*] Building dex dataset ...')

dex_CASE        = ddply(ex_CASE, .(as_id), get_stats_ex_boolean_v2, .progress = 'text')
dex_CASE$BF     = p.adjust(dex_CASE$w, 'bonferroni', n=sum(!is.na(dex_CASE$w)))
dex_CASE$BF.ks  = p.adjust(dex_CASE$k, 'bonferroni', n=sum(!is.na(dex_CASE$k)))
dex_CASE$FDR    = p.adjust(dex_CASE$w, 'fdr', n=sum(!is.na(dex_CASE$w)))
dex_CASE$FDR.ks = p.adjust(dex_CASE$k, 'fdr', n=sum(!is.na(dex_CASE$k)))

yyy = ddply(ex_CASE, .(as_id), summarise, sd_case=sd(value[type]), sd_cntr=sd(value[!type]), .progress = 'text')
dex_CASE_flig   = ddply(subset(ex_CASE,!as_id%in%subset(yyy,sd_case==0 & sd_cntr==0)$as_id), .(as_id), summarise, fligner=fligner.test(list(value[type],value[!type]))$p.value, .progress = 'text')

dex_CASE        = cbind(dex_CASE,dex_CASE_flig$fligner[match(dex_CASE$as_id,dex_CASE_flig$as_id)])
colnames(dex_CASE)[length(colnames(dex_CASE))] = 'fligner'
dex_CASE$fligner[which(is.na(dex_CASE$fligner))] = 1

saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))

for(i in c(3:5,7:9)) dex_CASE[,i]=dex_CASE[,i]*100

dex_CASE$delta.median           = with(dex_CASE, case.median - cntr.median)
dex_CASE$delta.mean             = with(dex_CASE, case.mean - cntr.mean)
dex_CASE$fc.median              = foldchange(dex_CASE$case.median,dex_CASE$cntr.median)
dex_CASE$fc.mean                = foldchange(dex_CASE$case.mean,dex_CASE$cntr.mean)

message('[*] Saving dex dataset ...')

saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))


# Add annotations to "dex" dataset ================

dex_CASE = readRDS(paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))
events_info = readRDS("Rdata/events_info.rds")

dex_CASE = cbind.data.frame( events_info[match(dex_CASE$as_id,events_info$event_id),1:6], dex_CASE[,2:ncol(dex_CASE)] )
rm(events_info)

 
dex_CASE$gene_id = sapply(strsplit(as.character(dex_CASE$gene_name),'\\.'),"[", 1 )


#  Add "overall.mean", "overall.sd" to "dex" dataset ==================

message('[*] Adding overall.mean and overall.sd to dex dataset and saving ...')

ex_CASE= readRDS(paste0(Data_out_dir,'/PRAD_selected_exons_',STRATA_code,'.rds'))

nas = ddply(ex_CASE, .(as_id), summarise
            , good = sum(!is.na(value)) # dovrebbe essere inutile a questo punto (abbiamo gia' fatto il filtro sui PSI)
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


# =  simulations for pemp and success rate ========


# LOADINGS =================

message('[*] Loading variables ...')

library(snowfall)
library(plyr)

dex = readRDS(paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))
ex  = readRDS(paste0(Data_out_dir,'/PRAD_selected_exons_',STRATA_code,'.rds'))
ex$type=ex$cases
ex$as_id=as.character(ex$as_id)
ex$splice_type=as.character(ex$splice_type)
ex$variable=as.character(ex$variable)

ex = subset(ex, as_id%in%subset(dex, w<=0.05 | k<=0.05 )$event_id )

tmp = subset(ex, as_id==ex$as_id[1])[,c('variable',"type")]
patients = unique(tmp)

message( '[*] Variables loaded, start bootstrapping ...')


#********************
set.seed(30580)
NCORES = 2
NSIM   = 1000
#********************

samples <- 1:NSIM
samples <- lapply(samples, function(x,y){ y$type=sample(y$type); return(y) } , y=patients )


ex$type        = NULL
ex$splice_type = NULL

message( '[*] Bootstrap P emp ...')

functionsToCluster = as.character(lsf.str(envir = .GlobalEnv))
sfInit(parallel=TRUE, cpus=NCORES, type="SOCK")
sfExport(list = c(functionsToCluster, "ex"))
sfLibrary(parallel)
sfLibrary(plyr)
sfClusterSetupRNG()
l <- sfLapply(samples, function(x){ tryCatch( execute_bootstrap(  sample=x
                                                                  , data = ex
                                                                  , type='p.empirical')
                                              , error = function(e) return(e)) })

w =  sfLapply(l, "[[", 2)

w = do.call(cbind.data.frame, w)

colnames(w) = NULL
rownames(w) = l[[1]][,1]

k =  sfLapply(l, "[[", 3)

k = do.call(cbind.data.frame, k)

colnames(k) = NULL
rownames(k) = l[[1]][,1]

f =  sfLapply(l, "[[", 4)

f = do.call(cbind.data.frame, f)

colnames(f) = NULL
rownames(f) = l[[1]][,1]

rm(ex)

sfStop()



saveRDS(samples, file = paste0(Data_out_dir,'/EX_bootstrap_pemp_samples_',STRATA_code,'.rds'))
saveRDS(w, file = paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_wilcoxon_',STRATA_code,'.rds'))
saveRDS(k, file = paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_kolgomorov_',STRATA_code,'.rds'))
saveRDS(f, file = paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_fligner_',STRATA_code,'.rds'))

rm(w)
rm(k)
rm(f)

message( '[*] Done')

# Sample size
pheno <- list()
pheno[[1]]  = as.character(subset(patients, type)$variable)
pheno[[2]]  = as.character(subset(patients, !type)$variable)
sample.size = length(pheno[[1]])

# set.seed(30580)
samples = lapply(1:NSIM, function(i) {
  sample_barcode = vector(mode = "numeric", length = 2*sample.size)
  sample_barcode[1:sample.size] = pheno[[1]]
  sample_barcode[(sample.size+1):length(sample_barcode)] = sample(pheno[[2]], sample.size)
  return(sample_barcode)
})

ex  = readRDS(paste0(Data_out_dir,'/PRAD_selected_exons_',STRATA_code,'.rds'))
ex$type=ex$cases

ex = subset(ex, as_id%in%subset(dex, w<=0.05 | k<=0.05 )$event_id )
ex$splice_type = NULL

message( '[*] Bootstrap success rate ...')
functionsToCluster = as.character(lsf.str(envir = .GlobalEnv))
sfInit(parallel=TRUE, cpus=NCORES, type="SOCK")
sfExport(list = c(functionsToCluster, "ex"))
sfLibrary(parallel)
sfLibrary(plyr)
sfClusterSetupRNG()
l <- sfLapply(samples, function(x){ tryCatch( execute_bootstrap(  sample=x
                                                                  , data = ex
                                                                  , type='sample.size')
                                              , error = function(e) return(e)) })


sfStop()

w =  lapply(l, "[[", 2)
w = do.call(cbind.data.frame, w)
colnames(w) = NULL
rownames(w) = l[[1]][,1]

k =  lapply(l, "[[", 3)
k = do.call(cbind.data.frame, k)
rownames(k) = l[[1]][,1]
colnames(k) = NULL

f =  lapply(l, "[[", 4)
f = do.call(cbind.data.frame, f)
rownames(f) = l[[1]][,1]
colnames(f) = NULL

saveRDS(samples, file = paste0(Data_out_dir,'/EX_bootstrap_SR_samples_',STRATA_code,'.rds'))
saveRDS( w, file = paste0(Data_out_dir,'/EX_bootstrap_SR_pv_wilcoxon_',STRATA_code,'.rds'))
saveRDS(k, file = paste0(Data_out_dir,'/EX_bootstrap_SR_pv_kolgomorov_',STRATA_code,'.rds'))
saveRDS(f, file = paste0(Data_out_dir,'/EX_bootstrap_SR_pv_fligner_',STRATA_code,'.rds'))

rm(w)
rm(k)
rm(f)

message( '[*] Done')
message('[*] Collecting bootstrapping results ...')


obs = readRDS(paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))

# P emp - Wilcoxon ----

EX_bootstrap_pemp_pv_wilcoxon <- readRDS(paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_wilcoxon_',STRATA_code,'.rds'))

ugo = melt(as.matrix(EX_bootstrap_pemp_pv_wilcoxon))
ugo = ugo[,c(1,3)]
colnames(ugo)[1]='event_id'
colnames(ugo)[2]='w'

l = ugo
l$obs.w <- obs$w[match(l$event_id,obs$as_id)]

l <- ddply(l, .(event_id), summarise
           
           , p.emp.w = (1+sum(w<obs.w, na.rm = T))/(1+length(na.omit(w)))
           , p.obs.w = sum(obs.w)/1000
)

l.w.pemp = l


# P emp - Kolmogorov ----

EX_bootstrap_pemp_pv_kolgomorov <- readRDS(paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_kolgomorov_',STRATA_code,'.rds'))

ugo = melt(as.matrix(EX_bootstrap_pemp_pv_kolgomorov))
ugo = ugo[,c(1,3)]
colnames(ugo)[1]='event_id'
colnames(ugo)[2]='k'

l = ugo
l$obs.k <- obs$k[match(l$event_id,obs$as_id)]

l <- ddply(l, .(event_id), summarise
           
           , p.emp.k = (1+sum(k<obs.k, na.rm = T))/(1+length(na.omit(k)))
           , p.obs.k = sum(obs.k)/1000
)

l.k.pemp = l


# P emp - Fligner ----

EX_bootstrap_pemp_pv_fligner <- readRDS(paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_fligner_',STRATA_code,'.rds'))

ugo = melt(as.matrix(EX_bootstrap_pemp_pv_fligner))
ugo = ugo[,c(1,3)]
colnames(ugo)[1]='event_id'
colnames(ugo)[2]='f'

l = ugo
l$obs.f <- obs$fligner[match(l$event_id,obs$as_id)]

l <- ddply(l, .(event_id), summarise
           
           , p.emp.f = (1+sum(f<obs.f, na.rm = T))/(1+length(na.omit(f)))
           , p.obs.f = sum(obs.f)/1000
)

l.f.pemp = l


# Sample size - Wilcoxon  ----

EX_bootstrap_SR_pv_wilcoxon <- readRDS(paste0(Data_out_dir,'/EX_bootstrap_SR_pv_wilcoxon_',STRATA_code,'.rds'))

ugo = melt(as.matrix(EX_bootstrap_SR_pv_wilcoxon))

ugo = ugo[,c(1,3)]
colnames(ugo)[1]='event_id'
colnames(ugo)[2]='w'

l = ugo

sig.threshold = 0.05

l <- ddply(l, .(event_id), summarise
           
           , success_rate.w = sum(w<sig.threshold, na.rm = T)/length(na.omit(w))
)

l.w.size = l


# Sample size - Kolmogorov  ----

EX_bootstrap_SR_pv_kolgomorov <- readRDS(paste0(Data_out_dir,'/EX_bootstrap_SR_pv_kolgomorov_',STRATA_code,'.rds'))

ugo = melt(as.matrix(EX_bootstrap_SR_pv_kolgomorov))

ugo = ugo[,c(1,3)]
colnames(ugo)[1]='event_id'
colnames(ugo)[2]='k'

l = ugo

sig.threshold = 0.05

l <- ddply(l, .(event_id), summarise
           
           , success_rate.k = sum(k<sig.threshold, na.rm = T)/length(na.omit(k))
)

l.k.size = l


# Sample size - Fligner  ----

EX_bootstrap_SR_pv_fligner <- readRDS(paste0(Data_out_dir,'/EX_bootstrap_SR_pv_fligner_',STRATA_code,'.rds'))

ugo = melt(as.matrix(EX_bootstrap_SR_pv_fligner))

ugo = ugo[,c(1,3)]
colnames(ugo)[1]='event_id'
colnames(ugo)[2]='f'

l = ugo

sig.threshold = 0.05

l <- ddply(l, .(event_id), summarise
           
           , success_rate.f = sum(f<sig.threshold, na.rm = T)/length(na.omit(f))
)

l.f.size = l


# Integration in the "dex" dataset ==================

message('[*] Integrating bootstrapping results into dex and saving ...')

dex_CASE = readRDS(file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))

dex_CASE$pemp.w = NA
dex_CASE$pemp.k = NA
dex_CASE$pemp.f = NA
dex_CASE$s.rate.w = NA
dex_CASE$s.rate.k = NA
dex_CASE$s.rate.f = NA

dex_CASE$pemp.w = l.w.pemp$p.emp.w[match(dex_CASE$event_id,l.w.pemp$event_id)]
dex_CASE$pemp.k = l.k.pemp$p.emp.k[match(dex_CASE$event_id,l.k.pemp$event_id)]
dex_CASE$pemp.f = l.f.pemp$p.emp.f[match(dex_CASE$event_id,l.f.pemp$event_id)]
dex_CASE$s.rate.w = l.w.size$success_rate.w[match(dex_CASE$event_id,l.w.size$event_id)]
dex_CASE$s.rate.k = l.k.size$success_rate.k[match(dex_CASE$event_id,l.k.size$event_id)]
dex_CASE$s.rate.f = l.f.size$success_rate.f[match(dex_CASE$event_id,l.f.size$event_id)]

saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'_bootstrap.rds'))


# Select significantly differentially included ASE ----

message('[*] Selecting svASE ...')

dex_CASE         = readRDS(paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'_bootstrap.rds'))
dex_CASE = subset(dex_CASE, overall.mean>0.01 & overall.mean<0.99)

dex_CASE$BF = NULL
dex_CASE$BF.ks = NULL
dex_CASE$BF.f = NULL
dex_CASE$FDR = NULL
dex_CASE$FDR.ks = NULL
dex_CASE$FDR.f = NULL

dex_CASE = subset(dex_CASE,!gene_id%in%names(table(dex_CASE$gene_id)[which(table(dex_CASE$gene_id)>500)]))
saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'_bootstrap.rds'))

# Add info about gene symbol and standard deviation ---


map_ENSnames_symbols <- readRDS("sources/map_ENSnames_symbols.rds")
dex_CASE$symbol   = NA
dex_CASE$symbol   = map_ENSnames_symbols$symbol[match(dex_CASE$gene_id,map_ENSnames_symbols$ensembl_gene_id)]
dex_CASE$case.sd  = dex_CASE$case.se*sqrt(dex_CASE$case.N)
dex_CASE$cntr.sd  = dex_CASE$cntr.se*sqrt(dex_CASE$cntr.N)
dex_CASE$fc.sd    = foldchange(dex_CASE$case.sd, dex_CASE$cntr.sd)
dex_CASE$delta.sd = dex_CASE$case.sd-dex_CASE$cntr.sd

qtd = quantile(dex_CASE$delta.mean, seq(0,1,0.01))
qts = quantile(dex_CASE$delta.sd, seq(0,1,0.01))

pdf(file=paste0(FIG_DIR2,'/Quantiles_delta_mean_standard_dev_',STRATA_code,'.pdf'), paper = 'a4', useDingbats = F)
par(mfrow=c(2,1))
plot(qtd, pch=1, col='grey25', main="Quantiles delta mean");abline(h=qtd['15%'], col='red');abline(h=qtd['85%'], col='red')
plot(qts, pch=1, col='grey25', main="Quantiles delta st.dev.");abline(h=qts['20%'], col='red');abline(h=qts['80%'], col='red')
dev.off()

dex_CASE$exon_type="invariant"
dex_CASE$exon_type[which(  (dex_CASE$delta.mean<=qtd["15%"] | dex_CASE$delta.mean>=qtd["85%"])  | 
                             (dex_CASE$delta.sd  <=qts['20%']  | dex_CASE$delta.sd>=qts["80%"])) ]="variant"

dex_CASE = subset(dex_CASE, exon_type=='variant')
message('[*] dex_CASE dim after subset on variant events: ')

dex_CASE$BF     = p.adjust(dex_CASE$w, 'bonferroni', n=sum(!is.na(dex_CASE$w))) 
dex_CASE$BF.ks  = p.adjust(dex_CASE$k, 'bonferroni', n=sum(!is.na(dex_CASE$k)))
dex_CASE$BF.f  = p.adjust(dex_CASE$fligner, 'bonferroni', n=sum(!is.na(dex_CASE$fligner)))
dex_CASE$FDR    = p.adjust(dex_CASE$w, 'fdr', n=sum(!is.na(dex_CASE$w)))
dex_CASE$FDR.ks = p.adjust(dex_CASE$k, 'fdr', n=sum(!is.na(dex_CASE$k)))
dex_CASE$FDR.f = p.adjust(dex_CASE$fligner, 'fdr', n=sum(!is.na(dex_CASE$fligner)))

message('[*] Saving dex dataset with Variant events only ...')

saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'_Variant.rds'))

dexSIGN_CASE = subset(dex_CASE, (FDR<0.05 & pemp.w<0.05 & s.rate.w>0.7) | (FDR.f<0.05 & pemp.f<0.05 & s.rate.f>0.7) & exon_type=='variant')
message('[*] dexSIGN dim with selection on wilcoxon and fligner tests: ')

message('[*] Saving dexSIGN dataset ...')

saveRDS(dexSIGN_CASE, file=paste0(Data_out_dir,'/dexSIGN_',STRATA_code,'.rds'))



# Scatter plots and histograms of the significant events ----

dexSIGN_CASE$delta.class.mean = NA
dexSIGN_CASE$delta.class.mean[which(dexSIGN_CASE$delta.mean>0)] = 'H mean'
dexSIGN_CASE$delta.class.mean[which(dexSIGN_CASE$delta.mean<0)] = 'L mean'
dexSIGN_CASE$delta.class.mean = as.factor(dexSIGN_CASE$delta.class.mean)
dexSIGN_CASE$delta.class.sd = NA
dexSIGN_CASE$delta.class.sd[which(dexSIGN_CASE$delta.sd>0)] = 'H sd'
dexSIGN_CASE$delta.class.sd[which(dexSIGN_CASE$delta.sd<0)] = 'L sd'
dexSIGN_CASE$delta.class.sd = as.factor(dexSIGN_CASE$delta.class.sd)


pdf(file=paste0(FIG_DIR2,'/Scatter_delta.mean_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
print(ggplot(dexSIGN_CASE, aes(x=overall.mean, y=overall.sd, color=delta.class.mean)) +
  geom_point(size=1) + ggtitle(paste0(STRATA_code,' svASEs (',nrow(dexSIGN_CASE),')')) +
  scale_color_manual(values=c(rgb(223/255,0,63/255),rgb(18/255,0,124/255))))# + stat_function(fun = y, color='black') + geom_vline(xintercept = 0.15) + geom_vline(xintercept = 0.85)
dev.off()

pdf(file=paste0(FIG_DIR2,'/Scatter_delta.sd_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
print(ggplot(dexSIGN_CASE, aes(x=overall.mean, y=overall.sd, color=delta.class.sd)) +
  geom_point(size=1) + ggtitle(paste0(STRATA_code,' svASEs (',nrow(dexSIGN_CASE),')')) +
  scale_color_manual(values=c(rgb(242/255,117/255,109/255),rgb(26/255,166/255,184/255))))# + stat_function(fun = y, color='black') + geom_vline(xintercept = 0.15) + geom_vline(xintercept = 0.85)
dev.off()


pdf(file=paste0(FIG_DIR2,'/Histogram_delta.mean_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
print(ggplot(dexSIGN_CASE,aes(x=overall.mean,col=delta.class.mean))+geom_histogram(alpha=0, position="identity", bins=50)+theme_classic() +
  scale_color_manual(values=c(rgb(223/255,0,63/255),rgb(18/255,0,124/255))))
dev.off()

pdf(file=paste0(FIG_DIR2,'/Histogram_delta.sd_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
print(ggplot(dexSIGN_CASE,aes(x=overall.mean,col=delta.class.sd))+geom_histogram(alpha=0, position="identity", bins=50)+theme_classic() +
  scale_color_manual(values=c(rgb(242/255,117/255,109/255),rgb(26/255,166/255,184/255))))
dev.off()



# Vector field plots ----



p1=ggplot(dexSIGN_CASE,aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("ALL")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p2=ggplot(subset(dexSIGN_CASE,event_type=='ES'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("ES")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p3=ggplot(subset(dexSIGN_CASE,event_type=='A3'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("A3")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p4=ggplot(subset(dexSIGN_CASE,event_type=='A5'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("A5")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p5=ggplot(subset(dexSIGN_CASE,event_type=='IR'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("IR")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p6=ggplot(subset(dexSIGN_CASE,event_type=='MEX'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("MEX")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
pdf(file=paste0(FIG_DIR2,'/Scatter_plot_vector_field_color_mean_ALL_EVENTS_AND_TYPES_dexSIGN_',STRATA_code,'.pdf'), width = 17, height = 6.5, useDingbats = F)
print(grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3))
dev.off()

p1=ggplot(dexSIGN_CASE,aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("ALL")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p2=ggplot(subset(dexSIGN_CASE,event_type=='ES'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("ES")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p3=ggplot(subset(dexSIGN_CASE,event_type=='A3'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("A3")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p4=ggplot(subset(dexSIGN_CASE,event_type=='A5'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("A5")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p5=ggplot(subset(dexSIGN_CASE,event_type=='IR'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("IR")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p6=ggplot(subset(dexSIGN_CASE,event_type=='MEX'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("MEX")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
pdf(file=paste0(FIG_DIR2,'/Scatter_plot_vector_field_color_sd_ALL_EVENTS_AND_TYPES_dexSIGN_',STRATA_code,'.pdf'), width = 17, height = 6.5, useDingbats = F)
print(grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3))
dev.off()



# Histograms in 4 classes ----


pvalue=binom.test(x=nrow(subset(dexSIGN_CASE,overall.mean<=0.15 & delta.mean<0)), n=nrow(subset(dexSIGN_CASE,overall.mean<=0.15)), alternative="two.sided")$p.val
if (pvalue<0.001) {
  pvalue = '***'
} else if (pvalue>0.001 & pvalue<0.01) {
  pvalue ='**'
} else if (pvalue>0.01 & pvalue<0.05) {
  pvalue ='*'
} else {pvalue='ns'}
p1=as.grob(~barplot(c(nrow(subset(dexSIGN_CASE,overall.mean<=0.15 & delta.mean<0)),nrow(subset(dexSIGN_CASE,overall.mean<=0.15 & delta.mean>0))),ylim = c(0,2000),xlab = '0-0.15',main=as.character(pvalue)))
pvalue=binom.test(x=nrow(subset(dexSIGN_CASE,overall.mean>0.15 & overall.mean<=0.5 & delta.mean<0)), n=nrow(subset(dexSIGN_CASE,overall.mean>0.15 & overall.mean<=0.5)), alternative="two.sided")$p.value
if (pvalue<0.001) {
  pvalue = '***'
} else if (pvalue>0.001 & pvalue<0.01) {
  pvalue ='**'
} else if (pvalue>0.01 & pvalue<0.05) {
  pvalue ='*'
} else {pvalue='ns'}
p2=as.grob(~barplot(c(nrow(subset(dexSIGN_CASE,overall.mean>0.15 & overall.mean<=0.5 & delta.mean<0)),nrow(subset(dexSIGN_CASE,overall.mean>0.15 & overall.mean<=0.5 & delta.mean>0))),ylim = c(0,2000),xlab = '0.15-0.5',main=as.character(pvalue)))
pvalue=binom.test(x=nrow(subset(dexSIGN_CASE,overall.mean>0.5 & overall.mean<=0.85 & delta.mean<0)), n=nrow(subset(dexSIGN_CASE,overall.mean>0.5 & overall.mean<=0.85)), alternative="two.sided")$p.value
if (pvalue<0.001) {
  pvalue = '***'
} else if (pvalue>0.001 & pvalue<0.01) {
  pvalue ='**'
} else if (pvalue>0.01 & pvalue<0.05) {
  pvalue ='*'
} else {pvalue='ns'}
p3=as.grob(~barplot(c(nrow(subset(dexSIGN_CASE,overall.mean>0.5 & overall.mean<=0.85 & delta.mean<0)),nrow(subset(dexSIGN_CASE,overall.mean>0.5 & overall.mean<=0.85 & delta.mean>0))),ylim = c(0,2000),xlab = '0.5-0.85',main=as.character(pvalue)))
pvalue=binom.test(x=nrow(subset(dexSIGN_CASE,overall.mean>0.85 & delta.mean<0)), n=nrow(subset(dexSIGN_CASE,overall.mean>0.85)), alternative="two.sided")$p.value
if (pvalue<0.001) {
  pvalue = '***'
} else if (pvalue>0.001 & pvalue<0.01) {
  pvalue ='**'
} else if (pvalue>0.01 & pvalue<0.05) {
  pvalue ='*'
} else {pvalue='ns'}
p4=as.grob(~barplot(c(nrow(subset(dexSIGN_CASE,overall.mean>0.85 & delta.mean<0)),nrow(subset(dexSIGN_CASE,overall.mean>0.85 & delta.mean>0))),ylim = c(0,2000),xlab = '0.85-1',main=as.character(pvalue)))
pdf(file=paste0(FIG_DIR2,'/Histogram_delta_mean_4_classes_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
print(grid.arrange(p1,p2,p3,p4,ncol=4))
dev.off()
pvalue=binom.test(x=nrow(subset(dexSIGN_CASE,overall.mean<=0.15 & delta.sd<0)), n=nrow(subset(dexSIGN_CASE,overall.mean<=0.15)), alternative="two.sided")$p.val
if (pvalue<0.001) {
  pvalue = '***'
} else if (pvalue>0.001 & pvalue<0.01) {
  pvalue ='**'
} else if (pvalue>0.01 & pvalue<0.05) {
  pvalue ='*'
} else {pvalue='ns'}
p1=as.grob(~barplot(c(nrow(subset(dexSIGN_CASE,overall.mean<=0.15 & delta.sd<0)),nrow(subset(dexSIGN_CASE,overall.mean<=0.15 & delta.sd>0))),ylim = c(0,2000),xlab = '0-0.15',main=as.character(pvalue)))
pvalue=binom.test(x=nrow(subset(dexSIGN_CASE,overall.mean>0.15 & overall.mean<=0.5 & delta.sd<0)), n=nrow(subset(dexSIGN_CASE,overall.mean>0.15 & overall.mean<=0.5)), alternative="two.sided")$p.value
if (pvalue<0.001) {
  pvalue = '***'
} else if (pvalue>0.001 & pvalue<0.01) {
  pvalue ='**'
} else if (pvalue>0.01 & pvalue<0.05) {
  pvalue ='*'
} else {pvalue='ns'}
p2=as.grob(~barplot(c(nrow(subset(dexSIGN_CASE,overall.mean>0.15 & overall.mean<=0.5 & delta.sd<0)),nrow(subset(dexSIGN_CASE,overall.mean>0.15 & overall.mean<=0.5 & delta.sd>0))),ylim = c(0,2000),xlab = '0.15-0.5',main=as.character(pvalue)))
pvalue=binom.test(x=nrow(subset(dexSIGN_CASE,overall.mean>0.5 & overall.mean<=0.85 & delta.sd<0)), n=nrow(subset(dexSIGN_CASE,overall.mean>0.5 & overall.mean<=0.85)), alternative="two.sided")$p.value
if (pvalue<0.001) {
  pvalue = '***'
} else if (pvalue>0.001 & pvalue<0.01) {
  pvalue ='**'
} else if (pvalue>0.01 & pvalue<0.05) {
  pvalue ='*'
} else {pvalue='ns'}
p3=as.grob(~barplot(c(nrow(subset(dexSIGN_CASE,overall.mean>0.5 & overall.mean<=0.85 & delta.sd<0)),nrow(subset(dexSIGN_CASE,overall.mean>0.5 & overall.mean<=0.85 & delta.sd>0))),ylim = c(0,2000),xlab = '0.5-0.85',main=as.character(pvalue)))
pvalue=binom.test(x=nrow(subset(dexSIGN_CASE,overall.mean>0.85 & delta.sd<0)), n=nrow(subset(dexSIGN_CASE,overall.mean>0.85)), alternative="two.sided")$p.value
if (pvalue<0.001) {
  pvalue = '***'
} else if (pvalue>0.001 & pvalue<0.01) {
  pvalue ='**'
} else if (pvalue>0.01 & pvalue<0.05) {
  pvalue ='*'
} else {pvalue='ns'}
p4=as.grob(~barplot(c(nrow(subset(dexSIGN_CASE,overall.mean>0.85 & delta.sd<0)),nrow(subset(dexSIGN_CASE,overall.mean>0.85 & delta.sd>0))),ylim = c(0,2000),xlab = '0.85-1',main=as.character(pvalue)))
pdf(file=paste0(FIG_DIR2,'/Histogram_delta_sd_4_classes_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
print(grid.arrange(p1,p2,p3,p4,ncol=4))
dev.off()



# Gene ontology ----

dexSIGN = readRDS(paste0(Data_out_dir,'/dexSIGN_',STRATA_code,'.rds'))
dexSIGN$ensembl = dexSIGN$gene_id
dexSIGN$Type = dexSIGN$event_type

cls = c('event_id','ensembl','Type')

ase = dexSIGN[,cls]



# % ORA with 186 KEGG ==================

pll = load('Rdata/KEGG_list.148.SPLICING_RELATED.Rdata')

kegg=read.csv("Tables/KEGG.csv")

k = msigdbr(species = "Homo sapiens", category = "C2", subcategory = 'KEGG')

m_t2g <-  k%>%  dplyr::select(gs_name, entrez_gene)


p = list()
for(j in c("all_events")){
  p[[paste0(j)]] = try(get_enrich(ase, m_t2g))
}

df = rbind.data.frame()
for(i in names(p)){  df = rbind.data.frame(df ,cbind.data.frame(cell=i, p[[i]]$go@result) )}
df$msigDB_158=df$Description
df$msigDB_158[df$msigDB_158=='SRG']='KEGG_SPLICEOSOME'

library(dplyr)

df = left_join(df, kegg, by='msigDB_158')
df$type=df$cell

df = ddply(df, .(type), mutate, rank= rank(p.adjust))
top_10 = unique(subset(df, rank<=10)$Description)
to_plot = subset(df, Description %in% top_10)

tmp = ddply(to_plot,.(ID),summarise,mean_padj=mean(p.adjust),mean_GR=mean(parse_ratio(GeneRatio)))
tmp = tmp[order(tmp$mean_GR,decreasing = T),]
to_plot$ID = factor(to_plot$ID, levels = rev(tmp$ID))


tmp_KEGG_186_TCGA=df
tmp_KEGG_186_TCGA$ID=as.character(tmp_KEGG_186_TCGA$ID)


pdf(file=paste0(FIG_DIR2,'/ORA_all_event_types_TCGA_KEGG_186_', STRATA_code, '.pdf'), height = unit(4,'cm'), width = unit(8, 'cm'), useDingbats = F)
to_plot=subset(to_plot,p.adjust<0.1)
print(ggplot( to_plot,
        aes(x= parse_ratio(GeneRatio), y= ID)) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=factor(p.adjust<=.25), fill=p.adjust, size = Count),shape=21) +
  geom_text(aes(label=rank)) +
  # scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=c('transparent','black'))+
  scale_fill_viridis_c(option = 'D',guide=guide_colorbar(reverse=T, draw.llim = T), direction = -1)+#,
  # , values  =scales::rescale(c(0,0.001,0.1,0.25,0.5,.75,1)), breaks=c(0,0.01,0.1,0.25,.5,1)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + xlab("Gene Ratio") +ylab(NULL)+
  facet_wrap(~type, nrow=1))
dev.off()

saveRDS(to_plot,paste0('Data_out_dir/ORA_all_event_types_TCGA_KEGG_186_', STRATA_code, '.rds'))



# Survival analysis on FLNA exon 30 ----

# FLNA exon 30 higher inclusion upon FOXA1 HE is significant only in the high purity group. We now want to check whether its inclusion is prognostic also in this set of samples.

library(survminer)
library(survival)

dexSIGN = readRDS(paste0(Data_out_dir,'/dexSIGN_',STRATA_code,'.rds'))

ex = readRDS(paste0(Data_out_dir,'/PRAD_selected_exons_',STRATA_code,'.rds'))
ex$as_id =  as.character(ex$as_id)
ex$splice_type =  as.character(ex$splice_type)
df = readRDS("Rdata/survival_data_20190304.rds")
df$first_event_months = df$first_event/(365/12) 

dexSIGN_pois_ess = subset(dexSIGN,event_id=='exon_skip_517627')

surv_summary_pois_ess = matrix(nrow = 5, ncol = nrow(dexSIGN_pois_ess))
dff = subset(df, first_event_months<110)
for (i in 1:nrow(dexSIGN_pois_ess)) {
  event = dexSIGN_pois_ess$event_id[i]
  incl = subset(ex,as_id==event)$value
  names(incl) = substr(subset(ex,as_id==event)$variable,1,12)
  
  qt=quantile(incl)
  incl2 = rep(NA, length(incl))
  names(incl2)=names(incl)
  incl2[names(incl[incl>=qt['75%']])] = 'HI'
  incl2[names(incl[incl<=qt['25%']])] = 'LI'
  
  tmp = rep(NA, nrow(dff))
  names(tmp) = dff$bcr_patient_barcode
  tmp = incl2[names(tmp)]
  
  dff = cbind(dff,tmp)
  dff$tmp = factor(dff$tmp,levels=c('LI','HI'))
  formula <- as.formula("Surv(time=first_event_months, event=DFS) ~ tmp")
  cox <- coxph(formula, dff)
  tmp = summary(cox)
  
  surv_summary_pois_ess[1,i] = tmp$sctest[3]
  surv_summary_pois_ess[2,i] = tmp$coefficients[2]
  tmp2 = cox.zph(cox)
  surv_summary_pois_ess[3,i] = tmp2$table[3]
  
  surv_summary_pois_ess[4,i] = tmp$conf.int[3]
  surv_summary_pois_ess[5,i] = tmp$conf.int[4]
  
  # dff$tmp = NULL
  
  print(i)
}
colnames(surv_summary_pois_ess) = dexSIGN_pois_ess$event_id
rownames(surv_summary_pois_ess) = c('logrank_pvalue','hazard_ratio','prop_ass_pvalue','lower_CI','upper_CI')


gfit <- survfit(formula = as.formula("Surv(time=first_event_months, event=DFS) ~ tmp"), data = subset(dff, first_event_months<110))
pdf(file=paste0(FIG_DIR2,'/kaplan_meier_FLNA_exon_30.pdf'), height = unit(5.5,'cm'), width = unit(5, 'cm'), useDingbats = F)
print(ggsurvplot(gfit, data = subset(dff, first_event_months<110), risk.table = TRUE, palette = c('green','red','darkgreen','orange'), pval = T, pval.method = T))
dev.off()




#  Cumulative distr. of events  ----
 


dexSIGN=readRDS(paste0(Data_out_dir,'/dexSIGN_',STRATA_code,'.rds'))


# MEAN ----
dexSIGN_sil=subset(dexSIGN,delta.mean<0)
dexSIGN_enh=subset(dexSIGN,delta.mean>0)

res = hist(dexSIGN_sil$overall.mean,nclass = 50)
res_df_sil=as.data.frame(cbind(res$mids,res$counts))
colnames(res_df_sil)=c('overall.mean.psi','counts')
res_df_sil$cum_counts=cumsum(res_df_sil$counts)
ggplot(res_df_sil,aes(x=overall.mean.psi,y=cum_counts))+geom_line()+theme_bw()

res = hist(dexSIGN_enh$overall.mean,nclass = 50)
res_df_enh=as.data.frame(cbind(res$mids,res$counts))
colnames(res_df_enh)=c('overall.mean.psi','counts')
res_df_enh$cum_counts=cumsum(res_df_enh$counts)
ggplot(res_df_enh,aes(x=overall.mean.psi,y=cum_counts))+geom_line()+theme_bw()

res_df_sil$type='sil'
res_df_sil$cum_counts_fraction=res_df_sil$cum_counts/res_df_sil$cum_counts[nrow(res_df_sil)]
res_df_enh$type='enh'
res_df_enh$cum_counts_fraction=res_df_enh$cum_counts/res_df_enh$cum_counts[nrow(res_df_enh)]

res_df=rbind(res_df_sil,res_df_enh)
p_mean=ggplot(res_df,aes(x=overall.mean.psi,y=cum_counts_fraction,col=type))+geom_line()+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed")
p_mean_counts=ggplot(res_df,aes(x=overall.mean.psi,y=cum_counts,col=type))+geom_line()+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed")
res_df_mean_true = res_df

res_df_2=cbind(res_df_sil,res_df_enh)
colnames(res_df_2)[1:5]=paste0(colnames(res_df_2)[1:5],'_SIL')
colnames(res_df_2)[6:10]=paste0(colnames(res_df_2)[6:10],'_ENH')
res_df_2$enh_sil_ratio=res_df_2$counts_ENH/(res_df_2$counts_ENH+res_df_2$counts_SIL)
p_ratio_mean=ggplot(res_df_2,aes(x=overall.mean.psi_ENH,y=enh_sil_ratio))+geom_line()+theme_bw()+geom_hline(yintercept = 0.5,linetype = "dashed")+ggtitle('Mean')+
  stat_smooth(aes(y=enh_sil_ratio, x=overall.mean.psi_ENH), method = lm, formula = y ~ poly(x, 10), se = FALSE, size=0.75)


# SD ----
dexSIGN_sil=subset(dexSIGN,delta.sd<0)
dexSIGN_enh=subset(dexSIGN,delta.sd>0)

res = hist(dexSIGN_sil$overall.mean,nclass = 50)
res_df_sil=as.data.frame(cbind(res$mids,res$counts))
colnames(res_df_sil)=c('overall.mean.psi','counts')
res_df_sil$cum_counts=cumsum(res_df_sil$counts)
ggplot(res_df_sil,aes(x=overall.mean.psi,y=cum_counts))+geom_line()+theme_bw()

res = hist(dexSIGN_enh$overall.mean,nclass = 50)
res_df_enh=as.data.frame(cbind(res$mids,res$counts))
colnames(res_df_enh)=c('overall.mean.psi','counts')
res_df_enh$cum_counts=cumsum(res_df_enh$counts)
ggplot(res_df_enh,aes(x=overall.mean.psi,y=cum_counts))+geom_line()+theme_bw()

res_df_sil$type='sil'
res_df_sil$cum_counts_fraction=res_df_sil$cum_counts/res_df_sil$cum_counts[nrow(res_df_sil)]
res_df_enh$type='enh'
res_df_enh$cum_counts_fraction=res_df_enh$cum_counts/res_df_enh$cum_counts[nrow(res_df_enh)]

res_df=rbind(res_df_sil,res_df_enh)
p_sd=ggplot(res_df,aes(x=overall.mean.psi,y=cum_counts_fraction,col=type))+geom_line()+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed")
p_sd_counts=ggplot(res_df,aes(x=overall.mean.psi,y=cum_counts,col=type))+geom_line()+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed")
res_df_sd_true = res_df

res_df_2=cbind(res_df_sil,res_df_enh)
colnames(res_df_2)[1:5]=paste0(colnames(res_df_2)[1:5],'_SIL')
colnames(res_df_2)[6:10]=paste0(colnames(res_df_2)[6:10],'_ENH')
res_df_2$enh_sil_ratio=res_df_2$counts_ENH/(res_df_2$counts_ENH+res_df_2$counts_SIL)
p_ratio_sd=ggplot(res_df_2,aes(x=overall.mean.psi_ENH,y=enh_sil_ratio))+geom_line()+theme_bw()+geom_hline(yintercept = 0.5,linetype = "dashed")+ggtitle('Standard deviation')+
  stat_smooth(aes(y=enh_sil_ratio, x=overall.mean.psi_ENH), method = lm, formula = y ~ poly(x, 10), se = FALSE, size=0.75)


# NULL MODEL ----

set.seed(8890)

for (i in 1:1000) {
  
  print(i)
  
  # MEAN ---
  
  dexSIGN$delta.mean_random = sample(x = c(-1,1),size = nrow(dexSIGN),replace = T)

  dexSIGN_sil=subset(dexSIGN,delta.mean_random<0)
  dexSIGN_enh=subset(dexSIGN,delta.mean_random>0)
  
  res = hist(dexSIGN_sil$overall.mean,nclass = 50)
  res_df_sil=as.data.frame(cbind(res$mids,res$counts))
  colnames(res_df_sil)=c('overall.mean.psi','counts')
  res_df_sil$cum_counts=cumsum(res_df_sil$counts)

  res = hist(dexSIGN_enh$overall.mean,nclass = 50)
  res_df_enh=as.data.frame(cbind(res$mids,res$counts))
  colnames(res_df_enh)=c('overall.mean.psi','counts')
  res_df_enh$cum_counts=cumsum(res_df_enh$counts)

  res_df_sil$type='sil'
  res_df_sil$cum_counts_fraction=res_df_sil$cum_counts/res_df_sil$cum_counts[nrow(res_df_sil)]
  res_df_enh$type='enh'
  res_df_enh$cum_counts_fraction=res_df_enh$cum_counts/res_df_enh$cum_counts[nrow(res_df_enh)]
  
  res_df=rbind(res_df_sil,res_df_enh)
  
  res_df_2=cbind(res_df_sil,res_df_enh)
  colnames(res_df_2)[1:5]=paste0(colnames(res_df_2)[1:5],'_SIL')
  colnames(res_df_2)[6:10]=paste0(colnames(res_df_2)[6:10],'_ENH')
  res_df_2$enh_sil_ratio=res_df_2$counts_ENH/(res_df_2$counts_ENH+res_df_2$counts_SIL)
  
  if (i==1) {
    res_df$rep=i
    res_df_mean = res_df
  } else {
    res_df$rep=i
    res_df_mean = rbind(res_df_mean,res_df)
  }
  
  # SD ---
  
  dexSIGN$delta.sd_random = sample(x = c(-1,1),size = nrow(dexSIGN),replace = T)

  dexSIGN_sil=subset(dexSIGN,delta.sd_random<0)
  dexSIGN_enh=subset(dexSIGN,delta.sd_random>0)
  
  res = hist(dexSIGN_sil$overall.mean,nclass = 50)
  res_df_sil=as.data.frame(cbind(res$mids,res$counts))
  colnames(res_df_sil)=c('overall.mean.psi','counts')
  res_df_sil$cum_counts=cumsum(res_df_sil$counts)

  res = hist(dexSIGN_enh$overall.mean,nclass = 50)
  res_df_enh=as.data.frame(cbind(res$mids,res$counts))
  colnames(res_df_enh)=c('overall.mean.psi','counts')
  res_df_enh$cum_counts=cumsum(res_df_enh$counts)

  res_df_sil$type='sil'
  res_df_sil$cum_counts_fraction=res_df_sil$cum_counts/res_df_sil$cum_counts[nrow(res_df_sil)]
  res_df_enh$type='enh'
  res_df_enh$cum_counts_fraction=res_df_enh$cum_counts/res_df_enh$cum_counts[nrow(res_df_enh)]
  
  res_df=rbind(res_df_sil,res_df_enh)
  
  res_df_2=cbind(res_df_sil,res_df_enh)
  colnames(res_df_2)[1:5]=paste0(colnames(res_df_2)[1:5],'_SIL')
  colnames(res_df_2)[6:10]=paste0(colnames(res_df_2)[6:10],'_ENH')
  res_df_2$enh_sil_ratio=res_df_2$counts_ENH/(res_df_2$counts_ENH+res_df_2$counts_SIL)
  
  if (i==1) {
    res_df$rep=i
    res_df_sd = res_df
  } else {
    res_df$rep=i
    res_df_sd = rbind(res_df_sd,res_df)
  }
  
}

res_df_sd$rep=as.character(res_df_sd$rep)
res_df_mean$rep=as.character(res_df_mean$rep)

if (purity_stratification=='high') {
  save(res_df_mean,res_df_sd,file='Rdata/null_model_splicing_prob_05_for_delta_pos_and_delta_neg__high_purity.Rdata')
} else if (purity_stratification=='low') {
  save(res_df_mean,res_df_sd,file='Rdata/null_model_splicing_prob_05_for_delta_pos_and_delta_neg__low_purity.Rdata')
}




p_sd_counts=
  ggplot(res_df_sd,aes(x=overall.mean.psi,y=cum_counts,col=type))+geom_line()+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed")

to_plot_enh=ddply(subset(res_df_mean,type=='enh'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts),median_c=median(cum_counts),c_95=quantile(cum_counts,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
to_plot_sil=ddply(subset(res_df_mean,type=='sil'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts),median_c=median(cum_counts),c_95=quantile(cum_counts,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
p_mean_compl = ggplot() +
  geom_line(data=res_df_mean_true,mapping=aes(x=overall.mean.psi,y=cum_counts,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


tmp=subset(res_df_mean_true,type=='sil')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts[which(tmp$overall.mean.psi>0.51)] = tmp$cum_counts[which(tmp$overall.mean.psi>0.51)]-tmp$cum_counts[which(tmp$overall.mean.psi==0.51)]
tmp_s=tmp
tmp=subset(res_df_mean_true,type=='enh')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts[which(tmp$overall.mean.psi>0.51)] = tmp$cum_counts[which(tmp$overall.mean.psi>0.51)]-tmp$cum_counts[which(tmp$overall.mean.psi==0.51)]
tmp_e=tmp
res_df_mean_true_TMP=rbind(tmp_s,tmp_e)

tmp=to_plot_enh
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi>0.51)] = tmp$mean_c[which(tmp$overall.mean.psi>0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_enh_TMP=tmp

tmp=to_plot_sil
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi>0.51)] = tmp$mean_c[which(tmp$overall.mean.psi>0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_sil_TMP=tmp

p_mean_compl_2 = ggplot() +
  geom_line(data=res_df_mean_true_TMP,mapping=aes(x=overall.mean.psi,y=cum_counts,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')




to_plot_enh=ddply(subset(res_df_sd,type=='enh'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts),median_c=median(cum_counts),c_95=quantile(cum_counts,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
to_plot_sil=ddply(subset(res_df_sd,type=='sil'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts),median_c=median(cum_counts),c_95=quantile(cum_counts,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
p_sd_compl = ggplot() +
  geom_line(data=res_df_sd_true,mapping=aes(x=overall.mean.psi,y=cum_counts,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta sd')



if (purity_stratification=='high') {
  pdf(paste0(FIG_DIR2, file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg__high_purity.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
  print(grid.arrange(p_mean_compl,p_sd_compl,ncol=2))
  dev.off()
} else if (purity_stratification=='low') {
  pdf(paste0(FIG_DIR2, file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg__low_purity.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
  print(grid.arrange(p_mean_compl,p_sd_compl,ncol=2))
  dev.off()
}




tmp=subset(res_df_sd_true,type=='sil')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts[which(tmp$overall.mean.psi>0.51)] = tmp$cum_counts[which(tmp$overall.mean.psi>0.51)]-tmp$cum_counts[which(tmp$overall.mean.psi==0.51)]
tmp_s=tmp
tmp=subset(res_df_sd_true,type=='enh')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts[which(tmp$overall.mean.psi>0.51)] = tmp$cum_counts[which(tmp$overall.mean.psi>0.51)]-tmp$cum_counts[which(tmp$overall.mean.psi==0.51)]
tmp_e=tmp
res_df_mean_true_TMP=rbind(tmp_s,tmp_e)

tmp=to_plot_enh
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi>0.51)] = tmp$mean_c[which(tmp$overall.mean.psi>0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_enh_TMP=tmp

tmp=to_plot_sil
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi>0.51)] = tmp$mean_c[which(tmp$overall.mean.psi>0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_sil_TMP=tmp

p_sd_compl_2 = ggplot() +
  geom_line(data=res_df_mean_true_TMP,mapping=aes(x=overall.mean.psi,y=cum_counts,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


if (purity_stratification=='high') {
  pdf(paste0(FIG_DIR2, file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_SPLITTED__high_purity.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
  print(grid.arrange(p_mean_compl_2,p_sd_compl_2,ncol=2))
  dev.off()
} else if (purity_stratification=='low') {
  pdf(paste0(FIG_DIR2, file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_SPLITTED__low_purity.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
  print(grid.arrange(p_mean_compl_2,p_sd_compl_2,ncol=2))
  dev.off()
}




# MEAN ----


tmp_e = subset(res_df_mean_true,type=='enh')
tmp_e=tmp_e[order(tmp_e$overall.mean.psi,decreasing = T),]
tmp_e$cum_counts_inverted=cumsum(tmp_e$counts)
tmp_s = subset(res_df_mean_true,type=='sil')
tmp_s=tmp_s[order(tmp_s$overall.mean.psi,decreasing = T),]
tmp_s$cum_counts_inverted=cumsum(tmp_s$counts)

res_df_mean_true_inverted=rbind(tmp_e,tmp_s)


for (i in unique(res_df_mean$rep)) {
  
  tmp = subset(res_df_mean,rep==i)
  tmp_e = subset(tmp,type=='enh')
  tmp_e=tmp_e[order(tmp_e$overall.mean.psi,decreasing = T),]
  tmp_e$cum_counts_inverted=cumsum(tmp_e$counts)
  tmp_s = subset(tmp,type=='sil')
  tmp_s=tmp_s[order(tmp_s$overall.mean.psi,decreasing = T),]
  tmp_s$cum_counts_inverted=cumsum(tmp_s$counts)
  
  tmp=rbind(tmp_e,tmp_s)
  
  if (i == unique(res_df_mean$rep)[1]) {
    res_df_mean_inverted = tmp
  } else {
    res_df_mean_inverted = rbind(res_df_mean_inverted,tmp)
  }
  
}


to_plot_enh=ddply(subset(res_df_mean_inverted,type=='enh'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts_inverted),median_c=median(cum_counts_inverted),c_95=quantile(cum_counts_inverted,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts_inverted,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
to_plot_sil=ddply(subset(res_df_mean_inverted,type=='sil'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts_inverted),median_c=median(cum_counts_inverted),c_95=quantile(cum_counts_inverted,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts_inverted,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
p_mean_compl_inv = ggplot() +
  geom_line(data=res_df_mean_true_inverted,mapping=aes(x=overall.mean.psi,y=cum_counts_inverted,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


tmp=subset(res_df_mean_true_inverted,type=='sil')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)] = tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)]-tmp$cum_counts_inverted[which(tmp$overall.mean.psi==0.51)]
tmp_s=tmp
tmp=subset(res_df_mean_true_inverted,type=='enh')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)] = tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)]-tmp$cum_counts_inverted[which(tmp$overall.mean.psi==0.51)]
tmp_e=tmp
res_df_mean_true_TMP=rbind(tmp_s,tmp_e)

tmp=to_plot_enh
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi<0.51)] = tmp$mean_c[which(tmp$overall.mean.psi<0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_enh_TMP=tmp

tmp=to_plot_sil
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi<0.51)] = tmp$mean_c[which(tmp$overall.mean.psi<0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_sil_TMP=tmp

p_mean_compl_2_inv = ggplot() +
  geom_line(data=res_df_mean_true_TMP,mapping=aes(x=overall.mean.psi,y=cum_counts_inverted,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')





# SD ----


tmp_e = subset(res_df_sd_true,type=='enh')
tmp_e=tmp_e[order(tmp_e$overall.mean.psi,decreasing = T),]
tmp_e$cum_counts_inverted=cumsum(tmp_e$counts)
tmp_s = subset(res_df_sd_true,type=='sil')
tmp_s=tmp_s[order(tmp_s$overall.mean.psi,decreasing = T),]
tmp_s$cum_counts_inverted=cumsum(tmp_s$counts)

res_df_sd_true_inverted=rbind(tmp_e,tmp_s)


for (i in unique(res_df_sd$rep)) {
  
  tmp = subset(res_df_sd,rep==i)
  tmp_e = subset(tmp,type=='enh')
  tmp_e=tmp_e[order(tmp_e$overall.mean.psi,decreasing = T),]
  tmp_e$cum_counts_inverted=cumsum(tmp_e$counts)
  tmp_s = subset(tmp,type=='sil')
  tmp_s=tmp_s[order(tmp_s$overall.mean.psi,decreasing = T),]
  tmp_s$cum_counts_inverted=cumsum(tmp_s$counts)
  
  tmp=rbind(tmp_e,tmp_s)
  
  if (i == unique(res_df_sd$rep)[1]) {
    res_df_sd_inverted = tmp
  } else {
    res_df_sd_inverted = rbind(res_df_sd_inverted,tmp)
  }
  
}


to_plot_enh=ddply(subset(res_df_sd_inverted,type=='enh'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts_inverted),median_c=median(cum_counts_inverted),c_95=quantile(cum_counts_inverted,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts_inverted,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
to_plot_sil=ddply(subset(res_df_sd_inverted,type=='sil'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts_inverted),median_c=median(cum_counts_inverted),c_95=quantile(cum_counts_inverted,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts_inverted,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
p_sd_compl_inv = ggplot() +
  geom_line(data=res_df_sd_true_inverted,mapping=aes(x=overall.mean.psi,y=cum_counts_inverted,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


tmp=subset(res_df_sd_true_inverted,type=='sil')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)] = tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)]-tmp$cum_counts_inverted[which(tmp$overall.mean.psi==0.51)]
tmp_s=tmp
tmp=subset(res_df_sd_true_inverted,type=='enh')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)] = tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)]-tmp$cum_counts_inverted[which(tmp$overall.mean.psi==0.51)]
tmp_e=tmp
res_df_sd_true_TMP=rbind(tmp_s,tmp_e)

tmp=to_plot_enh
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi<0.51)] = tmp$mean_c[which(tmp$overall.mean.psi<0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_enh_TMP=tmp

tmp=to_plot_sil
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi<0.51)] = tmp$mean_c[which(tmp$overall.mean.psi<0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_sil_TMP=tmp

p_sd_compl_2_inv = ggplot() +
  geom_line(data=res_df_sd_true_TMP,mapping=aes(x=overall.mean.psi,y=cum_counts_inverted,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')



if (purity_stratification=='high') {
  pdf(paste0(FIG_DIR2, file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_inverted__high_purity.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
  print(grid.arrange(p_mean_compl_inv,p_sd_compl_inv,ncol=2))
  dev.off()
  
  pdf(paste0(FIG_DIR2, file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_SPLITTED_inverted__high_purity.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
  print(grid.arrange(p_mean_compl_2_inv,p_sd_compl_2_inv,ncol=2))
  dev.off()
} else if (purity_stratification=='low') {
  pdf(paste0(FIG_DIR2, file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_inverted__low_purity.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
  print(grid.arrange(p_mean_compl_inv,p_sd_compl_inv,ncol=2))
  dev.off()
  
  pdf(paste0(FIG_DIR2, file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_SPLITTED_inverted__low_purity.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
  print(grid.arrange(p_mean_compl_2_inv,p_sd_compl_2_inv,ncol=2))
  dev.off()
}





 
# 6. ANALYSIS USING NORMALIZED RNASEQDB DATA ----
 

# Load and normalize datasets to counts per million reads ----

prad=read.delim2('Tables/RNAseqDB/prad-rsem-count-tcga-t.txt',stringsAsFactors = F)
rownames_prad=prad$Hugo_Symbol
prad$Hugo_Symbol=NULL
prad$Entrez_Gene_Id=NULL
prad=apply(prad, 2, as.numeric)
rownames(prad)=rownames_prad
prad=apply(prad, 2, function(x) x/sum(x)*1000000)
colnames(prad)=gsub(colnames(prad),pattern = '[.]',replacement = '-')

prostate=read.delim2('Tables/RNAseqDB/prostate-rsem-count-gtex.txt',stringsAsFactors = F)
rownames_prostate=prostate$Hugo_Symbol
prostate$Hugo_Symbol=NULL
prostate$Entrez_Gene_Id=NULL
prostate=apply(prostate, 2, as.numeric)
rownames(prostate)=rownames_prostate
prostate=apply(prostate, 2, function(x) x/sum(x)*1000000)


# Remove from normal samples outliers with high xCell ImmuneScore ----

source("sources/utils_RNA.R")
c=readRDS('Rdata/gene_length.rds')

library(xCell)

prostate_length=c$Length[match(rownames(prostate),c$gene_name)]
prostate_tpms = countToTpm2(prostate, prostate_length)
prostate_xcell=xCellAnalysis(prostate_tpms)
prostate_xcell=as.data.frame(t(prostate_xcell))
ggplot(prostate_xcell,aes(x=MicroenvironmentScore,y=ImmuneScore))+geom_point()+geom_smooth(method = 'lm')+stat_cor()+theme_bw()
pg=ggplot(prostate_xcell,aes(x=1,y=ImmuneScore))+geom_boxplot(notch = T)+geom_smooth(method = 'lm')+stat_cor()+theme_bw()+geom_jitter(size=0.5)+ggtitle('GTEX - Immune')
pgg=ggplot(prostate_xcell,aes(x=1,y=StromaScore))+geom_boxplot(notch = T)+geom_smooth(method = 'lm')+stat_cor()+theme_bw()+geom_jitter(size=0.5)+ggtitle('GTEX - Stroma')
prostate = prostate[,colnames(prostate)%in%rownames(subset(prostate_xcell,ImmuneScore<0.035))]
saveRDS(prostate_xcell,'Rdata/xCell_GTEX.rds')

prad_length=c$Length[match(rownames(prad),c$gene_name)]
prad_tpms = countToTpm2(prad, prad_length)
prad_xcell=xCellAnalysis(prad_tpms)
prad_xcell=as.data.frame(t(prad_xcell))
saveRDS(prad_xcell,'Rdata/xCell_primary_PC_RNAseqDB.rds')


#  Select highly pure TCGA tumor samples ----

prad_pure = prad[,unique(subset(rna,purity>=0.9)$Tumor_Sample_Barcode)[which(unique(subset(rna,purity>=0.9)$Tumor_Sample_Barcode)%in%colnames(prad))]]


#  GLMs and Differential expression with bootstrapping ----


c=readRDS("Rdata/gene_length.rds")
c$gene_name[which(c$gene_name=='SNU13')]='NHP2L1'

summarySE <-
  function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
           conf.interval=.95, .drop=TRUE) {
    library(plyr)
    
    # New version of length which can handle NA's: if na.rm==T, don't count them
    length2 <- function (x, na.rm=FALSE) {
      if (na.rm) sum(!is.na(x))
      else       length(x)
    }
    
    # This does the summary. For each group's data frame, return a vector with
    # N, mean, and sd
    datac <- ddply(data, groupvars, .drop=.drop,
                   .fun = function(xx, col) {
                     c(N    = length2(xx[[col]], na.rm=na.rm),
                       mean = mean   (xx[[col]], na.rm=na.rm),
                       sd   = sd     (xx[[col]], na.rm=na.rm)
                     )
                   },
                   measurevar
    )
    
    # Rename the "mean" column
    colnames(datac)[which(colnames(datac)=='mean')]=measurevar
    # datac <- rename(datac, c("mean" = measurevar))
    
    datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
    
    # Confidence interval multiplier for standard error
    # Calculate t-statistic for confidence interval:
    # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
    ciMult <- qt(conf.interval/2 + .5, datac$N-1)
    datac$ci <- datac$se * ciMult
    
    return(datac)
  }

ugo=load(paste0('Rdata/RBP_simulations.montecarlo_TCGA__purity_strata','high','.Rdata'))
rm(samples_boot,r2)
SRPs[which(SRPs=='HNRPLL')]='HNRNPLL'


relaimpo_tcga_list = list()
relaimpo_gtex_list = list()
coe_tcga_list = list()
coe_gtex_list = list()

DE_SRGs_tcga_list = list()
DE_SRGs_gtex_list = list()

set.seed(8890)
for (n in 1:100) {
  
  print(n)
  
  relaimpo_tcga=matrix(NA,nrow = 10,ncol = 4)
  colnames(relaimpo_tcga)=c('FOXA1.lmg','MYC.lmg','AR.lmg','ERG.lmg')
  rownames(relaimpo_tcga)=seq(0,0.9,0.1)
  relaimpo_gtex=matrix(NA,nrow = 10,ncol = 4)
  colnames(relaimpo_gtex)=c('FOXA1.lmg','MYC.lmg','AR.lmg','ERG.lmg')
  rownames(relaimpo_gtex)=seq(0,0.9,0.1)
  
  coe_tcga=matrix(NA,nrow = 10,ncol = 4)
  colnames(coe_tcga)=c('FOXA1','MYC','AR','ERG')
  rownames(coe_tcga)=seq(0,0.9,0.1)
  coe_gtex=matrix(NA,nrow = 10,ncol = 4)
  colnames(coe_gtex)=c('FOXA1','MYC','AR','ERG')
  rownames(coe_gtex)=seq(0,0.9,0.1)
  
  for (cont_perc in seq(0,0.9,0.1)) {
    
    print(cont_perc)
    
    simulated_tumors_tcga = prad_pure
    simulated_tumors_gtex = prad_pure
    
    simulated_tumors_gtex = simulated_tumors_gtex*(1-cont_perc)
    for (i in 1:ncol(simulated_tumors_gtex)) {
      simulated_tumors_gtex[,i]=simulated_tumors_gtex[,i]+prostate[,sample(c(1:ncol(prostate)),size = 1,replace = T)]*cont_perc
    }
    
    simulated_tumors_gtex_length=c$Length[match(rownames(simulated_tumors_gtex),c$gene_name)]
    simulated_tumors_gtex_tpms = countToTpm2(simulated_tumors_gtex, simulated_tumors_gtex_length)
    
    for (i in c('gtex')) {
      
      pl = readRDS("Rdata/KEGG_Genetic_information_process.GSECA.200105.rds")
      
      if (i=='tcga') {
        rna_sub=simulated_tumors_tcga_tpms
      } else if (i=='gtex') {
        rna_sub=simulated_tumors_gtex_tpms
      }
      
      rna_sub=melt(rna_sub)
      colnames(rna_sub)=c('symbol','sample','TPM')
      
      srp.p = ddply(subset(rna_sub, symbol%in%pl$SPLICING_RELATED), .(sample),summarise, tpm = sum(TPM))
      
      p = dcast(subset(rna_sub, symbol %in%c("FOXA1",'AR',"ERG","MYC")), sample~symbol, value.var = 'TPM')
      
      library(caret)
      p$SRP = srp.p$tpm[match(p$sample,srp.p$sample)]
      rownames(p) = p$sample; p = p[,-1]
      
      x   = predict(preProcess(p, method = c("center", "scale", "YeoJohnson", "nzv"))
                    , newdata = p)
      
      library(relaimpo)
      
      get_coef = function(f){
        coe = coef(f); coe  = coe[-1]; coe = sort(coe)
        pvalue = coef(summary(f))[,4]
        coe = data.frame(g = names(coe), coe = coe, pval = pvalue[names(coe)])
        coe$g = factor(coe$g, levels = coe$g)
        coe$pval[coe$pval>0.05] = 1
        coe
      }
      
      # fit.p = glm(data=x, SRP~FOXA1+MYC+AR+ERG ); coe1 = get_coef(fit.p)
      
      
      if (i=='tcga') {
        coe_tcga[as.character(cont_perc),] = coe1[c('FOXA1','MYC','AR','ERG'),'coe']
      } else if (i=='gtex') {
        coe_gtex[as.character(cont_perc),] = coe1[c('FOXA1','MYC','AR','ERG'),'coe']
      }
      
      
      boot <- boot.relimp(fit.p, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
      b = booteval.relimp(boot,sort=TRUE) # print result
      
      # plot(b)
      
      if (i=='tcga') {
        relaimpo_tcga[as.character(cont_perc),] = b@est[c('FOXA1.lmg','MYC.lmg','AR.lmg','ERG.lmg')]
      } else if (i=='gtex') {
        relaimpo_gtex[as.character(cont_perc),] = b@est[c('FOXA1.lmg','MYC.lmg','AR.lmg','ERG.lmg')]
      }
      
      
      # Differential SRG expression ---
      
      
      if (i=='tcga') {
        r2 = melt(simulated_tumors_tcga_tpms)
      } else if (i=='gtex') {
        r2 = melt(simulated_tumors_gtex_tpms)
      }
      
      colnames(r2)=c('symbol','Tumor_Sample_Barcode','TPM')
      tmp = subset(r2,symbol=='FOXA1')
      qt=quantile(tmp$TPM)
      tmp$FOXA1_overexpression=F
      tmp$FOXA1_overexpression[which(tmp$TPM>=qt['75%'])]=T
      r2$FOXA1_overexpression = tmp$FOXA1_overexpression[match(r2$Tumor_Sample_Barcode,tmp$Tumor_Sample_Barcode)]
      rm(tmp)
      x = subset(r2, symbol %in% SRPs)
      
      t = summarySE(x, measurevar = "TPM", groupvars = c('symbol','FOXA1_overexpression'), na.rm = T)
      
      t1 = ddply(x, .(symbol),
                 summarise,
                 median.TPM    = median(TPM[FOXA1_overexpression], na.rm=T),
                 median.TPM.wt = median(TPM[!FOXA1_overexpression], na.rm=T),
                 ks=ks.test(    TPM[FOXA1_overexpression], TPM[!FOXA1_overexpression], alternative = 'two.sided' )$p.value # %%%%%%%%%%%%%%%%%%%%%% CAMBIATO
      )
      
      t1$BF.ks  = p.adjust(t1$ks, 'bonferroni', n=sum(!is.na(t1$ks)))
      t1$FDR.ks = p.adjust(t1$ks, 'fdr', n=sum(!is.na(t1$ks)))
      
      t1$emp=0
      
      t = cbind(t, t1[match(t$symbol, t1$symbol),])
      tc = subset(t,FOXA1_overexpression)
      tc$TPM.wt = subset(t, !FOXA1_overexpression)$TPM
      tc$se.wt  = subset(t, !FOXA1_overexpression)$se
      tc$delta = with(tc, TPM-TPM.wt)
      tc$L2R = log2(tc$median.TPM)-log2(tc$median.TPM.wt)
      tc$FC  = tc$median.TPM/tc$median.TPM.wt
      
      tc$DEG.sign_FDR = with(tc,
                             FDR.ks<=0.01 &
                               emp<=0.01 & 
                               abs(L2R) >= 0.2
      )
      
      tc$direction = "DOWN"
      tc$direction[which(tc$delta>0)] = "UP"
      
      tc$cont_perc = cont_perc
      tc$cont_perc=as.character(tc$cont_perc)
      
      if (i=='tcga' & cont_perc==0) {
        tc_tot_tcga = tc
      } else if (i=='tcga' & !cont_perc==0) {
        tc_tot_tcga = rbind(tc_tot_tcga,tc)
      } else if (i=='gtex' & cont_perc==0) {
        tc_tot_gtex = tc
      } else if (i=='gtex' & !cont_perc==0) {
        tc_tot_gtex = rbind(tc_tot_gtex,tc)
      }
      
    }
    
  }
  
  relaimpo_gtex_list[[n]] = relaimpo_gtex
  coe_gtex_list[[n]] = coe_gtex
  
  DE_SRGs_gtex_list[[n]]=tc_tot_gtex
  
  rm(tc_tot_tcga,tc_tot_gtex)
  
}



relaimpo_gtex = Reduce(`+`, relaimpo_gtex_list) / length(relaimpo_gtex_list)
coe_gtex = Reduce(`+`, coe_gtex_list) / length(coe_gtex_list)

save(relaimpo_gtex,coe_gtex,file='Rdata/in_silico_contamination_with_100_bootstrapping.Rdata')
save(relaimpo_gtex_list,coe_gtex_list,file='Rdata/in_silico_contamination_with_100_bootstrapping_Lists.Rdata')

save(DE_SRGs_gtex_list,file='Rdata/in_silico_contamination_SRG_DIFERENTIAL_EXPRESSION_with_100_bootstrapping_Lists.Rdata')



# Summary of GLMs and TFs expression using line plots ====

library(RColorBrewer)
library(wesanderson)
library(ggsci)

ugo=load('Rdata/in_silico_contamination_with_100_bootstrapping_Lists.Rdata')

tmp=do.call(rbind,relaimpo_gtex_list)
for(i in 1:100) {
  v=rep(i,10)
  if (i==1) vrep=v else vrep=c(vrep,v)
}
tmp=as.data.frame(tmp)
tmp$rep=vrep
tmp$cont_perc=rep(c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9),10)
tmp=ddply(tmp,.(cont_perc),summarise,mean_FOXA1=mean(FOXA1.lmg),mean_MYC=mean(MYC.lmg),mean_AR=mean(AR.lmg),mean_ERG=mean(ERG.lmg),
          sd_FOXA1=sd(FOXA1.lmg),sd_MYC=sd(MYC.lmg),sd_AR=sd(AR.lmg),sd_ERG=sd(ERG.lmg))
to_plot=tmp[,1:5]
to_plot=melt(to_plot,id.vars = 'cont_perc')
tmp=tmp[,c(1,6:9)]
tmp=melt(tmp,id.vars = 'cont_perc')
to_plot$sd=tmp$value
to_plot$variable=gsub(to_plot$variable,pattern = 'mean_',replacement = '')
to_plot$variable=factor(to_plot$variable,levels=c("MYC","AR","ERG","FOXA1"))
p2=ggplot(to_plot,aes(x=(1-cont_perc),y=value,col=variable))+geom_line()+theme_bw()+geom_point()+geom_errorbar(mapping = aes(ymin = (value-sd/10), ymax = (value+sd/10)), width=0.05)+ggtitle('GTEX')+ylab('% of R2')+xlab('purity')+
  scale_color_jco()
p22=ggplot(to_plot,aes(x=(1-cont_perc),y=value,col=variable))+geom_line()+theme_bw()+geom_point()+geom_errorbar(mapping = aes(ymin = (value-sd), ymax = (value+sd)), width=0.05)+ggtitle('GTEX')+ylab('% of R2')+xlab('purity')+
  scale_color_jco()

pdf(paste0(FIG_DIR2, file="summary_glms_with_100_bootstrap.pdf"), height=unit(5,'cm'), width=unit(12,'cm'), useDingbats = F)
print(grid.arrange(p2,p22,ncol=2))
dev.off()





# HEATMAP OF DE SRGs  with GTEX normal contamination ----

library(ComplexHeatmap)

tt=load('Rdata/in_silico_contamination_SRG_DIFERENTIAL_EXPRESSION_with_100_bootstrapping_Lists.Rdata')

prova <- lapply(DE_SRGs_gtex_list, function(x) ddply(x, .(cont_perc), summarise, sum_de=sum(DEG.sign_FDR)))

prova = lapply(prova, function(x) {rownames(x)=x$cont_perc; return(x)})

prova = lapply(prova, function(x) {y=as.data.frame(x$sum_de); rownames(y)=rownames(x); return(y)})

prova = do.call(cbind, prova)

colnames(prova)=seq(1,100)
prova$cont_perc=rownames(prova)

prova2 = as.data.frame(melt(prova,id.vars = 'cont_perc'))


# ROW ANNO ---

prova3 = do.call(rbind, DE_SRGs_gtex_list)

prova3 = ddply(prova3,.(symbol,cont_perc),summarise,n_de=sum(DEG.sign_FDR))

prova4 = as.data.frame(dcast(prova3,formula = symbol~cont_perc,value.var = 'n_de'))
rownames(prova4)=prova4$symbol
prova4$symbol=NULL


# COLUMN ANNO ---

prova5 <- lapply(DE_SRGs_gtex_list, function(x) ddply(x, .(symbol), summarise, sum_de=sum(DEG.sign_FDR)))

prova5 = lapply(prova5, function(x) {rownames(x)=x$symbol; return(x)})

prova5 = lapply(prova5, function(x) {y=as.data.frame(x$sum_de); rownames(y)=rownames(x); return(y)})

prova5 = do.call(cbind, prova5)

colnames(prova5)=seq(1,100)
prova5$symbol=rownames(prova5)

# prova6 = as.data.frame(melt(prova5,id.vars = 'symbol'))


# PLOT HEATMAP ---

m = as.matrix(prova[,!colnames(prova)=='cont_perc'])
haB = rowAnnotation(number_DE_SRGs = anno_barplot(rowMeans(m), height = unit(4, "cm")))

library(genefilter)
n_DE_SRGs_means=rowMeans(m)
n_DE_SRGs_sds=rowSds(m)
to_plot=as.data.frame(cbind(n_DE_SRGs_means,n_DE_SRGs_sds))
colnames(to_plot)=c('mean','sd')
to_plot$cont_perc=as.numeric(rownames(to_plot))
pdf(paste0(FIG_DIR2, file="barplot_N_DE_SRGs_annotation_GTEX_with_100_bootstrap_with_errorbar.pdf"), height=unit(4,'cm'), width=unit(5,'cm'), useDingbats = F)
print(ggplot(to_plot,aes(x=(1-cont_perc),y=mean))+geom_col()+theme_bw()+geom_errorbar(mapping = aes(ymin = (mean-sd/10), ymax = (mean+sd/10)), width=0.05))
dev.off()
pdf(paste0(FIG_DIR2, file="barplot_N_DE_SRGs_annotation_GTEX_with_100_bootstrap_with_errorbar_SD.pdf"), height=unit(4,'cm'), width=unit(5,'cm'), useDingbats = F)
print(ggplot(to_plot,aes(x=(1-cont_perc),y=mean))+geom_col()+theme_bw()+geom_errorbar(mapping = aes(ymin = (mean-sd), ymax = (mean+sd)), width=0.05))
dev.off()


col_fun = circlize::colorRamp2(c(0, 50, 63, 83, 100), c("black", "red", "orange", "gold", "yellow"))

heat_red = prova4[rowMeans(prova4[,!colnames(prova4)=='symbol'])>5,]
prova5_red = prova5[rownames(heat_red),]
m = as.matrix(t(prova5_red[,!colnames(prova5_red)=='symbol']))
ha2_redB = HeatmapAnnotation(number_of_times_DE = anno_barplot(colMeans(m), height = unit(4, "cm")))

pdf(paste0(FIG_DIR2, file="DE_SRGs_heatmap_GTEX_with_100_bootstrap__anno_barplots.pdf"), height=unit(8,'cm'), width=unit(26,'cm'), useDingbats = F)
ht=draw(Heatmap(t(heat_red),cluster_rows = F,right_annotation = haB,top_annotation = ha2_redB, col = col_fun, rect_gp = gpar(col = "darkgrey", lwd = 0.5)))
print(ht)
dev.off()



library(genefilter)
n_times_DE_means=rowMeans(t(m))
n_times_DE_sds=rowSds(t(m))
to_plot=as.data.frame(cbind(n_times_DE_means,n_times_DE_sds))
colnames(to_plot)=c('mean','sd')
to_plot$SRG=rownames(to_plot)

ord=column_order(ht)
to_plot=to_plot[rev(ord),]
to_plot$SRG=factor(to_plot$SRG,levels=to_plot$SRG)

pdf(paste0(FIG_DIR2, file="barplot_N_times_DE_annotation_GTEX_with_100_bootstrap_with_errorbar.pdf"), height=unit(12,'cm'), width=unit(4,'cm'), useDingbats = F)
print(ggplot(to_plot,aes(y=(SRG),x=mean))+geom_col()+theme_bw()+geom_errorbar(mapping = aes(xmin = (mean-sd/10), xmax = (mean+sd/10)), height=0.05))
dev.off()
pdf(paste0(FIG_DIR2, file="barplot_N_times_DE_annotation_GTEX_with_100_bootstrap_with_errorbar_SD.pdf"), height=unit(12,'cm'), width=unit(4,'cm'), useDingbats = F)
print(ggplot(to_plot,aes(y=(SRG),x=mean))+geom_col()+theme_bw()+geom_errorbar(mapping = aes(xmin = (mean-sd), xmax = (mean+sd)), height=0.05))
dev.off()





# ******************** ----
# ******************** xCell scores in TCGA primary PC from RNAseqDB ----
# ******************** ----

source("sources/utils_RNA.R")
c=readRDS("Rdata/gene_length.rds")

rna = readRDS('Rdata/RNA_TPM_tumor.rds')

prad=read.delim2('Tables/RNAseqDB/prad-rsem-count-tcga-t.txt',stringsAsFactors = F)
rownames_prad=prad$Hugo_Symbol
prad$Hugo_Symbol=NULL
prad$Entrez_Gene_Id=NULL
prad=apply(prad, 2, as.numeric)
rownames(prad)=rownames_prad
prad=apply(prad, 2, function(x) x/sum(x)*1000000)
colnames(prad)=gsub(colnames(prad),pattern = '[.]',replacement = '-')
prad_length=c$Length[match(rownames(prad),c$gene_name)]
prad_tpms = countToTpm2(prad, prad_length)
prad_tpms=reshape2::melt(prad_tpms)
colnames(prad_tpms)=c('symbol','sample','TPM')

prad_tpms=subset(prad_tpms,sample%in%rna$Tumor_Sample_Barcode)

xcell_primary_PC_RNAseqDB_original=readRDS('Rdata/xCell_primary_PC_RNAseqDB.rds')

prad_tpms$ImmuneScore=xcell_primary_PC_RNAseqDB_original[prad_tpms$sample,'ImmuneScore']
prad_tpms$StromaScore=xcell_primary_PC_RNAseqDB_original[prad_tpms$sample,'StromaScore']
prad_tpms$MicroenvironmentScore=xcell_primary_PC_RNAseqDB_original[prad_tpms$sample,'MicroenvironmentScore']
prad_tpms$EpithelialCells=xcell_primary_PC_RNAseqDB_original[prad_tpms$sample,'Epithelial cells']
prad_tpms$SmoothMuscle=xcell_primary_PC_RNAseqDB_original[prad_tpms$sample,'Smooth muscle']

tmp=subset(prad_tpms,symbol=='FOXA1')
qt=quantile(tmp$TPM)
tmp$FOXA1_overexpression = cut(tmp$TPM,breaks = c(qt['0%'],qt['75%'],qt['100%']),include.lowest = T,labels = c('F','T'))
prad_tpms$FOXA1_overexpression = tmp$FOXA1_overexpression[match(prad_tpms$sample,tmp$sample)]

p1=ggplot()+geom_point(subset(prad_tpms,symbol=='FOXA1'),mapping=aes(x=TPM,y=ImmuneScore))+geom_smooth(subset(prad_tpms,symbol=='FOXA1'),mapping=aes(x=TPM,y=ImmuneScore),method = 'lm')+stat_cor(subset(prad_tpms,symbol=='FOXA1'),mapping=aes(x=TPM,y=ImmuneScore))+theme_bw()+xlab('FOXA1 (TPM)')+ggtitle('Immune score')
p2=ggplot()+geom_point(subset(prad_tpms,symbol=='FOXA1'),mapping=aes(x=TPM,y=StromaScore))+geom_smooth(subset(prad_tpms,symbol=='FOXA1'),mapping=aes(x=TPM,y=StromaScore),method = 'lm')+stat_cor(subset(prad_tpms,symbol=='FOXA1'),mapping=aes(x=TPM,y=StromaScore))+theme_bw()+xlab('FOXA1 (TPM)')+ggtitle('Stroma score')
p3=ggplot()+geom_point(subset(prad_tpms,symbol=='FOXA1'),mapping=aes(x=TPM,y=EpithelialCells))+geom_smooth(subset(prad_tpms,symbol=='FOXA1'),mapping=aes(x=TPM,y=EpithelialCells),method = 'lm')+stat_cor(subset(prad_tpms,symbol=='FOXA1'),mapping=aes(x=TPM,y=EpithelialCells))+theme_bw()+xlab('FOXA1 (TPM)')+ggtitle('EpithelialCells')
p4=ggplot(subset(prad_tpms,symbol=='FOXA1'),aes(x=FOXA1_overexpression,y=ImmuneScore))+geom_boxplot(notch = T)+stat_compare_means()+theme_bw()
p5=ggplot(subset(prad_tpms,symbol=='FOXA1'),aes(x=FOXA1_overexpression,y=StromaScore))+geom_boxplot(notch = T)+stat_compare_means()+theme_bw()
p6=ggplot(subset(prad_tpms,symbol=='FOXA1'),aes(x=FOXA1_overexpression,y=EpithelialCells))+geom_boxplot(notch = T)+stat_compare_means()+theme_bw()
pdf(paste0(FIG_DIR, file="XCELL_TCGA_scatter_and_boxplot_scores_vs_FOXA1__from_RNAseqDB.pdf"), height=unit(8,'cm'), width=unit(12,'cm'), useDingbats = F)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
dev.off()




# GLM STRATA ON XCELL SCORES  SELECTING THE MOST "CONTAMINATED" SAMPLES ----

pl = readRDS("Rdata/KEGG_Genetic_information_process.GSECA.200105.rds")

qt=quantile(subset(prad_tpms,symbol=='FOXA1')$StromaScore,seq(0,1,0.05))
prad_tpms$StromaScore_HighRest = 'rest'
prad_tpms$StromaScore_HighRest[which(prad_tpms$StromaScore>=qt['75%'])] = 'high'

qt=quantile(subset(prad_tpms,symbol=='FOXA1')$EpithelialCells,seq(0,1,0.05))
prad_tpms$EpithelialCells_HighRest = 'rest'
prad_tpms$EpithelialCells_HighRest[which(prad_tpms$EpithelialCells>=qt['75%'])] = 'high'

qt=quantile(subset(prad_tpms,symbol=='FOXA1')$ImmuneScore,seq(0,1,0.05))
prad_tpms$ImmuneScore_HighRest = 'rest'
prad_tpms$ImmuneScore_HighRest[which(prad_tpms$ImmuneScore>=qt['75%'])] = 'high'



# StromaScore ----

for(i in c('rest','high')) {
  
  if (i=='all') {
    rna_sub=prad_tpms
  } else if (i=='high' | i=='rest') {
    rna_sub=subset(prad_tpms,StromaScore_HighRest==i) } else if (i=='not_cont_by_stroma') {
      rna_sub=subset(prad_tpms,StromaScore_pvalue>0.05)
    } else if (i=='cont_by_stroma') {
      rna_sub=subset(prad_tpms,StromaScore_pvalue<0.05)
    }
  
  srp.p = ddply(subset(rna_sub, symbol%in%pl$SPLICING_RELATED), .(sample),summarise, tpm = sum(TPM))
  
  p = dcast(subset(rna_sub, symbol %in%c("FOXA1",'AR',"ERG","MYC")), sample~symbol, value.var = 'TPM')
  
  library(caret)
  p$SRP = srp.p$tpm[match(p$sample,srp.p$sample)]
  rownames(p) = p$sample; p = p[,-1]
  
  x   = predict(preProcess(p, method = c("center", "scale", "YeoJohnson", "nzv"))
                , newdata = p)
  
  library(relaimpo)
  
  get_coef = function(f){
    coe = coef(f); coe  = coe[-1]; coe = sort(coe)
    pvalue = coef(summary(f))[,4]
    coe = data.frame(g = names(coe), coe = coe, pval = pvalue[names(coe)])
    coe$g = factor(coe$g, levels = coe$g)
    coe$pval[coe$pval>0.05] = 1
    print(ggplot(coe, aes(x=g, y=coe, fill = coe, size=pval<0.05) )+
            geom_bar(stat='identity',color='black')+theme_bw()+coord_flip()+scale_fill_gradientn(limits=c(0,0.3),colours = c('white','red'))+#scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red")+
            scale_size_manual(values = c(0,1)))
    coe
  }
  
  set.seed(30580)
  
  pdf(paste0(FIG_DIR, file="SRP_TF_primary_glm_coeff_NEW__StromaScore_",i,"__from_RNAseqDB.pdf"), height=unit(3,'cm'), width=unit(3, 'cm'), useDingbats = F )
  fit.p = glm(data=x, SRP~FOXA1+MYC+AR+ERG ); coe1 = get_coef(fit.p)
  dev.off()
  
  set.seed(30580)
  boot <- boot.relimp(fit.p, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
  b = booteval.relimp(boot,sort=TRUE) # print result
  
  pdf(paste0(FIG_DIR, file="RelImpo_glm_primary__StromaScore_",i,"__from_RNAseqDB.pdf"), height=unit(5,'cm'), width=unit(5,'cm'), useDingbats = F)
  par(las=2)
  plot(b)
  dev.off()
}


# EpithelialCells ----

for(i in c('rest','high')) {
  
  if (i=='all') {
    rna_sub=prad_tpms
  } else if (i=='high' | i=='rest') {
    rna_sub=subset(prad_tpms,EpithelialCells_HighRest==i) } else if (i=='not_cont_by_epithelial') {
      rna_sub=subset(prad_tpms,EpithelialCells_pvalue>0.05)
    } else if (i=='cont_by_epithelial') {
      rna_sub=subset(prad_tpms,EpithelialCells_pvalue<0.05)
    }
  
  srp.p = ddply(subset(rna_sub, symbol%in%pl$SPLICING_RELATED), .(sample),summarise, tpm = sum(TPM))
  
  p = dcast(subset(rna_sub, symbol %in%c("FOXA1",'AR',"ERG","MYC")), sample~symbol, value.var = 'TPM')
  
  library(caret)
  p$SRP = srp.p$tpm[match(p$sample,srp.p$sample)]
  rownames(p) = p$sample; p = p[,-1]
  
  x   = predict(preProcess(p, method = c("center", "scale", "YeoJohnson", "nzv"))
                , newdata = p)
  
  library(relaimpo)
  
  get_coef = function(f){
    coe = coef(f); coe  = coe[-1]; coe = sort(coe)
    pvalue = coef(summary(f))[,4]
    coe = data.frame(g = names(coe), coe = coe, pval = pvalue[names(coe)])
    coe$g = factor(coe$g, levels = coe$g)
    coe$pval[coe$pval>0.05] = 1
    print(ggplot(coe, aes(x=g, y=coe, fill = coe, size=pval<0.05) )+
            geom_bar(stat='identity',color='black')+theme_bw()+coord_flip()+scale_fill_gradientn(limits=c(0,0.3),colours = c('white','red'))+#scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red")+
            scale_size_manual(values = c(0,1)))
    coe
  }
  
  set.seed(30580)
  
  pdf(paste0(FIG_DIR, file="SRP_TF_primary_glm_coeff_NEW__EpithelialCellsScore_",i,"__from_RNAseqDB.pdf"), height=unit(3,'cm'), width=unit(3, 'cm'), useDingbats = F )
  fit.p = glm(data=x, SRP~FOXA1+MYC+AR+ERG ); coe1 = get_coef(fit.p)
  dev.off()
  
  set.seed(30580)
  boot <- boot.relimp(fit.p, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
  b = booteval.relimp(boot,sort=TRUE) # print result
  
  pdf(paste0(FIG_DIR, file="RelImpo_glm_primary__EpithelialCellsScore_",i,"__from_RNAseqDB.pdf"), height=unit(5,'cm'), width=unit(5,'cm'), useDingbats = F)
  par(las=2)
  plot(b)
  dev.off()
}


# ImmuneScore ----

for(i in c('rest','high')) {
  
  if (i=='all') {
    rna_sub=prad_tpms
  } else if (i=='high' | i=='rest') {
    rna_sub=subset(prad_tpms,ImmuneScore_HighRest==i) } else if (i=='not_cont_by_immune') {
      rna_sub=subset(prad_tpms,ImmuneScore_pvalue>0.05)
    } else if (i=='cont_by_immune') {
      rna_sub=subset(prad_tpms,ImmuneScore_pvalue<0.05)
    }
  
  srp.p = ddply(subset(rna_sub, symbol%in%pl$SPLICING_RELATED), .(sample),summarise, tpm = sum(TPM))
  
  p = dcast(subset(rna_sub, symbol %in%c("FOXA1",'AR',"ERG","MYC")), sample~symbol, value.var = 'TPM')
  
  library(caret)
  p$SRP = srp.p$tpm[match(p$sample,srp.p$sample)]
  rownames(p) = p$sample; p = p[,-1]
  
  x   = predict(preProcess(p, method = c("center", "scale", "YeoJohnson", "nzv"))
                , newdata = p)
  
  library(relaimpo)
  
  get_coef = function(f){
    coe = coef(f); coe  = coe[-1]; coe = sort(coe)
    pvalue = coef(summary(f))[,4]
    coe = data.frame(g = names(coe), coe = coe, pval = pvalue[names(coe)])
    coe$g = factor(coe$g, levels = coe$g)
    coe$pval[coe$pval>0.05] = 1
    print(ggplot(coe, aes(x=g, y=coe, fill = coe, size=pval<0.05) )+
            geom_bar(stat='identity',color='black')+theme_bw()+coord_flip()+scale_fill_gradientn(limits=c(0,0.3),colours = c('white','red'))+#scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red")+
            scale_size_manual(values = c(0,1)))
    coe
  }
  
  set.seed(30580)
  
  pdf(file=paste0(FIG_DIR, "SRP_TF_primary_glm_coeff_NEW__ImmuneScore_",i,"__from_RNAseqDB.pdf"), height=unit(3,'cm'), width=unit(3, 'cm'), useDingbats = F )
  fit.p = glm(data=x, SRP~FOXA1+MYC+AR+ERG ); coe1 = get_coef(fit.p)
  dev.off()
  
  set.seed(30580)
  boot <- boot.relimp(fit.p, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
  b = booteval.relimp(boot,sort=TRUE) # print result
  
  pdf(file=paste0(FIG_DIR, "RelImpo_glm_primary__ImmuneScore_",i,"__from_RNAseqDB.pdf"), height=unit(5,'cm'), width=unit(5,'cm'), useDingbats = F)
  par(las=2)
  plot(b)
  dev.off()
}






# Merge high and low purity ORA results ----



to_plot=readRDS(paste0('Rdata/ORA_all_event_types_TCGA_KEGG_186_', 'high', '.rds'))

to_plot_high=to_plot
to_plot_high$purity='high'

to_plot=readRDS(paste0('Rdata/ORA_all_event_types_TCGA_KEGG_186_', 'low', '.rds'))

to_plot_low=to_plot
to_plot_low$purity='low'
to_plot=rbind(to_plot_high,to_plot_low)

to_plot=subset(to_plot,p.adjust<0.1)

pdf(paste0(FIG_DIR, file='ORA_all_event_types_TCGA_KEGG_186__KEGG_186__HIGH_AND_LOW_PURITY_MERGED.pdf'), height = unit(4,'cm'), width = unit(8, 'cm'), useDingbats = F)
to_plot2=to_plot[order(to_plot$p.adjust,decreasing = T),]
to_plot2$ID=factor(to_plot2$ID,levels=unique(to_plot2$ID))
ggplot( to_plot2,
        aes(x= parse_ratio(GeneRatio), y= ID)) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=factor(p.adjust<=.25), fill=p.adjust, size = Count, shape=purity)) +
  geom_text(aes(label=rank)) +
  scale_shape_manual(values=c(21,24))+
  scale_color_manual(values=c('transparent','black'))+
  scale_fill_viridis_c(option = 'D',guide=guide_colorbar(reverse=T, draw.llim = T), direction = -1)+#,
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + xlab("Gene Ratio") +ylab(NULL)+
  facet_wrap(~type, nrow=1)
dev.off()


