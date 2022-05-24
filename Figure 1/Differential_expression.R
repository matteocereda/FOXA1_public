# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Differential expression analysis ----


# ******************** ----
# ******************** DEA in primary and metastatic samples ----
# ******************** ----


source("Scripts/config.R")

# ∞∞∞ EDGR ======================

# query <- GDCquery(project = "TCGA-PRAD",
#                   data.category = "Gene expression",
#                   data.type = "Gene expression quantification",
#                   experimental.strategy = "RNA-Seq",
#                   platform = "Illumina HiSeq",
#                   file.type = "results",
#                   barcode = unique(rna$Tumor_Sample_Barcode),
#                   legacy = TRUE)
# # GDCdownload(query)
# prad <- GDCprepare(query,save = TRUE, save.filename = "Rdata/PRAD.Rdata")

load("Rdata/PRAD.Rdata")# load("Rdata/PRAD.Rdata")
rna=readRDS("Rdata/RNA_TPM_tumor.rds")# load("Prostate/RNA_TPM_tumor.Rdata")
load("Rdata/KEGG_list.148.SPLICING_RELATED.Rdata") 

library(TCGAbiolinks)
library(SummarizedExperiment)


tmp=subset(rna,symbol=='FOXA1')
qt=quantile(tmp$TPM)
tmp$FOXA1_overexpression=F
tmp$FOXA1_overexpression[which(tmp$TPM>=qt['75%'])]=T
rna$FOXA1_overexpression=tmp$FOXA1_overexpression[match(rna$Tumor_Sample_Barcode,tmp$Tumor_Sample_Barcode)]

rc <- assay(data, "raw_count")
rc.norm <- TCGAanalyze_Normalization(tabDF = rc, geneInfo =  TCGAbiolinks::geneInfo)

# quantile filter of genes
rc.norm.filt <- TCGAanalyze_Filtering(tabDF = rc.norm,
                                      method = "quantile",
                                      qnt.cut =  0.1)

# selection of normal samples "NT"
cntr <- unique( subset(rna, !FOXA1_overexpression)$Tumor_Sample_Barcode)

# selection of tumor samples "TP"
case <- unique( subset(rna, FOXA1_overexpression)$Tumor_Sample_Barcode)

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

save(rc.norm, rc.norm.filt, dataDEGs, file="Rdata/PRAD_EDGR_FOXA1_75th.Rdata")



# ∞∞∞ DESEQ ======================

get_DESEQ = function(counts, samples, cond, ctrl){ 
  
  # Perform DESeq analysis
  require(DESeq2)
  require(BiocParallel)
  
  tmp=counts[,samples]
  colData <- data.frame(condition=factor(cond))
  colData$condition <- relevel(colData$condition, ref = ctrl) 
  dds <- DESeqDataSetFromMatrix(tmp, colData, formula(~condition))
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds, parallel = T, BPPARAM = MulticoreParam(workers = 2))
  res <- results(dds)
  res
  
}

# 1. Load and prepare PRAD raw counts ====

library(TCGAbiolinks)

ttt=load('Rdata/PRAD.Rdata')
prad=assays(data)$raw_count

library(reshape2)
prad <- melt(prad)
colnames(prad) <- c("symbol","Tumor_Sample_Barcode","value")

rna = readRDS('Rdata/RNA_TPM_tumor.rds')
xxx=subset(rna,symbol=='FOXA1')
rm(rna)

qt=quantile(xxx$TPM)
xxx$FOXA1_overexpression=F
xxx$FOXA1_overexpression[which(xxx$TPM>=qt['75%'])]=T

prad$type = "CNTR"
prad$type[which(prad$Tumor_Sample_Barcode%in%subset(xxx,FOXA1_overexpression)$Tumor_Sample_Barcode)] = "CASE"

prad$value=round(prad$value)

prad.rwc <- dcast(prad, symbol ~ Tumor_Sample_Barcode, value.var = "value", fun.aggregate = mean)
rownames(prad.rwc) <- as.character(prad.rwc$symbol)
prad.rwc$symbol <- NULL

case = unique(subset(prad,type=='CASE')$Tumor_Sample_Barcode)
cntr = unique(subset(prad,type=='CNTR')$Tumor_Sample_Barcode)

pheno = c("CASE","CNTR")
pheno_nn = c(length(case),length(cntr))
prad.rwc = prad.rwc[,c(case,cntr)]

# 2. Run DESeq ====

myc.deout2 <- get_DESEQ(counts = prad.rwc, samples = c(case,cntr), 
                        cond = c(rep(pheno[1],pheno_nn[1]),rep(pheno[2],pheno_nn[2])),
                        ctrl = pheno[2]) # , tlen = ptlen

save(myc.deout2, file="Rdata/PRAD_DESEQ_FOXA1_75th.Rdata")




# SUBSET ON SRPs ========

library(clusterProfiler)
library(org.Hs.eg.db)

pl = readRDS('Rdata/KEGG_Genetic_information_process.GSECA.200105.rds')


# EDGER ----

edger=load('Rdata/PRAD_EDGR_FOXA1_75th.Rdata')

pl$SPLICING_RELATED[!pl$SPLICING_RELATED%in%rownames(dataDEGs)] # "PABPC1" "ZNF638" "QKI"    "SRSF8" 
rownames(dataDEGs)[which(rownames(dataDEGs)=='SFRS8')]='SRSF8'
pl$SPLICING_RELATED[!pl$SPLICING_RELATED%in%rownames(dataDEGs)] #"PABPC1" "ZNF638" "QKI"  
rownames(dataDEGs)[which(rownames(dataDEGs)=='QK')]='QKI'
pl$SPLICING_RELATED[!pl$SPLICING_RELATED%in%rownames(dataDEGs)] #"PABPC1" "ZNF638"
rownames(dataDEGs)[which(rownames(dataDEGs)=='PABPC2')]='PABPC1'
pl$SPLICING_RELATED[!pl$SPLICING_RELATED%in%rownames(dataDEGs)] #"ZNF638"
rownames(dataDEGs)[which(rownames(dataDEGs)=='ZFP638')]='ZNF638'
pl$SPLICING_RELATED[!pl$SPLICING_RELATED%in%rownames(dataDEGs)] #


dataDEGs = subset(dataDEGs,rownames(dataDEGs)%in%pl$SPLICING_RELATED)


saveRDS(dataDEGs,'Rdata/PRAD_EDGR_FOXA1_75th_148_SRPs.rds')


# DESEQ ----

deseq=load('Rdata/PRAD_DESEQ_FOXA1_75th.Rdata')

dataDEGs = as.data.frame(myc.deout2)
dataDEGs$symbol = sapply(strsplit(rownames(dataDEGs),'[|]'),`[`,1)

pl$SPLICING_RELATED[!pl$SPLICING_RELATED%in%dataDEGs$symbol] # 
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


dataDEGs = subset(dataDEGs,symbol%in%pl$SPLICING_RELATED)

saveRDS(dataDEGs,'Rdata/PRAD_DESEQ_FOXA1_75th_148_SRPs.rds')






# ∞∞∞ TCGA - load datasets and libraries ----

rna = readRDS('Rdata/RNA_TPM_tumor.rds')
pl = readRDS('Rdata/KEGG_Genetic_information_process.GSECA.200105.rds')
deseq=readRDS('Rdata/PRAD_DESEQ_FOXA1_75th_148_SRPs.rds')
dataDEGs=readRDS('Rdata/PRAD_EDGR_FOXA1_75th_148_SRPs.rds')

r2 = rna

library(plyr)
library(reshape2)

# ∞∞∞ TCGA - prepare simulations for empirical pvalue ----

SRPs = pl[['SPLICING_RELATED']]
SRPs%in%unique(r2$symbol)
sum(SRPs%in%unique(r2$symbol))
SRPs[!SRPs%in%unique(r2$symbol)]

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
r2$symbol[r2$symbol=='HNRPLL'] = 'HNRNPLL'
r2$symbol[r2$symbol=='NHP2L1'] = 'SNU13'

SRPs%in%unique(r2$symbol)
sum(SRPs%in%unique(r2$symbol))
SRPs[!SRPs%in%unique(r2$symbol)] 


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=subset(r2,symbol=='FOXA1')
qt=quantile(tmp$TPM)
tmp$FOXA1_overexpression=F
tmp$FOXA1_overexpression[which(tmp$TPM>=qt['75%'])]=T
r2$FOXA1_overexpression=tmp$FOXA1_overexpression[match(r2$Tumor_Sample_Barcode,tmp$Tumor_Sample_Barcode)]
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


# !!! DO NOT RUN EVERY TIME !!!! Per le simulazioni --------------

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

save(samples_boot, r2, SRPs, file="Rdata/RBP_simulations.montecarlo_TCGA__FOXA1_HE.Rdata") 


# ∞∞∞ TCGA - simulations for empirical pvalue ----

library(snow, verbose=F)
library(plyr)

load('Rdata/RBP_simulations.montecarlo_TCGA__FOXA1_HE.Rdata')

simulazione=function(samples, ri, pl){
  ri$type="WT"
  ri$type[which(ri$sample%in%samples)]="HE"
  x = subset(ri, symbol %in% pl)
  t1 = ddply(x, .(symbol),
             summarise, 
             ks=ks.test(TPM[type=="HE"], TPM[type=="WT"], alternative = 'two.sided' )$p.value 
  )
  t1$BF.ks = p.adjust(t1$ks, 'bonferroni', n=sum(!is.na(t1$ks)))
  
  return(t1[,c('symbol','ks','BF.ks')])
}


nclust = 4
print(nclust)

clus <- makeCluster(nclust)
clusterEvalQ(clus, library(plyr) ) 

print("EXPORTING")
clusterExport(cl=clus, c("r2", "SRPs", "simulazione", "samples_boot"))
print("RUNNING")
sim = parLapply(clus, samples_boot, simulazione,   ri=r2, pl=SRPs)

stopCluster(clus)

print("SAVING RESULTS")

save(sim,file='Rdata/20211130_RBP_sim_1000.montecarlo.ks.test_TCGA_TWO_SIDED__FOXA1_HE.Rdata')

print("DONE")


# ∞∞∞ TCGA - compute differentially expressed SRPs ----

load('Rdata/RBP_simulations.montecarlo_TCGA__FOXA1_HE.Rdata')

SRPs%in%unique(r2$symbol)
sum(SRPs%in%unique(r2$symbol))
SRPs[!SRPs%in%unique(r2$symbol)] 

r2$FOXA1_overexpression = r2$type

x = subset(r2, symbol %in% SRPs)

summarySE<-function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                    conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=T,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  colnames(datac)[colnames(datac)=='mean']=measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
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

load('Rdata/20211130_RBP_sim_1000.montecarlo.ks.test_TCGA_TWO_SIDED__FOXA1_HE.Rdata')
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


rownames(deseq)=deseq$symbol

colnames(deseq) = paste0("deseq.", colnames(deseq))

tc=cbind(tc, deseq[tc$symbol,])


tc$DEG.sign = with(tc, ( (
  abs(edger.logFC) >= 0.2 &
    edger.FDR<=0.01) |
    (
      abs(deseq.log2FoldChange)>=0.2 &
        deseq.padj<=0.01))
  & (BF.ks<=0.05) &
    emp<=0.05
)
sum(tc$DEG.sign) # 54

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


saveRDS(tc,'Rdata/final_dataset_148_SRPs_TCGA__FOXA1_HE.rds')


two_tailed=subset(tc,DEG.sign)
two_tailed_FDR=subset(tc,DEG.sign_FDR)


# ∞∞∞ SU2C - load datasets and libraries ----

pl = readRDS('Rdata/KEGG_Genetic_information_process.GSECA.200105.rds')

r2 = readRDS(file="Rdata/prad_su2c_2015_200716.rds")

r2$sample=as.character(r2$sample)
r2$symbol=as.character(r2$symbol)


# ∞∞∞ SU2C - prepare simulations for empirical pvalue ----

SRPs = pl[['SPLICING_RELATED']]
SRPs%in%unique(r2$symbol)
sum(SRPs%in%unique(r2$symbol))
SRPs[!SRPs%in%unique(r2$symbol)]

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
r2$symbol[r2$symbol=='HNRPLL'] = 'HNRNPLL'
r2$symbol[r2$symbol=='NHP2L1'] = 'SNU13'


SRPs%in%unique(r2$symbol)
sum(SRPs%in%unique(r2$symbol))
SRPs[!SRPs%in%unique(r2$symbol)] # 148



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp=subset(r2,symbol=='FOXA1')
qt=quantile(tmp$TPM)
tmp$FOXA1_overexpression=F
tmp$FOXA1_overexpression[which(tmp$TPM>=qt['75%'])]=T
r2$FOXA1_overexpression=tmp$FOXA1_overexpression[match(r2$sample,tmp$sample)]
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# !!! DO NOT EVERY TIME RUN !!!! Per le simulazioni --------------

library(snow, verbose=F)
library(plyr)

r2 = r2[, c("sample", "symbol", "TPM", "FOXA1_overexpression")]


r2$type = r2$FOXA1_overexpression
r2$FOXA1_overexpression = NULL

pat  = unique(r2[, c("sample", "type")])

n_foxa1_he = sum(pat$type) # = 30

samples_boot = list()
set.seed(30580) # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i in 1:10000) samples_boot[[i]]=sample(pat$sample, n_foxa1_he) 

save(samples_boot, r2, SRPs, file="Rdata/RBP_simulations.montecarlo_SU2C__FOXA1_HE.Rdata") # %%%%%%%%%%%%%%%%%%%%%%


# ∞∞∞ SU2C - simulations for empirical pvalue ----

library(snow, verbose=F)
library(plyr)

load('Rdata/RBP_simulations.montecarlo_SU2C__FOXA1_HE.Rdata')

simulazione=function(samples, ri, pl){
  ri$type="WT"
  ri$type[which(ri$sample%in%samples)]="HE"
  x = subset(ri, symbol %in% pl)
  t1 = ddply(x, .(symbol),
             summarise, 
             ks=ks.test(TPM[type=="HE"], TPM[type=="WT"], alternative = 'two.sided' )$p.value
  )
  t1$BF.ks = p.adjust(t1$ks, 'bonferroni', n=sum(!is.na(t1$ks)))
  
  return(t1[,c('symbol','ks','BF.ks')])
}


nclust = 4
print(nclust)

clus <- makeCluster(nclust)
clusterEvalQ(clus, library(plyr) ) 

print("EXPORTING")
clusterExport(cl=clus, c("r2", "SRPs", "simulazione", "samples_boot"))
print("RUNNING")
sim = parLapply(clus, samples_boot, simulazione,   ri=r2, pl=SRPs)

stopCluster(clus)

print("SAVING RESULTS")

save(sim,file='Rdata/20211201_RBP_sim_1000.montecarlo.ks.test_SU2C_TWO_SIDED__FOXA1_HE.Rdata')

print("DONE")



# ∞∞∞ SU2C - compute differentially expressed SRPs ----

load('Rdata/RBP_simulations.montecarlo_SU2C__FOXA1_HE.Rdata')

SRPs%in%unique(r2$symbol)
sum(SRPs%in%unique(r2$symbol))
SRPs[!SRPs%in%unique(r2$symbol)] 

r2$FOXA1_overexpression = r2$type

x = subset(r2, symbol %in% SRPs)

summarySE<-function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                    conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=T,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column
  colnames(datac)[colnames(datac)=='mean']=measurevar
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}

t = summarySE(x, measurevar = "TPM", groupvars = c('symbol','FOXA1_overexpression'), na.rm = T)

t1 = ddply(x, .(symbol),
           summarise,
           median.TPM    = median(TPM[FOXA1_overexpression], na.rm=T),
           median.TPM.wt = median(TPM[!FOXA1_overexpression], na.rm=T),
           ks=ks.test(    TPM[FOXA1_overexpression], TPM[!FOXA1_overexpression], alternative = 'two.sided' )$p.value
)

t1$BF.ks  = p.adjust(t1$ks, 'bonferroni', n=sum(!is.na(t1$ks)))
t1$FDR.ks = p.adjust(t1$ks, 'fdr', n=sum(!is.na(t1$ks)))

# simulazioni
load('Rdata/20211201_RBP_sim_1000.montecarlo.ks.test_SU2C_TWO_SIDED__FOXA1_HE.Rdata')
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


tc$DEG.sign = with(tc, ( (
  abs(L2R) >= 0.2 &   
    BF.ks<=0.05 &
    emp<=0.05) ) )

tc$DEG.sign_FDR = with(tc, ( (
  abs(L2R) >= 0.2 &   
    FDR.ks<=0.01 &
    emp<=0.01) ) )

tc$direction = "DOWN"
tc$direction[which(tc$delta>0)] = "UP"


saveRDS(tc,'Rdata/final_dataset_148_SRPs_SU2C__FOXA1_HE.rds')

two_tailed_FDR_met=subset(tc,DEG.sign_FDR)



# ∞∞∞∞ Save lists of SRPs ----

two_tailed_FDR_list = two_tailed_FDR$symbol
two_tailed_FDR_met_list = two_tailed_FDR_met$symbol
overlap_TCGA_met = subset(two_tailed_FDR,two_tailed_FDR$symbol%in%two_tailed_FDR_met$symbol)$symbol
save(two_tailed_FDR,two_tailed_FDR_met,two_tailed_FDR_list,two_tailed_FDR_met_list,overlap_TCGA_met,file='Rdata/NEW_differentially_expressed_lists.Rdata')

DEG_SRPs_TCGA = two_tailed_FDR
DEG_SRPs_SU2C = two_tailed_FDR_met
DEG_SRPs_18_common = subset(two_tailed_FDR,symbol%in%two_tailed_FDR_met$symbol)
save(DEG_SRPs_TCGA,DEG_SRPs_SU2C,DEG_SRPs_18_common,file='Rdata/NEW_differentially_expressed_SRPs.Rdata')





# ******************** ----
# ******************** DEA in VCaP and PC3 cells RNAseq data ----
# ******************** ----


# ******* VCaP ----


c=readRDS('Rdata/Counts_wt_si_FOXA1_VCAP.rds')

# Init ####
source("sources/utils_RNA.R")
library(tidyverse)
PCA <- function(expr, expr_cutoff=0.1, scale=T, center=T) {
  rna = t(expr)
  # remove lowly expressed
  rna = rna[,colSums(rna)>expr_cutoff] 
  # remove constant values
  rna_sd = apply(rna, 2, sd) 
  rna = rna[,rna_sd>0]
  # remove NA values
  rna_na = apply(rna, 2, function(x) sum(is.na(x)))  
  rna = rna[,rna_na==0]
  RNA_pca = prcomp(rna, scale=scale, center=center) 
  pca_summary = summary(RNA_pca)$importance
  pca_plot = as.data.frame(RNA_pca$x[,c("PC1", "PC2", "PC3")])
  list(pca_plot,pca_summary,dim(rna))
}
plot_PCA = function(x){
  y = as.data.frame(p[[1]]) 
  y$sample = rownames(y)
  y$group = "NSI"
  y$group[grep("siF",y$sample)]="siFOXA1"
  ggplot(y, aes(x=PC1, y=PC2, col=group, shape=group)) + geom_point(size=4) + 
    ggrepel::geom_text_repel(aes(label=sample),point.padding = 1) +
    xlab(paste0("PC1 (",round(p[[2]][2,1]*100,1),"%)")) + ylab(paste0("PC2 (",round(p[[2]][2,2]*100,1),"%)")) +
    ggpubr::theme_pubr() + ggsci::scale_color_d3() + 
    ggtitle(paste0("PCA (",p[[3]][2], ") TPMs"))+
    theme(legend.position = 'bottom' , legend.box = "vertical", panel.grid.minor = element_blank())#+
  
}
get_DESEQ = function(counts, samples, cond, ctrl, tlen, cores=10){
  
  # Perform DESeq analysis
  require(DESeq2)
  require(BiocParallel)
  
  tmp= floor(as.matrix(counts[,samples]))
  colData <- data.frame(condition=factor(cond))
  colData$condition <- relevel(colData$condition, ref = ctrl) 
  dds <- DESeqDataSetFromMatrix(tmp, colData, formula(~condition))
  mcols(dds)$basepairs=tlen
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds, parallel = T, BPPARAM = MulticoreParam(workers = cores))
  dds
}

# TPM #######
tpms = countToTpm2(c[,5:10], c$Length)
colnames(tpms) = gsub('_post_NGmerge','',colnames(tpms))
t = cbind.data.frame(c[,1:4],tpms)
saveRDS(t, "Rdata/TPMs_wt_si_FOXA1_VCAP.rds")

t2=readRDS("Rdata/TPMs_wt_si_FOXA1_VCAP.rds")


# deseq =========


rownames(c) = paste0(c$gene_name,'_', rownames(c))
deg = get_DESEQ(c
                , samples=c('nsi.R1','nsi.R2','nsi.R4', 'siFOXA1.R1','siFOXA1.R2', 'siFOXA1.R4')
                , cond = c(rep('NSI',3),rep('siFOXA1',3))
                , ctrl='siFOXA1'
                , c$Length)

res <- results(deg, name="condition_NSI_vs_siFOXA1")
res <- res[order(res$pvalue),]
saveRDS(res, file='Rdata/RNAseq_VCAP_deseq.rds')
res=readRDS('Rdata/RNAseq_VCAP_deseq.rds')

# PCA DESEQ ===========
vsd <- vst(deg, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")
library(ggplot2)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)#paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pm = pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=colors)

pdf(file='dist_VCAP.pdf', height = unit(3,'cm'), width = unit(3.5,'cm'), useDingbats = F )
print(pm)
dev.off()

pcaData <- DESeq2::plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file='deseq2_pca_VCAP.pdf', height = unit(3,'cm'), width = unit(3.5,'cm'), useDingbats = F )
ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) + ggrepel::geom_label_repel(aes(label=name))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+coord_fixed()+ggpubr::theme_pubr() 
dev.off()

res=readRDS('Rdata/RNAseq_VCAP_deseq.rds')
a=load('Rdata/NEW_differentially_expressed_SRPs.Rdata')

res$gene_name = sapply(strsplit(rownames(res),'\\_'),'[[',1)

saveRDS(res, file='Rdata/RNAseq_VCAP_deseq.rds')

as.data.frame(subset(res, gene_name %in% DEG_SRPs_18_common$symbol &  padj<0.1))

DEG_SRPs_18_common = cbind.data.frame(DEG_SRPs_18_common, res[match(DEG_SRPs_18_common$symbol, res$gene_name),1:6])

res$gene_name[which(res$gene_name=='HNRNPLL')]='HNRPLL'
DEG_SRPs_TCGA = cbind.data.frame(DEG_SRPs_TCGA, res[match(DEG_SRPs_TCGA$symbol, res$gene_name),1:6])
DEG_SRPs_SU2C = cbind.data.frame(DEG_SRPs_SU2C, res[match(DEG_SRPs_SU2C$symbol, res$gene_name),1:6])

save(DEG_SRPs_18_common, DEG_SRPs_SU2C, DEG_SRPs_TCGA, file='Rdata/NEW_differentially_expressed_SRPs_deseq2VCAP.Rdata')



# ******* PC3 ----



c=readRDS('Rdata/Counts_wt_si_FOXA1_PC3_nextseq_novaseq.rds')

count_matrix<-as.matrix(c[, 5:ncol(c)])

saveRDS(count_matrix, "Rdata/Counts_wt_si_FOXA1_PC3_nextseq_novaseq_matrix.rds")

# Init ####
source("Rdata/utils_RNA.R")
library(tidyverse)
PCA <- function(expr, expr_cutoff=0.1, scale=T, center=T) {
  rna = t(expr)
  # remove lowly expressed
  rna = rna[,colSums(rna)>expr_cutoff] 
  # remove constant values
  rna_sd = apply(rna, 2, sd) 
  rna = rna[,rna_sd>0]
  # remove NA values
  rna_na = apply(rna, 2, function(x) sum(is.na(x)))  
  rna = rna[,rna_na==0]
  RNA_pca = prcomp(rna, scale=scale, center=center) 
  pca_summary = summary(RNA_pca)$importance
  pca_plot = as.data.frame(RNA_pca$x[,c("PC1", "PC2", "PC3")])
  list(pca_plot,pca_summary,dim(rna))
}
plot_PCA = function(x){
  y = as.data.frame(p[[1]]) 
  y$sample = rownames(y)
  y$group = "NSI"
  y$group[grep("siF",y$sample)]="siFOXA1"
  ggplot(y, aes(x=PC1, y=PC2, col=group, shape=group)) + geom_point(size=4) + 
    ggrepel::geom_text_repel(aes(label=sample),point.padding = 1) +
    xlab(paste0("PC1 (",round(p[[2]][2,1]*100,1),"%)")) + ylab(paste0("PC2 (",round(p[[2]][2,2]*100,1),"%)")) +
    ggpubr::theme_pubr() + ggsci::scale_color_d3() + 
    ggtitle(paste0("PCA (",p[[3]][2], ") TPMs"))+
    theme(legend.position = 'bottom' , legend.box = "vertical", panel.grid.minor = element_blank())#+
  
}
get_DESEQ = function(counts, samples, cond, ctrl, tlen, cores=1){
  
  # Perform DESeq analysis
  require(DESeq2)
  require(BiocParallel)
  
  tmp= floor(as.matrix(counts[,samples]))
  colData <- data.frame(condition=factor(cond))
  colData$condition <- relevel(colData$condition, ref = ctrl) 
  dds <- DESeqDataSetFromMatrix(tmp, colData, formula(~condition))
  mcols(dds)$basepairs=tlen
  dds <- dds[ rowSums(counts(dds)) > 1, ]
  dds <- DESeq(dds, parallel = T, BPPARAM = MulticoreParam(workers = cores))
  dds
}

# TPM #######

tpms = countToTpm2(c[,5:10], c$Length)
colnames(tpms) = gsub('pc3.','',colnames(tpms))
t = cbind.data.frame(c[,1:4],tpms)

saveRDS(t, "Rdata/TPM_wt_si_FOXA1_PC3_nextseq_novaseq.rds")


t$mean_nsi=rowMeans(t[, c("nsi.r6", "nsi.r7", "nsi.r8")])
t$mean_si=rowMeans(t[, c("siFOXA1.r6", "siFOXA1.r7", "siFOXA1.r8")])


t$delta_mean=t$mean_nsi-t$mean_si

t$fc=t$mean_nsi/t$mean_si

saveRDS(t, "Rdata/TPM_wt_si_FOXA1_PC3_nextseq_novaseq.rds")


# deseq =========


c=readRDS("Rdata/Counts_wt_si_FOXA1_PC3_nextseq_novaseq.rds")

colnames(c) = gsub('pc3.','',colnames(c))
rownames(c) = paste0(c$gene_name,'_', rownames(c))


deg = get_DESEQ(c
                , samples=c("nsi.r6","nsi.r7","nsi.r8", "siFOXA1.r6","siFOXA1.r7","siFOXA1.r8")
                , cond = c(rep('NSI',3),rep('siFOXA1',3))
                , ctrl='siFOXA1'
                , c$Length)

res <- results(deg, name="condition_NSI_vs_siFOXA1")
res <- res[order(res$pvalue),]

res$gene_name=c$gene_name[match(rownames(res), rownames(c))]
res$gene_id=c$Geneid[match(rownames(res), rownames(c))]

saveRDS(res, "Rdata/RNAseq_PC3_nextseq_novaseq_deseq.rds")


# PCA DESEQ ===========
vsd <- vst(deg, blind=FALSE)
sampleDists <- dist(t(assay(vsd)))
library("RColorBrewer")
library("pheatmap")
library(ggplot2)
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)#paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pm = pheatmap(sampleDistMatrix,
              clustering_distance_rows=sampleDists,
              clustering_distance_cols=sampleDists,
              col=colors)

pdf(file='dist_PC3_next_nova.pdf', height = unit(3,'cm'), width = unit(3.5,'cm'), useDingbats = F )

print(pm)
dev.off()


pcaData <- DESeq2::plotPCA(vsd, intgroup=c("condition"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(file='deseq2_pca_PC3_next_nova.pdf', height = unit(3,'cm'), width = unit(3.5,'cm'), useDingbats = F )


ggplot(pcaData, aes(PC1, PC2, color=condition)) +
  geom_point(size=3) + ggrepel::geom_label_repel(aes(label=name))+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance"))+coord_fixed()+ggpubr::theme_pubr() 
dev.off()



res=readRDS('Rdata/RNAseq_PC3_nextseq_novaseq_deseq.rds')
a=load('Rdata/NEW_differentially_expressed_SRPs.Rdata')

res$gene_name = sapply(strsplit(rownames(res),'\\_'),'[[',1)

saveRDS(res, file='Rdata/RNAseq_PC3_nextseq_novaseq_deseq.rds')

as.data.frame(subset(res, gene_name %in% DEG_SRPs_18_common$symbol &  padj<0.1))

DEG_SRPs_18_common = cbind.data.frame(DEG_SRPs_18_common, res[match(DEG_SRPs_18_common$symbol, res$gene_name),1:6])

res$gene_name[which(res$gene_name=='HNRNPLL')]='HNRPLL'
DEG_SRPs_TCGA = cbind.data.frame(DEG_SRPs_TCGA, res[match(DEG_SRPs_TCGA$symbol, res$gene_name),1:6])
DEG_SRPs_SU2C = cbind.data.frame(DEG_SRPs_SU2C, res[match(DEG_SRPs_SU2C$symbol, res$gene_name),1:6])

save(DEG_SRPs_18_common, DEG_SRPs_SU2C, DEG_SRPs_TCGA, file='Rdata/NEW_differentially_expressed_SRPs_deseq2_PC3_nextseq_novaseq.Rdata')



