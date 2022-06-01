# ••••••••••••••••••••••••••••••••••••••••• ----
# *** GLMs to model cumulative SRGs expression as a function of the expression of FOXA1, AR, MYC and ERG ----

library(plyr)
library(reshape2)
library(maditr)
library(caret)
library(relaimpo)
library(ggpubr)
library(gridExtra)


FIGURE_DIR = # set figure directory

pl = readRDS("Rdata/KEGG_Genetic_information_process.GSECA.200105.rds")

rna=readRDS("Rdata/RNA_TPM_tumor.rds")
suc= readRDS("Rdata/prad_su2c_2015_200716.rds")

nepc = read.delim('Tables/nepc_wcm_2016/data_RNA_Seq_expression_median.txt',sep = '\t')
clinical = read.delim('Tables/nepc_wcm_2016/data_clinical_sample.txt',sep='\t',skip = 4)
nepc = cbind(nepc[,1:2],nepc[colnames(nepc)%in%subset(clinical,DISEASE_CODE=='CRPC-NE')$SAMPLE_ID])
nepc = melt(nepc,id.vars = c('Hugo_Symbol','Entrez_Gene_Id'))
nepc$Hugo_Symbol[which(nepc$Entrez_Gene_Id=='4809')] = 'SNU13'


srp.p = ddply(subset(rna, symbol%in%pl$SPLICING_RELATED), .(sample),summarise, tpm = sum(TPM))
srp.m = ddply(subset(suc, symbol%in%pl$SPLICING_RELATED), .(sample),summarise, tpm = sum(TPM))
srp.n = ddply(subset(nepc, Hugo_Symbol%in%pl$SPLICING_RELATED), .(variable),summarise, tpm = sum(value))

p = dcast(subset(rna, symbol %in%c("FOXA1",'AR',"ERG","MYC")), sample~symbol, value.var = 'TPM')
m = dcast(subset(suc, symbol %in%c("FOXA1",'AR',"ERG","MYC")), sample~symbol, value.var = 'TPM')
n = dcast(subset(nepc, Hugo_Symbol %in%c("FOXA1",'AR',"ERG","MYC")), variable~Hugo_Symbol, value.var = 'value')


p$SRP = srp.p$tpm[match(p$sample,srp.p$sample)]
rownames(p) = p$sample; p = p[,-1]
m$SRP = srp.m$tpm[match(m$sample,srp.m$sample)]
rownames(m) = m$sample; m = m[,-1]
n$SRP = srp.n$tpm[match(n$variable,srp.n$variable)]
rownames(n) = n$variable; n = n[,-1]

x   = predict(preProcess(p, method = c("center", "scale", "YeoJohnson", "nzv"))
              , newdata = p)
y   = predict(preProcess(m, method = c("center", "scale", "YeoJohnson", "nzv"))
              , newdata = m)
z   = predict(preProcess(n, method = c("center", "scale", "YeoJohnson", "nzv"))
              , newdata = n)


get_coef = function(f){
  coe = coef(f); coe  = coe[-1]; coe = sort(coe)
  pvalue = coef(summary(f))[,4]
  coe = data.frame(g = names(coe), coe = coe, pval = pvalue[names(coe)])
  coe$g = factor(coe$g, levels = coe$g)
  coe$pval[coe$pval>0.05] = 1
  print(ggplot(coe, aes(x=g, y=coe, fill = coe, size=pval<0.05) )+
          geom_bar(stat='identity',color='black')+theme_bw()+coord_flip()+scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red")+
          scale_size_manual(values = c(0,1)))
  coe
}

set.seed(30580)

pdf(file=paste0(FIGURE_DIR, "SRP_TF_primary_glm_coeff_NEW.pdf"), height=unit(3,'cm'), width=unit(3, 'cm'), useDingbats = F )
fit.p = glm(data=x, SRP~FOXA1+MYC+AR+ERG ); coe1 = get_coef(fit.p)
dev.off()
pdf(file=paste0(FIGURE_DIR, "SRP_TF_metastatic_glm_coeff.pdf"), height=unit(3,'cm'), width=unit(3, 'cm'), useDingbats = F )
fit.m = glm(data=y, SRP~FOXA1+MYC+AR+ERG ); coe2 = get_coef(fit.m)
dev.off()
pdf(file=paste0(FIGURE_DIR, "SRP_TF_neuroendocrine_glm_coeff__NO_HOXB13.pdf"), height=unit(3,'cm'), width=unit(3, 'cm'), useDingbats = F )
fit.n = glm(data=z, SRP~FOXA1+MYC+AR+ERG ); coe3 = get_coef(fit.n)
dev.off()


set.seed(30580)
boot <- boot.relimp(fit.p, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
b = booteval.relimp(boot,sort=TRUE) # print result
boot <- boot.relimp(fit.m, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
b2 = booteval.relimp(boot,sort=TRUE) # print result
boot <- boot.relimp(fit.n, b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
b3 = booteval.relimp(boot,sort=TRUE) # print result

pdf(file=paste0(FIGURE_DIR, "RelImpo_glm_primary.pdf"), height=unit(5,'cm'), width=unit(5,'cm'), useDingbats = F)
plot(b)
dev.off()
pdf(file=paste0(FIGURE_DIR, "RelImpo_glm_metastatic.pdf"), height=unit(5,'cm'), width=unit(5,'cm'), useDingbats = F)
plot(b2)
dev.off()
pdf(file=paste0(FIGURE_DIR, "RelImpo_glm_neuroendocrine.pdf"), height=unit(5,'cm'), width=unit(5,'cm'), useDingbats = F)
plot(b3)
dev.off()



