# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Survival analysis ----



# ******************** ----
# ******************** Survival analysis on "cumulative" inclusion of silenced and enhanced poison and essential exons ----
# ******************** ----

library(survminer)
library(survival)
library(plyr)

df = readRDS("Rdata/survival_data_20190304.rds") 
df$first_event_months = df$first_event/(365/12)
dexSIGN = readRDS("Rdata/dexSIGN_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds")
ex = readRDS("Rdata/PRAD_selected_exons.rds")
ex$as_id =  as.character(ex$as_id)
ex$splice_type =  as.character(ex$splice_type)


# Silenced poison exons ----

ex_poi = subset(ex,as_id%in%subset(dexSIGN,POISON & delta.mean<0)$event_id)

ex_poi = subset(ex_poi,patient%in%df$bcr_patient_barcode)

tmp = ddply(ex_poi,.(as_id),mutate,qt25=quantile(value)['25%'],qt75=quantile(value)['75%'])

tmp = ddply(tmp,.(patient),summarise,n25=sum(value<=qt25),n75=sum(value>=qt75))

tmp$ratio_n25_n75 = tmp$n25/tmp$n75

qt=quantile(tmp$ratio_n25_n75)

tmp$pois_ratio_lev=cut(tmp$ratio_n25_n75,breaks = c(0,qt['25%'],qt['50%'],qt['75%'],qt['100%']),labels = c('(0,25]','none','none','(75,100]'))
tmp = subset(tmp,!pois_ratio_lev=='none')
tmp$pois_ratio_lev = as.factor(as.character(tmp$pois_ratio_lev))
df$pois_ratio_lev = tmp$pois_ratio_lev[match(df$bcr_patient_barcode,tmp$patient)]
formula <- as.formula("Surv(time=first_event_months, event=DFS) ~ pois_ratio_lev") 
cox <- coxph(formula, df)
summary(cox)
cox.zph(cox)
pdf(file='Silenced_poison_svES_extreme_quartiles_forest_plot.pdf', width = 8, height = 6, useDingbats = F, onefile=F)
ggforest(cox)
dev.off()
gfit <- survfit(formula = as.formula("Surv(time = first_event_months, event = DFS) ~ pois_ratio_lev"), data = subset(df, first_event_months<110))
pdf(file='Silenced_poison_svES_extreme_quartiles_kaplan_meier_Red_high_incl_Blue_low_incl.pdf', width = 9, height = 9, useDingbats = F, onefile=F)
ggsurvplot(gfit, data = subset(df, first_event_months<110), risk.table = TRUE) 
dev.off()


# Enhanced poison exons ----

ex_poi = subset(ex,as_id%in%subset(dexSIGN,POISON & delta.mean>0)$event_id)

ex_poi = subset(ex_poi,patient%in%df$bcr_patient_barcode)

tmp = ddply(ex_poi,.(as_id),mutate,qt25=quantile(value)['25%'],qt75=quantile(value)['75%'])

tmp = ddply(tmp,.(patient),summarise,n25=sum(value<=qt25),n75=sum(value>=qt75))

tmp$ratio_n25_n75 = tmp$n25/tmp$n75

qt=quantile(tmp$ratio_n25_n75)

tmp$pois_ratio_lev=cut(tmp$ratio_n25_n75,breaks = c(0,qt['25%'],qt['50%'],qt['75%'],qt['100%']),labels = c('(0,25]','none','none','(75,100]'))
tmp = subset(tmp,!pois_ratio_lev=='none')
tmp$pois_ratio_lev = as.factor(as.character(tmp$pois_ratio_lev))
df$pois_ratio_lev = tmp$pois_ratio_lev[match(df$bcr_patient_barcode,tmp$patient)]
formula <- as.formula("Surv(time=first_event_months, event=DFS) ~ pois_ratio_lev") 
cox <- coxph(formula, df)
summary(cox)
cox.zph(cox)
pdf(file='Enhanced_poison_svES_extreme_quartiles_forest_plot.pdf', width = 8, height = 6, useDingbats = F, onefile=F)
ggforest(cox)
dev.off()
gfit <- survfit(formula = as.formula("Surv(time = first_event_months, event = DFS) ~ pois_ratio_lev"), data = subset(df, first_event_months<110))
pdf(file='Enhanced_poison_svES_extreme_quartiles_kaplan_meier_Red_high_incl_Blue_low_incl.pdf', width = 9, height = 9, useDingbats = F, onefile=F)
ggsurvplot(gfit, data = subset(df, first_event_months<110), risk.table = TRUE) 
dev.off()


# Silenced essential exons ----

ex_poi = subset(ex,as_id%in%subset(dexSIGN,ESSENTIAL & delta.mean<0)$event_id)

ex_poi = subset(ex_poi,patient%in%df$bcr_patient_barcode)

tmp = ddply(ex_poi,.(as_id),mutate,qt25=quantile(value)['25%'],qt75=quantile(value)['75%'])

tmp = ddply(tmp,.(patient),summarise,n25=sum(value<=qt25),n75=sum(value>=qt75))

tmp$ratio_n25_n75 = tmp$n25/tmp$n75

qt=quantile(tmp$ratio_n25_n75)

tmp$ess_ratio_lev=cut(tmp$ratio_n25_n75,breaks = c(0,qt['25%'],qt['50%'],qt['75%'],qt['100%']),labels = c('(0,25]','none','none','(75,100]'))
tmp = subset(tmp,!ess_ratio_lev=='none')
tmp$ess_ratio_lev = as.factor(as.character(tmp$ess_ratio_lev))
df$ess_ratio_lev = tmp$ess_ratio_lev[match(df$bcr_patient_barcode,tmp$patient)]
formula <- as.formula("Surv(time=first_event_months, event=DFS) ~ ess_ratio_lev")
cox <- coxph(formula, df)
summary(cox)
cox.zph(cox)
pdf(file='Silenced_essential_svES_extreme_quartiles_forest_plot.pdf', width = 8, height = 6, useDingbats = F, onefile=F)
ggforest(cox)
dev.off()
gfit <- survfit(formula = as.formula("Surv(time = first_event_months, event = DFS) ~ ess_ratio_lev"), data = subset(df, first_event_months<110))
pdf(file='Silenced_essential_svES_extreme_quartiles_kaplan_meier_Red_high_incl_Blue_low_incl.pdf', width = 9, height = 9, useDingbats = F, onefile=F)
ggsurvplot(gfit, data = subset(df, first_event_months<110), risk.table = TRUE)
dev.off()


# Enhanced essential exons ----

ex_poi = subset(ex,as_id%in%subset(dexSIGN,ESSENTIAL & delta.mean>0)$event_id)

ex_poi = subset(ex_poi,patient%in%df$bcr_patient_barcode)

tmp = ddply(ex_poi,.(as_id),mutate,qt25=quantile(value)['25%'],qt75=quantile(value)['75%'])

tmp = ddply(tmp,.(patient),summarise,n25=sum(value<=qt25),n75=sum(value>=qt75))

tmp$ratio_n25_n75 = tmp$n25/tmp$n75

qt=quantile(tmp$ratio_n25_n75)

tmp$ess_ratio_lev=cut(tmp$ratio_n25_n75,breaks = c(0,qt['25%'],qt['50%'],qt['75%'],qt['100%']),labels = c('(0,25]','none','none','(75,100]'))
tmp = subset(tmp,!ess_ratio_lev=='none')
tmp$ess_ratio_lev = as.factor(as.character(tmp$ess_ratio_lev))
df$ess_ratio_lev = tmp$ess_ratio_lev[match(df$bcr_patient_barcode,tmp$patient)]
formula <- as.formula("Surv(time=first_event_months, event=DFS) ~ ess_ratio_lev")
cox <- coxph(formula, df)
summary(cox)
cox.zph(cox)
pdf(file='Enhanced_essential_svES_extreme_quartiles_forest_plot.pdf', width = 8, height = 6, useDingbats = F, onefile=F)
ggforest(cox)
dev.off()
gfit <- survfit(formula = as.formula("Surv(time = first_event_months, event = DFS) ~ ess_ratio_lev"), data = subset(df, first_event_months<110))
pdf(file='Enhanced_essential_svES_extreme_quartiles_kaplan_meier_Red_high_incl_Blue_low_incl.pdf', width = 9, height = 9, useDingbats = F, onefile=F)
ggsurvplot(gfit, data = subset(df, first_event_months<110), risk.table = TRUE)
dev.off()








# ******************** ----
# ******************** Survival analysis for each poison or essential exon ----
# ******************** ----

library(survminer)
library(survival)

dexSIGN = readRDS("Rdata/dexSIGN_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds")

ex = readRDS("Rdata/PRAD_selected_exons.rds")
ex$as_id =  as.character(ex$as_id)
ex$splice_type =  as.character(ex$splice_type)
df = readRDS("Rdata/survival_data_20190304.rds")
df$first_event_months = df$first_event/(365/12) 

dexSIGN_pois_ess = subset(dexSIGN, (POISON & !ESSENTIAL) |
                            (ESSENTIAL & !POISON) )

surv_summary_pois_ess = matrix(nrow = 5, ncol = nrow(dexSIGN_pois_ess))
dff = df
for (i in 1:nrow(dexSIGN_pois_ess)) {
  event = dexSIGN_pois_ess$event_id[i]
  incl = subset(ex,as_id==event)$value
  names(incl) = substr(subset(ex,as_id==event)$variable,1,12)
  
  qt=quantile(incl)
  incl2 = rep(NA, length(incl))
  names(incl2)=names(incl)
  incl2[names(incl[incl>=qt['75%']])] = 'HI'
  incl2[names(incl[incl<=qt['25%']])] = 'LI'
  
  tmp = rep(NA, nrow(df))
  names(tmp) = df$bcr_patient_barcode
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
  
  dff$tmp = NULL
  
  print(i)
}
colnames(surv_summary_pois_ess) = dexSIGN_pois_ess$event_id
rownames(surv_summary_pois_ess) = c('logrank_pvalue','hazard_ratio','prop_ass_pvalue','lower_CI','upper_CI')





# ******************** ----
# ******************** GLM ----
# ******************** ----


library(survminer)
library(survival)
library(My.stepwise)
library(reshape2)
library(gridExtra)
library(ggrepel)

dexSIGN = readRDS("Rdata/dexSIGN_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds")
ex = readRDS("Rdata/PRAD_selected_exons.rds")
ex$as_id =  as.character(ex$as_id)
ex$splice_type =  as.character(ex$splice_type)
df = readRDS("Rdata/survival_data_20190304.rds")
df$first_event_months = df$first_event/(365/12) 

dexSIGN_pois_sil=subset(dexSIGN,(POISON & !ESSENTIAL) & delta.mean<0)
dexSIGN_pois_enh=subset(dexSIGN,(POISON & !ESSENTIAL) & delta.mean>0)
dexSIGN_ess_sil=subset(dexSIGN,(ESSENTIAL & !POISON) & delta.mean<0)
dexSIGN_ess_enh=subset(dexSIGN,(ESSENTIAL & !POISON) & delta.mean>0)

dexSIGN_pois_ess = rbind(dexSIGN_pois_sil,dexSIGN_pois_enh,dexSIGN_ess_sil,dexSIGN_ess_enh)

dexSIGN_pois_ess$logrank_FDR = p.adjust(dexSIGN_pois_ess$logrank_pval,method = 'fdr')


sel = as.character(subset(dexSIGN_pois_ess,logrank_FDR<0.05 & hazard_ratio>1)$event_id)

formula <- as.formula(paste("FOXA1~",
                            paste(sel, collapse = "+")))

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


ex$variable = as.character(ex$variable)
DFF = subset(df,df$FOXA1_levels%in%c('(0-25]','(75-100]'))
DFF$splic_event = NULL
for (i in 1:length(sel)) {
  event = sel[i]
  incl = subset(ex,as_id==event)$value
  names(incl) = substr(subset(ex,as_id==event)$variable,1,12)
  tmp = rep(NA, nrow(DFF))
  names(tmp) = DFF$bcr_patient_barcode
  tmp = incl[names(tmp)] # usare "scale()" non cambia i risultati %%%%%%%%%%%
  DFF = cbind(DFF,tmp)
  print(i)
}
colnames(DFF)[(ncol(DFF)-length(sel)+1):ncol(DFF)] = sel  

fit = glm(data=DFF, formula = formula ); coe = get_coef(fit) 

coe$symbol = dexSIGN_pois_ess$symbol[match(coe$g,dexSIGN_pois_ess$event_id)]
coe$code = paste0(sapply(strsplit(as.character(coe$g),'_'),'[',3),'_',coe$symbol)
coe$code = factor(coe$code,levels = coe$code)
pdf(file='FOXA1_surv_sign_FDR_pois_ess_glm_coeff__FOXA1_25_75__Harmful_events_only.pdf', width = 5, height = 6, useDingbats = F, onefile=F)
ggplot(coe, aes(x=g, y=coe, fill = coe, size=pval<0.05) )+
  geom_bar(stat='identity',color='black')+theme_bw()+coord_flip()+scale_fill_gradient2(midpoint=0, low="blue", mid="white", high="red")+
  scale_size_manual(values = c(0,1))
dev.off()

set.seed(30580)
boot <- boot.relimp(fit
                    , b = 1000, type = c("lmg"), rank = TRUE, diff = TRUE, rela = TRUE)
b = booteval.relimp(boot,sort=TRUE) # print result

pdf(file='RelImpo_glm_FOXA1_surv_sign_FDR_pois_ess__FOXA1_25_75__Harmful_events_only.pdf', width = 5, height = 5, useDingbats = F, onefile=F)
par(las=2)
plot(b)
dev.off()








# ******************** ----
# ******************** Univariable survival analysis on harmful exons ----
# ******************** ----




ex <- readRDS("Rdata/PRAD_selected_exons.rds")
ex$cases=NULL

library(survival)
library(survminer)

outlierMatrix_survSIGN = matrix(nrow=length(sel), ncol=length(unique(ex$variable)))
rownames(outlierMatrix_survSIGN) = sel
colnames(outlierMatrix_survSIGN) = unique(ex$variable)

# Create "outierMatrix" ---
n = 1
for (i in sel) {
  tmp = subset(ex, as_id==i)
  qt = quantile(tmp$value, seq(0,1,0.25))
  up = qt['75%']
  low = qt['25%']
  for (j in 1:ncol(outlierMatrix_survSIGN)) {
    if (tmp$value[j]<=low) outlierMatrix_survSIGN[n,j] = 'L'
    if (tmp$value[j]>=up) outlierMatrix_survSIGN[n,j] = 'U'
    if (!tmp$value[j]<=low & !tmp$value[j]>=up) outlierMatrix_survSIGN[n,j] = NA
  }
  n = n+1
  print(i)
}

colnames(outlierMatrix_survSIGN) = substr(colnames(outlierMatrix_survSIGN),1,12)


df = readRDS("Rdata/survival_data_20190304.rds")
df$first_event_months = df$first_event/(365/12) # problema delgi NA risolto

df2 = df

event = 'exon_skip_517627'

df2$splic_event = NA
for (j in 1:ncol(outlierMatrix_survSIGN)) {
  df2$splic_event[match(colnames(outlierMatrix_survSIGN)[j],df2$bcr_patient_barcode)] = outlierMatrix_survSIGN[event,j]
}
formula <- as.formula("Surv(time=first_event_months, event=DFS) ~ splic_event") 
cox <- coxph(formula, df2)
tmp = summary(cox)
tmp$sctest[3]
tmp$coefficients[2]
tmp = cox.zph(cox)
tmp$table[3]
ggforest(cox)
gfit <- survfit(formula = as.formula("Surv(time = first_event_months, event = DFS) ~ splic_event"), data = subset(df2, first_event_months<110))
ggsurvplot(gfit, data = subset(df2, first_event_months<110), risk.table = TRUE)


gfit_list = list()
n = 1
for (event in sel) {
  df2$splic_event = NA
  for (j in 1:ncol(outlierMatrix_survSIGN)) {
    df2$splic_event[match(colnames(outlierMatrix_survSIGN)[j],df2$bcr_patient_barcode)] = outlierMatrix_survSIGN[event,j]
  }
  gfit_list[[n]] <- survfit(formula = as.formula("Surv(time = first_event_months, event = DFS) ~ splic_event"), data = subset(df2, first_event_months<110))
  n = n+1
}
names(gfit_list) = sel

p1=ggsurvplot(gfit_list[[1]], data = subset(df2, first_event_months<110), risk.table = TRUE, title=names(gfit_list)[[1]])
p2=ggsurvplot(gfit_list[[2]], data = subset(df2, first_event_months<110), risk.table = TRUE, title=names(gfit_list)[[2]])
p3=ggsurvplot(gfit_list[[3]], data = subset(df2, first_event_months<110), risk.table = TRUE, title=names(gfit_list)[[3]])
p4=ggsurvplot(gfit_list[[4]], data = subset(df2, first_event_months<110), risk.table = TRUE, title=names(gfit_list)[[4]])
p5=ggsurvplot(gfit_list[[5]], data = subset(df2, first_event_months<110), risk.table = TRUE, title=names(gfit_list)[[5]])
p6=ggsurvplot(gfit_list[[6]], data = subset(df2, first_event_months<110), risk.table = TRUE, title=names(gfit_list)[[6]])

ggsurvlist <- list(
  p1,p2,p3,p4,p5,p6
)

# Arrange and save into pdf file
res <- arrange_ggsurvplots(ggsurvlist, print = FALSE, ncol = 2, nrow = 3)
ggsave("Figures/kaplan_meyer_all_harmful_events.pdf", res, height = 30, width = 21)






# Optimal cut point for FLNA exon 30 inclusion ----

df2=df
df2$FLNA_ex30 = subset(ex,as_id=='exon_skip_517627')$value[match(df2$bcr_patient_barcode,subset(ex,as_id=='exon_skip_517627')$patient)]

cutpoint=surv_cutpoint(df2,time='first_event_months',event='DFS',variables = c('FLNA_ex30'))

df2$FLNA_ex30_discr='L'
df2$FLNA_ex30_discr[which(df2$FLNA_ex30>=cutpoint$cutpoint$cutpoint)]='U'

formula <- as.formula("Surv(time=first_event_months, event=DFS) ~ FLNA_ex30_discr") 
cox <- coxph(formula, df2)
tmp = summary(cox)
tmp$sctest[3]
tmp$coefficients[2]
tmp = cox.zph(cox)
tmp$table[3]
pdf(file = '~/Dropbox (HuGeF)/Prostate_splicing/FOXA1/Paper/CellReports_211220/REVISION/original_figures/FLNA_survival_optimal_cutpoint_strata.pdf',width = unit(6,'cm'),height = unit(4,'cm'))
ggforest(cox)
dev.off()
gfit <- survfit(formula = as.formula("Surv(time = first_event_months, event = DFS) ~ FLNA_ex30_discr"), data = subset(df2, first_event_months<110))
pdf(file = '~/Dropbox (HuGeF)/Prostate_splicing/FOXA1/Paper/CellReports_211220/REVISION/original_figures/FLNA_survival_optimal_cutpoint_strata_Kaplan.pdf',width = unit(6,'cm'),height = unit(6,'cm'))
ggsurvplot(gfit, data = subset(df2, first_event_months<110), risk.table = TRUE)
dev.off()



rna <- readRDS("Rdata/RNA_TPM_tumor.rds")

df2$FOXA1_expr=subset(rna,symbol=='FOXA1')$TPM[match(df2$bcr_patient_barcode,subset(rna,symbol=='FOXA1')$patient)]

p1=ggplot(df2,aes(x=FOXA1_expr,y=FLNA_ex30,col=FLNA_ex30_discr,shape=FOXA1_overexpression))+geom_point()+theme_bw()+geom_vline(xintercept = min(subset(rna,symbol=='FOXA1' & FOXA1_overexpression)$TPM),linetype='dashed')+geom_hline(yintercept = cutpoint$cutpoint$cutpoint,linetype='dashed')

p2=ggplot(df2,aes(x=FLNA_ex30_discr,fill=FOXA1_overexpression))+geom_bar(position='fill')+theme_bw()+
  scale_y_continuous(labels = scales::percent)+
  scale_fill_manual(values = c('grey','darkgreen'))

to_test=table(df2$FOXA1_overexpression,df2$FLNA_ex30_discr)
to_test[2,]/colSums(to_test) # Proportion of FOXA1_overexpression patients in FLNA stratification according to survival optimal cutpoint
# L         U 
# 0.2086331 0.4155844 
fisher.test(to_test)
# Fisher's Exact Test for Count Data
# data:  to_test
# p-value = 0.0003662
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.511946 4.764315
# sample estimates:
# odds ratio 
#   2.688679


p3=ggplot(df2,aes(fill=FLNA_ex30_discr,x=FOXA1_overexpression))+geom_bar(position='fill')+theme_bw()+
  scale_y_continuous(labels = scales::percent)

to_test=table(df2$FOXA1_overexpression,df2$FLNA_ex30_discr)
to_test[,2]/rowSums(to_test) # Proportion of FLNA stratification with respect to FOXA1 stratification 
# 0         1 
# 0.1698113 0.3555556
fisher.test(to_test)
# Fisher's Exact Test for Count Data
# data:  to_test
# p-value = 0.0003662
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#  1.511946 4.764315
# sample estimates:
# odds ratio 
#   2.688679


pdf(file = '~/Dropbox (HuGeF)/Prostate_splicing/FOXA1/Paper/CellReports_211220/REVISION/original_figures/FLNA_inclusion_VS_FOXA1_expression__optimal_cutpoint.pdf',width = unit(24,'cm'),height = unit(6,'cm'))
grid.arrange(p1,p2,p3,ncol=3)
dev.off()


df2$key=paste0('FOXA1_',df2$FOXA1_overexpression,'_FLNA_',df2$FLNA_ex30_discr)

formula <- as.formula("Surv(time=first_event_months, event=DFS) ~ key") 
cox <- coxph(formula, df2)
tmp = summary(cox)
tmp$sctest[3]
tmp$coefficients[2]
tmp = cox.zph(cox)
tmp$table[3]
pdf(file = '~/Dropbox (HuGeF)/Prostate_splicing/FOXA1/Paper/CellReports_211220/REVISION/original_figures/FLNA_survival_optimal_cutpoint_strata_VS_FOXA1.pdf',width = unit(6,'cm'),height = unit(4,'cm'))
ggforest(cox)
dev.off()
gfit <- survfit(formula = as.formula("Surv(time = first_event_months, event = DFS) ~ key"), data = subset(df2, first_event_months<110))
pdf(file = '~/Dropbox (HuGeF)/Prostate_splicing/FOXA1/Paper/CellReports_211220/REVISION/original_figures/FLNA_survival_optimal_cutpoint_strata_VS_FOXA1_Kaplan.pdf',width = unit(6,'cm'),height = unit(6,'cm'))
ggsurvplot(gfit, data = subset(df2, first_event_months<110), risk.table = TRUE)
dev.off()




df2$SRSF1_expr=subset(rna,symbol=='SFRS1')$TPM[match(df2$bcr_patient_barcode,subset(rna,symbol=='SFRS1')$patient)]
my_comparisons=list(c('FOXA1_0_FLNA_L','FOXA1_0_FLNA_U'),c('FOXA1_1_FLNA_L','FOXA1_1_FLNA_U'),c('FOXA1_0_FLNA_L','FOXA1_1_FLNA_L'),c('FOXA1_0_FLNA_U','FOXA1_1_FLNA_U'),c('FOXA1_0_FLNA_L','FOXA1_1_FLNA_U'))
pdf(file = '~/Dropbox (HuGeF)/Prostate_splicing/FOXA1/Paper/CellReports_211220/REVISION/original_figures/Boxplot_SRSF1_expr_VS_FOXA1_and_FLNAex30.pdf',width = unit(8,'cm'),height = unit(6,'cm'))
ggplot(df2,aes(x=key,y=SRSF1_expr,fill=FLNA_ex30_discr))+geom_boxplot(notch = T)+geom_jitter()+theme_bw()+stat_compare_means(comparisons = my_comparisons)
dev.off()

pdf(file = '~/Dropbox (HuGeF)/Prostate_splicing/FOXA1/Paper/CellReports_211220/REVISION/original_figures/Boxplot_SRSF1_expr_VS_FOXA1_and_FLNAex30__no_Jitter.pdf',width = unit(8,'cm'),height = unit(6,'cm'))
ggplot(df2,aes(x=key,y=SRSF1_expr,fill=FLNA_ex30_discr))+geom_boxplot(notch = T)+theme_bw()+stat_compare_means(comparisons = my_comparisons)
dev.off()

