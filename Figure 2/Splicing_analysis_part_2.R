#!/usr/bin/env Rscript --vanilla

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%  Analysis of splicing variability in Prostate Cancer TCGA patients %%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PART 2 # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message('[****] SCRIPT: Splicing_analysis_part_2.R')

print('[****] SCRIPT: Splicing_analysis_part_2.R')


message('[*] Setting working path and creating directories ...')

gene = 'FOXA1'

print(paste0('STRATA: ',gene))

work_dir     = 'Automatic_splicing_analysis_results'
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
library(ggplotify)
library(ggpubr)

# XXX. Collect bootstrapping results ----

message('[*] Collecting bootstrapping results ...')

library(reshape2)

dex_CASE = readRDS(paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))
obs=dex_CASE[,c(1,15:21)]

# P emp - Wilcoxon ----

EX_bootstrap_pemp_pv_wilcoxon <- readRDS(paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_wilcoxon_',STRATA_code,'.rds'))

ugo = melt(as.matrix(EX_bootstrap_pemp_pv_wilcoxon))
ugo = ugo[,c(1,3)]
colnames(ugo)[1]='event_id'
colnames(ugo)[2]='w'

l = ugo
l$obs.w <- obs$w[match(l$event_id,obs$event_id)]

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
l$obs.k <- obs$k[match(l$event_id,obs$event_id)]

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
l$obs.f <- obs$f[match(l$event_id,obs$event_id)]

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


# Integration in the "dex" dataset ----

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


# XXX. Select significantly differentially included ASE ----

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
dim(dex_CASE)

saveRDS(dex_CASE, file=paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'_bootstrap.rds'))


# Add info about gene symbol and standard deviation ---

map_ENSnames_symbols <- readRDS("Rdata/map_ENSnames_symbols.rds")
dex_CASE$symbol   = NA
dex_CASE$symbol   = map_ENSnames_symbols$symbol[match(dex_CASE$gene_id,map_ENSnames_symbols$ensembl_gene_id)]
dex_CASE$case.sd  = dex_CASE$case.se*sqrt(dex_CASE$case.N)
dex_CASE$cntr.sd  = dex_CASE$cntr.se*sqrt(dex_CASE$cntr.N)
dex_CASE$fc.sd    = foldchange(dex_CASE$case.sd, dex_CASE$cntr.sd)
dex_CASE$delta.sd = dex_CASE$case.sd-dex_CASE$cntr.sd

qtd = quantile(dex_CASE$delta.mean, seq(0,1,0.01))
qts = quantile(dex_CASE$delta.sd, seq(0,1,0.01))

pdf(file=paste0(Fig_out_dir,'/Quantiles_delta_mean_standard_dev_',STRATA_code,'.pdf'), paper = 'a4', useDingbats = F)
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

# XXX. Characterization of significan ASE ----

message('[*] Saving piechart of dexSIGN ...')

# Piechart event type
bp <- ggplot(dexSIGN_CASE, aes(x='', fill=event_type)) + geom_bar()
pie <- bp + coord_polar("y", start=0) 
pdf(file=paste0(Fig_out_dir,'/Piechart_dexSIGN_',STRATA_code,'.pdf'), paper = 'a4', useDingbats = F)
pie
dev.off()



# Scatter plots and histograms of the significant events ----

dexSIGN_CASE$delta.class.mean = NA
dexSIGN_CASE$delta.class.mean[which(dexSIGN_CASE$delta.mean>0)] = 'H mean'
dexSIGN_CASE$delta.class.mean[which(dexSIGN_CASE$delta.mean<0)] = 'L mean'
dexSIGN_CASE$delta.class.mean = as.factor(dexSIGN_CASE$delta.class.mean)
dexSIGN_CASE$delta.class.sd = NA
dexSIGN_CASE$delta.class.sd[which(dexSIGN_CASE$delta.sd>0)] = 'H sd'
dexSIGN_CASE$delta.class.sd[which(dexSIGN_CASE$delta.sd<0)] = 'L sd'
dexSIGN_CASE$delta.class.sd = as.factor(dexSIGN_CASE$delta.class.sd)

library(ggplot2)

pdf(file=paste0(Fig_out_dir,'/Scatter_delta.mean_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
ggplot(dexSIGN_CASE, aes(x=overall.mean, y=overall.sd, color=delta.class.mean)) +
  geom_point(size=1) + ggtitle(paste0(STRATA_code,' svASEs (',nrow(dexSIGN_CASE),')')) +
  scale_color_manual(values=c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))# + stat_function(fun = y, color='black') + geom_vline(xintercept = 0.15) + geom_vline(xintercept = 0.85)
dev.off()

pdf(file=paste0(Fig_out_dir,'/Scatter_delta.sd_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
ggplot(dexSIGN_CASE, aes(x=overall.mean, y=overall.sd, color=delta.class.sd)) +
  geom_point(size=1) + ggtitle(paste0(STRATA_code,' svASEs (',nrow(dexSIGN_CASE),')')) +
  scale_color_manual(values=c(rgb(242/255,117/255,109/255),rgb(26/255,166/255,184/255)))# + stat_function(fun = y, color='black') + geom_vline(xintercept = 0.15) + geom_vline(xintercept = 0.85)
dev.off()


pdf(file=paste0(Fig_out_dir,'/Histogram_delta.mean_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
ggplot(dexSIGN_CASE,aes(x=overall.mean,col=delta.class.mean))+geom_histogram(alpha=0, position="identity", bins=50)+theme_classic() +
  scale_color_manual(values=c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))
dev.off()

pdf(file=paste0(Fig_out_dir,'/Histogram_delta.sd_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
ggplot(dexSIGN_CASE,aes(x=overall.mean,col=delta.class.sd))+geom_histogram(alpha=0, position="identity", bins=50)+theme_classic() +
  scale_color_manual(values=c(rgb(242/255,117/255,109/255),rgb(26/255,166/255,184/255)))
dev.off()



# Vector field plots ----

p1=ggplot(dexSIGN_CASE,aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("ALL")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p2=ggplot(subset(dexSIGN_CASE,event_type=='ES'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("ES")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p3=ggplot(subset(dexSIGN_CASE,event_type=='A3'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("A3")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p4=ggplot(subset(dexSIGN_CASE,event_type=='A5'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("A5")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p5=ggplot(subset(dexSIGN_CASE,event_type=='IR'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("IR")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p6=ggplot(subset(dexSIGN_CASE,event_type=='MEX'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.mean)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(18/255,0,124/255), high=rgb(223/255,0,63/255))+ggtitle("MEX")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
pdf(file=paste0(Fig_out_dir,'/Scatter_plot_vector_field_color_mean_ALL_EVENTS_AND_TYPES_dexSIGN_',STRATA_code,'.pdf'), width = 17, height = 6.5, useDingbats = F)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
dev.off()

p1=ggplot(dexSIGN_CASE,aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("ALL")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p2=ggplot(subset(dexSIGN_CASE,event_type=='ES'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("ES")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p3=ggplot(subset(dexSIGN_CASE,event_type=='A3'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("A3")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p4=ggplot(subset(dexSIGN_CASE,event_type=='A5'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("A5")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p5=ggplot(subset(dexSIGN_CASE,event_type=='IR'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("IR")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
p6=ggplot(subset(dexSIGN_CASE,event_type=='MEX'),aes(x=cntr.mean,y=cntr.sd,xend=case.mean,yend=case.sd,col=sign(delta.sd)))+geom_segment(arrow = arrow(length = unit(0.1,"cm")))+scale_color_gradient(low=rgb(26/255,166/255,184/255), high=rgb(242/255,117/255,109/255))+ggtitle("MEX")+theme_bw()+ylim(0,42) # aes(x = x1, y = y1, xend = x2, yend = y2, colour = "segment"), data = df
pdf(file=paste0(Fig_out_dir,'/Scatter_plot_vector_field_color_sd_ALL_EVENTS_AND_TYPES_dexSIGN_',STRATA_code,'.pdf'), width = 17, height = 6.5, useDingbats = F)
grid.arrange(p1,p2,p3,p4,p5,p6,ncol=3)
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
pdf(file=paste0(Fig_out_dir,'/Histogram_delta_mean_4_classes_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
grid.arrange(p1,p2,p3,p4,ncol=4)
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
pdf(file=paste0(Fig_out_dir,'/Histogram_delta_sd_4_classes_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
grid.arrange(p1,p2,p3,p4,ncol=4)
dev.off()



# Box plots and Wilcoxon of the significant events ----

x = melt(dexSIGN_CASE, measure.vars = c("case.mean","cntr.mean"))
x$trend=x$value>50 
pdf(file=paste0(Fig_out_dir,'/Boxplots_mean_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
ggplot(x, aes(y=value, x=variable, alpha=0.5, fill=variable))+geom_boxplot(notch=T)+facet_grid(event_type~trend, scales = "free")+coord_flip()+
  stat_compare_means(label = "p.signif",label.y = 40,label.x = 1.5)
dev.off()

y=ddply(x, .(event_type, trend), summarise, pv=wilcox.test(value[variable=="case.mean"], value[variable=='cntr.mean'])$p.value)
subset(y, pv<0.05)

x = melt(dexSIGN_CASE, measure.vars = c("case.sd","cntr.sd"))
pdf(file=paste0(Fig_out_dir,'/Boxplots_sd_dexSIGN_',STRATA_code,'.pdf'), width = 10, height = 7, useDingbats = F)
ggplot(x, aes(y=value, x=variable, alpha=0.5, fill=variable))+geom_boxplot(notch=T)+facet_grid(. ~ x$event_type, scales = "free")+
  stat_compare_means(label = "p.signif",label.y = 30,label.x = 1.3)
dev.off()

y=ddply(x, .(event_type), summarise, pv=wilcox.test(value[variable=="case.sd"], value[variable=='cntr.sd'])$p.value)
subset(y, pv<0.05)



# Box plots and Wilcoxon of the significant events 4 CLASSES----

x = melt(dexSIGN_CASE, measure.vars = c("case.mean","cntr.mean"))

x$trend=NA
x$trend[which(x$overall.mean<=0.15)] = 1 
x$trend[which(x$overall.mean>0.15 & x$overall.mean<=0.5)] = 2
x$trend[which(x$overall.mean>0.5 & x$overall.mean<=0.85)] = 3
x$trend[which(x$overall.mean>0.85)] = 4
pdf(file=paste0(Fig_out_dir,'/Boxplots_mean_dexSIGN_',STRATA_code,'_4_classes_015_05_085.pdf'), width = 10, height = 7, useDingbats = F)
ggplot(x, aes(y=value, x=variable, alpha=0.5, fill=variable))+geom_boxplot(notch=T,outlier.shape = NA)+facet_grid(event_type~trend, scales = "free")+coord_flip()+theme_test()+
  stat_compare_means(label = "p.signif",label.y = (-Inf),label.x = 1.5)
dev.off()

y=ddply(x, .(event_type, trend), summarise, pv=wilcox.test(value[variable=="case.mean"], value[variable=='cntr.mean'])$p.value)
yy=ddply(x, .(event_type, trend), summarise, m.mean.case=mean(value[variable=="case.mean"]), m.mean.cntr=mean(value[variable=='cntr.mean']))
y$delta.mean.mean = yy$m.mean.case-yy$m.mean.cntr
subset(y, pv<0.05)

x = melt(dexSIGN_CASE, measure.vars = c("case.sd","cntr.sd"))
pdf(file=paste0(Fig_out_dir,'/Boxplots_sd_dexSIGN_',STRATA_code,'_noOutliers.pdf'), width = 10, height = 7, useDingbats = F)
ggplot(x, aes(y=value, x=variable, alpha=0.5, fill=variable))+geom_boxplot(notch=T,outlier.shape = NA)+facet_grid(. ~ x$event_type, scales = "free")+theme_test()
dev.off()

y=ddply(x, .(event_type), summarise, pv=wilcox.test(value[variable=="case.sd"], value[variable=='cntr.sd'])$p.value)
subset(y, pv<0.05)


# sd in 4 classes---
x = melt(dexSIGN_CASE, measure.vars = c("case.sd","cntr.sd"))
x$trend=NA
x$trend[which(x$overall.mean<=0.15)] = 1 
x$trend[which(x$overall.mean>0.15 & x$overall.mean<=0.5)] = 2
x$trend[which(x$overall.mean>0.5 & x$overall.mean<=0.85)] = 3
x$trend[which(x$overall.mean>0.85)] = 4
pdf(file=paste0(Fig_out_dir,'/Boxplots_sd_4_classes_dexSIGN_',STRATA_code,'_noOutliers.pdf'), width = 10, height = 7, useDingbats = F)
ggplot(x, aes(y=value, x=variable, alpha=0.5, fill=variable))+geom_boxplot(notch=T,outlier.shape = NA)+facet_grid(event_type~trend, scales = "free")+coord_flip()+theme_test()+ylim(0,20)+
  stat_compare_means(label = "p.signif",label.y = (-Inf),label.x = 1.5)
dev.off()

y=ddply(x, .(event_type, trend), summarise, pv=wilcox.test(value[variable=="case.sd"], value[variable=='cntr.sd'])$p.value)
yy=ddply(x, .(event_type, trend), summarise, m.mean.case=mean(value[variable=="case.sd"]), m.mean.cntr=mean(value[variable=='cntr.sd']))
y$delta.mean.mean = yy$m.mean.case-yy$m.mean.cntr
subset(y, pv<0.05)





# Subset of ES events ----


# dexSIGN ---- 

message('[*] Loading dexSIGN dataset...')

dexSIGN <- readRDS(paste0(Data_out_dir,'/dexSIGN_',STRATA_code,'.rds'))

dexSIGN_ES = subset(dexSIGN, event_type=='ES')

dexSIGN_ES$event_chr = paste0('chr',dexSIGN_ES$event_chr)

dexSIGN_ES$alt_coord_1 = sapply(strsplit(as.character(dexSIGN_ES$alt_region_coordinates),"[:]"), `[`, 1)
dexSIGN_ES$alt_coord_2 = sapply(strsplit(as.character(dexSIGN_ES$alt_region_coordinates),"[:]"), `[`, 2)

dexSIGN_ES$ensemble = sapply(strsplit(as.character(dexSIGN_ES$gene_name),"[.]"), `[`, 1)
dexSIGN_ES$strand = coord$strand[match(dexSIGN_ES$ensemble,coord$GeneName)] 

saveRDS(dexSIGN_ES, file=paste0(Data_out_dir,'/dexSIGN_EXON_SKIPPING_only_',STRATA_code,'.rds'))



# dex ---- 

message('[*] Loading dexVAR dataset...')

dexVAR <- readRDS(paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'_Variant.rds'))

dexVAR_ES = subset(dexVAR, event_type=='ES')

dexVAR_ES$event_chr = paste0('chr',dexVAR_ES$event_chr)

dexVAR_ES$alt_coord_1 = sapply(strsplit(as.character(dexVAR_ES$alt_region_coordinates),"[:]"), `[`, 1)
dexVAR_ES$alt_coord_2 = sapply(strsplit(as.character(dexVAR_ES$alt_region_coordinates),"[:]"), `[`, 2)

dexVAR_ES$ensemble = sapply(strsplit(as.character(dexVAR_ES$gene_name),"[.]"), `[`, 1)
dexVAR_ES$strand = coord$strand[match(dexVAR_ES$ensemble,coord$GeneName)] 

saveRDS(dexVAR_ES, file=paste0(Data_out_dir,'/dexVAR_EXON_SKIPPING_only_',STRATA_code,'.rds'))






# Gene ontology ----

library(DOSE)
library(igraph)
library(org.Hs.eg.db)
library(tidyverse)
library(magrittr)
library(plyr)
library(clusterProfiler)
library(msigdbr)

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
  ego <- enricher(unique(ids$ENTREZID)  
                  , TERM2GENE = TERM2GENE
                  , pAdjustMethod = pAdjustMethod
                  , pvalueCutoff=pvalueCutoff
                  , qvalueCutoff=qvalueCutoff
                  , minGSSize = minGSSize )
  ego = format_enrichr(ego, org=OrgDb, key=typeFormat)
  list('id'=ids, 'go'=ego)
}

kegg = read.csv('Tables/KEGG.csv')

dexSIGN = readRDS(paste0(Data_out_dir,'/dexSIGN_',STRATA_code,'.rds'))
dexSIGN$ensembl = dexSIGN$gene_id
dexSIGN$Type = dexSIGN$event_type

cls = c('event_id','ensembl','Type')

ase = dexSIGN[,cls]



# % ORA with 186 KEGG ==================

pll = load('Rdata/KEGG_list.148.SPLICING_RELATED.Rdata')

library(msigdbr)
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


pdf(file=paste0(Fig_out_dir,'/ORA_all_event_types_TCGA_KEGG_186.pdf'), height = unit(4,'cm'), width = unit(8, 'cm'), useDingbats = F)
to_plot=subset(to_plot,p.adjust<0.1)
ggplot( to_plot,
        aes(x= parse_ratio(GeneRatio), y= ID)) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=factor(p.adjust<=.25), fill=p.adjust, size = Count),shape=21) +
  geom_text(aes(label=rank)) +
  scale_color_manual(values=c('transparent','black'))+
  scale_fill_viridis_c(option = 'D',guide=guide_colorbar(reverse=T, draw.llim = T), direction = -1)+#,
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + xlab("Gene Ratio") +ylab(NULL)+
  facet_wrap(~type, nrow=1)
dev.off()



pdf(file=paste0(Fig_out_dir,'/ORA_all_event_types_TCGA_KEGG_186__richFactor.pdf'), height = unit(4,'cm'), width = unit(8, 'cm'), useDingbats = F)
to_plot=subset(to_plot,p.adjust<0.1)
ggplot( to_plot,
        aes(x= richFactor, y= ID)) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=factor(p.adjust<=.25), fill=p.adjust, size = Count),shape=21) +
  geom_text(aes(label=rank)) +
  scale_color_manual(values=c('transparent','black'))+
  scale_fill_viridis_c(option = 'D',guide=guide_colorbar(reverse=T, draw.llim = T), direction = -1)+#,
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + xlab("richFactor") +ylab(NULL)+
  facet_wrap(~type, nrow=1)
dev.off()




# ******************** ----
# ******************** Cumulative distr. of events ----
# ******************** ----


library(ggplot2)
library(gridExtra)
library(ggplotify)

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



# grid.arrange(p_mean,p_sd,p_mean_counts,p_sd_counts,ncol=2)
# 
# grid.arrange(p_ratio_mean,p_ratio_sd,ncol=2)



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

save(res_df_mean,res_df_sd,file='Rdata/null_model_splicing_prob_05_for_delta_pos_and_delta_neg.Rdata')




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



pdf(file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg.pdf',width = unit(10,'cm'),height = unit(3.5,'cm'))
grid.arrange(p_mean_compl,p_sd_compl,ncol=2)
dev.off()




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


pdf(file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_SPLITTED.pdf',width = unit(10,'cm'),height = unit(3.5,'cm'))
grid.arrange(p_mean_compl_2,p_sd_compl_2,ncol=2)
dev.off()



# %%%% INVERTED ----


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



pdf(file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_inverted.pdf',width = unit(10,'cm'),height = unit(3.5,'cm'))
grid.arrange(p_mean_compl_inv,p_sd_compl_inv,ncol=2)
dev.off()
  
pdf(file = 'TCGA_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_SPLITTED_inverted.pdf',width = unit(10,'cm'),height = unit(3.5,'cm'))
grid.arrange(p_mean_compl_2_inv,p_sd_compl_2_inv,ncol=2)
dev.off()





