# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Splcing analysis of RNAseq data of PC3 cells ----


options(stringsAsFactors = F)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gtools)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(igraph)
library(org.Hs.eg.db)
library(gridExtra)
library(rtracklayer)

setwd("sources/Whippet_raw_files/PC3")

anno=readRDS('../../../Rdata/gencode.v28.annotation.rds')

# PSI ====

f  <- list.files(path=".",
                 recursive=F,
                 pattern="psi"
                 ,full.names=T)

a=read.delim2(f[1], header=T)
name="Psi_nsi-r6"
a$key=paste0(a$Gene, "_", a$Node, "_", a$Type)
a$Psi=as.numeric(a$Psi)

a=a[, c("key", "Psi")]

colnames(a)[ncol(a)]=name

for(i in f[2:length(f)]){
  
  name=gsub("./pc3-", "", i)
  name=gsub(".psi.gz", "", name)
  
  b=read.delim2(i, header=T)
  
  b$key=paste0(b$Gene, "_", b$Node, "_", b$Type)
  b$Psi=as.numeric(b$Psi)
  
  a=cbind(a, b$Psi[match(a$key, b$key)])
  colnames(a)[ncol(a)]=paste0("Psi_",name)
  
}


saveRDS(a,"../../../Rdata/Psi_all_events_all_samples_PC3_nextseq_nova.rds")



# CI_width ==========


f  <- list.files(path=".",
                 recursive=F,
                 pattern="psi"
                 ,full.names=T)

a=read.delim2(f[1], header=T)
name="CI_Width_nsi-r6"
a$key=paste0(a$Gene, "_", a$Node, "_", a$Type)
a$CI_Width=as.numeric(a$CI_Width)

a=a[, c("key", "CI_Width")]

colnames(a)[ncol(a)]=name

for(i in f[2:length(f)]){
  
  name=gsub("./pc3-", "", i)
  name=gsub(".psi.gz", "", name)
  
  b=read.delim2(i, header=T)
  
  b$key=paste0(b$Gene, "_", b$Node, "_", b$Type)
  b$CI_Width=as.numeric(b$CI_Width)
  
  a=cbind(a, b$CI_Width[match(a$key, b$key)])
  colnames(a)[ncol(a)]=paste0("CI_Width_",name)
  
}


saveRDS(a,"../../../Rdata/CI_width_all_events_all_samples_PC3_nextseq_nova.rds")



# Build final dataframe ----

prova = read.table('Comparisons_nsi_siFOXA1.diff.gz',sep = '\t',header = T,row.names = NULL)
dim(prova)
tmp = colnames(prova)
# 326373
prova$Entropy = NULL
colnames(prova) = tmp[2:length(tmp)]
prova$chr = sapply(strsplit(as.character(prova$Coord),':'),'[',1)
prova$coord_range = sapply(strsplit(as.character(prova$Coord),':'),'[',2)
prova$start = sapply(strsplit(as.character(prova$coord_range),'-'),'[',1)
prova$end = sapply(strsplit(as.character(prova$coord_range),'-'),'[',2)
prova$start = as.numeric(prova$start)
prova$end = as.numeric(prova$end)

# Filter events with CI_width>0.2 ---
delta = prova
delta$key = paste0(delta$Gene,'_',delta$Node,'_',delta$Type)

saveRDS(delta, "../../../Rdata/Whippet_PC3_nextseq_nova.rds")


all_CI_width=readRDS("../../../Rdata/CI_width_all_events_all_samples_PC3_nextseq_nova.rds")
CI_ok=subset(all_CI_width, `CI_Width_nsi-r6`<0.2 &
               `CI_Width_nsi-r7`<0.2 & 
               `CI_Width_nsi-r8`<0.2 & 
               `CI_Width_siFOXA1-r6`<0.2 &
               `CI_Width_siFOXA1-r7`<0.2 &
               `CI_Width_siFOXA1-r8`<0.2)


delta = subset(delta,key%in%CI_ok$key)

saveRDS(delta, "../../../Rdata/Whippet_PC3_nextseq_novaseq_filter_0_2.rds")

# Select significantly differentially included events ---
delta$STRINGENT = F
delta$STRINGENT[(abs(delta$DeltaPsi)>0.05 & delta$Probability>0.9 & !delta$Complexity=='K0')] = T

saveRDS(delta, "../../../Rdata/Whippet_PC3_nextseq_novaseq_filter_0_2.rds")


# *** ANALYSIS ----

options(stringsAsFactors = F)
library(plyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gtools)
library(clusterProfiler)
library(ReactomePA)
library(DOSE)
library(igraph)
library(org.Hs.eg.db)
library(gridExtra)
library(rtracklayer)
library(ggExtra)

setwd("code")

delta=readRDS("Rdata/Whippet_PC3_nextseq_novaseq_filter_0_2.rds")

delta=subset(delta, !Type%in%c("TE","TS"))

saveRDS(delta, "Rdata/Whippet_PC3_nextseq_novaseq_filter_0_2.rds")

delta=readRDS("Rdata/Whippet_PC3_nextseq_novaseq_filter_0_2.rds")

delta$overall.mean = with(delta, Psi_A + Psi_B)/2

prova_sign = subset(delta,STRINGENT & !Complexity=='K0')


# Piechart event type:
bp <- ggplot(prova_sign, aes(x='', fill=Type)) + geom_bar()

stats=as.data.frame(table(prova_sign$Type))
stats$distr="events"


pdf(file='piechart_event_types_STRINGENT.pdf', paper = 'a4', useDingbats = F)

ggplot(stats, aes(x = distr, y = Freq, fill =Var1)) +theme_bw()+
  geom_bar(stat="identity", width = 0.7)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_text(aes(label=Freq),position = position_stack(0.6) )+
  scale_fill_manual(values= c("RI"="darkorchid1",
                              "CE"="royalblue1",
                              "AL"="turquoise4",
                              "AF"="seagreen1",
                              "AA"="aquamarine",
                              "AD"="aquamarine3"
  ))+coord_polar("y", start=0)


dev.off()


# Histogram colored with DeltaPsi:

prova_sign$overall.mean = (prova_sign$Psi_A+prova_sign$Psi_B)/2
prova_sign$class_delta = '<0'
prova_sign$class_delta[which(sign(prova_sign$DeltaPsi)>0)] = '>0' 
pdf(file='histogram_colored_DeltaPsi_CIwidth_less_02_STRINGENT.pdf', height = unit(3,'cm'), width = unit(5.5,'cm'), useDingbats = F )
ggplot(prova_sign,aes(x=overall.mean,color=class_delta))+geom_histogram(alpha=0, position="identity", bins=50)+theme_classic() +
  scale_color_manual(values=c(rgb(18/255,0,124/255),rgb(223/255,0,63/255)))
dev.off()

pdf(file='histogram_deltaPSI_single_event_types_STRINGENT.pdf', height = unit(8,'cm'), width = unit(6,'cm'), useDingbats = F )
ggplot(prova_sign,aes(x=overall.mean,color=class_delta))+geom_histogram(alpha=0, position="identity", bins=50)+theme_classic() +
  scale_color_manual(values=c(rgb(18/255,0,124/255),rgb(223/255,0,63/255))) + facet_wrap(~ Type, ncol=2)
dev.off()



# Boxplots DeltaPsi:

x = melt(prova_sign, measure.vars = c("Psi_B","Psi_A"))
x$trend=x$overall.mean>0.5  
library(ggpubr)
pdf(file='boxplot_meanPSI_single_event_types_CIwidth_less_02_STRINGENT.pdf', height = unit(7,'cm'), width = unit(10,'cm'), useDingbats = F )
ggplot(x, aes(y=value, x=variable, alpha=0.5, fill=variable))+geom_boxplot(notch=T)+facet_grid(Type~trend, scales = "free")+coord_flip()+theme_test()+
  stat_compare_means(label = "p.signif",label.y = 0.4,label.x = 1.5)
dev.off()
pdf(file='boxplot_meanPSI_All_event_types_CIwidth_less_02_STRINGENT.pdf', height = unit(3,'cm'), width = unit(10,'cm'), useDingbats = F )
ggplot(x, aes(y=value, x=variable, alpha=0.5, fill=variable))+geom_boxplot(notch=T)+facet_grid(~trend, scales = "free")+coord_flip()+theme_test()+
  stat_compare_means(label = "p.signif",label.y = 0.4,label.x = 1.5)
dev.off()

x = melt(prova_sign, measure.vars = c("Psi_B","Psi_A"))
x$trend=NA
x$trend[which(x$overall.mean<=0.15)] = 1 # %%%%%%%%%%%%%%%%%%%%%%% 
x$trend[which(x$overall.mean>0.15 & x$overall.mean<=0.5)] = 2
x$trend[which(x$overall.mean>0.5 & x$overall.mean<=0.85)] = 3
x$trend[which(x$overall.mean>0.85)] = 4
pdf(file='boxplot_meanPSI_single_event_types_4_classes_015_05_085_CIwidth_less_02_STRINGENT.pdf', height = unit(7,'cm'), width = unit(10,'cm'), useDingbats = F )
ggplot(x, aes(y=value, x=variable, alpha=0.5, fill=variable))+geom_boxplot(notch=T,outlier.shape = NA)+facet_grid(Type~trend, scales = "free")+coord_flip()+theme_test()+
  stat_compare_means(label = "p.signif",label.y = (-Inf),label.x = 1.5)
dev.off()
pdf(file='boxplot_meanPSI_All_event_types_4_classes_015_05_085_CIwidth_less_02_STRINGENT.pdf', height = unit(3,'cm'), width = unit(10,'cm'), useDingbats = F )
ggplot(x, aes(y=value, x=variable, alpha=0.5, fill=variable))+geom_boxplot(notch=T,outlier.shape = NA)+facet_grid(~trend, scales = "free")+coord_flip()+theme_test()+
  stat_compare_means(label = "p.signif",label.y = (-Inf),label.x = 1.5)
dev.off()


