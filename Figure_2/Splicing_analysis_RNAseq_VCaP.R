# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Splcing analysis of RNAseq data of VCaP cells ----

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

repo = '/Volumes/Prule/repo/FOXA1_public/'
FIG_DIR = '~/Downloads/FOXA1/Figure2/'
setwd(paste0(repo, "sources/Whippet_raw_files/VCaP"))


anno=readRDS('../../../Rdata/gencode.v28.annotation.rds')

# PSI ==== 

# NSI ----

prova = read.table('nsi-R1.psi',sep = '\t',header = T)
prova$chr = sapply(strsplit(as.character(prova$Coord),':'),'[',1)
prova$coord_range = sapply(strsplit(as.character(prova$Coord),':'),'[',2)
prova$start = sapply(strsplit(as.character(prova$coord_range),'-'),'[',1)
prova$end = sapply(strsplit(as.character(prova$coord_range),'-'),'[',2)
prova$start = as.numeric(prova$start)
prova$end = as.numeric(prova$end)

nsi_1_psi = prova
colnames(nsi_1_psi)[6]="Psi_nsi_r1"
nsi_1_psi$gene_name=anno$gene_name[match(nsi_1_psi$Gene, anno$gene_id)]


prova = read.table('nsi-R2.psi',sep = '\t',header = T)
prova$chr = sapply(strsplit(as.character(prova$Coord),':'),'[',1)
prova$coord_range = sapply(strsplit(as.character(prova$Coord),':'),'[',2)
prova$start = sapply(strsplit(as.character(prova$coord_range),'-'),'[',1)
prova$end = sapply(strsplit(as.character(prova$coord_range),'-'),'[',2)
prova$start = as.numeric(prova$start)
prova$end = as.numeric(prova$end)

nsi_2_psi = prova
colnames(nsi_2_psi)[6]="Psi_nsi_r2"


prova = read.table('nsi-R4.psi',sep = '\t',header = T)
prova$chr = sapply(strsplit(as.character(prova$Coord),':'),'[',1)
prova$coord_range = sapply(strsplit(as.character(prova$Coord),':'),'[',2)
prova$start = sapply(strsplit(as.character(prova$coord_range),'-'),'[',1)
prova$end = sapply(strsplit(as.character(prova$coord_range),'-'),'[',2)
prova$start = as.numeric(prova$start)
prova$end = as.numeric(prova$end)

nsi_4_psi = prova
colnames(nsi_4_psi)[6]="Psi_nsi_r4"


# siFOXA1 ----

prova = read.table('siFOXA1-R1.psi',sep = '\t',header = T)
prova$chr = sapply(strsplit(as.character(prova$Coord),':'),'[',1)
prova$coord_range = sapply(strsplit(as.character(prova$Coord),':'),'[',2)
prova$start = sapply(strsplit(as.character(prova$coord_range),'-'),'[',1)
prova$end = sapply(strsplit(as.character(prova$coord_range),'-'),'[',2)
prova$start = as.numeric(prova$start)
prova$end = as.numeric(prova$end)

si_1_psi = prova
colnames(si_1_psi)[6]="Psi_siFOXA1_r1"


prova = read.table('siFOXA1-R2.psi',sep = '\t',header = T)
prova$chr = sapply(strsplit(as.character(prova$Coord),':'),'[',1)
prova$coord_range = sapply(strsplit(as.character(prova$Coord),':'),'[',2)
prova$start = sapply(strsplit(as.character(prova$coord_range),'-'),'[',1)
prova$end = sapply(strsplit(as.character(prova$coord_range),'-'),'[',2)
prova$start = as.numeric(prova$start)
prova$end = as.numeric(prova$end)

si_2_psi = prova
colnames(si_2_psi)[6]="Psi_siFOXA1_r2"


prova = read.table('siFOXA1-R4.psi',sep = '\t',header = T)
prova$chr = sapply(strsplit(as.character(prova$Coord),':'),'[',1)
prova$coord_range = sapply(strsplit(as.character(prova$Coord),':'),'[',2)
prova$start = sapply(strsplit(as.character(prova$coord_range),'-'),'[',1)
prova$end = sapply(strsplit(as.character(prova$coord_range),'-'),'[',2)
prova$start = as.numeric(prova$start)
prova$end = as.numeric(prova$end)

si_4_psi = prova
colnames(si_4_psi)[6]="Psi_siFOXA1_r4"

all_psi=cbind.data.frame(nsi_1_psi[,c("Gene", "gene_name", "Node","Coord","Strand","Type", "Psi_nsi_r1")],nsi_2_psi$Psi_nsi_r2, nsi_4_psi$Psi_nsi_r4,
                         si_1_psi$Psi_siFOXA1_r1, si_2_psi$Psi_siFOXA1_r2, si_4_psi$Psi_siFOXA1_r4)

colnames(all_psi)[8:12]=c("Psi_nsi_r2", "Psi_nsi_r4", "Psi_siFOXA1_r1", "Psi_siFOXA1_r2", "Psi_siFOXA1_r4")
saveRDS(all_psi,"../../../Rdata/Psi_all_events_all_samples_VCaP.rds")



# CI_width ==========

prova = read.table('nsi-R1.psi',sep = '\t',header = T)

nsi_1_psi = prova
colnames(nsi_1_psi)[7]="CI_width_nsi_r1"


prova = read.table('nsi-R2.psi',sep = '\t',header = T)
nsi_2_psi = prova
colnames(nsi_2_psi)[7]="CI_width_nsi_r2"



prova = read.table('nsi-R4.psi',sep = '\t',header = T)
nsi_4_psi = prova
colnames(nsi_4_psi)[7]="CI_width_nsi_r4"



prova = read.table('siFOXA1-R1.psi',sep = '\t',header = T)
si_1_psi = prova
colnames(si_1_psi)[7]="CI_width_siFOXA1_r1"


prova = read.table('siFOXA1-R2.psi',sep = '\t',header = T)
si_2_psi = prova
colnames(si_2_psi)[7]="CI_width_siFOXA1_r2"


prova = read.table('siFOXA1-R4.psi',sep = '\t',header = T)
si_4_psi = prova
colnames(si_4_psi)[7]="CI_width_siFOXA1_r4"


all_CI_width=cbind.data.frame(nsi_1_psi[,c("Gene", "Node","Coord","Strand","Type", "CI_width_nsi_r1")],nsi_2_psi$CI_width_nsi_r2, nsi_4_psi$CI_width_nsi_r4,
                              si_1_psi$CI_width_siFOXA1_r1, si_2_psi$CI_width_siFOXA1_r2, si_4_psi$CI_width_siFOXA1_r4)

colnames(all_CI_width)[7:11]=c("CI_width_nsi_r2", "CI_width_nsi_r4", "CI_width_siFOXA1_r1", "CI_width_siFOXA1_r2", "CI_width_siFOXA1_r4")
saveRDS(all_CI_width,"../../../Rdata/CI_width_all_events_all_samples_VCaP.rds")




# Build final dataframe ----

prova = read.table('Comparisons_nsi_siFOXA1.diff',sep = '\t',header = T,row.names = NULL)
dim(prova)
tmp = colnames(prova)
# 256295
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

saveRDS(delta, "../../../Rdata/Whippet_VCAP_index_BAM.rds")

all_CI_width=readRDS("../../../Rdata/CI_width_all_events_all_samples_VCaP.rds")
all_CI_width$key=paste0(all_CI_width$Gene,'_',all_CI_width$Node,'_',all_CI_width$Type)
CI_ok=subset(all_CI_width, CI_width_nsi_r1<0.2 &
               CI_width_nsi_r2<0.2 & 
               CI_width_nsi_r4<0.2 & 
               CI_width_siFOXA1_r1<0.2 &
               CI_width_siFOXA1_r2<0.2 &
               CI_width_siFOXA1_r4<0.2)$key


delta = subset(delta,key%in%CI_ok)

saveRDS(delta, "../../../Rdata/Whippet_VCAP_filter_0_2.rds")

# Select significantly differentially included events ---
delta$STRINGENT = F
delta$STRINGENT[(abs(delta$DeltaPsi)>0.05 & delta$Probability>0.9 & !delta$Complexity=='K0')] = T

saveRDS(delta, "../../../Rdata/Whippet_VCAP_filter_0_2.rds")



# *** ANALYSIS --------

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

setwd(repo)

delta=readRDS("Rdata/Whippet_VCAP_filter_0_2.rds")

delta=subset(delta, !Type%in%c("TE","TS"))

saveRDS(delta,"Rdata/Whippet_VCAP_filter_0_2_no_TS_TE.rds")

delta=readRDS("Rdata/Whippet_VCAP_filter_0_2_no_TS_TE.rds")

delta$overall.mean = with(delta, Psi_A + Psi_B)/2

prova_sign = subset(delta,STRINGENT & !Complexity=='K0')


# Piechart event type:
bp <- ggplot(prova_sign, aes(x='', fill=Type)) + geom_bar()

stats=as.data.frame(table(prova_sign$Type))
stats$distr="events"


pdf(file=paste0(FIG_DIR, 'piechart_event_types_STRINGENT.pdf'), paper = 'a4', useDingbats = F)

ggplot(stats, aes(x = distr, y = Freq, fill =Var1)) +theme_bw()+
  geom_bar(stat="identity", width = 0.7)+theme(axis.text.x = element_text(angle = 90, hjust = 1))+geom_text(aes(label=Freq),position = position_stack(0.6) )+
  scale_fill_manual(values= c("RI"="darkorchid1",
                              "CE"="royalblue1",
                              "TE"="plum2",
                              "TS"="turquoise3",
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
pdf(file=paste0(FIG_DIR, 'histogram_colored_DeltaPsi_CIwidth_less_02_STRINGENT.pdf'), height = unit(3,'cm'), width = unit(5.5,'cm'), useDingbats = F )
ggplot(prova_sign,aes(x=overall.mean,color=class_delta))+geom_histogram(alpha=0, position="identity", bins=50)+theme_classic() +
  scale_color_manual(values=c(rgb(18/255,0,124/255),rgb(223/255,0,63/255)))
dev.off()

pdf(file=paste0(FIG_DIR, 'histogram_deltaPSI_single_event_types_STRINGENT.pdf'), height = unit(8,'cm'), width = unit(6,'cm'), useDingbats = F )
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










