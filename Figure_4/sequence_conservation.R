# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Sequence conservation ----


library(tidyverse)
library(dplyr)
library(plyr)
library(ggpubr)
library(tidyr)

# FIG_DIR = '/Volumes/Prule/remove_me_Figures_FOXA1/'
# system(paste0('mkdir ', FIG_DIR))

# ex = readRDS('Rdata/variable_exons_with_splice_sites_strength_annotation.rds')
ex = readRDS('variable_exons_with_splice_sites_strength_annotation.rds')
ex$upstream_5ss_donor_motif = ex$left_donor_5ss_motif
ex$upstream_5ss_donor_motif[ex$strand=='-'] = ex$right_donor_5ss_motif[ex$strand=='-']
ex$upstream_3ss_acceptor_motif = ex$left_acceptor_3ss_motif
ex$upstream_3ss_acceptor_motif[ex$strand=='-'] = ex$right_acceptor_3ss_motif[ex$strand=='-']

ex$downstream_5ss_donor_motif = ex$right_donor_5ss_motif
ex$downstream_5ss_donor_motif[ex$strand=='-'] = ex$left_donor_5ss_motif[ex$strand=='-']
ex$downstream_3ss_acceptor_motif = ex$right_acceptor_3ss_motif
ex$downstream_3ss_acceptor_motif[ex$strand=='-'] = ex$left_acceptor_3ss_motif[ex$strand=='-']

# Conservation ======
# library(GenomicScores)

cns = ex %>% subset(event_type=='ES' ) %>% 
  dplyr::select(event_chr,event_coordinates,gene_id, strand,
                gene_id,symbol,event_id, delta.mean, exon_type, FOXA1_sign,
                pois_ess, exon_length, intron_length) %>% 
  separate(event_coordinates,c('up3ss','up5ss','3ss','5ss','dw3ss','dw5ss'), convert=T) 

cns$ASE_set = 'vASE'
cns$ASE_set[cns$FOXA1_sign] = "FOXA1_regulated"
cns$ASE_set=factor(cns$ASE_set, levels=c('FOXA1_regulated', 'vASE'))
cns$type = c('Sil','Enh')[(sign(cns$delta.mean)>0)+1]

cns = cns[order(cns$event_chr, cns$`3ss`),]


cns$r3s = ifelse(cns$strand=='+', cns$`3ss`-150, cns$`5ss`-25)
cns$r3e = ifelse(cns$strand=='+', cns$`3ss`+25,  cns$`5ss`+150)

cns$r4s = ifelse(cns$strand=='+', cns$`5ss`-25, cns$`3ss`-150)
cns$r4e = ifelse(cns$strand=='+', cns$`5ss`+150, cns$`3ss`+25)


library(GenomicRanges)
library(rtracklayer)
r3 = GRanges( seqnames = cns$event_chr, IRanges(start = cns$r3s, end = cns$r3e)
              , strand = cns$strand, event_id=cns$event_id)


# DOWNLOAD CONSERVATION DATA FROM http://hgdownload.cse.ucsc.edu/goldenPath/hg19/phyloP100way/README.txt 


library(doParallel)
# ON CLUSTER ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
registerDoParallel(cores=28) ## each process may need up to 20Gb of RAM
ll=list()
ll = foreach (chr=levels(seqnames(r3))) %dopar% {
  rawscores=import.bw(BigWigFile(
    file.path("/hpcnfs/data/reference/ucsc/hg19/hg19.100way.phyloP/hg19.100way.phyloP100way/"
              , sprintf("%s.phyloP100way.bw", chr))))
  ex = subset(r3, seqnames==chr)
  ov = findOverlaps(ex, rawscores)
  res=as.data.frame(ex[queryHits(ov),])[,c('event_id','strand')]
  res$score=score(rawscores[subjectHits(ov),])
  res$event_id = as.character(res$event_id)
  res$event_id = factor(res$event_id,unique(res$event_id))
  res$strand = as.character(res$strand)
  res$strand = factor(res$strand,unique(res$strand))
  res = ddply(res, .(event_id), mutate, pos = ifelse(strand=='+', 1:176, 176:1), .progress = 'text')
  ll[[chr]]=res
}
# ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞


r3 = do.call(rbind.data.frame, ll)
r3 = left_join(r3, cns[, c('ASE_set','pois_ess','type','event_id')], key='event_id')
r3$ss = "3'ss"

# saveRDS(r3,'Rdata/phyloP_TCGA_3ss_ALL_EVENTS.rds')
saveRDS(r3,'phyloP_TCGA_3ss_ALL_EVENTS.rds')


subset(r3, !is.na(pois_ess)) %>% ggplot( aes(pos, score, group=ASE_set, color=ASE_set))+
  geom_smooth(method = 'gam', se = T)+ggpubr::theme_pubr()+facet_wrap(pois_ess~type)

ss = ddply(r3, .(ASE_set, pois_ess, type, pos),
           summarise, m = mean(score), md=median(score), sd=sd(score))


x11()
ggplot(ss, aes(pos, md, group=ASE_set, color=ASE_set))+
  geom_line()+
  ggpubr::theme_pubr()+facet_wrap(pois_ess~type)

ddply(r3, .(ASE_set, pos),
      summarise, m = mean(score), md=median(score), sd=sd(score)) %>%
  ggplot(aes(pos, md, group=ASE_set, color=ASE_set))+
  geom_errorbar(aes(ymin=md-sd, ymax=md+sd))+
  geom_line()+
  ggpubr::theme_pubr()
+facet_wrap(~type)


r4 = GRanges( seqnames = cns$event_chr, IRanges(start = cns$r4s, end = cns$r4e)
              , strand = cns$strand, event_id=cns$event_id)
# r4 = left_join(r4, cns[, c('ASE_set','pois_ess','type','event_id')], key='event_id')
# Error in UseMethod("left_join") : 
  # no applicable method for 'left_join' applied to an object of class "c('GRanges', 'GenomicRanges', 'Ranges', 'GenomicRanges_OR_missing', 'GenomicRanges_OR_GenomicRangesList', 'GenomicRanges_OR_GRangesList', 'List', 'Vector', 'list_OR_List', 'Annotated', 'vector_OR_Vector')"
# r4 = merge(r4, cns[, c('ASE_set','pois_ess','type','event_id')], by = 'event_id', all.x = T)

# 
# r3$ss = "3'ss"
# r4$ss = "5'ss"
# 
# prof = rbind.data.frame(r3,r4)
# 
# ttest = ddply(prof, .(ss, ASE_set, pos),
#               summarise, m = mean(score), md=median(score), sd=sd(score)) %>%
#   ddply( .(ss), function(x) with(x,broom::tidy(t.test(md[ASE_set!='vASE'], md[ASE_set=='vASE'], paired = T))))
# 
# ttest2=ddply(prof, .(ss, pois_ess, ASE_set, pos),
#              summarise, m = mean(score), md=median(score), sd=sd(score)) %>%
#   ddply( .(pois_ess,ss), function(x) with(x,broom::tidy(t.test(md[ASE_set!='vASE'], md[ASE_set=='vASE'], paired = T))))
# 
# ttest3=ddply(prof, .(ss, pois_ess, type, ASE_set, pos),
#              summarise, m = mean(score), md=median(score), sd=sd(score)) %>%
#   ddply( .(pois_ess,type,ss), function(x) with(x,broom::tidy(t.test(md[ASE_set!='vASE'], md[ASE_set=='vASE'], paired = T))))%>%
#   arrange(dplyr::desc(p.value))
# 
# # saveRDS(prof, file='Rdata/phylop_3ss_5ss_TCGA_ALL_EVENTS.rds')
# saveRDS(prof, file='phylop_3ss_5ss_TCGA_ALL_EVENTS.rds')
# 
# options(bitmapType="cairo")
# subset(prof, !is.na(pois_ess)) %>% ggplot( aes(pos, score, group=ASE_set, color=ASE_set))+
#   geom_smooth(method = 'gam', se = T)+ggpubr::theme_pubr()+facet_wrap(pois_ess~type+ss, ncol=2)


# ON CLUSTER ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞
registerDoParallel(cores=28) ## each process may need up to 20Gb of RAM
ll = list()
ll = foreach (chr=levels(seqnames(r4))) %dopar% {
  
  rawscores=import.bw(BigWigFile(
    file.path("/hpcnfs/data/reference/ucsc/hg19/hg19.100way.phyloP/hg19.100way.phyloP100way"
              , sprintf("%s.phyloP100way.bw", chr))))
  ex = subset(r4, seqnames==chr)
  ov = findOverlaps(ex, rawscores)
  res=as.data.frame(ex[queryHits(ov),])[,c('event_id','strand')]
  res$score=score(rawscores[subjectHits(ov),])
  res$event_id = as.character(res$event_id)
  res$event_id = factor(res$event_id,unique(res$event_id))
  res$strand = as.character(res$strand)
  res$strand = factor(res$strand,unique(res$strand))
  res = ddply(res, .(event_id), mutate, pos = ifelse(strand=='+', 1:176, 176:1), .progress = 'text')
  ll[[chr]]=res
}
# ∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞∞

r4 = do.call(rbind.data.frame, ll)
r4$ss = "5'ss"
r4 = left_join(r4, cns[, c('ASE_set','pois_ess','type','event_id')], key='event_id')

# saveRDS(r4,'Rdata/phyloP_TCGA_5ss_ALL_EVENTS.rds')
saveRDS(r4,'phyloP_TCGA_5ss_ALL_EVENTS.rds')

ddply(r3, .(ASE_set, pos),
      summarise, m = mean(score), md=median(score), sd=sd(score)) %>%
  ggplot(aes(pos, md, group=ASE_set, color=ASE_set))+
  geom_errorbar(aes(ymin=md-sd, ymax=md+sd))+
  geom_line()+
  ggpubr::theme_pubr()
+facet_wrap(~type)




# r3 = readRDS('Rdata/phyloP_TCGA_3ss_ALL_EVENTS.rds')
# r4 = readRDS('Rdata/phyloP_TCGA_5ss_ALL_EVENTS.rds')

r3 = readRDS('phyloP_TCGA_3ss_ALL_EVENTS.rds')
r4 = readRDS('phyloP_TCGA_5ss_ALL_EVENTS.rds')

prof = rbind.data.frame(r3,r4)

ttest = ddply(prof, .(ss, ASE_set, pos),
              summarise, m = mean(score), md=median(score), sd=sd(score)) %>%
  ddply( .(ss), function(x) with(x,broom::tidy(t.test(md[ASE_set!='vASE'], md[ASE_set=='vASE'], paired = T))))

ttest2=ddply(prof, .(ss, pois_ess, ASE_set, pos),
             summarise, m = mean(score), md=median(score), sd=sd(score)) %>%
  ddply( .(pois_ess,ss), function(x) with(x,broom::tidy(t.test(md[ASE_set!='vASE'], md[ASE_set=='vASE'], paired = T))))

ttest3=ddply(prof, .(ss, pois_ess, type, ASE_set, pos),
             summarise, m = mean(score), md=median(score), sd=sd(score)) %>%
  ddply( .(pois_ess,type,ss), function(x) with(x,broom::tidy(t.test(md[ASE_set!='vASE'], md[ASE_set=='vASE'], paired = T))))%>%
  arrange(dplyr::desc(p.value))

# dexSIGN=readRDS('Rdata/dexSIGN_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds')
# dexVAR=readRDS('Rdata/dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds')

dexSIGN=readRDS('dexSIGN_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds')
dexVAR=readRDS('dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds')


prof$ASE_set = 'rest'
prof$ASE_set[which(prof$event_id%in%dexVAR$event_id)] = 'vASE'
prof$ASE_set[which(prof$event_id%in%dexSIGN$event_id)] = 'FOXA1_regulated'

# saveRDS(prof, file='Rdata/phylop_3ss_5ss_TCGA_ALL_EVENTS.rds')
saveRDS(prof, file='phylop_3ss_5ss_TCGA_ALL_EVENTS.rds')


# prof=readRDS('Rdata/phylop_3ss_5ss_TCGA_ALL_EVENTS.rds')
prof=readRDS('phylop_3ss_5ss_TCGA_ALL_EVENTS.rds')

options(bitmapType="cairo")


pdf(file="conservation_gam_MDG__all_and_delta_strata_ALL_EVENTS.pdf", h=10, w=10, useDingbats = F)
prof %>% ggplot( aes(pos, score, group=ASE_set, color=ASE_set))+
  geom_smooth(method = 'gam', se = T)+ggpubr::theme_pubr()+facet_wrap(~ss, ncol=2)
dev.off()


