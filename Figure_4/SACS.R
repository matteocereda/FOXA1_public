# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Enrichment of SACS-marked exons ----


library(readxl)
library(tidyr)
library(GenomicRanges)
library(GenomicFeatures)
library(plyr)
library(ggplot2)


dexVAR_gene <- readRDS('Rdata/dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds')
FIG_DIR = '/Volumes/Prule/remove_me_Figures_FOXA1/'
system(paste0('mkdir ', FIG_DIR))

# ∞∞∞ Intersection with SACS-marked exons ----

# SACS ========
scs =list()
for(i in 1:7) {
  scs[[i]] =  read_xlsx('Tables/SACS/41467_2021_20979_MOESM6_ESM.xlsx',sheet = i+1)
  scs[[i]]$sacs = i
  scs[[i]] = scs[[i]][,c(1,ncol(scs[[i]]))]
}
scs = do.call(rbind, scs)

scs = scs %>% separate('ALL SPLICING EVENTS (SACS-MARKED EXONS)'
                       , into=c('chr','up', 'start','end',"dw", 'strand','symbol','geneID')
                       , sep='\\:', convert=T) 


ex = dexVAR_gene
es = ex %>% subset(event_type=='ES') %>%
  dplyr::select(event_id, event_chr, event_coordinates, strand, gene_id,  delta.mean, FOXA1_sign,pois_ess) # ,exon_length, intron_length
es$ASE_set = 'vASE'
es$ASE_set[es$FOXA1_sign] = "FOXA1_regulated"
es$FOXA1_sign = NULL
es = es %>% separate(event_coordinates, into=c("u.start","u.end",'start','end','d.start','d.end'), convert=T)

es$e.length  = es$end-es$start
es$u.length  = with(es, u.end   - u.start )
es$d.length  = with(es, d.end   - d.start )
es$iu.length = with(es, start   - u.end  )
es$id.length = with(es, d.start - end )



q = with(es, GRanges(seqnames = event_chr, IRanges(start = start, end=end)))
t = with(scs, GRanges(seqnames = chr, IRanges(start = start, end=end)))
o = findOverlaps(q, t)
overlaps <- pintersect(q[queryHits(o)], t[subjectHits(o)])
percentOverlap <- width(overlaps) / width(t[subjectHits(o)])
o <- o[percentOverlap > 0.9]

es$sacs = 'REST'
es$sacs[unique(queryHits(o))] = 'SACS'
es$sacs = factor(es$sacs, levels=c('REST','SACS'))

scs = ddply(scs, .(chr, start, end), summarise, sacs=paste(unique(sacs),collapse=';') )
t = with(scs, GRanges(seqnames = chr, IRanges(start = start, end=end)))
o = findOverlaps(q, t)
overlaps <- pintersect(q[queryHits(o)], t[subjectHits(o)])
percentOverlap <- width(overlaps) / width(t[subjectHits(o)])
o <- o[percentOverlap > 0.9]

es$sacs_type = NA
es$sacs_type[queryHits(o)] = scs$sacs[subjectHits(o)]

es$sacs = factor(es$sacs, levels=c('REST','SACS'))

es = es %>% as_tibble()


# ∞∞∞ Enrichment of SACS ----

to_fisher = table(es$sacs,es$ASE_set)

fisher.test(to_fisher)

df_sacs=as.data.frame(table(es$sacs_type,es$ASE_set))
df_sacs$prop = NA
df_sacs$prop[which(df_sacs$Var2=='FOXA1_regulated')]=df_sacs$Freq[which(df_sacs$Var2=='FOXA1_regulated')]/table(es$ASE_set)[1]
df_sacs$prop[which(df_sacs$Var2=='vASE')]=df_sacs$Freq[which(df_sacs$Var2=='vASE')]/table(es$ASE_set)[2]

pdf(file=paste0(FIG_DIR, "SACS_proportion.pdf"), h=5, w=5, useDingbats = F)
ggplot(df_sacs,aes(x=Var2,y=prop,fill=Var1))+geom_col()+theme_bw()
dev.off()
