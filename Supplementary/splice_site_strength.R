# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Splice site strength analysis ----


library(tidyverse)
library(magrittr)
library(org.Hs.eg.db)
library(clusterProfiler)
library(ggpubr)

ex = readRDS('Rdata/variable_exons_with_splice_sites_strength_annotation.rds')

es = ex %>% subset(event_type=='ES')
es=es[, c("gene_id", "event_id", "delta.mean", "overall.mean", "overall.sd", "exon_type", "FOXA1_sign",
          "pois_ess", 
          "upstream_5ss_donor_motif_STRENGTH", "upstream_3ss_acceptor_motif_STRENGTH",
          "downstream_5ss_donor_motif_STRENGTH", "downstream_3ss_acceptor_motif_STRENGTH"
          ,"exon_length", "intron_length")]

es=reshape2::melt(es, measure.vars=c('upstream_5ss_donor_motif_STRENGTH', 'upstream_3ss_acceptor_motif_STRENGTH',
                                'downstream_5ss_donor_motif_STRENGTH', 'downstream_3ss_acceptor_motif_STRENGTH'))
es$variable = sub('_motif_STRENGTH','', es$variable)
es$variable = sub('_donor|_acceptor','', es$variable)
es$variable = factor(es$variable, levels= c('upstream_5ss', 'upstream_3ss','downstream_5ss','downstream_3ss' ) )
es$ASE_set = 'vASE'
es$ASE_set[es$FOXA1_sign] = "FOXA1_regulated"
es$ASE_set=factor(es$ASE_set, levels=c('FOXA1_regulated', 'vASE'))


p1=es %>% ggboxplot( x = "variable", y = "value", fill = "ASE_set",
                     notch = T )+ scale_fill_manual(values=c('green','grey'))+
  stat_compare_means( aes(group = ASE_set), label = "p.format", method = 'wilcox.test'
  )
p2=  es %>% ggboxplot( x = "variable", y = "value", fill = "ASE_set",
                       notch = T )+ scale_fill_manual(values=c('green','grey'))+
  scale_y_sqrt(breaks=(c(1,5,10,15)),labels=as.character(c(1,5,10,15)))

pdf(file='splice_site_strength.pdf', height = unit(8,'cm'), width = unit(5, 'cm'), useDingbats = F)
gridExtra::grid.arrange(p1,p2)
dev.off()
