# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Exon-intron length ----


es = readRDS('Rdata/dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds')


es$u_intron_start = as.numeric(sapply(strsplit(as.character(es$event_coordinates),':'),`[`,2))
es$u_intron_end = as.numeric(sapply(strsplit(as.character(es$event_coordinates),':'),`[`,3))

es$d_intron_start = as.numeric(sapply(strsplit(as.character(es$event_coordinates),':'),`[`,4))
es$d_intron_end = as.numeric(sapply(strsplit(as.character(es$event_coordinates),':'),`[`,5))

es$exon_start = as.numeric(sapply(strsplit(as.character(es$event_coordinates),':'),`[`,3))
es$exon_end = as.numeric(sapply(strsplit(as.character(es$event_coordinates),':'),`[`,4))


es$e.length = es$exon_end-es$exon_start
es$iu.length = es$u_intron_end-es$u_intron_start
es$id.length = es$d_intron_end-es$d_intron_start

es$ASE_set = 'vASE'
es$ASE_set[which(es$FOXA1_sign)] = 'FOXA1_regulated'


to_plot=reshape2::melt(es[,63:66])
to_plot$variable=factor(to_plot$variable,levels = c('iu.length','e.length','id.length'))

pdf(file="exon_intron_length_boxplots.pdf", h=8, w=8, useDingbats = F)
ggplot(to_plot,aes(x=ASE_set,fill=ASE_set,y= log10(value)))+geom_boxplot(notch = T)+stat_compare_means()+theme_bw()+facet_wrap(~variable)
dev.off()
