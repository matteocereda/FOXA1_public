# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Gene ontology on genes harbouring FOXA1-regulated events from RNAseq data of VCaP and PC3 cells ----


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
  ego <- enricher(unique(ids$ENTREZID)  #, universe = unique(bg$ENTREZID)
                  , TERM2GENE = TERM2GENE
                  , pAdjustMethod = pAdjustMethod
                  , pvalueCutoff=pvalueCutoff
                  , qvalueCutoff=qvalueCutoff
                  , minGSSize = minGSSize )
  ego = format_enrichr(ego, org=OrgDb, key=typeFormat)
  list('id'=ids, 'go'=ego)
}

repo = '/Volumes/Prule/repo/FOXA1_public/'
FIG_DIR = '~/Downloads/FOXA1/Figure2/'
setwd(repo)

vcap = readRDS('Rdata/Whippet_VCAP_filter_0_2.rds')
vcap = subset(vcap,!Type%in%c('TE','TS'))
vcap$STRINGENT = F
vcap$STRINGENT[which(abs(vcap$DeltaPsi)>0.05 & vcap$Probability>0.9 & !vcap$Complexity=='K0')] = T
vcap$Type[which(vcap$Type=='AF')] = 'AA'
vcap$Type[which(vcap$Type=='AL')] = 'AD'
vcap$ensembl = sapply(strsplit(vcap$Gene, '\\.'),'[[',1)

pc3 = readRDS('Rdata/Whippet_PC3_nextseq_novaseq_filter_0_2.rds')
pc3 = subset(pc3,!Type%in%c('TE','TS'))
pc3$STRINGENT = F
pc3$STRINGENT[which(abs(pc3$DeltaPsi)>0.05 & pc3$Probability>0.9 & !pc3$Complexity=='K0')] = T
pc3$Type[which(pc3$Type=='AF')] = 'AA'
pc3$Type[which(pc3$Type=='AL')] = 'AD'
pc3$ensembl = sapply(strsplit(pc3$Gene, '\\.'),'[[',1)

vcap$cell='VCaP'
pc3$cell='PC3'

cls = c('cell','ensembl','Type','STRINGENT')
ase = rbind.data.frame(pc3[,cls],vcap[,cls])
# ase_p = subset(ase, SENSITIVE); ase_p = ase_p[,1:3]; ase_p=unique(ase_p); ase_p$set='permissive'
ase_s = subset(ase, STRINGENT); ase_s = ase_s[,1:3]; ase_s=unique(ase_s); ase_s$set='stringent'

ase = ase_s
# ase = rbind.data.frame(ase_p, ase_s)



# % ORA with 186 KEGG single event types =======

pll = load('Rdata/KEGG_list.148.SPLICING_RELATED.Rdata')

library(msigdbr)
k = msigdbr(species = "Homo sapiens", category = "C2", subcategory = 'KEGG')

m_t2g <-  k%>%  dplyr::select(gs_name, entrez_gene)


p = list()
for(s in c("stringent") ) for(i in c('VCaP', 'PC3')) for(j in c("CE", "AD", "AA", "RI")){
  p[[paste0(s, '_',i, '_',j)]] = try(get_enrich(subset(ase, set==s & cell==i & Type==j), m_t2g))
}

df = rbind.data.frame()
for(i in names(p)){  df = rbind.data.frame(df ,cbind.data.frame(cell=i, p[[i]]$go@result) )}
df = df %>% separate(cell,into=c("set", 'cell', 'type'))

df = ddply(df, .(set, cell,type), mutate, rank= rank(p.adjust))

top_10_2 = unique(subset(df, rank<=10 & set=='stringent')$Description)

to_plot2 = subset(df, Description %in% top_10_2 & set=='stringent')

tmp2 = ddply(to_plot2,.(ID),summarise,mean_padj=mean(p.adjust),mean_GR=mean(parse_ratio(GeneRatio)))
tmp2 = tmp2[order(tmp2$mean_GR,decreasing = T),]
to_plot2$ID = factor(to_plot2$ID, levels = rev(tmp2$ID))

pdf(file = paste0(FIG_DIR, 'ORA_single_event_types_Cells_KEGG_186.pdf'), height=unit(8,'cm'), width=unit(22,'cm') ,useDingbats = F)
to_plot2=subset(to_plot2,p.adjust<0.25)
ggplot( to_plot2,
           aes(x= parse_ratio(GeneRatio), y= ID)) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=factor(p.adjust<=.1), fill=p.adjust, size = Count, shape=cell)) +
  geom_text(aes(label=rank)) +
  scale_shape_manual(values=c(24,21))+
  scale_color_manual(values=c('transparent','black'))+
  scale_fill_viridis_c(option = 'D',guide=guide_colorbar(reverse=T, draw.llim = T), direction = -1)+#,
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + xlab("Gene Ratio") +ylab(NULL)+
  facet_wrap(~type, nrow=1, scales='free_x')
dev.off()



# % ORA with 186 KEGG all event types together =======

pll = load('Rdata/KEGG_list.148.SPLICING_RELATED.Rdata')

library(msigdbr)
k = msigdbr(species = "Homo sapiens", category = "C2", subcategory = 'KEGG')

m_t2g <-  k%>%  dplyr::select(gs_name, entrez_gene)

p = list()
for(s in c("stringent") ) for(i in c('VCaP', 'PC3')){
  p[[paste0(s, '_',i)]] = try(get_enrich(subset(ase, set==s & cell==i), m_t2g))
}

df = rbind.data.frame()
for(i in names(p)){  df = rbind.data.frame(df ,cbind.data.frame(cell=i, p[[i]]$go@result) )}
df = df %>% separate(cell,into=c("set", 'cell'))

df = ddply(df, .(set, cell), mutate, rank= rank(p.adjust))

top_10_2 = unique(subset(df, rank<=10 & set=='stringent')$Description)

to_plot2 = subset(df, Description %in% top_10_2 & set=='stringent')

tmp2 = ddply(to_plot2,.(ID),summarise,mean_padj=mean(p.adjust),mean_GR=mean(parse_ratio(GeneRatio)))
tmp2 = tmp2[order(tmp2$mean_GR,decreasing = T),]
to_plot2$ID = factor(to_plot2$ID, levels = rev(tmp2$ID))

pdf(file=paste0(FIG_DIR, 'ORA_all_event_types_Cells_KEGG_186.pdf'), height = unit(8,'cm'), width = unit(22, 'cm'), useDingbats = F)
to_plot2=subset(to_plot2,p.adjust<0.1)
ggplot( to_plot2,
           aes(x= parse_ratio(GeneRatio), y= ID)) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=factor(p.adjust<=.25), fill=p.adjust, size = Count, shape=cell)) +
  geom_text(aes(label=rank)) +
  scale_shape_manual(values=c(24,21))+
  scale_color_manual(values=c('transparent','black'))+
  scale_fill_viridis_c(option = 'D',guide=guide_colorbar(reverse=T, draw.llim = T), direction = -1)+#,
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + xlab("Gene Ratio") +ylab(NULL)+
  facet_wrap(~set, ncol=1, scales='free_x')
dev.off()


pdf(file=paste0(FIG_DIR, 'ORA_all_event_types_Cells_KEGG_186__richFactor.pdf'), height = unit(4,'cm'), width = unit(8, 'cm'), useDingbats = F)
to_plot2=subset(to_plot2,p.adjust<0.1)
ggplot( to_plot2,
        aes(x= richFactor, y= ID)) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=factor(p.adjust<=.25), fill=p.adjust, size = Count, shape=cell)) +
  geom_text(aes(label=rank)) +
  scale_shape_manual(values=c(24,21))+
  scale_color_manual(values=c('transparent','black'))+
  scale_fill_viridis_c(option = 'D',guide=guide_colorbar(reverse=T, draw.llim = T), direction = -1)+#,
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + xlab("richFactor") +ylab(NULL)+
  facet_wrap(~set, ncol=1, scales='free_x')
dev.off()

