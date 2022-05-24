# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Integrated analysis of TF binding in prostate cancer ----

# • Download selected datasets from ChIP-Atlas (See Supplementary Table 2)  ----
# One folder for each TF and cell line

# • Merge and sort datasets ----

options(stringsAsFactors = F)

library(tidyr)
library(plyr)
library(ggplot2)



# FOXA1  LNCaP ============
# bash ============
# in LNCaP_FOXA1/

# cat * > FOXA1_selected_lncap_hg19_chip_atlas_05.bed
# sort -k1,1 -k2,2n FOXA1_selected_lncap_hg19_chip_atlas_05.bed > FOXA1_selected_lncap_hg19_chip_atlas_05_sorted.bed
# mergeBed -i FOXA1_selected_lncap_hg19_chip_atlas_05_sorted.bed > FOXA1_selected_lncap_hg19_chip_atlas_05_sorted_merged.bed

# with atac =============

setwd('LNCaP_FOXA1/')

system("a=FOXA1_selected_lncap_hg19_chip_atlas_05_sorted_merged.bed;
b=Tables/chip_atac/ATAC_hg19_PRAD.bed;
intersectBed -a $a -b $b -u > FOXA1_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed")


# Promoters LNCAP intersection =================

inp="Tables/chip_atac/lncap_promoters_Rhie_2000_2000_in_ATAC.bed"
out="FOXA1_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed"


system(sprintf("a=FOXA1_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))


# Analysis only enhancers in 1 Mb ===============

inp="Tables/chip_atac/LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"
out="Intersection_FOXA1_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"


system(sprintf("a=FOXA1_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))




# FOXA1  VCaP ============
# bash ============
# in VCaP_FOXA1/

# cat * > FOXA1_selected_vcap_hg19_chip_atlas_05.bed
# sort -k1,1 -k2,2n FOXA1_selected_vcap_hg19_chip_atlas_05.bed > FOXA1_selected_vcap_hg19_chip_atlas_05_sorted.bed
# mergeBed -i FOXA1_selected_vcap_hg19_chip_atlas_05_sorted.bed > FOXA1_selected_vcap_hg19_chip_atlas_05_sorted_merged.bed


# with atac =============

setwd('VCaP_FOXA1/')

system("a=FOXA1_selected_vcap_hg19_chip_atlas_05_sorted_merged.bed;
b=Tables/chip_atac/ATAC_hg19_PRAD.bed;
intersectBed -a $a -b $b -u > FOXA1_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed")


# Promoters LNCAP intersection =================

inp="Tables/chip_atac/lncap_promoters_Rhie_2000_2000_in_ATAC.bed"
out="FOXA1_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed"


system(sprintf("a=FOXA1_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))


# Analysis only enhancers in 1 Mb ===============

inp="Tables/chip_atac/VCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"
out="Intersection_FOXA1_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_VCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"


system(sprintf("a=FOXA1_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))





# ERG  LNCaP ============
# bash ============
# in LNCaP_ERG/

# cat * > ERG_selected_lncap_hg19_chip_atlas_05.bed
# sort -k1,1 -k2,2n ERG_selected_lncap_hg19_chip_atlas_05.bed > ERG_selected_lncap_hg19_chip_atlas_05_sorted.bed
# mergeBed -i ERG_selected_lncap_hg19_chip_atlas_05_sorted.bed > ERG_selected_lncap_hg19_chip_atlas_05_sorted_merged.bed

# with atac =============

setwd('LNCaP_ERG/')

system("a=ERG_selected_lncap_hg19_chip_atlas_05_sorted_merged.bed;
b=Tables/chip_atac/ATAC_hg19_PRAD.bed;
intersectBed -a $a -b $b -u > ERG_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed")


# Promoters LNCAP intersection =================

inp="Tables/chip_atac/lncap_promoters_Rhie_2000_2000_in_ATAC.bed"
out="ERG_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed"


system(sprintf("a=ERG_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))


# Analysis only enhancers in 1 Mb ===============

inp="Tables/chip_atac/LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"
out="Intersection_ERG_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"


system(sprintf("a=ERG_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))




# ERG  VCaP ============
# bash ============
# in VCaP_ERG/

# cat * > ERG_selected_vcap_hg19_chip_atlas_05.bed
# sort -k1,1 -k2,2n ERG_selected_vcap_hg19_chip_atlas_05.bed > ERG_selected_vcap_hg19_chip_atlas_05_sorted.bed
# mergeBed -i ERG_selected_vcap_hg19_chip_atlas_05_sorted.bed > ERG_selected_vcap_hg19_chip_atlas_05_sorted_merged.bed


# with atac =============

setwd('VCaP_ERG/')

system("a=ERG_selected_vcap_hg19_chip_atlas_05_sorted_merged.bed;
b=Tables/chip_atac/ATAC_hg19_PRAD.bed;
intersectBed -a $a -b $b -u > ERG_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed")


# Promoters LNCAP intersection =================

inp="Tables/chip_atac/lncap_promoters_Rhie_2000_2000_in_ATAC.bed"
out="ERG_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed"


system(sprintf("a=ERG_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))


# Analysis only enhancers in 1 Mb ===============

inp="Tables/chip_atac/VCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"
out="Intersection_ERG_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_VCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"


system(sprintf("a=ERG_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))



# AR  LNCaP ============
# bash ============
# in LNCaP_AR/

# cat * > AR_selected_lncap_hg19_chip_atlas_05.bed
# sort -k1,1 -k2,2n AR_selected_lncap_hg19_chip_atlas_05.bed > AR_selected_lncap_hg19_chip_atlas_05_sorted.bed
# mergeBed -i AR_selected_lncap_hg19_chip_atlas_05_sorted.bed > AR_selected_lncap_hg19_chip_atlas_05_sorted_merged.bed

# with atac =============

setwd('LNCaP_AR/')

system("a=AR_selected_lncap_hg19_chip_atlas_05_sorted_merged.bed;
b=Tables/chip_atac/ATAC_hg19_PRAD.bed;
intersectBed -a $a -b $b -u > AR_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed")


# Promoters LNCAP intersection =================

inp="Tables/chip_atac/lncap_promoters_Rhie_2000_2000_in_ATAC.bed"
out="AR_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed"


system(sprintf("a=AR_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))


# Analysis only enhancers in 1 Mb ===============

inp="Tables/chip_atac/LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"
out="Intersection_AR_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"


system(sprintf("a=AR_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))



# AR  VCaP ============
# bash ============
# in VCaP_AR/

# cat *.bed > AR_selected_vcap_hg19_chip_atlas_05.bed
# sort -k1,1 -k2,2n AR_selected_vcap_hg19_chip_atlas_05.bed > AR_selected_vcap_hg19_chip_atlas_05_sorted.bed
# mergeBed -i AR_selected_vcap_hg19_chip_atlas_05_sorted.bed > AR_selected_vcap_hg19_chip_atlas_05_sorted_merged.bed


# with atac =============

setwd('VCaP_AR/')

system("a=AR_selected_vcap_hg19_chip_atlas_05_sorted_merged.bed;
b=Tables/chip_atac/ATAC_hg19_PRAD.bed;
intersectBed -a $a -b $b -u > AR_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed")


# Promoters LNCAP intersection =================

inp="Tables/chip_atac/lncap_promoters_Rhie_2000_2000_in_ATAC.bed"
out="AR_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed"


system(sprintf("a=AR_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))


# Analysis only enhancers in 1 Mb ===============

inp="Tables/chip_atac/VCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"
out="Intersection_AR_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_VCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"


system(sprintf("a=AR_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))




# MYC  LNCaP ============
# bash ============
# in LNCaP_MYC/

# cat *.bed > MYC_selected_lncap_hg19_chip_atlas_05.bed
# sort -k1,1 -k2,2n MYC_selected_lncap_hg19_chip_atlas_05.bed > MYC_selected_lncap_hg19_chip_atlas_05_sorted.bed
# mergeBed -i MYC_selected_lncap_hg19_chip_atlas_05_sorted.bed > MYC_selected_lncap_hg19_chip_atlas_05_sorted_merged.bed

# with atac =============

setwd('LNCaP_MYC/')

system("a=MYC_selected_lncap_hg19_chip_atlas_05_sorted_merged.bed;
b=Tables/chip_atac/ATAC_hg19_PRAD.bed;
intersectBed -a $a -b $b -u > MYC_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed")


# Promoters LNCAP intersection =================

inp="Tables/chip_atac/lncap_promoters_Rhie_2000_2000_in_ATAC.bed"
out="MYC_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed"


system(sprintf("a=MYC_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))


# Analysis only enhancers in 1 Mb ===============

inp="Tables/chip_atac/LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"
out="Intersection_MYC_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed"


system(sprintf("a=MYC_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC.bed;
       b=%s;
       intersectBed -a $a -b $b -wo > %s", inp, out))




# • Integrate analysis ----

library(clusterProfiler)
library(DOSE)
library(igraph)
library(org.Hs.eg.db)
library(tidyverse)
library(magrittr)


format_enrichr = function(x,  org=org.Hs.eg.db, key="ENTREZID"){
  y <- setReadable(x, OrgDb = org, keyType=key)
  y <- mutate(y, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  y <- mutate(y, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  y
}

format_enrichKEGG = function(x, kegg, org=org.Hs.eg.db, key="UNIPROT"){
  y <- setReadable(x, OrgDb = org, keyType=key)
  y <- mutate(y, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  y <- mutate(y, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  y@result$kegg_id=gsub('hsa0','',y$ID)
  y@result = left_join(y@result, kegg, by='kegg_id')
  y
}

format_enrichr_msig = function(x, kegg, org=org.Hs.eg.db, key="ENTREZID"){
  y <- setReadable(x, OrgDb = org, keyType=key)
  y <- clusterProfiler::mutate(y, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
  y <- clusterProfiler::mutate(y, FoldEnrichment = parse_ratio(GeneRatio) / parse_ratio(BgRatio))
  y@result$msigDB_158=y@result$Description
  y@result$msigDB_158[y@result$msigDB_158=='SRG']='KEGG_SPLICEOSOME'
  y@result = left_join(y@result, kegg, by='msigDB_158')
  y
}

# LOADING DATA ===========

kegg = read.csv('Tables/KEGG.csv')


# Choose a TF and cell line on which the analysis will be performed: ----

# FOXA1 ----

# LNCaP %%%%%%%%%%
prom=read.delim2("LNCaP_FOXA1/FOXA1_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed",header = F)
prom$which = 'promoter'
prom = prom[,c('which','V10')]
colnames(prom)[2] = 'gene_name'
prom = unique(prom)
enh=read.delim2("LNCaP_FOXA1/Intersection_FOXA1_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed",header = F)
enh$which = 'enhancer'
enh = enh[,c('which','V7')]
colnames(enh)[2] = 'gene_name'
enh = unique(enh)

# VCaP %%%%%%%%%%
prom=read.delim2("VCaP_FOXA1/FOXA1_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed",header = F)
prom$which = 'promoter'
prom = prom[,c('which','V10')]
colnames(prom)[2] = 'gene_name'
prom = unique(prom)
enh=read.delim2("VCaP_FOXA1/Intersection_FOXA1_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_VCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed",header = F)
enh$which = 'enhancer'
enh = enh[,c('which','V7')]
colnames(enh)[2] = 'gene_name'
enh = unique(enh)

# MYC ----

# LNCaP %%%%%%%%%%
prom=read.delim2("LNCaP_MYC/MYC_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed",header = F)
prom$which = 'promoter'
prom = prom[,c('which','V10')]
colnames(prom)[2] = 'gene_name'
prom = unique(prom)
enh=read.delim2("LNCaP_MYC/Intersection_MYC_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed",header = F)
enh$which = 'enhancer'
enh = enh[,c('which','V7')]
colnames(enh)[2] = 'gene_name'
enh = unique(enh)

# VCaP %%%%%%%%%%
# NO DATA


# AR ----

# LNCaP %%%%%%%%%%
prom=read.delim2("/LNCaP_AR/AR_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed",header = F)
prom$which = 'promoter'
prom = prom[,c('which','V10')]
colnames(prom)[2] = 'gene_name'
prom = unique(prom)
enh=read.delim2("LNCaP_AR/Intersection_AR_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed",header = F)
enh$which = 'enhancer'
enh = enh[,c('which','V7')]
colnames(enh)[2] = 'gene_name'
enh = unique(enh)

# VCaP %%%%%%%%%%
prom=read.delim2("/VCaP_AR/AR_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed",header = F)
prom$which = 'promoter'
prom = prom[,c('which','V10')]
colnames(prom)[2] = 'gene_name'
prom = unique(prom)
enh=read.delim2("VCaP_AR/Intersection_AR_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_VCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed",header = F)
enh$which = 'enhancer'
enh = enh[,c('which','V7')]
colnames(enh)[2] = 'gene_name'
enh = unique(enh)


# ERG ----

# LNCaP %%%%%%%%%%
prom=read.delim2("LNCaP_ERG/ERG_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed",header = F)
prom$which = 'promoter'
prom = prom[,c('which','V10')]
colnames(prom)[2] = 'gene_name'
prom = unique(prom)
enh=read.delim2("LNCaP_ERG/Intersection_ERG_selected_lncap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_LNCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed",header = F)
enh$which = 'enhancer'
enh = enh[,c('which','V7')]
colnames(enh)[2] = 'gene_name'
enh = unique(enh)

# VCaP %%%%%%%%%%
prom=read.delim2("VCaP_ERG/ERG_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_promoters_2000_2000_in_ATAC.bed",header = F)
prom$which = 'promoter'
prom = prom[,c('which','V10')]
colnames(prom)[2] = 'gene_name'
prom = unique(prom)
enh=read.delim2("VCaP_ERG/Intersection_ERG_selected_vcap_hg19_chip_atlas_05_sorted_merged_in_ATAC_Ramanand_VCAP_enhancers_associated_to_genes_in_ATAC_not_on_promoters_of_the_same_gene_only_active_genes_in_1Mb.bed",header = F)
enh$which = 'enhancer'
enh = enh[,c('which','V7')]
colnames(enh)[2] = 'gene_name'
enh = unique(enh)



df = rbind.data.frame(cbind.data.frame('type'='P', SYMBOL=unique(prom$gene_name))
                      ,cbind.data.frame('type'='E', SYMBOL=unique(enh$gene_name)) )


ids <- bitr(unique(df$SYMBOL), fromType="SYMBOL", toType=c("UNIPROT", "ENSEMBL", 'ENTREZID','ALIAS'), OrgDb="org.Hs.eg.db")
df = left_join(df, ids, by='SYMBOL')



# ORA with 186 KEGG ==================

library(msigdbr)
k = msigdbr(species = "Homo sapiens", category = "C2", subcategory = 'KEGG')

m_t2g <-  k%>%  dplyr::select(gs_name, entrez_gene)


eprom <- enricher(unique(subset(df, type=='P')$ENTREZID)
                  , TERM2GENE = m_t2g
                  , pAdjustMethod = "BH", pvalueCutoff=1, qvalueCutoff=1 )

eprom = format_enrichr_msig(eprom, kegg )

eenh <-  enricher(unique(subset(df, type=='E')$ENTREZID)
                  , TERM2GENE = m_t2g
                  , pAdjustMethod = "BH", pvalueCutoff=1, qvalueCutoff=1 )

eenh = format_enrichr_msig(eenh, kegg )

# manual 
p =eprom@result
p$re = 'P'
p$Description=fct_reorder(p$Description, as.numeric(sub("/\\d+", "", p$GeneRatio)), .desc = F)
p$rank = rank(p$p.adjust)

e = eenh@result
e$re = 'E'
e$Description=factor(e$Description, levels = levels(p$Description))
e$rank = rank(e$p.adjust)

tmp = rbind.data.frame(p,e)


tmp$type='KEGG_186'

# Choose ===
# LNCaP
tmp$cell_line='LNCaP'
tmp_KEGG_186_LNCaP=tmp
tmp_KEGG_186_LNCaP$region='Promoter'
tmp_KEGG_186_LNCaP$region[which(tmp_KEGG_186_LNCaP$re=='E')]='Enhancer'
# write.csv(tmp_KEGG_186_LNCaP[,c('ID','GeneRatio','BgRatio','pvalue','p.adjust','geneID','region')],'ORA_LNCAP_enhancer_promoter_FOXA1.csv')
# write.csv(tmp_KEGG_186_LNCaP[,c('ID','GeneRatio','BgRatio','pvalue','p.adjust','geneID','region')],'ORA_LNCAP_enhancer_promoter_MYC.csv')
# write.csv(tmp_KEGG_186_LNCaP[,c('ID','GeneRatio','BgRatio','pvalue','p.adjust','geneID','region')],'ORA_LNCAP_enhancer_promoter_AR.csv')
write.csv(tmp_KEGG_186_LNCaP[,c('ID','GeneRatio','BgRatio','pvalue','p.adjust','geneID','region')],'ORA_LNCAP_enhancer_promoter_ERG.csv')
# VCaP
tmp$cell_line='VCaP'
tmp_KEGG_186_VCaP=tmp
tmp_KEGG_186_VCaP$region='Promoter'
tmp_KEGG_186_VCaP$region[which(tmp_KEGG_186_VCaP$re=='E')]='Enhancer'
# write.csv(tmp_KEGG_186_VCaP[,c('ID','GeneRatio','BgRatio','pvalue','p.adjust','geneID','region')],'ORA_VCAP_enhancer_promoter_FOXA1.csv')
# write.csv(tmp_KEGG_186_VCaP[,c('ID','GeneRatio','BgRatio','pvalue','p.adjust','geneID','region')],'ORA_VCAP_enhancer_promoter_AR.csv')
write.csv(tmp_KEGG_186_VCaP[,c('ID','GeneRatio','BgRatio','pvalue','p.adjust','geneID','region')],'ORA_VCAP_enhancer_promoter_ERG.csv')



# pdf(file='MYC_LNCAP_RHIE_promoters_ramanand_enhancers_ORA_in_ATAC_not_on_promoters_of_their_gene_only_active_genes_limited_exp_hallmark_1Mb__KEGG_186.pdf', height = unit(4,'cm'), width = unit(10, 'cm'), useDingbats = F)
# pdf(file='AR_LNCAP_RHIE_promoters_ramanand_enhancers_ORA_in_ATAC_not_on_promoters_of_their_gene_only_active_genes_limited_exp_hallmark_1Mb__KEGG_186.pdf', height = unit(4,'cm'), width = unit(10, 'cm'), useDingbats = F)
# pdf(file='AR_VCAP_RHIE_promoters_ramanand_enhancers_ORA_in_ATAC_not_on_promoters_of_their_gene_only_active_genes_limited_exp_hallmark_1Mb__KEGG_186.pdf', height = unit(4,'cm'), width = unit(10, 'cm'), useDingbats = F)
# pdf(file='ERG_LNCAP_RHIE_promoters_ramanand_enhancers_ORA_in_ATAC_not_on_promoters_of_their_gene_only_active_genes_limited_exp_hallmark_1Mb__KEGG_186.pdf', height = unit(4,'cm'), width = unit(10, 'cm'), useDingbats = F)
pdf(file='ERG_VCAP_RHIE_promoters_ramanand_enhancers_ORA_in_ATAC_not_on_promoters_of_their_gene_only_active_genes_limited_exp_hallmark_1Mb__KEGG_186.pdf', height = unit(4,'cm'), width = unit(10, 'cm'), useDingbats = F)
# pdf(file='FOXA1_VCAP_RHIE_promoters_ramanand_enhancers_ORA_in_ATAC_not_on_promoters_of_their_gene_only_active_genes_limited_exp_hallmark_1Mb__KEGG_186.pdf', height = unit(4,'cm'), width = unit(10, 'cm'), useDingbats = F)
# pdf(file='FOXA1_LNCAP_RHIE_promoters_ramanand_enhancers_ORA_in_ATAC_not_on_promoters_of_their_gene_only_active_genes_limited_exp_hallmark_1Mb__KEGG_186.pdf', height = unit(4,'cm'), width = unit(10, 'cm'), useDingbats = F)


tmp=subset(tmp,p.adjust<0.1)
tmp=ddply(tmp,.(ID),mutate,mean_enr=mean(parse_ratio(GeneRatio)))
tmp=tmp[order(tmp$mean_enr,decreasing = F),]
tmp$ID=factor(tmp$ID,levels=unique(tmp$ID))
ggplot(subset(tmp, ID %in% subset(tmp, rank<=10)$ID)
       ,aes(
         x=  parse_ratio(GeneRatio)
         ,y = ID)
) +
  geom_segment(aes(xend=0, yend = ID)) +
  geom_point(aes(color=factor(p.adjust<=0.1), fill=p.adjust, size = Count, shape=re)) +
  scale_shape_manual(values=c(24,21))+
  scale_color_manual(values=c('transparent','black'))+
  scale_fill_viridis_c(option = 'D',guide=guide_colorbar(reverse=T, draw.llim = T), direction = -1,
                       # values  =scales::rescale(c(0,0.005,0.01,0.02,0.03,0.04)), breaks=c(0,0.005,0.01,0.02,0.03,0.04))+
                       breaks=c(0,0.001,0.01,0.02,0.05,0.08,0.09,0.1),limits = c(0, 0.1))+
  scale_size_continuous(range=c(2, 10),breaks = c(22,30,40,50,60,70,80,90)) +
  theme_minimal() +
  xlab("Gene Ratio") +
  ylab(NULL) +
  theme(legend.position="right", legend.direction="horizontal")


dev.off()





# MERGE RESULTS ----

files = list.files(pattern = 'ORA_',full.names = T)

for (i in 1:length(files)) {
  
  ego = read.csv(files[[i]])
  cell_line=strsplit(files[[i]],split = '_')[[1]][4]
  TF = strsplit(files[[i]],split = '_')[[1]][7]
  TF = strsplit(TF,split = '[.]')[[1]][1]
  ego$TF = TF
  ego$cell_line = cell_line
  
  if (i==1) {
    ego_df = ego
  } else {
    ego_df = rbind(ego_df,ego)
  }
  
}

subset(ego_df,ID=='KEGG_SPLICEOSOME' & cell_line=='LNCAP' & region=='Promoter')[,c(2:6,8:10)]
subset(ego_df,ID=='KEGG_SPLICEOSOME' & cell_line=='LNCAP' & region=='Enhancer')[,c(2:6,8:10)]

subset(ego_df,ID=='KEGG_SPLICEOSOME' & cell_line=='VCAP' & region=='Promoter')[,c(2:6,8:10)]
subset(ego_df,ID=='KEGG_SPLICEOSOME' & cell_line=='VCAP' & region=='Enhancer')[,c(2:6,8:10)]

saveRDS(ego_df,'ORA_df_all_TFs.rds')


