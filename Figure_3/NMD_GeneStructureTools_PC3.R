# Load packages
library(GeneStructureTools)
library(GenomicRanges)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg19)
library(Gviz)
library(rtracklayer)


# list of files in the whippet directory
whippet_file_directory <- "DIR_WITH_WHIPPET_DATA"

# read in files as a whippetDataSet
wds <- readWhippetDataSet(whippet_file_directory)

tmp = slot(wds, "coordinates")
tmp = keepStandardChromosomes(tmp, pruning.mode = 'coarse')
slot(wds, "coordinates") = tmp

tmp = wds@diffSplicingResults
tmp$chromosome = sapply(strsplit(tmp$coord,':'),`[`,1)
chromosomes = c('chr1','chr10','chr11','chr12' ,'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrX','chrY')
tmp = subset(tmp,chromosome%in%chromosomes)
tmp$chromosome = NULL
wds@diffSplicingResults = tmp


message('[*] Reading gtf ...')

# read in gtf annotation
gtf <- rtracklayer::import("gencode.v28lift37.annotation.gtf")
chromosomes = c('chr1','chr10','chr11','chr12' ,'chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr2','chr20','chr21','chr22','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chrX','chrY')
gtf = gtf[seqnames(gtf)%in%chromosomes]
gtf2=as.data.frame(gtf)
gtf2$seqnames=as.character(gtf2$seqnames)
gtf2$seqnames=factor(gtf2$seqnames)
gtf2=makeGRangesFromDataFrame(gtf2,keep.extra.columns=T)
exons <- gtf2[gtf2$type=="exon"]
transcripts <- gtf2[gtf2$type=="transcript"]

# add first/last annotation (speeds up later steps)
if(!("first_last" %in% colnames(mcols(exons)))){
  t <- as.data.frame(table(exons$transcript_id))
  exons$first_last <- NA
  exons$first_last[exons$exon_number == 1] <- "first"
  exons$first_last[exons$exon_number == 
                     t$Freq[match(exons$transcript_id, t$Var1)]] <- "last"
}

# specify the BSGenome annotation
g <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19


message('[*] Start NMD analysis ...')

delta = readRDS('Rdata/Whippet_PC3_nextseq_novaseq_filter_0_2.rds')
delta$key = paste0(delta$Gene,'_',delta$Node,'_',delta$Type)

library(notNMD)
wds2 = wds
tmp = wds2@diffSplicingResults
tmp$key = paste0(tmp$gene,'_',tmp$node,'_',tmp$type)
tmp = subset(tmp,key%in%subset(delta,Type=='CE' & !Complexity=='K0')$key) #  & abs(DeltaPsi)>0.01
wds2@diffSplicingResults = tmp
tmp = wds2@readCounts
tmp$key = paste0(tmp$Gene,'_',tmp$Node,'_',tmp$Type)
tmp = subset(tmp,key%in%subset(delta,Type=='CE' & !Complexity=='K0')$key) # & abs(DeltaPsi)>0.01
wds2@readCounts = tmp
wds.ce <- filterWhippetEvents(wds2, eventTypes = 'CE', psiDelta=0,probability=0)

# find exons in the gtf that overlap the skipped exon event
exons.ce <- findExonContainingTranscripts(wds.ce, 
                                          exons = exons,
                                          transcripts = transcripts)

# make skipped and included exon transcripts
# removes the skipped exon from all transcripts which contain it
skippedExonTranscripts <- skipExonInTranscript(skippedExons = exons.ce,
                                               exons=exons, 
                                               match="exact",
                                               whippetDataSet=wds.ce)
# skippedExonTranscripts_OLD = skippedExonTranscripts
# make non-skipped exon transcripts
normalTranscripts <- exons[exons$transcript_id %in% exons.ce$transcript_id]

# get ORF details for each set of transcripts
# note that notNMD requires upstream orf annotations
message('[*] orfs_normal')
orfs_normal <- getOrfs(normalTranscripts, BSgenome = g,
                       returnLongestOnly = FALSE, allFrames = TRUE, uORFs=TRUE)
save(orfs_normal,file='Rdata/NMD_analysis_GeneStructureTools_PC3.Rdata')
message('[*] orfs_skipped')
orfs_skipped <- getOrfs(skippedExonTranscripts[skippedExonTranscripts$set ==
                                                 "skipped_exon"],
                        BSgenome = g,
                        returnLongestOnly = FALSE, allFrames = TRUE, uORFs=TRUE)
save(orfs_normal,orfs_skipped,file='Rdata/NMD_analysis_GeneStructureTools_PC3.Rdata')
message('[*] orfs_included')

orfs_included <- getOrfs(skippedExonTranscripts[skippedExonTranscripts$set == 
                                                  "included_exon"],
                         BSgenome = g,
                         returnLongestOnly = FALSE, allFrames = TRUE, uORFs=TRUE)
save(orfs_normal,orfs_skipped,orfs_included,file='Rdata/NMD_analysis_GeneStructureTools_PC3.Rdata')


# calculate NMD probability
# --- note that if you have a different method for assessing NMD potential, you may substitute the values here
message('[*] Calculate NMD probability')
orfs_normal$nmd_prob <- notNMD::predictNMD(orfs_normal, "prob")
orfs_normal$nmd_class <- notNMD::predictNMD(orfs_normal)
orfs_skipped$nmd_prob <- notNMD::predictNMD(orfs_skipped, "prob")
orfs_skipped$nmd_class <- notNMD::predictNMD(orfs_skipped)
orfs_included$nmd_prob <- notNMD::predictNMD(orfs_included, "prob")
orfs_included$nmd_class <- notNMD::predictNMD(orfs_included)

orfs_normal <- orfs_normal[which(!is.na(orfs_normal$orf_length)),]
orfs_skipped <- orfs_skipped[which(!is.na(orfs_skipped$orf_length)),]
orfs_included <- orfs_included[which(!is.na(orfs_included$orf_length)),]

save(orfs_normal,orfs_skipped,orfs_included,file='Rdata/NMD_analysis_GeneStructureTools_PC3.Rdata')


# compare normal and skipped isoforms
# this time setting filterNMD to TRUE, which removes NMD targeted frames/isoforms where possible
message('[*] orfChange')
orfChange <- orfDiff(orfsX = orfs_included, 
                     orfsY = orfs_skipped, 
                     filterNMD = TRUE,
                     geneSimilarity = TRUE,
                     compareUTR=TRUE,
                     allORFs = orfs_normal)
save(orfChange,orfs_normal,orfs_skipped,orfs_included,file='Rdata/NMD_analysis_GeneStructureTools_PC3.Rdata')
message('[*] nmdChange')
nmdChange <- attrChangeAltSpliced(orfs_included,orfs_skipped,
                                  attribute="nmd_prob",
                                  compareBy="gene",
                                  useMax=FALSE)
m <- match(orfChange$id, nmdChange$id)
orfChange <- cbind(orfChange, nmdChange[m,-1])

orfChange$delta_nmd_prob = orfChange$nmd_prob_bygroup_x - orfChange$nmd_prob_bygroup_y

save(orfChange,nmdChange,orfs_normal,orfs_skipped,orfs_included,file='Rdata/NMD_analysis_GeneStructureTools_PC3.Rdata')
