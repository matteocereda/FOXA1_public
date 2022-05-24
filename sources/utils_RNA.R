# Tools for RNA-seq data analysis

# Load and Prepare Data ====
loadCounts <- function(countsDIR, pattern=NULL) {

  f  <- list.files(countsDIR        ,
                   pattern    = pattern,
                   full.names = T)

  counts <- do.call(cbind, lapply(f, read.delim, header=F, row.names=(1)))
  names(counts) <- gsub(".txt", "", basename(f))

  return(counts)


}

loadTPM <- function(countsDIR) {

  tpmData <- list.files(countsDIR,
                        pattern    = "TPM",
                        full.names = T)

  tpm <- read.csv2(tpmData,
                   stringsAsFactors = F,
                   header = T)

  x <- lapply(tpm[,-1], as.numeric)
  x <- as.matrix(as.data.frame(x))
  rownames(x) <- tpm$Gene

  return(x)
}

loadFunc <- function(type="TPM", ...) {

  if( type=="COUNTS" ) {
    return(loadCounts(...))
  } else if ( type=="TPM" ) {
    return(loadTPM(...))
  }
}

# Quality Controls ====

GrepStats <- function(f, pipeline="Hisat2") {

  # Pipeline Hisat2
  if( pipeline=="Hisat2" ) {
    grp <- paste0("grep 'reads\\|alignment' ",f,"  | cut -d ' ' -f1")
  } else {
    stop(message("[-] Please, provide a valid RNA-seq pipeline."))
  }

}

AlignStat <- function(STATDIR, ...) {

  fl <- list.files(STATDIR
                   , pattern = "txt"
                   , full.names = T)

  grps <- vapply(fl, GrepStats, FUN.VALUE = character(1), ...)

  tmp <- lapply(grps, system, intern=T)
  names(tmp) <- basename(fl)
  tmp <- do.call(rbind, tmp)

  nm <- rownames(tmp)

  tmp <- gsub("%","",tmp)
  tmp <- as.data.frame(apply(tmp, 2, as.numeric))
  colnames(tmp) <- c("library_size","mapping_rate")
  tmp$sample <- gsub("\\.txt","",nm)

  return(tmp)

}

# Normalization ====
countToTpm <- function(counts, effLen)
{
  # Transcripts Per Million (TPM)
  tpms <- apply(counts, 2,
                function(x) {
                  rate <- log(x) - log(effLen)
                  denom <- log(sum(exp(rate)))
                  c <- exp(rate - denom + log(1e6))
                  return(c)
                }
  )
  return(tpms)
}
#' @title Compute TPM for a read count matrix
#' @param counts A matrix of read counts with samples (columns) and genes (rows).
#' @param effLen A vector of gene cds length equal to number of rows of dfr.
#'
countToTpm2 = function(counts, effLen){
  tpms <- apply(counts, 2,
                function(x,y) {
                  rate <- log(x) - log(y)
                  denom <- log(sum(exp(rate), na.rm = T))
                  c <- exp(rate - denom + log(1e6))
                  return(c)
                }, y=effLen
  )
  return(tpms)
}

countToFpkm <- function(counts, effLen, pcg)
{
  # Fragments Per Kilobases per Million (FPKM)
  effLen <- effLen[rownames(counts)]
  pcg    <- intersect(pcg, rownames(counts))

  fpkms <- apply(counts, 2,
                 function(x) {
                   N <- sum(x[pcg])
                   c <- exp( log(x) + log(1e9) - log(effLen) - log(N) )
                   return(c)
                 }
  )
  return(fpkms)
}

# Principal Component Analysis ====

PCA <- function(expr, expr_cutoff=0.1, scale=T, center=T) {
  rna = t(expr)
  # remove lowly expressed
  rna = rna[,colSums(rna)>expr_cutoff]
  # remove constant values
  rna_sd = apply(rna, 2, sd)
  rna = rna[,rna_sd>0]
  # remove NA values
  rna_na = apply(rna, 2, function(x) sum(is.na(x)))
  rna = rna[,rna_na==0]
  RNA_pca = prcomp(rna, scale=scale, center=center)
  pca_summary = summary(RNA_pca)$importance
  pca_plot = as.data.frame(RNA_pca$x[,c("PC1", "PC2", "PC3")])
  list(pca_plot,pca_summary,dim(rna))
}

getPCA <- function(x) {

  rna <- t(x)
  rna <- rna[,colSums(rna)>1]
  RNA_pca <- prcomp(rna, scale. = T, center = T)
  pca_summary <- summary(RNA_pca)$importance

  return(list("RNA_pca"     = RNA_pca,
              "pca_summary" = pca_summary))
}

plotPCA <- function(x
                    , labels  = T
                    , groupBy = NULL
                    , pal     = NULL) {

  require(ggplot2)
  require(ggrepel)
  require(reshape2)

  my_theme <- theme(legend.position = "bottom"
                  , line = element_line(size=0.5)
                  , legend.title   = element_text(color="black",size=8)#rel(1.5)
                  , legend.text    = element_text(color="black",size=8)
                  , axis.title     = element_text(color="black",size=8)
                  , axis.text      = element_text(color="black",size=8)
                  , axis.line      = element_line(size=0.25)
                  , panel.grid.major = element_line(size=0.25)
                  , panel.grid.minor = element_blank()
                  , strip.text       = element_text(size=8)
                  , strip.background = element_rect(fill = NA))

  rna <- t(x)
  rna <- rna[,colSums(rna)>1]
  RNA_pca <- prcomp(rna, scale. = T, center = T)
  pca_summary <- summary(RNA_pca)$importance

  pca_plot <- as.data.frame(RNA_pca$x[,c("PC1", "PC2", "PC3", "PC4")])

  if( !is.null(groupBy) ) {
    pca_plot$sample <- groupBy[rownames(pca_plot)]
  } else {
    pca_plot$sample <- rownames(pca_plot)
  }

  if(is.null(pal)) {
    pal <- function(x) (viridis::viridis(x))
  }

  p <- ggplot(pca_plot, aes(x=PC1, y=PC2, col=sample)) + geom_point(size=3) +
    xlab(paste0("PC1 (",round(pca_summary[2,1]*100,1),"%)")) +
    ylab(paste0("PC2 (",round(pca_summary[2,2]*100,1),"%)")) +
    theme_bw() + ggtitle("RNA-seq PCA") + my_theme +
    theme(panel.grid.minor = element_blank()
          , plot.title = element_text(face="bold", hjust = 0.5, size=10)
          , aspect.ratio = 1) +
    scale_color_manual(values=pal(length(unique(pca_plot$sample))))


  if(labels) {
    p <- p + geom_label_repel(aes(label = rownames(pca_plot)),
                              fontface = 'bold', color = 'black', size=3,
                              box.padding = 0.35, point.padding = 0.5,
                              segment.color = 'grey50')
  }
  return(p)
}


# Correlation ====
plotCorrelation <- function(x, method="pearson", ...) {

  require(ComplexHeatmap)
  require(RColorBrewer)
  mcor <- cor(x, method=method, ...)

  myPalette <- colorRampPalette(brewer.pal(9, "YlOrRd"))
  hname <- paste0(method, " correlation")
  cHM <- Heatmap(mcor
                 , col  = myPalette(6)
                 , cell_fun = function(j, i, x, y, w, h, col) {
                         grid.text(round(mcor[i, j], digits = 2), x, y,gp = gpar(col='black', fontsize=6))
                       }
                 , name = hname
                 , row_names_gp = gpar(fontsize = 8)
                 , column_names_gp = gpar(fontsize = 8)
                 , heatmap_legend_param = list(title_position = "topcenter",
                                               # legend_width  = unit(4, "cm"),
                                               # legend_height = unit(0.5, "mm"),
                                               values_gp     = gpar(fontsize=8),
                                               legend_direction = "horizontal"))
  return(draw(cHM, heatmap_legend_side = "bottom"))
}

# Differential Expression Analysis ====
calculateDiffExpr <- function(m
                              , group  = NULL
                              , filter = T
                              , filter.cpm.th = 1
                              , filter.sample.th = 2
                              , method = "exact"
                              , reference = NULL # Control group
                              , design = NULL # Custom design matrix
                              , cf = NULL # Testing coefficient (default: all columns in design matrix against control)
                              , contrast = NULL # Contrast matrix
                              , return.y  = F) {

  # Compute differential expression with EdgeR

  require(edgeR)
  message("[*] Run EdgeR for Differential Expression Analysis")

  if(is.null(group)) {
    group <- as.factor(colnames(m))
  } else {
    m <- m[,names(group)]
    if (!is.factor(group))   group <- as.factor(group)
    if (!is.null(reference)) group <- relevel(group, reference)
  }

  message(" -- Condition: ", paste0(levels(group), collapse = "-"))

  y <- DGEList(counts=m, genes=rownames(m), group = group)
  y <- calcNormFactors(y)

  # Clean environment
  rm(m)
  gc(verbose = F)

  if(filter) {
    message(" -- Filtering lowly expressed genes")
    message("    -- Threshods:")
    message("     * CPM = "     , filter.cpm.th)
    message("     * Samples = " , filter.sample.th)

    if(filter.sample.th>=1) {
      # Number of samples
      fs <- filter.sample.th
    } else {
      # Percentage of samples
      fs <- floor(ncol(m)*filter.sample.th)
    }

    keep <- rowSums(cpm(y)>filter.cpm.th) >= fs
    y <- y[keep, , keep.lib.sizes=FALSE]
  }

  if(is.null(design)) {
    # Standard design matrix
    # First column is control ( = reference in group)
    message(" -- Standard design matrix")
    design <- model.matrix(~group)

  } else {
    # Customized design matrix
    message(" -- Custom model matrix")
  }

  if(is.null(cf)) {
    # Test all columns in design matrix
    # against control group
    cf  <- 2:ncol(design)
  }

  rownames(design) <- colnames(y)
  y <- estimateDisp(y, design, robust=TRUE)

  message(" -- Testing differential expression, method: ", method)

  if(method=="exact") {
    # Exact test
    de <- exactTest(y)
    de <- topTags(de, n = Inf)

  } else if(method=="lrt") {
    # Likelihood-ratio test
    fit <- glmFit(y, design)
    lrt <- glmLRT(fit, coef=cf)
    de  <- topTags(lrt, n = Inf)

  } else if(method=="qlf") {
    # Quasi-likelihood F-test
    fit <- glmQLFit(y, design)
    if(is.null(contrast)){
      # Standard comparison
      qlf <- glmQLFTest(fit, coef=cf)
    } else {
      # Anova-like for multiple group comparison
      message(" -- Anova-like for multiple group comparison")
      qlf <- glmQLFTest(fit, contrast = contrast)
    }
    de  <- topTags(qlf, n = Inf)

  } else {
    stop(message("[-] Method not available. Please, provide a valid one ( exact / lrt / qlf )."))
  }

  if(return.y) {

    res <- list("y"  = y,
                "de" = de)
    return(res)

  } else {

    return(de)
  }
}

saveXLSresEdgeR <- function(res, outfile) {

  require(xlsx)

  message("[+] Saving results to file: ", outfile)

  if( is.list(res) & !is.data.frame(res)) {
    for( i in 1:length(res) ) {
      tox <- res[[i]]
      ap  <- F
      if (i>1) ap <- T
      write.xlsx2(tox, file = outfile, sheetName = names(res)[i], append = ap, row.names = F)
    }
  } else {
    if(!is.data.frame(res)) res <- as.data.frame(res)
    tox <- res
    write.xlsx2(tox, file = outfile, sheetName = "edgeR", append = ap, row.names = F)
  }

}

getDEres <- function(x, genes) {
  # Get DE results for specific genes
  if(!is.data.frame(x)) x <- x$table
  x[intersect(genes, rownames(x)),]
}

getDEgenes <- function(x, fdrTh=0.1, fcTh=0.5) {
  # Get DE genes satisfying thresholds
  message("[+] Get differentially expressed genes, thresholds:")
  message(" - logFC = ", fcTh)
  message(" - FDR = ", fdrTh)

  if(!is.data.frame(x)) x <- x$table

  if( !is.null(fdrTh) ) {
    x <- x[ x[,'FDR'] <= fdrTh , ]
  }

  if(nrow(x) == 0) {
    stop(message("[-] No differentially expressed genes with provided FDR"))
  }

  if( !is.null(fcTh) ) {
    fcidx <- grep("logFC", colnames(x))

    if( length(fcidx)==1 ) {
      x <- x[ abs(x[,fcidx]) >= fcTh , ]
    } else {
      x <- x[ rowSums( abs(x[,fcidx]) >= fcTh) >= 1,]
    }

    if(nrow(x) == 0) {
      stop(message("[-] No differentially exxpressed genes with provided logFC"))
    }
  }

  return(x)

}

getDEgsigned <- function(deg, signed=(-1))
{
  fcidx <- grep("logFC", colnames(deg))
  x     <- deg[,fcidx]

  max.idx <- apply(x, 1, function(m) which.max(abs(m)))

  y <- vector(mode = 'numeric', length = nrow(x))
  names(y) <- names(max.idx)

  for(i in names(max.idx)) {
    y[i] <- x[i, max.idx[i]]
  }

  gene.idx <- names(y[which(sign(y)==signed)])

  return(x[gene.idx,])
}

plotRNAVolcanos <- function(de, lfcTh=1, pvTh=0.05, top=5)
{
  require(ggplot2)
  require(RColorBrewer)
  require(ggrepel)

  pal <- brewer.pal(9, "Set1")

  fcIdx <- grep("logFC|log2FoldChange", colnames(de), value = T)
  pvIdx <- grep("adj.P.Val|FDR", colnames(de), value = T)
  tmp   <- de[,c(fcIdx,pvIdx)]
  colnames(tmp) <- c("lfc","padj")
  tmp$status <- "none"
  tmp$status[which(tmp$lfc>=lfcTh & tmp$padj<=pvTh)] <- "Up-regulated"
  tmp$status[which(tmp$lfc<=(-lfcTh) & tmp$padj<=pvTh)] <- "Down-regulated"
  tmp$status <- factor(tmp$status)

  topUP <- tmp[tmp[,'status']=="Up-regulated",]
  topUP <- rownames(topUP[1:top,])
  topDW <- tmp[tmp[,'status']=="Down-regulated",]
  topDW <- rownames(topDW[1:top,])

  tmp$lab <- rownames(tmp)
  tmp$lab[-which(rownames(tmp)%in%c(topDW, topUP))] <- ""

  p <- ggplot(tmp, aes(x=lfc, y=-log10(padj),col=status, label=lab)) + geom_point(size=2) + theme_bw() +
    geom_hline(yintercept = -log10(pvTh), linetype = 'dashed') +
    geom_vline(xintercept = lfcTh, linetype = 'dashed') +
    geom_vline(xintercept = -lfcTh, linetype = 'dashed') +
    xlab("logFC") + ylab("adjusted P-value") +
    theme_bw() + theme(panel.grid = element_blank(),legend.position = "right",
                       aspect.ratio = 0.8, plot.title = element_text(face="bold", hjust = 0.5)) +
    scale_color_manual(values = c(pal[2],'black',pal[1])) +
    geom_text_repel()

  return(p)
}
# Clustering ====
findClustSSE <- function(scaledata, nKM=20)
{
  wss <- (nrow(scaledata)-1)*sum(apply(scaledata,2,var))
  for (i in 2:nKM) wss[i] <- sum(kmeans(scaledata,
                                       centers=i)$withinss)
  plot(1:nKM, wss, type="b", xlab="Number of Clusters",
       ylab="Within groups sum of squares")
}

findClustASW <- function(scaledata, nKM=20, ret=F)
{
  require(cluster)
  sil <- rep(0, nKM)
  #repeat k-means for 1:n and extract silhouette:
  for(i in 2:nKM){
    k1to20 <- kmeans(scaledata, centers = i, nstart = 25, iter.max = 20)
    ss <- silhouette(k1to20$cluster, dist(scaledata))
    sil[i] <- mean(ss[, 3])
  }

  # Plot the  average silhouette width
  plot(1:nKM, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
  abline(v = which.max(sil), lty = 2)

  sil.best <- as.numeric(which.max(sil))
  cat("Average silhouette width optimal number of clusters:", sil.best, "\n")

  if(ret) return(sil.best)
}

findClustCal <- function(scaledata, nKM=20, ret=F)
{
  # Calinski-Harabasz index
  require(vegan)
  fit <- cascadeKM(scaledata, 1, nKM, iter = 100)
  plot(fit, sortg = TRUE, grpmts.plot = TRUE)
  calinski.best <- as.numeric(which.max(fit$results[2,]))
  cat("Calinski criterion optimal number of clusters:", calinski.best, "\n")

  if(ret) return(calinski.best)
}

setKMClusters <- function(m, method="silhouette", ...)
{
  scaledata <- t(scale(t(m)))

  if( method=="sse" ) {
    return(findClustSSE(scaledata, ...))
  }

  if( method=="silhouette") {
    return(findClustASW(scaledata, ...))
  }

  if( method=="calinski") {
    return(findClustCal(scaledata, ...))
  }
}
# Heatmap ====
get_heatmap3 <- function(m
                         , annotDF  = NULL
                         , annotCol = NULL
                         , fig_out  = NULL
                         , fig_h = unit(4,'cm')
                         , fig_w = unit(4,'cm')
                         , retHm    = F
                         , ...){

  require(ComplexHeatmap)
  require(circlize)
  require(RColorBrewer)

  base_mean <- rowMeans(m)
  m_scaled <- t(apply(m, 1, scale))
  colnames(m_scaled) <- colnames(m)

  bPalette <- colorRampPalette(brewer.pal(9, "Reds"))
  myPalette <- c("blue","black","red")
  ramp <- colorRamp2(c(-2, 0, 2), myPalette)


  if (!is.null(annotDF)) {
    if (!is.null(annotCol)) {
      ha_column <- HeatmapAnnotation(df  = annotDF,
                                     col = annotCol,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8))
                                     ,annotation_height = unit(4, "mm"))
    } else {
      ha_column <- HeatmapAnnotation(df  = annotDF,
                                     annotation_legend_param = list(title_gp  = gpar(fontsize=8),
                                                                    values_gp = gpar(fontsize=8))
                                     , annotation_height = unit(4, "mm"))
    }
  } else {
    ha_column <- new("HeatmapAnnotation")
  }

  hm <- Heatmap(m_scaled,
                col = ramp,
                show_row_dend = T,
                row_names_side = "left",
                row_names_gp = gpar(fontsize=8),
                column_names_gp = gpar(fontsize=8),
                column_title_gp = gpar(fontsize=10, fontface="bold"),
                heatmap_legend_param = list(title = "row Z-score",
                                            title_gp = gpar(fontsize=8),
                                            title_position = "topcenter",
                                           # legend_width  = unit(4, "cm"),
                                            # legend_height = unit(0.5, "mm"),
                                            values_gp     = gpar(fontsize=8),
                                            legend_direction = "vertical")
                , top_annotation = ha_column
                , ...)

  bmscale <- summary(base_mean)
  bmramp <- colorRamp2(c(bmscale[1],bmscale[3],bmscale[5]), bPalette(3))
  bmh <- Heatmap(base_mean
                 , name = "Mean Expression"
                 , column_names_gp = gpar(fontsize=8)
                 , show_row_names = FALSE
                 , width = unit(3, "mm")
                 , col = bmramp
                 , heatmap_legend_param = list(title = "Base Mean",title_gp = gpar(fontsize=8)))

  hmOut <- hm + bmh

  if(!is.null(fig_out)){
    pdf(file = fig_out, useDingbats = F, h=fig_h, w=fig_w, paper = "a4")
    draw(hmOut, heatmap_legend_side = "right")
    dev.off()
  } else{
    draw(hmOut, heatmap_legend_side = "right")
  }

  if(retHm) return(hmOut)
}

get_clusters <- function(m, hm){
  # Retrieve clusters from K-means
  clusters <- lapply(row_order(hm),
                     function(x){
                       rownames(m[x,])
                     }
  )
  return(clusters)
}

# Gene Ontology ====
getGO <- function(x, reg="UP", ORGANISM = "mmusculus", gl_input=T){

  require(gProfileR)
  # Perform GO term enrichment using gProfileR
  if(gl_input) {
    gene_list <- x
  } else {
    if ( reg!="both") {
      gene_list <- subset(x, regulation==reg)$symbol
    } else {
      gene_list <- x$symbol
    }
  }

  gprofiler( gene_list,
             organism = ORGANISM,
             # custom_bg = background,
             max_p_value = 0.05,
             ordered_query =T,
             significant = T,
             correction_method='gSCS',
             min_set_size = 10,
             # evcodes=T,
             # exclude_iea = T,
             hier_filtering ='moderate',
             # src_filter = collections
             # png_fn = "Prostate/Figures/Goprofiler.AS.png"
             # include_graph=T
             src_filter = c("GO",'KEGG',"REAC","CORUM","OMIM",  "HP")
  )

}

processGO <- function(x, pvTh=0.05, cut=F){

  x <- subset(x, p.value<=0.05)
  x <- x[,c("term.name", "p.value")]
  if( any(duplicated(x$term.name))) x <- x[-which(duplicated(x$term.name)),]
  x$term.name <- factor(x$term.name,
                        levels = x$term.name[order(x$p.value, decreasing = F)])
  if(cut){
    if(nrow(x)>25) x <- x[1:25,]
  }
  x <- x[with(x, order(p.value, decreasing = F)),]
  return(x)
}

plotGObars <- function(GO, tool='clusterProfiler'){

  source("/sto1/andrea/Rtools/Plots/theme_setting.R")

  if( tool == 'gprofiler') {
    ggplot(GO, aes(x=term.name, y=-log10(p.value))) + geom_col(position = "dodge") +
      coord_flip() + theme_bw() + my_theme + theme(panel.grid = element_blank())
  } else if( tool == 'clusterProfiler') {
    ggplot(GO, aes(x=Description, y=-log10(qvalue))) + geom_col(position = "dodge") +
      coord_flip() + theme_bw() + my_theme + theme(panel.grid = element_blank())
  }
}

getGO_v2 <- function(geneList
                     , species      = 'mm'
                     , custom       = F
                     , qvalueCutoff = 0.05
                     , ont          = "BP"
                     , simplify     = T
                     , ...)
{
  require(clusterProfiler)

  if(species=='mm') {
    require(org.Mm.eg.db)
    db <- org.Mm.eg.db

  } else if(species=='hg') {
    require(org.Hs.eg.db)
    db <- org.Hs.eg.db

  }

  if(is.null(ont)) ont <- "BP"

  if(!custom){
    ego <- enrichGO(   gene         = geneList
                       , OrgDb         = db
                       , keyType       = 'SYMBOL'
                       , ont           = ont
                       , pAdjustMethod = "BH"
                       , pvalueCutoff  = 0.01
                       , qvalueCutoff  = qvalueCutoff )
  } else {
    ego <- enrichGO(  gene         = geneList
                      , OrgDb        = get(db)
                      , ... )
  }

  if(simplify) {
    ego <- simplify(ego)
  }

  return(ego)
}

# Modelling ====
find.variable.genes <- function(m
                                , n=1000
                                , x0 = c(-0.5, 0.5)
                                , ret.plot=T)
{

  log2_cv2 <- log2(apply(m, 1, function(r) (sd(r)/mean(r))**2))
  log2_mn  <- log2(apply(m, 1, function(r) mean(r)))

  idx      <- names(which(!is.na(log2_cv2)))

  log2_cv2 <- log2_cv2[idx]
  log2_mn  <- log2_mn[idx]

  noise.model <- function(log2_mn, log2_cv2)
  {
    function(x) sum(abs(log2((2 ** log2_mn) ** x[1] + x[2]) - log2_cv2))
  }

  xopt <- optim(x0, noise.model(log2_mn, log2_cv2))
  log2_cv2_fit  <- log2(( 2 ** log2_mn) ** xopt$par[1] + xopt$par[2])
  log2_cv2_diff <- log2_cv2 - log2_cv2_fit

  idx.ord <- order(log2_cv2_diff, decreasing = T)

  log2_cv2_diff[idx.ord][1:n] <- 'TRUE'
  log2_cv2_diff[idx.ord][(n+1):length(log2_cv2_diff)] <- 'FALSE'

  y <- cbind.data.frame(    'log2_cv2'      = log2_cv2[idx.ord]
                          , 'log2_mean'     = log2_mn[idx.ord]
                          , 'log2_cv2_fit'  = log2_cv2_fit[idx.ord]
                          , 'High'          = log2_cv2_diff[idx.ord] )

  genes <- rownames(y)[y$High=='TRUE']

  if(ret.plot) {

    source("/sto1/andrea/Rtools/Plots/theme_setting.R")

    p <- ggplot(y, aes(x=log2_mean, y=log2_cv2, col=High)) +
      geom_point(size=1) +
      geom_line(aes(y=log2_cv2_fit), col='#CC0000') +
      theme_bw() + my_theme + scale_color_manual(values=c('black', '#0066CC')) +
      xlab("log2(mean)") + ylab("log2(cv^2)")

    return(list("genes" = genes,
                "plot" = p))

  } else {
    return(genes)
  }
}
