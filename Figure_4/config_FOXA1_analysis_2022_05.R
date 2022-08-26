# Libraries ----

library(circlize)
library(ComplexHeatmap)
library(dplyr)
library(snow)
library(viridis)
library(tidyr)
library(plyr)
library(pbapply)


# Functions ----

p   = function(...){paste0(...)}

cha = function(...){as.character(...)}  

num = function(...){as.numeric(...)}

get_exons_with_tetramer_region1=function(tet, RESULTS,SPLICING_FILE){
  sp = read.table(SPLICING_FILE, sep=";", h=F)
  f = list.files(RESULTS, tet, full.names = T, recursive = T)
  rc  = read.delim(f[grep("region_count.tsv",f)])
  bed = read.table(f[grep("bed",f)], sep="\t", skip=1)
  sel = subset(rc, type!="0" & (hits_region1!=0))
  sp.sel = subset(sp, V1%in%sel$myRID)
  sel = cbind(sp.sel[,3:8],sel[match(sp.sel$V1,sel$myRID),])
  if(length(sel$V3)==0){
    return(data.frame(key=character()))
  }
  sel$key = with(sel,p(V3,":",V6,"-",V7))
  return(sel)
}

get_exons_with_tetramer_region2=function(tet, RESULTS,SPLICING_FILE){
  sp = read.table(SPLICING_FILE, sep=";", h=F)
  f = list.files(RESULTS, tet, full.names = T, recursive = T)
  rc  = read.delim(f[grep("region_count.tsv",f)])
  bed = read.table(f[grep("bed",f)], sep="\t", skip=1)
  sel = subset(rc, type!="0" & (hits_region2!=0))
  sp.sel = subset(sp, V1%in%sel$myRID)
  sel = cbind(sp.sel[,3:8],sel[match(sp.sel$V1,sel$myRID),])
  if(length(sel$V3)==0){
    return(data.frame(key=character()))
  }
  sel$key = with(sel,p(V3,":",V6,"-",V7))
  return(sel)
}

get_exons_with_tetramer_region3=function(tet, RESULTS,SPLICING_FILE){
  sp = read.table(SPLICING_FILE, sep=";", h=F)
  f = list.files(RESULTS, tet, full.names = T, recursive = T)
  rc  = read.delim(f[grep("region_count.tsv",f)])
  bed = read.table(f[grep("bed",f)], sep="\t", skip=1)
  sel = subset(rc, type!="0" & (hits_region3!=0))
  sp.sel = subset(sp, V1%in%sel$myRID)
  sel = cbind(sp.sel[,3:8],sel[match(sp.sel$V1,sel$myRID),])
  if(length(sel$V3)==0){
    return(data.frame(key=character()))
  }
  sel$key = with(sel,p(V3,":",V6,"-",V7))
  return(sel)
}

lmb.cluster.fisher.test = function(d, tot.d, contr, tot.contr ){
  pvd = matrix(0,nrow=nrow(d),ncol=ncol(d))
  colnames(pvd) = colnames(d)
  for(j in 1:ncol(pvd))
    for(i in 1:nrow(pvd)){
      test = matrix( data=c( d[i,j], contr[i,j], tot.d - d[i,j] , tot.contr - contr[i,j] ), nrow=2 )
      pvd[i,j] = fisher.test(test, alternative = "greater")$p.value
    }
  pvd
}

BC = function(x,y){ 
  x = x/sum(x)
  y = y/sum(y)
  return(sum(sqrt(x*y)))
}

changeParam = function(inIntronOld,inExonOld,inIntronNew,inExonNew){
  pos_adj = c(1 + inExonOld   + (-inExonNew:inIntronNew),
              1 + inExonOld   + 2*inIntronOld + (-inIntronNew:inExonNew),
              1 + 3*inExonOld + 2*inIntronOld + (-inExonNew:inIntronNew),
              1 + 3*inExonOld + 4*inIntronOld + (-inIntronNew:inExonNew))
  return(pos_adj)
}

# 
# SEeCLIPpeaks = function(datExpr,peaks,inIntron, inExon, n_cores){
#   # datExpr è del tipo 1;10046;chr6;+;43746655;43749692;43749824;43752277;1
#   # peaks è del tipo chr1	17496	17498
#   
#   # primo return: posizioni assolute genomiche delle regioni di interesse
#   # secondo return: matrice di 1 e 0: 1 se overlap con eclip, 0 elsewise
#   
#   # inExon = 58
#   # inIntron = 208
#   # matrice di (nrow) X (267). ogni riga va da -58 a 208
#   m1 = matrix(rep(-inExon:inIntron,nrow(datExpr)),
#               nrow  = nrow(datExpr),
#               ncol  = length(-inExon:inIntron),
#               byrow = T)
#   
#   m2 = matrix(rep(-inIntron:inExon,nrow(datExpr)),
#               nrow  = nrow(datExpr),
#               ncol  = length(-inIntron:inExon),
#               byrow = T)
#   
#   pos = cbind(m1 + datExpr$V5,
#               m2 + datExpr$V6,
#               m1 + datExpr$V7,
#               m2 + datExpr$V8)
#   #library(Rmpi) 
#   clus = snow::makeCluster(n_cores, type = "SOCK")
#   snow::clusterExport(clus,list = c("datExpr","peaks","inIntron","inExon"))
#   values_t = snow::parSapply(clus, 1:nrow(datExpr), function(x){
#     
#     m1 = matrix(rep(-inExon:inIntron,nrow(datExpr)),
#                 nrow  = nrow(datExpr),
#                 ncol  = length(-inExon:inIntron),
#                 byrow = T)
#     
#     m2 = matrix(rep(-inIntron:inExon,nrow(datExpr)),
#                 nrow  = nrow(datExpr),
#                 ncol  = length(-inIntron:inExon),
#                 byrow = T)
#     
#     pos = cbind(m1 + datExpr$V5,
#                 m2 + datExpr$V6,
#                 m1 + datExpr$V7,
#                 m2 + datExpr$V8)
#     
#     tmp0 = subset(peaks, V1==datExpr$V3[x])
#     tmp  = subset(tmp0, V2%in%pos[x,] | V3%in%pos[x,])
#     
#     if(nrow(tmp)==0){
#       return(pos[x,]*0)
#     }
#     
#     tmp         = do.call(rbind, lapply(1:nrow(tmp),function(y){cbind(seq(as.numeric(tmp$V2[y]),as.numeric(tmp$V3[y])),1)}))
#     a           = tmp[match(pos[x,],tmp[,1]),2]
#     a[is.na(a)] = 0
#     
#     if(datExpr$V4[x]=="-"){a = a[length(a):1]}
#     
#     return(a)
#   }) 
#   
#   stopCluster(clus)
#   
#   res = list(pos = pos, values = t(values_t))
#   
#   return(res)
# }


SEeCLIPpeaks = function(datExpr,peaks,inIntron, inExon){

  m1 = matrix(rep(-inExon:inIntron,nrow(datExpr)),
              nrow  = nrow(datExpr),
              ncol  = length(-inExon:inIntron),
              byrow = T)

  m2 = matrix(rep(-inIntron:inExon,nrow(datExpr)),
              nrow  = nrow(datExpr),
              ncol  = length(-inIntron:inExon),
              byrow = T)

  pos = cbind(m1 + datExpr$V5,
              m2 + datExpr$V6,
              m1 + datExpr$V7,
              m2 + datExpr$V8)

  values_t = pbapply::pbsapply(1:nrow(datExpr),function(x){
    tmp0 = subset(peaks, V1==datExpr$V3[x])
    tmp  = subset(tmp0,V2%in%pos[x,] | V3%in%pos[x,])

    if(nrow(tmp)==0){
      return(pos[x,]*0)
    }

    tmp         = do.call(rbind,lapply(1:nrow(tmp),function(y){cbind(seq(num(tmp$V2[y]),num(tmp$V3[y])),1)}))
    a           = tmp[match(pos[x,],tmp[,1]),2]
    a[is.na(a)] = 0

    if(datExpr$V4[x]=="-"){a = a[length(a):1]}

    return(a)
  })

  res = list(pos = pos, values = t(values_t))

  return(res)

}

anno_reg = function(v, col1 = "red", col2 = "blue", col3 = "yellow", h = 4){
  anno_simple(gp     = gpar(col="white"),
              x      = v, 
              col    = c("white"  = "white",
                         "red"    = col1,
                         "blue"   = col2,
                         "yellow" = col3),
              width  = unit(3,"mm"),
              height = unit(h,"mm"))
}

anno_mt_row = function(n){
  anno_simple(which = "row", 
              gp    = gpar(col="white"),
              x     = rep("a", n), 
              col   = c("a" = "white"),
              width = unit(1,"mm"))
}

anno_sign = function(v1, v2, col){
  anno_simple(which = "row",
              gp    = gpar(col="white"),
              x     = (sapply(strsplit(v1,"[.]"),head,1) %in% v2) %>% cha, 
              col   = c("TRUE" = col, "FALSE" = "white"),
              width = unit(3,"mm"))
}