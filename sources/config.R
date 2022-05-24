# library(TCGAbiolinks)
library(SummarizedExperiment)
library(data.table)
library(ggrepel)
# SELECTION OF PTEN LOSS SAMPLE

library(tidyr)
library(grid)
library(gridExtra)
library(statmod)
library(R.utils)

library(RColorBrewer)

plot_GSECA <- function(df, FDR=.1) {
  sel = unique(subset(df, pv2!="" & fdr<FDR )$path)
  toplot = subset(df, path%in%sel)[,c('path','variable','type','perc','pv','pv2','fdr')]
  tpm    = ddply(toplot, .(path, variable), summarise
                 , L2R =  log2(perc[type=="PTEN_LOSS"]/perc[type!="PTEN_LOSS"]) )
  toplot= subset(toplot, type=="PTEN_LOSS")
  toplot$id = paste0(toplot$path,".",toplot$variable)
  tpm$id = paste0(tpm$path,".",tpm$variable)
  toplot$L2R = tpm[match(toplot$id,tpm$id),'L2R']
  toplot$path =fix_names(toplot$path)
  toplot$variable=factor(as.character(toplot$variable), levels=c('not_expr','le','me','he') )
  toplot
}

myPalette <- colorRampPalette(rev(brewer.pal(9, "OrRd")), space="Lab")

mapply_pb = function(FUN, X, Y,  ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- mapply(wrapper, X, Y, ...)
  close(pb)
  res
}

lapply_pb = function(X, FUN, ...){
  env <- environment()
  pb_Total <- length(X)
  counter <- 0
  pb <- txtProgressBar(min = 0, max = pb_Total, style = 3)

  # wrapper around FUN
  wrapper <- function(...){
    curVal <- get("counter", envir = env)
    assign("counter", curVal +1 ,envir=env)
    setTxtProgressBar(get("pb", envir=env), curVal +1)
    FUN(...)
  }
  res <- lapply(X, wrapper, ...)
  close(pb)
  res
}
fix_names = function(p) {
  return(capitalize(gsub(pattern = "_", replacement = " ",fixed = T,x = tolower(gsub(pattern = "KEGG_",replacement = "",fixed = T,x = as.character(as.character(p)))))))
}

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )

  # Rename the "mean" column
  datac <- rename(datac, c("mean" = measurevar))

  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean

  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval:
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult

  return(datac)
}
summarySE_v2=function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)

  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }

  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     median = median   (xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     se   = sd(xx[[col]], na.rm=na.rm)/ sqrt(length2(xx[[col]], na.rm=na.rm))
                   )
                 },
                 measurevar
  )
  return(datac)
}


col.new = c('#af8dc3','#f7f7f7','#7fbf7b',"#00BFC4","darkblue") #,'#f5f5f5'

get_fisher2 = function (u){
  u$pv = NA
  u$or = NA
  u$pw = NA
  for(i in 1:nrow(u)){
    f = fisher.test( matrix( unlist(u[i,c("a","ra","b","rb")]), nrow=2 ) ) #, alt="greater"
    u$pw[i] = NA
    # u$pw[i] = power.fisher.test(u[i,"a"]/(u[i,"a"]+u[i,"ra"]), u[i,"b"]/(u[i,"b"]+u[i,"rb"]),
                                # u[i,"a"]+u[i,"ra"], u[i,"b"]+u[i,"rb"], alpha=.05, nsim=10, alternative="two.sided")
    u$pv[i] = f$p
    u$or[i] = f$estimate
  }
  return(u)
}

make_pathway_list = function(pathMsigDbFile) {
  inputFile <- pathMsigDbFile
  con  <- file(inputFile, open = "r")

  c = 1
  pathway.list <- vector(mode="list",length=0)
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0) {
    print(c)
    myVector <- do.call("rbind",strsplit(oneLine, "\t"))
    t = vector(mode="list",length=1)
    t[[1]] = myVector[3:length(myVector)]
    names(t) = myVector[1]
    pathway.list = c(pathway.list,t)
    c = c+1
  }

  close(con)
  return(pathway.list)
}

get_summary_expr_class=function(r, pl, lvls=c("he","me","not_expr","le"), cl="PTEN_LOSS",
                                cols = c('PTEN_LOSS','rest.PTEN_LOSS','WT','rest.WT','pv','OR','pw')){
  # create progress bar
  pb <- txtProgressBar(min = 0, max = length(pl), style = 3)
  res = list()
  ii = 0
  for(i in 1:length(pl)){
    setTxtProgressBar(pb, ii)
    ii= ii+1
    sel = cbind(subset(r, symbol%in%pl[[i]]), path = names(pl)[i])


    ho =  ddply(sel, .(symbol), summarise,
                a   = sum(expr_class[type==cl]=="HE"),  ra  = sum(expr_class[type==cl]!="HE")
                , b   = sum(expr_class[type!=cl]=="HE"),  rb  = sum(expr_class[type!=cl]!="HE")
    )

    mo =    ddply(sel, .(symbol), summarise,
                  a   = sum(expr_class[type==cl]=="ME"),  ra  = sum(expr_class[type==cl]!="ME")
                  , b   = sum(expr_class[type!=cl]=="ME"),  rb  = sum(expr_class[type!=cl]!="ME")
    )

    lu =   ddply(sel, .(symbol), summarise,
                 a   = sum(expr_class[type==cl]=="LE"),  ra  = sum(expr_class[type==cl]!="LE")
                 , b   = sum(expr_class[type!=cl]=="LE"),  rb  = sum(expr_class[type!=cl]!="LE")
    )

    ne =   ddply(sel, .(symbol), summarise,
                 a   = sum(expr_class[type==cl]=="not Expr"),  ra  = sum(expr_class[type==cl]!="not Expr")
                 , b   = sum(expr_class[type!=cl]=="not Expr"),  rb  = sum(expr_class[type!=cl]!="not Expr")
    )
    # nen =  ddply(subset(sel, cancer=="COAD"), c("cancer","gene"), summarise, a   = sum(expr_class[sync=="SYNC"]=="not Expr NT"),  ra  = sum(expr_class[sync=="SYNC"]!="not Expr NT"),  b   = sum(expr_class[sync!="SYNC"]=="not Expr NT"), rb  = sum(expr_class[sync!="SYNC"]!="not Expr NT"))

    ho = get_fisher2(ho);
    mo = get_fisher2(mo);
    lu = get_fisher2(lu);
    ne = get_fisher2(ne)
    # nen = get_fisher2(nen)
    colnames(ho)[2:8]= paste0("he.",cols);
    colnames(mo)[2:8]= paste0("me.",cols);
    colnames(lu)[2:8]= paste0("le.",cols);
    colnames(ne)[2:8]= paste0("not_expr.",cols);
    # colnames(nen)[3:9]= paste("not_expr_nt.",cols, sep="", coll="");
    st = cbind(ho,
               mo[,2:8],
               lu[,2:8],
               ne[,2:8]
               # , nen[,3:9]
    )

    st$sig.in.one.cat = with(st, (he.pv<.05 | me.pv<.05 | le.pv<.05 | not_expr.pv<.05  )) # | not_expr_nt.pv<.05

    #     pxx = melt(st, id.vars=c('cancer','gene'), measure.vars = c('he.pv','me.pv','le.pv','not_expr.pv','not_expr_nt.pv'))
    #     pxx2 = melt(st, id.vars=c('cancer','gene'), measure.vars = c('he.OR','me.OR','le.pv','not_expr.OR','not_expr_nt.OR'))
    #     pxx = cbind(pxx, 'OR'=pxx2[,'value'])
    #     pxx$variable = sapply(strsplit(as.character(pxx$variable),"\\."), function(x) x[1])
    #     pxx$variable=factor(pxx$variable, levels=c("not Expr NT","not Expr","LE","ME","HE"),ordered=T)
    #     pxx$sig =pxx$value<0.05
    res[[i]] = st
  }
  close(pb)
  res
}

global_summary_expr_class = function(res,pl, lvls=c("he","me","le","not_expr"), cl=c("PTEN_LOSS","WT")){
  dfx = rbind()
  for(i in names(res)){
    # print(i)
    a = res[[i]]
    sa = c(apply(a[,paste0(lvls,".",cl[1])],2,sum),'type'=cl[1]); names(sa)[1:length(lvls)]=lvls
    sb = c(apply(a[,paste0(lvls,".",cl[2])],2,sum),'type'=cl[2]); names(sb)[1:length(lvls)]=lvls
    px = rbind(sa,sb); px=as.data.frame(px)

    for(z in 1:length(lvls)) px[,z] = as.numeric(as.character(px[,z]))
    px$n = c( sum( px[1,1:length(lvls)] ),sum( px[2,1:length(lvls)] ))
    pxx = melt(px, id.vars=c("type","n"))
    pxx$perc = (pxx$value/pxx$n)*100
    pxx$variable=factor(pxx$variable, levels=lvls, ordered=T)
    pxx$type = factor(pxx$type, levels=cl )
    pxx$rest = pxx$n - pxx$value
    pxx$pv =NA
    pxx$pw =NA
    for( j in pxx$variable ){
      u = subset(pxx, variable==j)[c("value","rest")]
      f = fisher.test(u) #, alt="greater"
      if(u[1,"value"]>0 & u[2,"value"]>0 ){
        pw = power.fisher.test(u[1,"value"]/(u[1,"value"]+u[1,"rest"]), u[2,"value"]/(u[2,"value"]+u[2,"value"]),
                               u[1,"value"]+u[1,"rest"], u[2,"value"]+u[2,"rest"], alpha=.05, nsim=10, alternative="two.sided")
        pxx$pw[which(pxx$variable==j)]=pw
      }
      pxx$pv[which(pxx$variable==j)]=f$p
    }

    pxx = ddply(pxx, .(variable), mutate, delta=perc[which(type==cl[1])]-perc[which(type==cl[2])])
    # pxx$delta=NA
    # for( j in pxx$variable ){
      # pxx$delta[which(pxx$variable==j)]=tmp$delta[which(tmp$variable==j)]
    # }
    testres = c("D","E")
    pxx$testres = testres[as.numeric(pxx$delta>0)+1]
    pxx$pv2 = paste(as.character(round(pxx$pv,2)),"(", pxx$testres ,")", sep="",coll="")
    pxx$pv2[which(pxx$pv>0.05 | pxx$type==cl[2])]=""
    pxx$path=i
    pxx=pxx[rev(order(pxx$variable)),]
    dfx = rbind(dfx,pxx)
  }
  df <- ddply(dfx, .(path,type), transform, cum.perc = Reduce('+', list(perc/2,cumsum(c(0,head(perc,-1))))))
  df[,c('type', 'variable','n', 'value', 'rest','perc','delta','testres','pv','pw','pv2','path','cum.perc')]
}


get_stats_ex = function(x){

  # n.p = sum(x$type=="PTEN_LOSS")
  na.p = sum(is.na(x$value[x$type=="PTEN_LOSS"]))
  s.p = s(x$value[x$type=="PTEN_LOSS"])[1:6]; names(s.p)=paste0("pt.",c("min","1qt","median","mean","3qt","max"))

  # n.w = sum(x$type=="WT")
  na.w = sum(is.na(x$value[x$type=="WT"]))
  s.w = s(x$value[x$type=="WT"])[1:6]; names(s.w)=paste0("wt.",c("min","1qt","median","mean","3qt","max"))

  w = (with(x,wilcox.test(value[type=="PTEN_LOSS"], value[type=="WT"])$p.value))
  k = with(x,ks.test(value[type=="PTEN_LOSS"], value[type=="WT"])$p.value)
  # t =  with(x,t.test(value[type=="PTEN_LOSS"], value[type=="WT"])$p.value)
  #
  return(c(
    # "N_pten" = n.p,
    "N_NA_pten" = na.p,
    s.p,
    # # "N_wt" = n.w,
    "N_NA_wt" = na.w,
    s.w,
    'w'=w,
    'k'=k
    # ,'t'=t
  ))
}

get_stats_ex_boolean = function(x){

  # n.p = sum(x$type=="PTEN_LOSS")
  na.p = sum(is.na(x$value[x$type]))
  s.p = s(x$value[x$type])[1:6]; names(s.p)=paste0("pt.",c("min","1qt","median","mean","3qt","max"))

  # n.w = sum(x$type=="WT")
  na.w = sum(is.na(x$value[!x$type]))
  s.w = s(x$value[!x$type])[1:6]; names(s.w)=paste0("wt.",c("min","1qt","median","mean","3qt","max"))

  w = (with(x,wilcox.test(value[type], value[!type])$p.value))
  k = with(x, ks.test(    value[type], value[!type])$p.value)
  # t =  with(x,t.test(value[type=="PTEN_LOSS"], value[type=="WT"])$p.value)
  #
  return(c(
    # "N_pten" = n.p,
    "N_NA_pten" = na.p,
    s.p,
    # # "N_wt" = n.w,
    "N_NA_wt" = na.w,
    s.w,
    'w'=w,
    'k'=k
    # ,'t'=t
  ))
}

get_stats_ex_boolean_v2 = function(x){
  na.p = sum(is.na(x$value[x$type]))
  na.w = sum(is.na(x$value[!x$type]))
  se = summarySE_v2(x, measurevar = 'value', groupvars = 'type', na.rm = T);
  s.p = se[which(se$type), 2:5]
  s.w = se[which(!se$type), 2:5]
  names(s.p)=paste0("case.",names(s.p))
  names(s.w)=paste0("cntr.",names(s.w))
  w = with(x,wilcox.test(value[type], value[!type])$p.value)
  k = with(x, ks.test(    value[type], value[!type])$p.value)
  return(unlist(c( s.p, s.w, 'w'=w, 'k'=k)) )
}

global_summary_expr_class.ONE_TAIL = function(res,pl, lvls=c("he","me","le","not_expr"), cl=c("PTEN_LOSS","WT")){
  dfx = rbind()
  for(i in names(res)){
    # print(i)
    a = res[[i]]
    sa = c(apply(a[,paste0(lvls,".",cl[1])],2,sum),'type'=cl[1]); names(sa)[1:length(lvls)]=lvls
    sb = c(apply(a[,paste0(lvls,".",cl[2])],2,sum),'type'=cl[2]); names(sb)[1:length(lvls)]=lvls

    px = rbind(sa,sb); px=as.data.frame(px)
    for(z in 1:length(lvls)) px[,z] = as.numeric(as.character(px[,z]))
    px$n = c( sum( px[1,1:length(lvls)] ),sum( px[2,1:length(lvls)] ))

    pxx = melt(px, id.vars=c("type","n"))
    pxx$perc     = (pxx$value/pxx$n)*100
    pxx$variable = factor(pxx$variable, levels=lvls, ordered=T)
    pxx$type     = factor(pxx$type, levels=cl )
    pxx$rest    = pxx$n - pxx$value
    pxx = ddply(pxx, .(variable), mutate, delta=perc[which(type==cl[1])]-perc[which(type==cl[2])])

    pxx$pv =1
    pxx$pw =1
    for( j in pxx$variable ){
      u = subset(pxx, variable==j)[c("value","rest")]

      if(unique(subset(pxx, variable==j)$delta)!=0){
        if(unique(subset(pxx, variable==j)$delta)>0) f = fisher.test(u, alt='greater') else f = fisher.test(u, alt='less')
        pxx$pv[which(pxx$variable==j)]=f$p
      }

      if(u[1,"value"]>0 & u[2,"value"]>0 ){
        pw = power.fisher.test(u[1,"value"]/(u[1,"value"]+u[1,"rest"]), u[2,"value"]/(u[2,"value"]+u[2,"value"]),
                               u[1,"value"]+u[1,"rest"], u[2,"value"]+u[2,"rest"],
                               alpha=.05, nsim=10, alternative="two.sided")
        pxx$pw[which(pxx$variable==j)]=pw
      }
    }

    testres = c("D","E")
    pxx$testres = testres[as.numeric(pxx$delta>0)+1]
    pxx$pv2 = paste(as.character(round(pxx$pv,2)),"(", pxx$testres ,")", sep="",coll="")
    pxx$pv2[which(pxx$pv>0.05 | pxx$type==cl[2])]=""
    pxx$path=i
    pxx=pxx[rev(order(pxx$variable)),]
    dfx = rbind(dfx,pxx)
  }
  df <- ddply(dfx, .(path,type), transform, cum.perc = Reduce('+', list(perc/2,cumsum(c(0,head(perc,-1))))))
  df[,c('type', 'variable','n', 'value', 'rest','perc','delta','testres','pv','pw','pv2','path','cum.perc')]
}

get_coord=function(x, cx){
  cx    = subset(cx,Symbol==x$symbol )
  xx    = as.numeric(unlist(strsplit(x$exons,":")))
  x$chr             = paste0('chr',cx$Chromosome[1])
  x$strand          = cx$Strand[1]
  x$exonStart_0base = paste(subset(cx, Exon%in%xx)$Chr_Start, collapse="|")
  x$exonEnd         = paste(subset(cx, Exon%in%xx)$Chr_Stop,  collapse="|")

  x$upstreamES      = ifelse(x$from_exon=='null', NA, subset(cx, Exon==x$from_exon)$Chr_Start)
  x$upstreamEE      = ifelse(x$from_exon=='null', NA, subset(cx, Exon==x$from_exon)$Chr_Stop)
  x$downstreamES    = ifelse(x$to_exon=='null', NA, subset(cx, Exon==x$to_exon)$Chr_Start)
  x$downstreamEE    = ifelse(x$to_exon=='null', NA, subset(cx, Exon==x$to_exon)$Chr_Stop)

  return(x)
}

get_coord_ME=function(x, cx){
  cx    = subset(cx,Symbol==x$symbol )
  xx    = as.numeric(unlist(strsplit(x$exons,"|")))
  x$chr             = paste0('chr',cx$Chromosome[1])
  x$strand          = cx$Strand[1]
  x$exonStart_0base = paste(subset(cx, Exon%in%xx)$Chr_Start, collapse="|")
  x$exonEnd         = paste(subset(cx, Exon%in%xx)$Chr_Stop,  collapse="|")

  x$upstreamES      = ifelse(x$from_exon=='null', NA, subset(cx, Exon==x$from_exon)$Chr_Start)
  x$upstreamEE      = ifelse(x$from_exon=='null', NA, subset(cx, Exon==x$from_exon)$Chr_Stop)
  x$downstreamES    = ifelse(x$to_exon=='null', NA, subset(cx, Exon==x$to_exon)$Chr_Start)
  x$downstreamEE    = ifelse(x$to_exon=='null', NA, subset(cx, Exon==x$to_exon)$Chr_Stop)

  return(x)
}
#
# get_coord_ES = function(x, cx){
#
#   res = rbind()
#   for(i in 1:nrow(x)){
#     print(i)
#     xx = unlist(strsplit(x$exons[i],"\\:"))
#     start=xx[1]
#     end=xx[length(xx)]
#
#     res = rbind(res, c(
#       'ID'	= x$as_id[i]
#       ,'GeneID'	= x$symbol[i]
#       ,'geneSymbol'= x$symbol[i]
#
#       ,'chr'	= paste0('chr',subset(cx, Symbol==x$symbol[i] & Exon==start)$Chromosome)
#       ,'strand'	= subset(cx, Symbol==x$symbol[i] & Exon==start)$Strand
#       ,'exonStart_0base'	= subset(cx, Symbol==x$symbol[i] & Exon==start)$Chr_Start
#       ,'exonEnd' = subset(cx, Symbol==x$symbol[i] & Exon==end)$Chr_Stop
#
#       ,'upstreamES' = ifelse(x$from_exon[i]=='null', NA, subset(cx, Symbol==x$symbol[i] & Exon==x$from_exon[i])$Chr_Start)
#       ,'upstreamEE' = ifelse(x$from_exon[i]=='null', NA, subset(cx, Symbol==x$symbol[i] & Exon==x$from_exon[i])$Chr_Stop)
#
#       ,'downstreamES' = ifelse(x$to_exon[i]=='null', NA, subset(cx, Symbol==x$symbol[i] & Exon==x$to_exon[i])$Chr_Start)
#       ,'downstreamEE' = ifelse(x$to_exon[i]=='null', NA, subset(cx, Symbol==x$symbol[i] & Exon==x$to_exon[i])$Chr_Stop)
#       # , x[i,]
#       # ,'ID'	=x$as_id[i]
#       # ,'IC_SAMPLE_1'= x$pt.median[i]
#       # ,'SC_SAMPLE_1'= x$wt.median[i]
#       # ,'IncFormLen' = NA
#       # ,'SkipFormLen' =NA
#       #
#       # ,'PValue' = x$k[i]
#       # ,'FDR' = x$FDR.ks[i]
#       # ,'IncLevel1' = x$pt.median[i]
#       # ,'IncLevel2'= x$wt.median[i]
#       ,'IncLevelDifference' = x$delta.median[i]
#     ))
#   }
#   return(res)
# }
#
# get_coord_ES_list = function(x, cx){
#     xx = unlist(strsplit(x$exons,"\\:"))
#     start=xx[1]
#     end=xx[length(xx)]
#
#     return( c(
#       'ID'	= x$as_id
#       ,'GeneID'	= x$symbol
#       ,'geneSymbol'= x$symbol
#
#       ,'chr'	= paste0('chr',subset(cx, Symbol==x$symbol & Exon==start)$Chromosome)
#       ,'strand'	= subset(cx, Symbol==x$symbol & Exon==start)$Strand
#       ,'exonStart_0base'	= subset(cx, Symbol==x$symbol & Exon==start)$Chr_Start
#       ,'exonEnd' = subset(cx, Symbol==x$symbol & Exon==end)$Chr_Stop
#
#       ,'upstreamES' = ifelse(x$from_exon=='null', NA, subset(cx, Symbol==x$symbol & Exon==x$from_exon)$Chr_Start)
#       ,'upstreamEE' = ifelse(x$from_exon=='null', NA, subset(cx, Symbol==x$symbol & Exon==x$from_exon)$Chr_Stop)
#
#       ,'downstreamES' = ifelse(x$to_exon=='null', NA, subset(cx, Symbol==x$symbol & Exon==x$to_exon)$Chr_Start)
#       ,'downstreamEE' = ifelse(x$to_exon=='null', NA, subset(cx, Symbol==x$symbol & Exon==x$to_exon)$Chr_Stop)
#       # , x[i,]
#       # ,'ID'	=x$as_id[i]
#       # ,'IC_SAMPLE_1'= x$pt.median[i]
#       # ,'SC_SAMPLE_1'= x$wt.median[i]
#       # ,'IncFormLen' = NA
#       # ,'SkipFormLen' =NA
#       #
#       # ,'PValue' = x$k[i]
#       # ,'FDR' = x$FDR.ks[i]
#       # ,'IncLevel1' = x$pt.median[i]
#       # ,'IncLevel2'= x$wt.median[i]
#        ,'IncLevelDifference' = x$delta.median
#     ))
# }
#

from_MATS_to_RNAmotifs = function(exons, fname){
  #row_id;second_id;chrom;strand;upstream_exon_end_position;exon_start_position;exon_end_position;dwstream_exon_start_position;diRank
  exons$DI = 0
  exons$DI[exons$type=="Silenced"]=(-1)
  exons$DI[exons$type=="Enhanced"]=1
  cn = c("ID",'chr','strand','upstreamEE','exonStart_0base','exonEnd','downstreamES',"DI")
  res = cbind(1:nrow(exons), exons[,cn])
  write.table(res, file=fname, col.names = F, row.names = F, quote = F, sep=";")
}

bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)

  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }

  if(length(df)<3){
    stop("File has less than 3 columns")
  }

  header <- c('chr','start','end','id','score','strand')
  names(df) <- header[1:length(names(df))]

  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }

  library("GenomicRanges")

  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}

PFAM_GO_binomial_test = function(x, case, alt ='g'){
  require(plyr)
  counts = table(case$pfam.acc)

  x$n_case = counts[x$pfam.acc]
  x$t_case = nrow(case)
  x        = x[which(!is.na(x$n_case)),]

  x = ddply(x, .(goterm, desc, go, level, t_exons, t_case), summarise
            , n_exons=sum(n_exons, na.rm = T)
            , n_case=sum(n_case, na.rm =T)
            , pfam.acc = paste(pfam.acc, collapse = "|") )

  x$exp     = with(x, n_exons/t_exons )
  x$obs    = with(x, n_case/t_case )
  # x        = x[order(x$n_exons, decreasing = T),]
  x$binom = apply(x[,c('n_case','t_case','exp')], 1, function(x) binom.test(x[1],x[2],p=x[3], alternative = alt)$p.value)
  x = dlply(x, ~go)

  x = lapply(x, function(y){
    y$BF  = p.adjust(y$binom, 'bonferroni')
    y$FDR = p.adjust(y$binom, 'fdr')
    return(y)
  })

  return(as.data.frame(do.call(rbind,x)))
}

read.gmt.file = function(pathMsigDbFile) {
  inputFile <- pathMsigDbFile
  con  <- file(inputFile, open = "r")

  c = 1
  pathway.list <- vector(mode="list",length=0)
  while (length(oneLine <- readLines(con, n = 1, warn = FALSE)) > 0)
  {
    myVector <- do.call("rbind",strsplit(oneLine, "\t"))
    t = vector(mode="list",length=1)
    t[[1]] = myVector[3:length(myVector)]
    names(t) = myVector[1]
    pathway.list = c(pathway.list,t)
    c = c+1
  }

  close(con)
  return(pathway.list)
}

GO_slimmer = function(x){
  require(GSEABase)
  fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
  slim <- getOBOCollection(fl)
  # BP = goSlim(GOCollection(na.omit(unique(x$pfam.GO_BP.id))), slim, "BP"); BP = subset(BP, Count!=0); BP = BP[order(BP$Percent, decreasing = T),]
  # MF = goSlim(GOCollection(na.omit(unique(x$pfam.GO_MF.id))), slim, "MF"); MF = subset(MF, Count!=0); MF = MF[order(MF$Percent, decreasing = T),]
  # CC = goSlim(GOCollection(na.omit(unique(x$pfam.GO_CC.id))), slim, "CC"); CC = subset(CC, Count!=0); CC = CC[order(CC$Percent, decreasing = T),]
  BP = goSlim(GOCollection(na.omit((x$pfam.GO_BP.id))), slim, "BP"); BP = subset(BP, Count!=0); BP = BP[order(BP$Percent, decreasing = T),]
  MF = goSlim(GOCollection(na.omit((x$pfam.GO_MF.id))), slim, "MF"); MF = subset(MF, Count!=0); MF = MF[order(MF$Percent, decreasing = T),]
  CC = goSlim(GOCollection(na.omit((x$pfam.GO_CC.id))), slim, "CC"); CC = subset(CC, Count!=0); CC = CC[order(CC$Percent, decreasing = T),]
  list(BP=BP,MF=BP,CC=BP)
}
# envirName <-paste("GO",onto,"ANCESTOR",sep="")
# envi<-eval(parse(text = envirName))
#
# le = lapply_pb(l, function(x,e) AnnotationDbi::mget(x,e,ifnotfound=NA), e=envi)
#
# AncestorsLst[[i]]<-unique(AncestorsLst[[i]]) [-1]}
#
#
# GO_slimmer=function(x){
#   require(GSEABase)
#   fl <- system.file("extdata", "goslim_plant.obo", package="GSEABase")
#   slim <- getOBOCollection(fl)
#   # BP = goSlim(GOCollection(na.omit(unique(x$pfam.GO_BP.id))), slim, "BP"); BP = subset(BP, Count!=0); BP = BP[order(BP$Percent, decreasing = T),]
#   # MF = goSlim(GOCollection(na.omit(unique(x$pfam.GO_MF.id))), slim, "MF"); MF = subset(MF, Count!=0); MF = MF[order(MF$Percent, decreasing = T),]
#   # CC = goSlim(GOCollection(na.omit(unique(x$pfam.GO_CC.id))), slim, "CC"); CC = subset(CC, Count!=0); CC = CC[order(CC$Percent, decreasing = T),]
#   BP = goSlim(GOCollection(na.omit((x$pfam.GO_BP.id))), slim, "BP"); BP = subset(BP, Count!=0); BP = BP[order(BP$Percent, decreasing = T),]
#   MF = goSlim(GOCollection(na.omit((x$pfam.GO_MF.id))), slim, "MF"); MF = subset(MF, Count!=0); MF = MF[order(MF$Percent, decreasing = T),]
#   CC = goSlim(GOCollection(na.omit((x$pfam.GO_CC.id))), slim, "CC"); CC = subset(CC, Count!=0); CC = CC[order(CC$Percent, decreasing = T),]
#   list(BP=BP,MF=BP,CC=BP)
# }
#
# x=goSlim(GOCollection(PGO.MF[,2]), slim, "MF")
