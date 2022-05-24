binomial_test_AS_category=function(case, control, fname=NULL){
  require(plyr)
  require(dplyr)
  spty = c('AA', 'AD', 'AP', 'AT', 'ES', 'ME', 'RI')
  t = as.data.frame(table(control$splice_type)[spty]); colnames(t)=c('splice_type','n_exons')
  t$splice_type=spty
  if(length(which(is.na(t$n_exons)))>0) t$n_exons[which(is.na(t$n_exons))]=0
  t$t_exons    = nrow(control)
  t$background = t$n_exons/t$t_exons

  t2 = as.data.frame(table(case$splice_type)[spty]); colnames(t2)=c('splice_type','n_cases')
  t2$splice_type=spty
  if(length(which(is.na(t2$n_cases)))>0) t2$n_cases[which(is.na(t2$n_cases))]=0
  t2$t_cases    = nrow(case)
  t2$p_cases = t2$n_cases/t2$t_cases
  if(length(which(is.na(t2$n_cases)))>0) t2$p_cases[which(is.na(t2$n_cases))]=0

  t = cbind(t, t2[,2:4])
  t$direction = "Enriched"
  t$direction[which(t$p<(t$background) )]="Depleted"
  t = ddply(t, .(splice_type), dplyr::mutate, binom=binom.test(n_cases, t_cases, background )$p.value)
  t$BF = p.adjust(t$binom, "bonferroni")
  t$FDR = p.adjust(t$binom, "fdr")

  m  = melt(t, id.vars = c('splice_type'), measure.vars =c('background','p_cases'))
  m2 = melt(t, id.vars = c('splice_type'), measure.vars =c('n_exons','n_cases'))
  m$label = m2$value

  if(!is.null(fname)) pdf(file=fname,paper='a4', useDingbats = F)
  p1 = ggplot(m, aes(x=splice_type, y=100*value, fill=variable))+geom_bar(stat="identity", position = 'dodge')+
    geom_text(aes(label=label), angle=90)+
    theme_classic()+xlab('')+ylab('Exons(%)')#+theme(legend.position = "none")
  print(p1)
  if(!is.null(fname)) dev.off()

  t
}

fisher_test_AS_category=function(case, control, fname=NULL){
  require(plyr)
  require(dplyr)
  spty = c('ES', 'A3', 'A5', 'IR', 'MEX')
  t = as.data.frame(table(control$event_type)[spty]); colnames(t)=c('splice_type','n_exons')
  t$splice_type=spty
  if(length(which(is.na(t$n_exons)))>0) t$n_exons[which(is.na(t$n_exons))]=0
  t$t_exons    = nrow(control)
  t$background = t$n_exons/t$t_exons

  t2 = as.data.frame(table(case$event_type)[spty]); colnames(t2)=c('splice_type','n_cases')
  t2$splice_type=spty
  if(length(which(is.na(t2$n_cases)))>0) t2$n_cases[which(is.na(t2$n_cases))]=0
  t2$t_cases    = nrow(case)
  t2$p_cases = t2$n_cases/t2$t_cases
  if(length(which(is.na(t2$n_cases)))>0) t2$p_cases[which(is.na(t2$n_cases))]=0

  t = cbind(t, t2[,2:4])
  t$direction = "Enriched"
  t$direction[which(t$p<(t$background) )]="Depleted"
  t = ddply(t, .(splice_type), dplyr::mutate, pv=fisher.test(matrix(c(n_cases, t_cases-n_cases, n_exons, t_exons-n_exons), nrow=2, byrow = T))$p.value)
  t$BF = p.adjust(t$pv, "bonferroni")
  t$FDR = p.adjust(t$pv, "fdr")

  m  = melt(t, id.vars = c('splice_type',"pv"), measure.vars =c('background','p_cases'))
  m2 = melt(t, id.vars = c('splice_type',"pv"), measure.vars =c('n_exons','n_cases'))
  m$label = m2$value

  if(!is.null(fname)) pdf(file=fname,paper='a4', useDingbats = F)
  p=ggplot(m, aes(x=splice_type, y=100*value, fill=variable))+geom_bar(stat="identity", position = 'dodge')+
          geom_text(aes(label=label),position = position_dodge(0.9), angle=45)+
          theme_classic()+xlab('')+ylab('Exons(%)')#+theme(legend.position = "none")
  tmp = subset(m, pv<=0.05 & variable=='p_cases' )
  p=p+geom_text(data=tmp, aes(y=50*value),label="*",position = position_dodge(0.9), size=5)
  print(p)
  if(!is.null(fname)) dev.off()

  t
}

GO_binomial_test = function(x, case, control=NULL, alt ='g'){
  require(plyr)
  if(!is.null(control)){
    # counts = table(control$symbol)
    #
    # x$n_cntrl = counts[match(x$symbol, names(counts))] #(m)
    # x$t_cntrl = sum(control$symbol%in%x$symbol) #(N)

    counts = table(control$symbol)

    x$n_cntrl = counts[match(x$symbol, names(counts))] #(m)
    # db$t_cntrl = sum(dex$symbol%in%db$symbol) #(N)
    # db$t_cntrl = nrow(coord)
    x$t_cntrl = nrow(control)

  }

  counts = table(case$symbol)

  x$n_case = counts[match(x$symbol, names(counts))] # (k)
  x$t_case = sum(case$symbol%in%x$symbol) #(n)

  x        = x[which(!is.na(x$n_case)),]

  x = ddply(x, .(goterm, term, go, level, t_case, t_cntrl), summarise
            , n_case=sum(n_case, na.rm =T)
            , n_cntrl=sum(n_cntrl, na.rm = T)
            , n_genes = length(unique(symbol))
            , genes = paste(unique(symbol), collapse = "|") )

  x$exp     = with(x, n_cntrl/t_cntrl ) #(m/N)
  x$obs    = with(x, n_case/t_case )



  x$binom = apply(x[,c('n_case','t_case','exp')], 1, function(y) binom.test(y[1],y[2],p=y[3], alternative = alt)$p.value)

  # x = dlply(x, ~go)

  # x = lapply(x, function(y){
  # y=ddply(y, .(level), mutate,
  #   BF  = p.adjust(binom, 'bonferroni', n=length(level)),
  #   FDR = p.adjust(binom, 'fdr', n=length(level))
  # )
  # y$BF  = p.adjust(y$binom, 'bonferroni')
  # y$FDR = p.adjust(y$binom, 'fdr')
  # # return(y)
  # })
  x$BF  = p.adjust(x$binom, 'bonferroni')
  x$FDR = p.adjust(x$binom, 'fdr')
  x = x[,c('goterm', 'term', 'go', 'level', 'n_case','t_case','n_cntrl', 't_cntrl', 'obs','exp', 'binom','BF','FDR','n_genes','genes')]
  return(x)
  # return(as.data.frame(do.call(rbind,x)))
}

KEGG_binomial_test = function(x, case, control=NULL, alt ='g'){
  require(plyr)
  if(!is.null(control)){
    # counts = table(control$symbol)
    #
    # x$n_cntrl = counts[match(x$symbol, names(counts))] #(m)
    # x$t_cntrl = sum(control$symbol%in%x$symbol) #(N)
    counts = table(control$symbol)

    x$n_cntrl = counts[match(x$symbol, names(counts))] #(m)
    # db$t_cntrl = sum(dex$symbol%in%db$symbol) #(N)
    # db$t_cntrl = nrow(coord)
    x$tmp = counts[match(x$alias, names(counts))] #(m)
    x$n_cntrl[which(!is.na(x$tmp))] = x$tmp[which(!is.na(x$tmp))]
    x$tmp = NULL

    x$t_cntrl = nrow(control) #(N)

  }

  counts = table(case$symbol)

  x$n_case = counts[match(x$symbol, names(counts))] # (k)
  x$t_case = sum(case$symbol%in%x$symbol) #(n)
  x        = x[which(!is.na(x$n_case)),]

  x = ddply(x, .(term, t_case, t_cntrl), summarise
            , n_case=sum(n_case, na.rm =T)
            , n_cntrl=sum(n_cntrl, na.rm = T)
            , n_genes = length(unique(symbol))
            , genes = paste(unique(symbol), collapse = "|") )

  x$exp     = with(x, n_cntrl/t_cntrl ) #(m/N)
  x$obs    = with(x, n_case/t_case )



  x$binom = apply(x[,c('n_case','t_case','exp')], 1, function(y) binom.test(y[1],y[2],p=y[3], alternative = alt)$p.value)

  x$BF  = p.adjust(x$binom, 'bonferroni')
  x$FDR = p.adjust(x$binom, 'fdr')
  x = x[,c('term', 'n_case','t_case','n_cntrl', 't_cntrl', 'obs','exp', 'binom','BF','FDR','n_genes','genes')]
  return(x)
  # return(as.data.frame(do.call(rbind,x)))
}

get_gprofiler =function( set, id, background="", collections=c("GO")){
  pro = gprofiler( set,
                   custom_bg = background,
                   max_p_value = 0.05,
                   ordered_query =T,
                   significant = T,
                   correction_method='gSCS',
                   min_set_size = 10,
                   # evcodes=T,
                   # exclude_iea = T,
                   hier_filtering ='moderate',
                   src_filter = collections
                   # png_fn = "Prostate/Figures/Goprofiler.AS.png"
                   # include_graph=T
                   # src_filter = c("GO",'KEGG',"REAC","CORUM","OMIM",  "HP")
  )

  write.csv(pro, file=paste0("Prostate/AS_events/GOprofiler.",id,".csv"), row.names = F)
  em = c("GO.ID","Description","p.Val","FDR","Phenotype","Genes")
  map = cbind( pro[,c('term.id','term.name','p.value','p.value')],1, pro$intersection)
  colnames(map) = em
  write.table(map, file=paste0("Prostate/AS_events/GOprofiler.",id,".map.txt"), row.names = F, col.names = F, quote = F, sep='\t')
  pro
}

