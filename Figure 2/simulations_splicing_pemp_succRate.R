#!/usr/bin/env Rscript --vanilla

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%% Splicing analysis : SIMULATIONS FOR pemp AND success rate %%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

message('[****] SCRIPT: simulations_splicing_pemp_succRate.R')

print('[****] SCRIPT: simulations_splicing_pemp_succRate.R')


# Set working path : ----

gene = 'FOXA1'

print(paste0('STRATA: ',gene))

work_dir     = 'Automatic_splicing_analysis_results'
STRATA_code  = paste0(gene,'_HE') 
work_dir     = paste0(work_dir,STRATA_code)
Fig_out_dir  = paste0(work_dir,'/Figures/splicing_analysis_',STRATA_code)
Data_out_dir = paste0(work_dir,'/Rdata/splicing_analysis_',STRATA_code)


# Functions =========
get_wk_test <- function(x){ 
  w = with(x, wilcox.test(value[type], value[!type])$p.value)
  k = with(x, ks.test(    value[type], value[!type])$p.value)
  f = with(x, fligner.test(list(value[type],value[!type]))$p.value)
  return( c('w'=w,'k'=k,'f'=f) )
}


execute_bootstrap <- function(samples, data, type='sample.size'){
  require(plyr)
  if(type=='p.empirical'){
    data$type <- samples$type[match(as.character(data$variable), as.character(samples$variable)) ]
  }else if(type=='sample.size'){
    data <- subset(data, as.character(variable)%in%samples)
  }
  rx  = ddply(data, .(as_id), get_wk_test)
  return(rx)
}



# LOADINGS =================

message('[*] Loading variables ...')

library(snowfall)
library(plyr)

dex = readRDS(paste0(Data_out_dir,'/PRAD_dex_',STRATA_code,'.rds'))
ex  = readRDS(paste0(Data_out_dir,'/PRAD_selected_exons_',STRATA_code,'.rds'))
ex$type=ex$cases

ex = subset(ex, as_id%in%subset(dex, w<=0.05 | k<=0.05 )$event_id )

tmp = subset(ex, as_id=='exon_skip_52')[,c('variable',"type")]
patients = unique(tmp)

message('[*] Variables loaded, start bootstrapping ...')

# P.EMPIRICAL and Sample Size ==========

#********************
set.seed(30580)
NCORES = 4
NSIM   = 1000
#********************

samples <- 1:NSIM
samples <- lapply(samples, function(x,y){ y$type=sample(y$type); return(y) } , y=patients )

ex$type        = NULL
ex$splice_type = NULL


functionsToCluster = as.character(lsf.str(envir = .GlobalEnv))
sfInit(parallel=TRUE, cpus=NCORES, type="SOCK")
sfExport(list = c(functionsToCluster, "ex"))
sfLibrary(parallel)
sfLibrary(plyr)
sfClusterSetupRNG()
l <- sfLapply(samples, function(x){ tryCatch( execute_bootstrap(  sample=x
                                                                  , data = ex
                                                                  , type='p.empirical')
                                              , error = function(e) return(e)) })

w =  sfLapply(l, "[[", 2)

w = do.call(cbind.data.frame, w)

colnames(w) = NULL
rownames(w) = l[[1]][,1]

k =  sfLapply(l, "[[", 3)

k = do.call(cbind.data.frame, k)

colnames(k) = NULL
rownames(k) = l[[1]][,1]

f =  sfLapply(l, "[[", 4)

f = do.call(cbind.data.frame, f)

colnames(f) = NULL
rownames(f) = l[[1]][,1]

rm(ex)

sfStop()

rm(ex)

saveRDS(samples, file = paste0(Data_out_dir,'/EX_bootstrap_pemp_samples_',STRATA_code,'.rds'))
saveRDS(w, file = paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_wilcoxon_',STRATA_code,'.rds'))
saveRDS(k, file = paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_kolgomorov_',STRATA_code,'.rds'))
saveRDS(f, file = paste0(Data_out_dir,'/EX_bootstrap_pemp_pv_fligner_',STRATA_code,'.rds'))

# Sample size 
pheno <- list()
pheno[[1]]  = as.character(subset(patients, type)$variable)
pheno[[2]]  = as.character(subset(patients, !type)$variable)
sample.size = length(pheno[[1]])

samples = lapply(1:NSIM, function(i) {
  sample_barcode = vector(mode = "numeric", length = 2*sample.size)
  sample_barcode[1:sample.size] = pheno[[1]]
  sample_barcode[(sample.size+1):length(sample_barcode)] = sample(pheno[[2]], sample.size)
  return(sample_barcode)
})

ex  = readRDS(paste0(Data_out_dir,'/PRAD_selected_exons_',STRATA_code,'.rds'))
ex$type=ex$cases

ex = subset(ex, as_id%in%subset(dex, w<=0.05 | k<=0.05 )$event_id )
ex$splice_type = NULL


functionsToCluster = as.character(lsf.str(envir = .GlobalEnv))
sfInit(parallel=TRUE, cpus=NCORES, type="SOCK")
sfExport(list = c(functionsToCluster, "ex"))
sfLibrary(parallel)
sfLibrary(plyr)
sfClusterSetupRNG()
l <- sfLapply(samples, function(x){ tryCatch( execute_bootstrap(  sample=x
                                                                  , data = ex
                                                                  , type='sample.size')
                                              , error = function(e) return(e)) })


sfStop()

w =  lapply(l, "[[", 2)
w = do.call(cbind.data.frame, w)
colnames(w) = NULL
rownames(w) = l[[1]][,1]

k =  lapply(l, "[[", 3)
k = do.call(cbind.data.frame, k)
rownames(k) = l[[1]][,1]
colnames(k) = NULL

f =  lapply(l, "[[", 4)
f = do.call(cbind.data.frame, f)
rownames(f) = l[[1]][,1]
colnames(f) = NULL

saveRDS(samples, file = paste0(Data_out_dir,'/EX_bootstrap_SR_samples_',STRATA_code,'.rds'))
saveRDS(w, file = paste0(Data_out_dir,'/EX_bootstrap_SR_pv_wilcoxon_',STRATA_code,'.rds'))
saveRDS(k, file = paste0(Data_out_dir,'/EX_bootstrap_SR_pv_kolgomorov_',STRATA_code,'.rds'))
saveRDS(f, file = paste0(Data_out_dir,'/EX_bootstrap_SR_pv_fligner_',STRATA_code,'.rds'))
