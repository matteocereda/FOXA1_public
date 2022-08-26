# RNAmotifs ====


# source("/Volumes/LaCie/ALBERTO/repo_FOXA1/Figure 4/config_FOXA1_analysis_2022_05.R")
repo = '/Volumes/Prule/repo/FOXA1_public/'
source('Figure_4/config_FOXA1_analysis_2022_05.R')

pr_folder = "input_FOXA1_regulated_vs_FOXA1_unregulated_mean_09_or_01"

# Paths ----

PATH_RNAMOTIFS     = 'sources/RNAmotifs/' #"/Volumes/LaCie/ALBERTO/repo_FOXA1/RNAmotifs"
PATH_TO_DB         = p(repo, 'sources/RNAmotifs/mCrossPWMs/') #"/Volumes/LaCie/ALBERTO/repo_FOXA1/Rdata/mCrossPWMs"        # PWMs are stored as [cell_line].[RBP].[PWM_version].txt 
PATH_TO_MACROAPE   = 'sources/RNAmotifs/Macro-ape/' # "/Volumes/LaCie/ALBERTO/repo_FOXA1/Macro-ape"
PATH_TO_PEAKS      = 'sources/RNAmotifs/Nobby/' #"/Volumes/LaCie/ALBERTO/repo_FOXA1/Rdata/Nobby/"
# PATH_ALL_DATASET   = '~/PRAD_dex_FOXA1_HE_ALL_ES_events_noMore500.rds'
PATH_DEX_DATASET   = 'Rdata/dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds'  #'/Volumes/LaCie/FOXA1_public/Rdata/dexVAR_EXON_SKIPPING_only_with_pfam_and_NMD_FOXA1_HE.rds'
PATH_DEX_FOXA1     = 'sources/RNAmotifs/NEW_differentially_expressed_lists.Rdata' #"/Volumes/LaCie/ALBERTO/repo_FOXA1/Rdata/NEW_differentially_expressed_lists.Rdata"
PATH_MCROSS_SCORES = 'sources/RNAmotifs/score_list_mCross.RData' #"/Volumes/LaCie/ALBERTO/repo_FOXA1/Rdata/score_list_mCross.RData"

PATH_RES_RNAMOTIFS = 'sources/RNAmotifs/results/' #p(PATH_RNAMOTIFS, "/results/",pr_folder,"/")
PATH_INP_RNAMOTIFS = p('sources/RNAmotifs/',pr_folder,'.txt') # p(PATH_RNAMOTIFS, "/input/",pr_folder,".txt") # _input_RNAmotifs
PATH_PR_FOLDER     = p(repo, 'sources/RNAmotifs/m3_light_results/',pr_folder,'/')  #p(PATH_RNAMOTIFS, "/m3_light/results/", pr_folder,"/")

# Parameters ----

ncpus_Macroape = 4
PWM_start      = 2         # start core positions (PWMs from mCROSS)
PWM_end        = 9         # end   core positions (PWMs from mCROSS)
ncol_PWMs      = 8         # number of columns of PWMs from RNAmotifs 
eCLIP_weights  = T
cl             = "HepG2"
percentile     = 75
pemp_thresh    = 5e-4

cell_line = cl

PATH_PWM_FOLDER    = p("PWM_", ncol_PWMs, "/")
PATH_RES           = p("./",pr_folder,"_",ncol_PWMs,"_columns/")
# PATH_RES           = p(repo, PATH_RNAMOTIFS,  pr_folder,"_",ncol_PWMs,"_columns/")

# Data preparation ----

d2            = readRDS(PATH_DEX_DATASET)
d2$delta.mean = d2$delta.mean/100


d2_reg            = d2 %>% subset(FOXA1_sign)

reg  = with(d2_reg, data.frame(event_id  = event_id %>% cha %>% strsplit("_") %>% sapply(tail,1) %>% num,
                           event_chr = event_chr,
                           strand    = strand %>% cha,
                           coord     = event_coordinates %>% cha %>%  strsplit(":") %>% sapply(function(y){num(y[2:5])}) %>% t,
                           diRank    = (delta.mean>0)*2-1))

d2_cntr            = d2 %>% subset(!FOXA1_sign & (overall.mean<=0.1 | overall.mean>=0.9))

cntr <- with(d2_cntr, data.frame(event_id  = event_id %>% cha %>% strsplit("_") %>% sapply(tail,1) %>% num,
                            event_chr = paste0("chr", event_chr),
                            strand    = strand %>% cha,
                            coord     = event_coordinates %>% cha %>%  strsplit(":") %>% sapply(function(y){num(y[2:5])}) %>% t,
                            diRank    = 0))

# enhanced/silenced (analogue of dIRank)
#  1 --> FOXA1_sign == TRUE & delta.mean > 0
#  0 --> FOXA1_sign == FALSE 
# -1 --> FOXA1_sign == TRUE & delta.mean < 0


# Each row must be formatted as:
# row_id; second_id; chrom; strand; up_exon_end; exon_start; exon_end; dw_exon_start; diRank

out = data.frame(id = 1:(nrow(reg)+nrow(cntr)), 
                 rbind(reg, cntr))

write.table(out, PATH_INP_RNAMOTIFS, quote = F, sep = ";", row.names = F, col.names = F)




# # RNAmotifs ----
# 
setwd(PATH_RNAMOTIFS) # on cluster
sysetm("./RNAmotifs.sh /hpcnfs/data/cgb/RNAmotifs_devel_Alberto input_FOXA1_regulated_vs_FOXA1_unregulated_mean_09_or_01 input_FOXA1_regulated_vs_FOXA1_unregulated_mean_09_or_01.txt hg19 15 0.5 200 10000 4 0.1 0.0005 1000 1000 > logs/input_FOXA1_regulated_vs_FOXA1_unregulated_mean_09_or_01.log &")

# Compute PWMs ----

# system(p("python ./extended_motifs.py ",pr_folder)) # FATTO SUL CLUSTER PERCHE' IN LOCALE NON ANDAVA


# MACRO-APE ----

clus <- makeCluster(ncpus_Macroape, type = "SOCK")

setwd(PATH_PR_FOLDER)
dir.create(PATH_PWM_FOLDER, showWarnings = F)

enr_tets = read.table(p(repo,PATH_RES_RNAMOTIFS,"enriched_tetramers.txt"), stringsAsFactors = F)$V1
motifs   = c()

for(i in enr_tets){
  temp_motif = list.files(PATH_PR_FOLDER, p(i,"_PWM.txt$"))
  motifs     = c(motifs,temp_motif)
  full_PWM   = read.table(p(PATH_PR_FOLDER, temp_motif),sep="\t")
  write.table(full_PWM[PWM_start:PWM_end,], p(PATH_PR_FOLDER, PATH_PWM_FOLDER, temp_motif), quote=F, sep="\t", row.names = F, col.names = F)
}

RNAmotifs = tail(strsplit(PATH_RNAMOTIFS,"/")[[1]],1)

clusterExport(clus,list("pr_folder", 'PATH_TO_MACROAPE', "ncol_PWMs", "RNAmotifs", "motifs","p", "PATH_RES", "PATH_PR_FOLDER","PATH_PWM_FOLDER", "PATH_TO_DB"))

setwd( PATH_TO_DB)
PWMs <- grep("HepG2",list.files(pattern=".txt$"), value = T)

setwd(p(repo, PATH_TO_MACROAPE))
dir.create(PATH_RES, showWarnings = F)

invisible(suppressWarnings(
  parSapply(clus, PWMs, function(x){
    print(x)
    name_PWM = substr(x,1,nchar(x)-4)
    
    for(m in motifs){
      print(m)
      tetr = substr(m,1,4)
      
      for(off in 0:3){
        tmp_file = p(tetr,"_",off,"_",name_PWM,".txt")
        write.table("", p(PATH_RES, tmp_file))
        system(p("java -cp ape.jar ru.autosome.macroape.EvalSimilarity ",  PATH_PR_FOLDER, PATH_PWM_FOLDER, m," ", PATH_TO_DB, "/", x," --position ", off-3,",direct --first-pcm --second-ppm > ", PATH_RES, tmp_file))
      }
      
      # Output example (without hashtags)
      #
      # S	6.426033208716181E-4
      # D	0.9993573966791284
      # L	16
      # SH	5
      # OR	direct
      # A1	>>>>>>>>>>......
      # A2	.....>>>>>>>>>>>
      # W	3265.0
      # W1	2768896.0
      # P1	6.44683837890625E-4
      # T1	6.3
      # W2	2315264.0
      # P2	5.390644073486328E-4
      # T2	6.3
      
      tmp_simil  <- 0
      tmp_offset <- 0
      
      for(off in 0:3){
        tmp_file = p(tetr,"_",off,"_",name_PWM,".txt")
        tmp <- read.table(p(PATH_RES,tmp_file), stringsAsFactors = F)
        
        if(as.numeric(tmp$V2[1]) > tmp_simil){
          tmp_simil  <- as.numeric(tmp$V2[1])
          tmp_offset <- off
        }
      }
      
      for(off in 0:3){
        if(off == tmp_offset){
          file.rename(p(PATH_RES,tetr,"_",tmp_offset,"_",name_PWM,".txt"),p(PATH_RES,tetr,"_",name_PWM,".txt"))
        } else {
          unlink(p(PATH_RES,tmp_file))
        }
      }
    }
  },simplify = T)))
stopCluster(clus)

setwd(p(pr_folder, "_", ncol_PWMs, "_columns"))   
files_mCross = list.files()                      # [tetr]_[cell_line].[RBP].[PWM_version].txt 

# Information to retrieve: tetr, RBM numerical ID, S, L, A1, A2, SH, OR, 

res_mCROSS <- data.frame(Tetramer     = character(),
                         Target_ID    = character(),
                         Similarity   = character(),
                         Total_length = character(),
                         Shift        = character(),
                         Orientation  = character(),
                         Alignment_1  = character(),
                         Alignment_2  = character())

for(i in 1:ncol(res_mCROSS)){
  res_mCROSS[,i] <- as.character(res_mCROSS[,i])
}

for(f in files_mCross){
  res_mCROSS[nrow(res_mCROSS) + 1,] <- c(substr(f, 1, 4), substr(f, start = 6, stop = nchar(f)-4), t(read.table(f, sep = "\t", stringsAsFactors = F)[c(1, 3, 4, 5, 6, 7), 2]))
}

res_mCROSS$RBP = sapply( strsplit(res_mCROSS$Target_ID,"[.]"),"[",2)
res_mCROSS     = res_mCROSS[order(res_mCROSS$Tetramer, -as.numeric(res_mCROSS$Similarity)),]

res            = res_mCROSS
res            = res[order(res$Tetramer, -as.numeric(res$Similarity)),]

setwd("..")
unlink(p(pr_folder,"_",ncol_PWMs,"_columns"),recursive = T)
write.csv(res, p(repo, PATH_TO_MACROAPE,"/",pr_folder,"_",ncol_PWMs,"_columns.csv"), row.names = F)


# Heatmap ----
setwd(repo)
ss = c()

input_df   = read.table(PATH_INP_RNAMOTIFS, header = F, sep=";", stringsAsFactors = F)
tbl_inp_df = table(input_df$V9)

load(PATH_DEX_FOXA1)
load(p(PATH_RES_RNAMOTIFS,"/bootstrap_10000.Rdata"))
df_enriched_tetramers <- read.csv(p(PATH_RES_RNAMOTIFS,"/",list.files(PATH_RES_RNAMOTIFS,".csv")[1]))
load(p(PATH_RES_RNAMOTIFS,"/",list.files(PATH_RES_RNAMOTIFS,"Rdata")[1]))
load(PATH_MCROSS_SCORES) 

p_enh = outRes %>% dplyr::select(p("r",rep(1:3,each=2),"enh_p",rep(c("Fis","Emp"),3)))
p_sil = outRes %>% dplyr::select(p("r",rep(1:3,each=2),"sil_p",rep(c("Fis","Emp"),3)))

p_cutoff =  p_enh %>% plyr::rename(c("r1enh_pFis" = "r1sil_pFis",
                               "r1enh_pEmp" = "r1sil_pEmp",
                               "r2enh_pFis" = "r2sil_pFis",
                               "r2enh_pEmp" = "r2sil_pEmp",
                               "r3enh_pFis" = "r3sil_pFis",
                               "r3enh_pEmp" = "r3sil_pEmp")) %>% 
  bind_rows(p_sil)                                 %>% 
  dplyr::select(p("r",1:3,"sil_pFis"))             %>% 
  plyr::rename(c("r1sil_pFis" = "r1_pFis_0.01",
           "r2sil_pFis" = "r2_pFis_0.01",
           "r3sil_pFis" = "r3_pFis_0.01"))         %>% 
  apply(2,quantile,0.01)


p_cutoff    = p_cutoff * (p_cutoff <= 0.05) + 0.05 * (p_cutoff > 0.05)
res         = outRes
res$is.sign = (
  (res$r1enh_pFis <= p_cutoff[1]  & res$r1enh_pEmp <= pemp_thresh) |
    (res$r2enh_pFis <= p_cutoff[2]  & res$r2enh_pEmp <= pemp_thresh) |
    (res$r3enh_pFis <= p_cutoff[3]  & res$r3enh_pEmp <= pemp_thresh) |
    
    (res$r1sil_pFis <= p_cutoff[1]  & res$r1sil_pEmp <= pemp_thresh) |
    (res$r2sil_pFis <= p_cutoff[2]  & res$r2sil_pEmp <= pemp_thresh) |
    (res$r3sil_pFis <= p_cutoff[3]  & res$r3sil_pEmp <= pemp_thresh) 
)

res = res %>% subset(is.sign) %>% as_tibble()


enh_ind = which(    (res$r1enh_pFis <= p_cutoff[1]  & res$r1enh_pEmp <= pemp_thresh) |
                      (res$r2enh_pFis <= p_cutoff[2]  & res$r2enh_pEmp <= pemp_thresh) |
                      (res$r3enh_pFis <= p_cutoff[3]  & res$r3enh_pEmp <= pemp_thresh))

sil_ind = which(    (res$r1sil_pFis <= p_cutoff[1]  & res$r1sil_pEmp <= pemp_thresh) |
                      (res$r2sil_pFis <= p_cutoff[2]  & res$r2sil_pEmp <= pemp_thresh) |
                      (res$r3sil_pFis <= p_cutoff[3]  & res$r3sil_pEmp <= pemp_thresh) )

res$exonType          = NA
res$exonType[enh_ind] = "enh"
res$exonType[sil_ind] = "sil"

if(length(intersect(enh_ind,sil_ind))>0){
  res$exonType[intersect(enh_ind,sil_ind)] = "both"
}

res$R1 = res$R2 = res$R3 = "0"

for(i in 1:3){
  
  res[,p("R",i)][  which( res[,p("r",i,"enh_pFis")] <= p_cutoff[i] & res[,p("r",i,"enh_pEmp")] <= pemp_thresh) , 1] = "Enhanced"
  res[,p("R",i)][  which( res[,p("r",i,"sil_pFis")] <= p_cutoff[i] & res[,p("r",i,"sil_pEmp")] <= pemp_thresh) , 1] = "Silenced"
  
  both_ind       = which((res[,p("r",i,"enh_pFis")] <= p_cutoff[i] & res[,p("r",i,"enh_pEmp")] <= pemp_thresh) & 
                           (res[,p("r",i,"sil_pFis")] <= p_cutoff[i] & res[,p("r",i,"sil_pEmp")] <= pemp_thresh))
  
  if(length(both_ind)>0){
    res[,p("R",i)][both_ind ,1] = "Both"
  }
}

x           = as.matrix(res[,c("R1","R2","R3")])
rownames(x) = res$tetramer

enr_tets_complete = df_enriched_tetramers %>% subset((r1_pf <= p_cutoff[1] & r1_pe <= pemp_thresh) |
                                                       (r2_pf <= p_cutoff[2] & r2_pe <= pemp_thresh) |
                                                       (r3_pf <= p_cutoff[3] & r3_pe <= pemp_thresh))

for(i in 1:3) enr_tets_complete[,p("R",i,"_sign")]  = enr_tets_complete[,p("r",i,"_pf")] <= p_cutoff[i] & enr_tets_complete[,p("r",i,"_pe")] <= pemp_thresh
for(i in 1:3) enr_tets_complete[,p("where.sign",i)] = enr_tets_complete[,c("R1_sign","R2_sign","R3_sign")] %>% apply(1,function(x){which(x)[i]})

enr_tets_complete$is.sign    = enr_tets_complete[,c("R1_sign","R2_sign","R3_sign")] %>% apply(1,function(x){sum(x)>0})
enr_tets_complete$where.sign = enr_tets_complete %>% with(paste(where.sign1,where.sign2,where.sign3,sep=","))

enr_tets = enr_tets_complete$tetramer %>% unique %>% cha


reg = x[enr_tets,] %>% t

if(nrow(reg)==1){reg = reg %>% t; colnames(reg) = enr_tets}

reg[reg == "0"]        = "white"
reg[reg == "Enhanced"] = "red"
reg[reg == "Silenced"] = "blue"
reg[reg == "Both"]     = "yellow"

names(reg) = enr_tets


pFis = matrix(0, 2, length(enr_tets), dimnames = list(c("enh","sil"), enr_tets))

Enh        = seq_along(enr_tets) %>% sapply(function(i){ifelse(sum(reg[,i] %in% c("red", "yellow"))>0, "red", "white")})
names(Enh) = colnames(reg)

Sil        = seq_along(enr_tets) %>% sapply(function(i){ifelse(sum(reg[,i] %in% c("blue", "yellow"))>0, "blue", "white")})
names(Sil) = colnames(reg)

num_exons  = matrix(0, 1, length(enr_tets), dimnames = list(c("num"),c(enr_tets)))



l1 = try({list_type_barplot1 = do.call(list,lapply(enr_tets, function(x){get_exons_with_tetramer_region1(x, paste0(repo, PATH_RES_RNAMOTIFS), PATH_INP_RNAMOTIFS)}))}, silent=T)
l2 = try({list_type_barplot2 = do.call(list,lapply(enr_tets, function(x){get_exons_with_tetramer_region2(x,PATH_RES_RNAMOTIFS, PATH_INP_RNAMOTIFS)}))}, silent=T)
l3 = try({list_type_barplot3 = do.call(list,lapply(enr_tets, function(x){get_exons_with_tetramer_region3(x,PATH_RES_RNAMOTIFS, PATH_INP_RNAMOTIFS)}))}, silent=T)

if(grepl("Error",l1[1])){list_type_barplot1 = lapply(rep(0,length(enr_tets)),function(x){x})}
if(grepl("Error",l2[1])){list_type_barplot2 = lapply(rep(0,length(enr_tets)),function(x){x})}
if(grepl("Error",l3[1])){list_type_barplot3 = lapply(rep(0,length(enr_tets)),function(x){x})}

names(list_type_barplot1) = names(list_type_barplot2) = names(list_type_barplot3) = enr_tets

t1 = try({type_barplot1    = do.call(rbind,lapply(enr_tets, function(x){table(list_type_barplot1[[x]]$type)}))},silent=T)
t2 = try({type_barplot2    = do.call(rbind,lapply(enr_tets, function(x){table(list_type_barplot2[[x]]$type)}))},silent=T)
t3 = try({type_barplot3    = do.call(rbind,lapply(enr_tets, function(x){table(list_type_barplot3[[x]]$type)}))},silent=T)

if(grepl("Error",t1[1])){type_barplot1 = matrix(0,length(enr_tets),2,dimnames=list(enr_tets,c(-1,1)))} else {rownames(type_barplot1) = enr_tets}
if(grepl("Error",t2[1])){type_barplot2 = matrix(0,length(enr_tets),2,dimnames=list(enr_tets,c(-1,1)))} else {rownames(type_barplot1) = enr_tets}
if(grepl("Error",t3[1])){type_barplot3 = matrix(0,length(enr_tets),2,dimnames=list(enr_tets,c(-1,1)))} else {rownames(type_barplot1) = enr_tets}

for(i in 1:3){
  if(ncol(get(p("type_barplot",i)))<2){
    tmp = matrix(0,length(enr_tets),2,dimnames = list(enr_tets,c(-1,1)))
    tmp[,colnames(get(p("type_barplot",i)))] = get(p("type_barplot",i))
    assign(p("type_barplot",i),tmp)
  }
}

type_barplot <- 0*type_barplot1

y = apply(x,1,function(z){z!="0"})

tetr_exons = list()

for(i in enr_tets){
  ind_sign = which(y[,i])
  exons = rbind()
  
  for(j in ind_sign){
    for(l in 1:3){
      
      if(j==l){
        tmp_df = get(p("list_type_barplot",l))[[i]]
        exons  = rbind(exons, tmp_df)
        tmp_df = tmp_df[which(!duplicated(tmp_df$key)),]
        nr_enh = tmp_df %>% subset(type == "1") %>% nrow
        nr_sil = tmp_df %>% subset(type == "-1") %>% nrow
        
        num_exons[1,i]   = num_exons[1,i] + nr_enh*(reg[j,i]=="red") + nr_sil*(reg[j,i]=="blue") + (nr_sil + nr_enh)*(reg[j,i]=="yellow")
      }
    }
  }
  
  exons              = exons[!duplicated(exons$key),]
  ln                 = length(tetr_exons)
  tetr_exons[[ln+1]] = exons
  type_barplot[i,]   = table(exons$type)
}

names(tetr_exons) = enr_tets
type_barplot      = type_barplot %>% t

alphabet = c("A","C","G","T")

res = read.csv(p(PATH_TO_MACROAPE,"/",pr_folder,"_",ncol_PWMs,"_columns.csv"),stringsAsFactors = F)
res = subset(res,Tetramer %in% enr_tets)

res$cell_line = sapply(res$Target_ID, function(x){strsplit(x,"[.]")[[1]][1]})
res$cell_line[!(res$cell_line %in% c("HepG2", "K562")) ] = "Unknown"

res2    = res
res     = res2[res2$cell_line == cl,] %>% as_tibble()
res$RBP = res$RBP %>% strsplit("[.]") %>% sapply(head,1)

unique_tets = res$Tetramer %>% unique
L_tets      = unique_tets  %>% length
RBPs        = res$RBP      %>% unique

M = matrix(NA, length(RBPs), L_tets, dimnames = list(RBPs,unique_tets))

p_file  = grep("^p[.]ena", list.files(PATH_RES_RNAMOTIFS), value=T) %>% tail(1)
p.ena   = read.csv(p(PATH_RES_RNAMOTIFS, p_file))
p_file  = grep("^p[.]sil", list.files(PATH_RES_RNAMOTIFS), value=T) %>% tail(1)
p.sil   = read.csv(p(PATH_RES_RNAMOTIFS, p_file))
datExpr = read.table(PATH_INP_RNAMOTIFS,sep=";", header=F)

peak_values  = list()

inInt = 200
inEx  = 50

for(rbp in RBPs){
  print(noquote(p(which(RBPs == rbp),"/",length(RBPs), "  ",rbp)))
  peaks1                      = read.delim(p(PATH_TO_PEAKS,grep("ordered_merged.bed$",list.files(pattern = rbp, path = PATH_TO_PEAKS),value = T)[1]), header=F)
  ln                          = length(peak_values)
  peak_values[[ln+1]]         = SEeCLIPpeaks(datExpr,peaks1,inInt,inEx)
  names(peak_values)[ln]      = rbp
}

if(nrow(p.ena)>1004){
  p.ena = p.ena[changeParam(1000,200,200,50),]
  p.sil = p.sil[changeParam(1000,200,200,50),]
}

ind_cntr = which(datExpr$V9==0)
for(i in seq_along(unique_tets)){
  temp = subset(res,Tetramer==unique_tets[i])
  for(j in unique(temp$RBP)){
    
    zz = 1
    if(eCLIP_weights){
      
      print(noquote(paste(i,"/",length(unique_tets),"  ",which(unique(temp$RBP)==j),"/",length(unique(temp$RBP)))))
      exon_ind = tetr_exons[[unique_tets[i]]]$myRID
      
      tmp_datExpr = datExpr[c(exon_ind,ind_cntr),]
      
      e1    = peak_values[[j]]$values[intersect(which(datExpr$V9>0), exon_ind),]
      s1    = peak_values[[j]]$values[intersect(which(datExpr$V9<0), exon_ind),]
      cntr1 = peak_values[[j]]$values[intersect(which(datExpr$V9==0),ind_cntr),]
      
      if(length(e1)>0){
        A1 = lmb.cluster.fisher.test( as.matrix(colSums(e1>0)),nrow(e1), as.matrix(colSums(cntr1>0)),nrow(cntr1))
      } else {
        A1 = rep(0,1004)
      }
      if(length(e1)>0){
        B1 = lmb.cluster.fisher.test( as.matrix(colSums(s1>0)),nrow(s1), as.matrix(colSums(cntr1>0)),nrow(cntr1))
      } else {
        B1 = rep(0,1004)
      }
      
      a1 = -2*log(A1)
      b1 = -2*log(B1)
      
      p_thr = -2*log(0.05)
      
      a1[a1 < p_thr] = 0
      b1[b1 < p_thr] = 0
      
      p.e = -2*log(p.ena)
      p.s = -2*log(p.sil)
      
      p.e[p.e < p_thr] = 0
      p.s[p.s < p_thr] = 0
      
      sc = c( BC(p.e, a1), 
              BC(p.s, b1))
      
      sc[is.na(sc)] = 0
      
      zz = mean(sc)
    }
    
    M[which(RBPs==j),i] = mean(temp$Similarity[which(temp$RBP==j)])*zz
    
  }
}

enriched_tetramers = enr_tets_complete$tetramer

M    = M[,enriched_tetramers]
orig = M

m = mean(M)
s = sd(M)

M = (M-m)/s

q_M = quantile(M, percentile/100)
q   = quantile(orig, c(0,percentile/100,0.95,1)) 
lq  = log(q + c(1e-4,rep(0,3)))

colfun3   = colorRamp2(breaks = (c((exp(seq(lq[2],lq[3],length.out=20))-m)/s,
                                  (exp(seq(lq[3],lq[4],length.out=50))-m)/s)+0.137),
                       colors = c(viridis::viridis(100)[31:100]))

# ord_dint  = c("WTCY","YTCY","TTTC","YTTS","WTTS","YTTY",
#               "WTTY","YATY","GTTT","YCTY","SGTW","WCCY",
#               "SCTY","CTTA","WGCW","YGCY","RTGY","WCAS",
#               "WCAY","YGCW","CTAA","ACAA","TTAC")

ord_dint  = c("RTGW","RTGY","YGTS","ATGT","YCTY","YTCY","YTTS","TGTG",
              "WACW","WCAS","YGCW","YGCY","CTAA")

saveRDS(M, p(PATH_TO_MACROAPE,"/",pr_folder,"_",ncol_PWMs,"_columns_MATRIX_FOR_HEATMAP.rds"))

N         = M[,ord_dint]

N[N<=q_M] = q_M - 1
N         = N[which(apply(N, 1, function(x){sum(x>q_M)>0})),]
N         = N[,which(apply(N,2, function(x){sum(x>q_M)>0}))]

ord_row = order(apply(N, 1, function(x){mean(x[x>q_M])}),decreasing=T)
N       = N[ord_row,]
M       = M[rownames(N),colnames(N)]

colfun_notsign = colorRamp2(breaks = (seq(q[1],q[2],length.out=60)-m)/s,colors = viridis::viridis(200)[1:60])

N = N[rownames(N) %in% c(two_tailed_FDR_list,two_tailed_FDR_met_list),]
M = M[rownames(N),]

if (min(N)<0) N = N-min(N) 

anno_PC    = anno_sign(rownames(N), two_tailed_FDR_list, "black")
anno_mCRPC = anno_sign(rownames(N), two_tailed_FDR_met_list, "orange")

anno_R1    = anno_reg(reg[1,colnames(N)])
anno_R2    = anno_reg(reg[2,colnames(N)])
anno_R3    = anno_reg(reg[3,colnames(N)])

anno_mt    = anno_reg(reg[1,colnames(N)], "white","white", "white", 1)
anno_mt_r  = anno_mt_row(nrow(N))


hmap = Heatmap(N,
               width               = ncol(N)*unit(4, "mm"), 
               height              = nrow(N)*unit(4, "mm"),
               border              = TRUE,
               cluster_rows        = FALSE,
               cluster_columns     = FALSE,
               col                 = colfun3,
               column_names_side   = "top",
               show_heatmap_legend = F,
               cell_fun            = function(j, i, x, y, width, height, fill) {
                 if(M[i,j]>=q_M){
                   grid.rect(x = x, y = y, width = width, height = height, gp = gpar( fill = NA))
                   grid.text(sprintf("%.3f", M[i, j]-min(M)), x, y, gp = gpar(col="white",fontsize = 3))
                 } else {
                   grid.rect(x = x, y = y, width = width, height = height, gp = gpar( fill = "white"))
                   grid.rect(x = x, y = y, width = width, height = height, gp = gpar(alpha = 0.5, fill = colfun_notsign(M[i,j])))
                   grid.text(sprintf("%.3f", M[i, j]-min(M)), x, y, gp = gpar(col="white",fontsize = 3))
                 }
               },
               top_annotation   = HeatmapAnnotation('N sil'       = anno_text(type_barplot["-1",colnames(N)],    gp = gpar(fontsize =9)),
                                                    'N Exons SIL' = anno_barplot(type_barplot["-1",colnames(N)], gp = gpar(fill = Sil[colnames(N)])),
                                                    '  '          = anno_mt,
                                                    'N enh'       = anno_text(type_barplot["1",colnames(N)],     gp = gpar(fontsize =9)),
                                                    'N Exons ENH' = anno_barplot(type_barplot["1",colnames(N)],  gp = gpar(fill = Enh[colnames(N)])),
                                                    '   '         = anno_mt,
                                                    'R1'          = anno_R1,
                                                    'R2'          = anno_R2,
                                                    'R3'          = anno_R3),
               right_annotation =     rowAnnotation('  '          = anno_mt_r,
                                                    '   '         = anno_mt_r,
                                                    'N'           = anno_barplot(apply(N,1,function(x){sum(x>q_M)})),
                                                    '    '        = anno_mt_r,
                                                    'Mean'        = anno_barplot(apply(N,1,function(x){mean(x[x>q_M])})),
                                                    '     '       = anno_mt_r),
               left_annotation  =     rowAnnotation('    '        = anno_mt_r,
                                                    'pri-PC'      = anno_PC, 
                                                    '  '          = anno_mt_r,
                                                    'mCRPC'       = anno_mCRPC,
                                                    ' '           = anno_mt_r,
                                                    annotation_name_side = "top")
)

hmap_lgd = Legend(col_fun = colfun3)

pdf(p(PATH_RES_RNAMOTIFS,"Heatmap_",format(Sys.time(),"%Y%m%d_%H_%M"),".pdf"), width = 14, height = 9)
print(hmap)
draw(hmap_lgd, x = unit(0.75,"npc"), y = unit(0.5,"npc"))
dev.off()
