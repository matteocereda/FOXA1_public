# ••••••••••••••••••••••••••••••••••••••••• ----
# *** Cumulative event distribution of FOXA1-regulated events from RNAseq data of VCaP and PC3 cells ----


library(ggplot2)
library(gridExtra)
repo = '/Volumes/Prule/repo/FOXA1_public/'
FIG_DIR = '~/Downloads/FOXA1/Figure2/'
setwd(repo)
system(paste0('mkdir ', FIG_DIR))



pc3=readRDS('Rdata/Whippet_PC3_nextseq_novaseq_filter_0_2.rds')
pc3=subset(pc3,!Type%in%c('TE','TS'))
pc3=subset(pc3,STRINGENT)
pc3=subset(pc3,!Complexity=='K0')

vcap=readRDS('Rdata/Whippet_VCAP_filter_0_2.rds')
vcap=subset(vcap,!Type%in%c('TE','TS'))
vcap=subset(vcap,STRINGENT)
vcap=subset(vcap,!Complexity=='K0')


# PC3 ----
dexSIGN=pc3

dexSIGN$overall.mean = (dexSIGN$Psi_A+dexSIGN$Psi_B)/2

dexSIGN_sil=subset(dexSIGN,DeltaPsi<0)
dexSIGN_enh=subset(dexSIGN,DeltaPsi>0)

res = hist(dexSIGN_sil$overall.mean,nclass = 50)
res_df_sil=as.data.frame(cbind(res$mids,res$counts))
colnames(res_df_sil)=c('overall.mean.psi','counts')
res_df_sil$cum_counts=cumsum(res_df_sil$counts)
ggplot(res_df_sil,aes(x=overall.mean.psi,y=cum_counts))+geom_line()+theme_bw()

res = hist(dexSIGN_enh$overall.mean,nclass = 50)
res_df_enh=as.data.frame(cbind(res$mids,res$counts))
colnames(res_df_enh)=c('overall.mean.psi','counts')
res_df_enh$cum_counts=cumsum(res_df_enh$counts)
ggplot(res_df_enh,aes(x=overall.mean.psi,y=cum_counts))+geom_line()+theme_bw()

res_df_sil$type='sil'
res_df_sil$cum_counts_fraction=res_df_sil$cum_counts/res_df_sil$cum_counts[nrow(res_df_sil)]
res_df_enh$type='enh'
res_df_enh$cum_counts_fraction=res_df_enh$cum_counts/res_df_enh$cum_counts[nrow(res_df_enh)]

res_df=rbind(res_df_sil,res_df_enh)
p_pc3=ggplot(res_df,aes(x=overall.mean.psi,y=cum_counts_fraction,col=type))+geom_line()+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed")+ggtitle('PC3')
p_pc3_counts=ggplot(res_df,aes(x=overall.mean.psi,y=cum_counts,col=type))+geom_line()+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed")
res_df_pc3_true = res_df


# VCaP ----

dexSIGN=vcap

dexSIGN$overall.mean = (dexSIGN$Psi_A+dexSIGN$Psi_B)/2

dexSIGN_sil=subset(dexSIGN,DeltaPsi<0)
dexSIGN_enh=subset(dexSIGN,DeltaPsi>0)

res = hist(dexSIGN_sil$overall.mean,nclass = 50)
res_df_sil=as.data.frame(cbind(res$mids,res$counts))
colnames(res_df_sil)=c('overall.mean.psi','counts')
res_df_sil$cum_counts=cumsum(res_df_sil$counts)
ggplot(res_df_sil,aes(x=overall.mean.psi,y=cum_counts))+geom_line()+theme_bw()

res = hist(dexSIGN_enh$overall.mean,nclass = 50)
res_df_enh=as.data.frame(cbind(res$mids,res$counts))
colnames(res_df_enh)=c('overall.mean.psi','counts')
res_df_enh$cum_counts=cumsum(res_df_enh$counts)
ggplot(res_df_enh,aes(x=overall.mean.psi,y=cum_counts))+geom_line()+theme_bw()

res_df_sil$type='sil'
res_df_sil$cum_counts_fraction=res_df_sil$cum_counts/res_df_sil$cum_counts[nrow(res_df_sil)]
res_df_enh$type='enh'
res_df_enh$cum_counts_fraction=res_df_enh$cum_counts/res_df_enh$cum_counts[nrow(res_df_enh)]

res_df=rbind(res_df_sil,res_df_enh)
p_vcap=ggplot(res_df,aes(x=overall.mean.psi,y=cum_counts_fraction,col=type))+geom_line()+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed")+ggtitle('VCaP')
p_vcap_counts=ggplot(res_df,aes(x=overall.mean.psi,y=cum_counts,col=type))+geom_line()+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed")
res_df_vcap_true = res_df


pdf(file = paste0(FIG_DIR, 'cumulative_plots_WHIPPET.pdf'),width = unit(10,'cm'),height = unit(6,'cm'))
grid.arrange(p_vcap,p_pc3,p_vcap_counts,p_pc3_counts,ncol=2)
dev.off()



# NULL MODEL ----
set.seed(8890)

for (i in 1:100) {
  
  print(i)
  
  # VCaP ---
  
  dexSIGN=vcap
  
  dexSIGN$overall.mean = (dexSIGN$Psi_A+dexSIGN$Psi_B)/2
  
  dexSIGN$delta.mean_random = sample(x = c(-1,1),size = nrow(dexSIGN),replace = T)

  dexSIGN_sil=subset(dexSIGN,delta.mean_random<0)
  dexSIGN_enh=subset(dexSIGN,delta.mean_random>0)
  
  res = hist(dexSIGN_sil$overall.mean,nclass = 50)
  res_df_sil=as.data.frame(cbind(res$mids,res$counts))
  colnames(res_df_sil)=c('overall.mean.psi','counts')
  res_df_sil$cum_counts=cumsum(res_df_sil$counts)

  res = hist(dexSIGN_enh$overall.mean,nclass = 50)
  res_df_enh=as.data.frame(cbind(res$mids,res$counts))
  colnames(res_df_enh)=c('overall.mean.psi','counts')
  res_df_enh$cum_counts=cumsum(res_df_enh$counts)

  res_df_sil$type='sil'
  res_df_sil$cum_counts_fraction=res_df_sil$cum_counts/res_df_sil$cum_counts[nrow(res_df_sil)]
  res_df_enh$type='enh'
  res_df_enh$cum_counts_fraction=res_df_enh$cum_counts/res_df_enh$cum_counts[nrow(res_df_enh)]
  
  res_df=rbind(res_df_sil,res_df_enh)
  
  res_df_2=cbind(res_df_sil,res_df_enh)
  colnames(res_df_2)[1:5]=paste0(colnames(res_df_2)[1:5],'_SIL')
  colnames(res_df_2)[6:10]=paste0(colnames(res_df_2)[6:10],'_ENH')
  res_df_2$enh_sil_ratio=res_df_2$counts_ENH/(res_df_2$counts_ENH+res_df_2$counts_SIL)
  
  if (i==1) {
    res_df$rep=i
    res_df_vcap = res_df
  } else {
    res_df$rep=i
    res_df_vcap = rbind(res_df_vcap,res_df)
  }
  
  # PC3 ---
  
  dexSIGN=pc3
  
  dexSIGN$overall.mean = (dexSIGN$Psi_A+dexSIGN$Psi_B)/2
  
  dexSIGN$delta.mean_random = sample(x = c(-1,1),size = nrow(dexSIGN),replace = T)

  dexSIGN_sil=subset(dexSIGN,delta.mean_random<0)
  dexSIGN_enh=subset(dexSIGN,delta.mean_random>0)
  
  res = hist(dexSIGN_sil$overall.mean,nclass = 50)
  res_df_sil=as.data.frame(cbind(res$mids,res$counts))
  colnames(res_df_sil)=c('overall.mean.psi','counts')
  res_df_sil$cum_counts=cumsum(res_df_sil$counts)

  res = hist(dexSIGN_enh$overall.mean,nclass = 50)
  res_df_enh=as.data.frame(cbind(res$mids,res$counts))
  colnames(res_df_enh)=c('overall.mean.psi','counts')
  res_df_enh$cum_counts=cumsum(res_df_enh$counts)

  res_df_sil$type='sil'
  res_df_sil$cum_counts_fraction=res_df_sil$cum_counts/res_df_sil$cum_counts[nrow(res_df_sil)]
  res_df_enh$type='enh'
  res_df_enh$cum_counts_fraction=res_df_enh$cum_counts/res_df_enh$cum_counts[nrow(res_df_enh)]
  
  res_df=rbind(res_df_sil,res_df_enh)
  
  res_df_2=cbind(res_df_sil,res_df_enh)
  colnames(res_df_2)[1:5]=paste0(colnames(res_df_2)[1:5],'_SIL')
  colnames(res_df_2)[6:10]=paste0(colnames(res_df_2)[6:10],'_ENH')
  res_df_2$enh_sil_ratio=res_df_2$counts_ENH/(res_df_2$counts_ENH+res_df_2$counts_SIL)
  
  if (i==1) {
    res_df$rep=i
    res_df_pc3 = res_df
  } else {
    res_df$rep=i
    res_df_pc3 = rbind(res_df_pc3,res_df)
  }
  
}


res_df_vcap$rep=as.character(res_df_vcap$rep)
res_df_pc3$rep=as.character(res_df_pc3$rep)

save(res_df_vcap,res_df_pc3,file='Rdata/null_model_splicing_RNAseq_prob_05_for_delta_pos_and_delta_neg.Rdata')


# vcap

to_plot_enh=ddply(subset(res_df_vcap,type=='enh'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts),median_c=median(cum_counts),c_95=quantile(cum_counts,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
to_plot_sil=ddply(subset(res_df_vcap,type=='sil'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts),median_c=median(cum_counts),c_95=quantile(cum_counts,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
p_vcap_compl = ggplot() +
  geom_line(data=res_df_vcap_true,mapping=aes(x=overall.mean.psi,y=cum_counts,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


tmp=subset(res_df_vcap_true,type=='sil')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts[which(tmp$overall.mean.psi>0.51)] = tmp$cum_counts[which(tmp$overall.mean.psi>0.51)]-tmp$cum_counts[which(tmp$overall.mean.psi==0.51)]
tmp_s=tmp
tmp=subset(res_df_vcap_true,type=='enh')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts[which(tmp$overall.mean.psi>0.51)] = tmp$cum_counts[which(tmp$overall.mean.psi>0.51)]-tmp$cum_counts[which(tmp$overall.mean.psi==0.51)]
tmp_e=tmp
res_df_vcap_true_TMP=rbind(tmp_s,tmp_e)

tmp=to_plot_enh
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi>0.51)] = tmp$mean_c[which(tmp$overall.mean.psi>0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_enh_TMP=tmp

tmp=to_plot_sil
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi>0.51)] = tmp$mean_c[which(tmp$overall.mean.psi>0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_sil_TMP=tmp

p_vcap_compl_2 = ggplot() +
  geom_line(data=res_df_vcap_true_TMP,mapping=aes(x=overall.mean.psi,y=cum_counts,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


# pc3

to_plot_enh=ddply(subset(res_df_pc3,type=='enh'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts),median_c=median(cum_counts),c_95=quantile(cum_counts,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
to_plot_sil=ddply(subset(res_df_pc3,type=='sil'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts),median_c=median(cum_counts),c_95=quantile(cum_counts,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
p_pc3_compl = ggplot() +
  geom_line(data=res_df_pc3_true,mapping=aes(x=overall.mean.psi,y=cum_counts,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta sd')


pdf(file = paste0(FIG_DIR, 'Whippet_RNAseq_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
grid.arrange(p_vcap_compl,p_pc3_compl,ncol=2)
dev.off()



tmp=subset(res_df_pc3_true,type=='sil')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts[which(tmp$overall.mean.psi>0.51)] = tmp$cum_counts[which(tmp$overall.mean.psi>0.51)]-tmp$cum_counts[which(tmp$overall.mean.psi==0.51)]
tmp_s=tmp
tmp=subset(res_df_pc3_true,type=='enh')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts[which(tmp$overall.mean.psi>0.51)] = tmp$cum_counts[which(tmp$overall.mean.psi>0.51)]-tmp$cum_counts[which(tmp$overall.mean.psi==0.51)]
tmp_e=tmp
res_df_pc3_true_TMP=rbind(tmp_s,tmp_e)

tmp=to_plot_enh
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi>0.51)] = tmp$mean_c[which(tmp$overall.mean.psi>0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_enh_TMP=tmp

tmp=to_plot_sil
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi>0.51)] = tmp$mean_c[which(tmp$overall.mean.psi>0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_sil_TMP=tmp

p_pc3_compl_2 = ggplot() +
  geom_line(data=res_df_pc3_true_TMP,mapping=aes(x=overall.mean.psi,y=cum_counts,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


pdf(file = paste0(FIG_DIR, 'Whippet_RNAseq_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_SPLITTED.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
grid.arrange(p_vcap_compl_2,p_pc3_compl_2,ncol=2)
dev.off()




# INVERTED ----

# VCaP

tmp_e = subset(res_df_vcap_true,type=='enh')
tmp_e=tmp_e[order(tmp_e$overall.mean.psi,decreasing = T),]
tmp_e$cum_counts_inverted=cumsum(tmp_e$counts)
tmp_s = subset(res_df_vcap_true,type=='sil')
tmp_s=tmp_s[order(tmp_s$overall.mean.psi,decreasing = T),]
tmp_s$cum_counts_inverted=cumsum(tmp_s$counts)

res_df_vcap_true_inverted=rbind(tmp_e,tmp_s)


for (i in unique(res_df_vcap$rep)) {
  
  tmp = subset(res_df_vcap,rep==i)
  tmp_e = subset(tmp,type=='enh')
  tmp_e=tmp_e[order(tmp_e$overall.mean.psi,decreasing = T),]
  tmp_e$cum_counts_inverted=cumsum(tmp_e$counts)
  tmp_s = subset(tmp,type=='sil')
  tmp_s=tmp_s[order(tmp_s$overall.mean.psi,decreasing = T),]
  tmp_s$cum_counts_inverted=cumsum(tmp_s$counts)
  
  tmp=rbind(tmp_e,tmp_s)
  
  if (i == unique(res_df_vcap$rep)[1]) {
    res_df_vcap_inverted = tmp
  } else {
    res_df_vcap_inverted = rbind(res_df_vcap_inverted,tmp)
  }
  
}


to_plot_enh=ddply(subset(res_df_vcap_inverted,type=='enh'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts_inverted),median_c=median(cum_counts_inverted),c_95=quantile(cum_counts_inverted,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts_inverted,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
to_plot_sil=ddply(subset(res_df_vcap_inverted,type=='sil'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts_inverted),median_c=median(cum_counts_inverted),c_95=quantile(cum_counts_inverted,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts_inverted,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
p_vcap_compl_inv = ggplot() +
  geom_line(data=res_df_vcap_true_inverted,mapping=aes(x=overall.mean.psi,y=cum_counts_inverted,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


tmp=subset(res_df_vcap_true_inverted,type=='sil')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)] = tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)]-tmp$cum_counts_inverted[which(tmp$overall.mean.psi==0.51)]
tmp_s=tmp
tmp=subset(res_df_vcap_true_inverted,type=='enh')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)] = tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)]-tmp$cum_counts_inverted[which(tmp$overall.mean.psi==0.51)]
tmp_e=tmp
res_df_vcap_true_TMP=rbind(tmp_s,tmp_e)

tmp=to_plot_enh
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi<0.51)] = tmp$mean_c[which(tmp$overall.mean.psi<0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_enh_TMP=tmp

tmp=to_plot_sil
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi<0.51)] = tmp$mean_c[which(tmp$overall.mean.psi<0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_sil_TMP=tmp

p_vcap_compl_2_inv = ggplot() +
  geom_line(data=res_df_vcap_true_TMP,mapping=aes(x=overall.mean.psi,y=cum_counts_inverted,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


# PC3

tmp_e = subset(res_df_pc3_true,type=='enh')
tmp_e=tmp_e[order(tmp_e$overall.mean.psi,decreasing = T),]
tmp_e$cum_counts_inverted=cumsum(tmp_e$counts)
tmp_s = subset(res_df_pc3_true,type=='sil')
tmp_s=tmp_s[order(tmp_s$overall.mean.psi,decreasing = T),]
tmp_s$cum_counts_inverted=cumsum(tmp_s$counts)

res_df_pc3_true_inverted=rbind(tmp_e,tmp_s)


for (i in unique(res_df_pc3$rep)) {
  
  tmp = subset(res_df_pc3,rep==i)
  tmp_e = subset(tmp,type=='enh')
  tmp_e=tmp_e[order(tmp_e$overall.mean.psi,decreasing = T),]
  tmp_e$cum_counts_inverted=cumsum(tmp_e$counts)
  tmp_s = subset(tmp,type=='sil')
  tmp_s=tmp_s[order(tmp_s$overall.mean.psi,decreasing = T),]
  tmp_s$cum_counts_inverted=cumsum(tmp_s$counts)
  
  tmp=rbind(tmp_e,tmp_s)
  
  if (i == unique(res_df_pc3$rep)[1]) {
    res_df_pc3_inverted = tmp
  } else {
    res_df_pc3_inverted = rbind(res_df_pc3_inverted,tmp)
  }
  
}


to_plot_enh=ddply(subset(res_df_pc3_inverted,type=='enh'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts_inverted),median_c=median(cum_counts_inverted),c_95=quantile(cum_counts_inverted,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts_inverted,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
to_plot_sil=ddply(subset(res_df_pc3_inverted,type=='sil'),.(overall.mean.psi),summarise,mean_c=mean(cum_counts_inverted),median_c=median(cum_counts_inverted),c_95=quantile(cum_counts_inverted,seq(0,1,0.01))['95%'],c_05=quantile(cum_counts_inverted,seq(0,1,0.01))['5%'],
                  mean_cf=mean(cum_counts_fraction),median_cf=median(cum_counts_fraction),cf_95=quantile(cum_counts_fraction,seq(0,1,0.01))['95%'],cf_05=quantile(cum_counts_fraction,seq(0,1,0.01))['5%'])
p_pc3_compl_inv = ggplot() +
  geom_line(data=res_df_pc3_true_inverted,mapping=aes(x=overall.mean.psi,y=cum_counts_inverted,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')


tmp=subset(res_df_pc3_true_inverted,type=='sil')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)] = tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)]-tmp$cum_counts_inverted[which(tmp$overall.mean.psi==0.51)]
tmp_s=tmp
tmp=subset(res_df_pc3_true_inverted,type=='enh')
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)] = tmp$cum_counts_inverted[which(tmp$overall.mean.psi<0.51)]-tmp$cum_counts_inverted[which(tmp$overall.mean.psi==0.51)]
tmp_e=tmp
res_df_pc3_true_TMP=rbind(tmp_s,tmp_e)

tmp=to_plot_enh
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi<0.51)] = tmp$mean_c[which(tmp$overall.mean.psi<0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_enh_TMP=tmp

tmp=to_plot_sil
tmp$overall.mean.psi=round(tmp$overall.mean.psi,digits = 2)
tmp$mean_c[which(tmp$overall.mean.psi<0.51)] = tmp$mean_c[which(tmp$overall.mean.psi<0.51)]-tmp$mean_c[which(tmp$overall.mean.psi==0.51)]
to_plot_sil_TMP=tmp

p_pc3_compl_2_inv = ggplot() +
  geom_line(data=res_df_pc3_true_TMP,mapping=aes(x=overall.mean.psi,y=cum_counts_inverted,col=type))+theme_bw()+scale_color_manual(values = c(rgb(223/255,0,63/255),rgb(18/255,0,124/255)))+
  geom_vline(xintercept = 0.15,linetype = "dashed")+geom_vline(xintercept = 0.5,linetype = "dashed")+geom_vline(xintercept = 0.85,linetype = "dashed") +
  geom_ribbon(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="red", alpha = .15) +
  geom_line(data=to_plot_enh_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(223/255,0,63/255)) + ylab("data")+theme_bw() +
  geom_ribbon(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c, ymin = c_05 , ymax = c_95), fill="blue", alpha = .15) +
  geom_line(data=to_plot_sil_TMP, mapping = aes(x=overall.mean.psi, y=mean_c), linetype = "dashed", col=rgb(18/255,0,124/255)) + ylab("data")+theme_bw() + ggtitle('Delta mean')






pdf(file = paste0(FIG_DIR, 'Whippet_RNAseq_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_inverted.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
grid.arrange(p_vcap_compl_inv,p_pc3_compl_inv,ncol=2)
dev.off()


pdf(file =paste0(FIG_DIR, 'Whippet_RNAseq_cumulative_with_null_model__prob_05_for_delta_pos_and_delta_neg_SPLITTED_inverted.pdf'),width = unit(10,'cm'),height = unit(3.5,'cm'))
grid.arrange(p_vcap_compl_2_inv,p_pc3_compl_2_inv,ncol=2)
dev.off()




