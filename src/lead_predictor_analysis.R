library(fBasics)
#library(GGally)
#library(future)
#library(future.apply)

setwd('z:/Synthetic-Forecast-v2-FIRO-DISES/')

source('z:/Synthetic-Forecast_Verification/src/forecast_verification_functions.R')

loc <- 'NHG'

#load in the prepared data
load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))

keysite <- which(site_names=="NHGC1")

#plan(multicore)
#ecrps_vec<-future_sapply(1:n_hefs,function(x){eCRPS(hefs_forward[keysite,,x,ld],obs_forward_all_leads_hind[keysite,x,ld])})

fit_knn_hefs_data <- hefs_forward_cumul[keysite,,]
#gen_knn_shefs_data <- cumul_gen[keysite,,]
fit_cumul_obs_fwd <- apply(obs_forward_all_leads_hind[keysite,,],1,sum)
#gen_cumul_obs_fwd <- apply(obs_forward_all_leads_gen[keysite,,],1,sum)
fit_cumul_hefs <- hefs_forward_cumul_ens_avg[keysite,]
hefs_forward_frac[is.na(hefs_forward_frac)==T]<-0

rnk_fun <- function(x){
  rnk_cdf <- rank(x,ties.method = 'random') / (length(x)+1)
  return(rnk_cdf[1])
}

frac_fun <- function(x){
  out<-x / sum(x)
  out[is.na(out)==T]<-0
  return(out)}



hefs_skew <- apply(fit_knn_hefs_data,2,skewness)
hefs_kurt <- apply(fit_knn_hefs_data,2,kurtosis)
hefs_sd <- apply(fit_knn_hefs_data,2,sd)
hefs_mn <- apply(fit_knn_hefs_data,2,mean)
rnk_cdf_hefs <- apply(rbind(fit_cumul_obs_fwd,fit_knn_hefs_data),2,rnk_fun)

dat = cbind(hefs_frac_mean,hefs_frac_med,ecrps_vec,hefs_mn,hefs_sd,rnk_cdf_hefs,hefs_kurt,hefs_skew)
colnames(dat)=c('ens_mn_frac','ens_med_frac','ecrps','ens_mean','ens_sdev','ens_obs_rank','ens_kurt','ens_skew')

#ggpairs(dat)


#---------------------------------------------------
#plots to look at fractionation in fitted data
#HEFS
ld = 1

srt_cumul_obs <- order(fit_cumul_obs_fwd,decreasing = T)
srt_obs <- order(obs_forward_all_leads_hind[keysite,,ld],decreasing=T)
obs_frac <- apply(obs_forward_all_leads_hind[keysite,,],1,frac_fun)
hefs_frac_mean <- apply(hefs_forward_frac[keysite,,,],c(2,3),function(x){mean(x,na.rm=T)})
hefs_frac_sd <- apply(hefs_forward_frac[keysite,,,],c(2,3),function(x){sd(x,na.rm=T)})
#hefs_frac_med <- apply(hefs_forward_frac[keysite,,,],c(2,3),function(x){median(x,na.rm=T)})

library(scales)

#my_power <- -1    #larger neg values = faster decay, more weight on early leads
#decay <- (1:30)^my_power / sum(((1:30)^my_power)) #k dim vector
my_power <- 0    #larger neg values = faster decay, more weight on early leads
decay <- (1:15)^my_power / sum(((1:15)^my_power)) #k dim vector

#////////////////////////////////////////////////////////////////////////////
#selecting by hefs frac
par(mfrow=c(4,3),mar=c(4,4,0,0),mgp=c(1.5,.5,0))
for(i in 1:12){
  plot(1:15,obs_frac[,srt_obs[i]],type='l',ylim=c(0,0.5),lwd=3,xlab='leads',ylab='fractionation')
  lines(1:15,hefs_frac_mean[srt_obs[i],],col='gray50',lwd=2)
  knn_dist <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist <- order(knn_dist)[1:10]
  #knn_dist <- sqrt(colSums(decay*(c(hefs_frac_mean[srt_obs[i],],hefs_frac_sd[srt_obs[i],]) - t(cbind(hefs_frac_mean,hefs_frac_sd)))^2))
  #ord_dist <- order(knn_dist)[1:10]
  #frac_vec <- c()
  #frac_vec[seq(1,30,2)]<-hefs_frac_mean[srt_obs[i],]
  #frac_vec[seq(2,30,2)]<-hefs_frac_sd[srt_obs[i],]
  #frac_arr<-array(NA,c(dim(hefs_frac_mean)[1],30))
  #frac_arr[,seq(1,30,2)]<-hefs_frac_mean
  #frac_arr[,seq(2,30,2)]<-hefs_frac_sd
  #knn_dist <- sqrt(colSums(decay*(frac_vec - t(frac_arr))^2))
  #ord_dist <- order(knn_dist)[1:10]
  for(k in 1:10){
    lines(1:15,hefs_frac_mean[ord_dist[k],],col=alpha(colour=k,alpha=0.5))
  }
}

par(mfrow=c(4,3),mar=c(4,4,0,0),mgp=c(1.5,.5,0))
for(i in 1:12){
  plot(1:15,hefs_frac_sd[srt_obs[i],],col='gray50',lwd=2)
  knn_dist <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist <- order(knn_dist)[1:10]
  #knn_dist <- sqrt(colSums(decay*(c(hefs_frac_mean[srt_obs[i],],hefs_frac_mean[srt_obs[i],]) - t(cbind(hefs_frac_mean,hefs_frac_sd)))^2))
  #ord_dist <- order(knn_dist)[1:10]
  #frac_vec <- c()
  #frac_vec[seq(1,30,2)]<-hefs_frac_mean[srt_obs[i],]
  #frac_vec[seq(2,30,2)]<-hefs_frac_sd[srt_obs[i],]
  #frac_arr<-array(NA,c(dim(hefs_frac_mean)[1],30))
  #frac_arr[,seq(1,30,2)]<-hefs_frac_mean
  #frac_arr[,seq(2,30,2)]<-hefs_frac_sd
  #knn_dist <- sqrt(colSums(decay*(frac_vec - t(frac_arr))^2))
  #ord_dist <- order(knn_dist)[1:10]
  for(k in 1:10){
    lines(1:15,hefs_frac_sd[ord_dist[k],],col=alpha(colour=k,alpha=0.5))
  }
}

#can you get skill back with same cumul ensemble and different fractionations?
#ensemble scatter plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist <- order(knn_dist)[1:10]
  #hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  plot(0,obs_forward_all_leads_hind[keysite,srt_obs[i],ld],xlim=c(0,10),ylim=c(0,20),xlab='NN #',ylab='ens flow',cex=3,pch=18)
  points(rep(0,42),hefs_forward[keysite,,srt_obs[i],ld],col=alpha('blue',alpha=0.25))
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    points(rep(k,42),hefs_ens_pred[,ld],col=alpha('orange',alpha=0.25))
  }
}

#eCRPS plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist <- order(knn_dist)[1:10]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  boxplot(samp_hefs_ecrps,ylab='eCRPS',border='orange',col='white')
  points(1,hefs_ecrps,col='blue',cex=3,pch=17)
}

#-----------------------
#does sorting matter
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist <- order(knn_dist)[1:10]
  #hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  plot(0,obs_forward_all_leads_hind[keysite,srt_obs[i],ld],xlim=c(0,10),ylim=c(0,20),xlab='NN #',ylab='ens flow',pch=18,cex=3)
  points(rep(0,42),hefs_forward[keysite,,srt_obs[i],ld],col=alpha('blue',alpha=0.25))
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    points(rep(k,42),hefs_ens_pred[,ld],col=alpha('orange',alpha=0.25))
  }
}

#eCRPS plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist <- order(knn_dist)[1:10]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  boxplot(samp_hefs_ecrps,ylab='eCRPS',border='orange',col='white')
  points(1,hefs_ecrps,col='blue',cex=3,pch=17)
}

#------------------------
#bias plots
#nosort
ecrps_err<-c()

for(i in 1:100){
  knn_dist <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist <- order(knn_dist)[1:30]
  #knn_dist <- sqrt(colSums(decay*(c(hefs_frac_mean[srt_obs[i],],hefs_frac_sd[srt_obs[i],]) - t(cbind(hefs_frac_mean,hefs_frac_sd)))^2))
  #frac_vec <- c()
  #frac_vec[seq(1,30,2)]<-hefs_frac_mean[srt_obs[i],]
  #frac_vec[seq(2,30,2)]<-hefs_frac_sd[srt_obs[i],]
  #frac_arr<-array(NA,c(dim(hefs_frac_mean)[1],30))
  #frac_arr[,seq(1,30,2)]<-hefs_frac_mean
  #frac_arr[,seq(2,30,2)]<-hefs_frac_sd
  #knn_dist <- sqrt(colSums(decay*(frac_vec - t(frac_arr))^2))
  #ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='HEFS-frac;no-sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

#sorted
#my_power <- -10    #larger neg values = faster decay, more weight on early leads
#decay <- (1:15)^my_power / sum(((1:15)^my_power)) #k dim vector

ecrps_err<-c()

for(i in 1:100){
  knn_dist <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist <- order(knn_dist)[1:30]
  #knn_dist <- sqrt(colSums(decay*(c(hefs_frac_mean[srt_obs[i],],hefs_frac_sd[srt_obs[i],]) - t(cbind(hefs_frac_mean,hefs_frac_sd)))^2))
  #ord_dist <- order(knn_dist)[1:30]
  #frac_vec <- c()
  #frac_vec[seq(1,30,2)]<-hefs_frac_mean[srt_obs[i],]
  #frac_vec[seq(2,30,2)]<-hefs_frac_sd[srt_obs[i],]
  #frac_arr<-array(NA,c(dim(hefs_frac_mean)[1],30))
  #frac_arr[,seq(1,30,2)]<-hefs_frac_mean
  #frac_arr[,seq(2,30,2)]<-hefs_frac_sd
  #knn_dist <- sqrt(colSums(decay*(frac_vec - t(frac_arr))^2))
  #ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='HEFS-frac;sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

#////////////////////////////////////////////////////////////////////////////
#selecting by obs fractionation
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  plot(1:15,obs_frac[,srt_obs[i]],type='l',ylim=c(0,0.5),lwd=3,xlab='leads',ylab='fractionation')
  lines(1:15,hefs_frac_mean[srt_obs[i],],col='gray50',lwd=2)
  knn_dist <- sqrt(colSums(decay*(obs_frac[,srt_obs[i]] - obs_frac)^2))
  ord_dist <- order(knn_dist)[1:10]
  for(k in 1:10){
    lines(1:15,hefs_frac_mean[ord_dist[k],],col=alpha(colour=k,alpha=0.5))
  }
}

#can you get skill back with same cumul ensemble and different fractionations?
#ensemble scatter plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(obs_frac[,srt_obs[i]] - obs_frac)^2))
  ord_dist <- order(knn_dist)[1:10]
  #hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  plot(0,obs_forward_all_leads_hind[keysite,srt_obs[i],ld],xlim=c(0,10),ylim=c(0,20),xlab='NN #',ylab='ens flow',cex=3,pch=18)
  points(rep(0,42),hefs_forward[keysite,,srt_obs[i],ld],col=alpha('blue',alpha=0.25))
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    points(rep(k,42),hefs_ens_pred[,ld],col=alpha('orange',alpha=0.25))
  }
}

#eCRPS plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(obs_frac[,srt_obs[i]] - obs_frac)^2))
  ord_dist <- order(knn_dist)[1:10]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  boxplot(samp_hefs_ecrps,ylab='eCRPS',border='orange',col='white')
  points(1,hefs_ecrps,col='blue',cex=3,pch=17)
}

#-----------------------
#does sorting matter
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(obs_frac[,srt_obs[i]] - obs_frac)^2))
  ord_dist <- order(knn_dist)[1:10]
  #hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  plot(0,obs_forward_all_leads_hind[keysite,srt_obs[i],ld],xlim=c(0,10),ylim=c(0,20),xlab='NN #',ylab='ens flow',pch=18,cex=3)
  points(rep(0,42),hefs_forward[keysite,,srt_obs[i],ld],col=alpha('blue',alpha=0.25))
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    points(rep(k,42),hefs_ens_pred[,ld],col=alpha('orange',alpha=0.25))
  }
}

#eCRPS plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(obs_frac[,srt_obs[i]] - obs_frac)^2))
  ord_dist <- order(knn_dist)[1:10]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  boxplot(samp_hefs_ecrps,ylab='eCRPS',border='orange',col='white')
  points(1,hefs_ecrps,col='blue',cex=3,pch=17)
}

#------------------------
#bias plots
#nosort
ecrps_err<-c()

for(i in 1:100){
  knn_dist <- sqrt(colSums(decay*(obs_frac[,srt_obs[i]] - obs_frac)^2))
  ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:30){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='Obs-frac;no-sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

#sorted
ecrps_err<-c()

for(i in 1:100){
  knn_dist <- sqrt(colSums(decay*(obs_frac[,srt_obs[i]] - obs_frac)^2))
  ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:30){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='Obs-frac;sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

#////////////////////////////////////////////////////////////////////////////
#selecting by obs fwd
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  plot(1:15,obs_frac[,srt_obs[i]],type='l',ylim=c(0,0.5),lwd=3,xlab='leads',ylab='fractionation')
  lines(1:15,hefs_frac_mean[srt_obs[i],],col='gray50',lwd=2)
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:10]
  for(k in 1:10){
    lines(1:15,hefs_frac_mean[ord_dist[k],],col=alpha(colour=k,alpha=0.5))
  }
}

#can you get skill back with same cumul ensemble and different fractionations?
#ensemble scatter plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:10]
  #hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  plot(0,obs_forward_all_leads_hind[keysite,srt_obs[i],ld],xlim=c(0,10),ylim=c(0,20),xlab='NN #',ylab='ens flow',cex=3,pch=18)
  points(rep(0,42),hefs_forward[keysite,,srt_obs[i],ld],col=alpha('blue',alpha=0.25))
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    points(rep(k,42),hefs_ens_pred[,ld],col=alpha('orange',alpha=0.25))
  }
}

#eCRPS plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:10]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  boxplot(samp_hefs_ecrps,ylab='eCRPS',border='orange',col='white')
  points(1,hefs_ecrps,col='blue',cex=3,pch=17)
}

#-----------------------
#does sorting matter
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:10]
  #hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  plot(0,obs_forward_all_leads_hind[keysite,srt_obs[i],ld],xlim=c(0,10),ylim=c(0,20),xlab='NN #',ylab='ens flow',pch=18,cex=3)
  points(rep(0,42),hefs_forward[keysite,,srt_obs[i],ld],col=alpha('blue',alpha=0.25))
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    points(rep(k,42),hefs_ens_pred[,ld],col=alpha('orange',alpha=0.25))
  }
}

#eCRPS plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:10]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  boxplot(samp_hefs_ecrps,ylab='eCRPS',border='orange',col='white')
  points(1,hefs_ecrps,col='blue',cex=3,pch=17)
}

#------------------------
#bias plots
#nosort
ecrps_err<-c()

for(i in 1:100){
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='Obs-fwd;no-sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

#sorted
ecrps_err<-c()

for(i in 1:100){
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='Obs-fwd;sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)






#//////////////////////////////////////////////////////////////////////////////////

#sorted w/sampled cumul
ecrps_err<-c()

for(i in 1:100){
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,ord_dist[k]]
    rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='Obs-fwd+samp_cumul;sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

#empirical frac w/sample cumul

ecrps_err<-c()

for(i in 1:100){
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,ord_dist[k]]
    rnk<-rank(hefs_forward_cumul[keysite,,srt_obs[i]])
    hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,srt_obs[i],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='Obs-fwd+emp frac;sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

#empirical cumul w/sample frac

ecrps_err<-c()

for(i in 1:100){
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    #rnk<-rank(hefs_forward_cumul[keysite,,srt_obs[i]])
    #hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='Obs-fwd+emp cumul;sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

#empirical cumul w/1-day frac samp

ecrps_err<-c()



for(i in 1:100){
  knn_dist <- sqrt(colSums((sort(hefs_forward_frac[keysite,,srt_obs[i],ld]) - apply(hefs_forward_frac[keysite,,,ld],2,sort))^2))
  ord_dist <- order(knn_dist)[1:30]
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    #rnk<-rank(hefs_forward_cumul[keysite,,ord_dist[k]])
    #hefs_cumul <- sort(hefs_cumul)[rnk]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='Obs-fwd+emp cumul+1dfrac;sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

##############################################END###############################
