library(fBasics)
#library(GGally)
#library(future)
#library(future.apply)
library(betareg)

data("GasolineYield",package='betareg')
data("FoodExpenditure",package="betareg")

estBetaParams <- function(mu, var) {
  alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
  beta <- alpha * (1 / mu - 1)
  return(params = list(alpha = alpha, beta = beta))
}

gy2 = betareg(yield ~ batch + temp | batch + temp, data=GasolineYield)

mu = predict(gy2,type='response')
var = predict(gy2,type='variance')

par = estBetaParams(mu,var)

par(mfrow=c(4,3))
for(i in 1:12){
  hist(rbeta(1000,shape1 = par$alpha[i],shape2 = par$beta[i]))
}


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
fit_cumul_hefs_sd <- apply(hefs_forward_cumul[keysite,,],2,sd)
hefs_forward_frac[is.na(hefs_forward_frac)==T]<-0

rnk_fun <- function(x){
  rnk_cdf <- rank(x,ties.method = 'random') / (length(x)+1)
  return(rnk_cdf[1])
}

frac_fun <- function(x){
  out<-x / sum(x)
  out[is.na(out)==T]<-0
  return(out)}


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
#beta regression
options(digits = 4)

hefs_frac_mean[hefs_frac_mean<=0]<-0.0000001
hefs_frac_mean[hefs_frac_mean>=1]<-0.999999999
beta_mat <- cbind(hefs_frac_mean[,1],obs_forward_all_leads_hind[keysite,,1],fit_cumul_obs_fwd,fit_cumul_hefs,fit_cumul_hefs_sd)
colnames(beta_mat) <- c('frac','obs_fwd1','cumul_obs','cumul_hefs','cumul_hefs_sd')
beta_mat <- as.data.frame(beta_mat)
beta_reg1 = betareg(frac ~ obs_fwd1 + cumul_obs + cumul_hefs + cumul_hefs_sd | obs_fwd1 + cumul_obs + cumul_hefs + cumul_hefs_sd, data=beta_mat)
beta_reg1 = betareg(frac ~ obs_fwd1 + cumul_obs | obs_fwd1 + cumul_obs , data=beta_mat)

mu = predict(beta_reg1,type='response')
var = predict(beta_reg1,type='variance')

par = estBetaParams(mu,var)

#ensemble scatter plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  plot(0,obs_forward_all_leads_hind[keysite,srt_obs[i],ld],xlim=c(0,10),ylim=c(0,20),xlab='NN #',ylab='ens flow',cex=3,pch=18)
  points(rep(0,42),hefs_forward[keysite,,srt_obs[i],ld],col=alpha('blue',alpha=0.25))
  knn_dist <- sqrt((obs_forward_all_leads_hind[keysite,srt_obs[i],ld] - obs_forward_all_leads_hind[keysite,,ld])^2)
  ord_dist <- order(knn_dist)[1:10]
  for(k in 1:10){
    hefs_mn <- mean(hefs_forward[keysite,,ord_dist[k],ld])
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-rbeta(42,shape1 = par$alpha[srt_obs[i]],shape2 = par$beta[srt_obs[i]])*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    pred_mn <- mean(hefs_ens_pred)
    mn_scale <- hefs_mn/obs_forward_all_leads_hind[keysite,ord_dist[k],ld] * obs_forward_all_leads_hind[keysite,srt_obs[i],ld] / pred_mn
    hefs_ens_pred <- hefs_ens_pred * mn_scale
    points(rep(k,42),hefs_ens_pred[,ld],col=alpha('orange',alpha=0.25))
  }
}

#eCRPS plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 89:100){
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  knn_dist <- sqrt((obs_forward_all_leads_hind[keysite,srt_obs[i],ld] - obs_forward_all_leads_hind[keysite,,ld])^2)
  ord_dist <- order(knn_dist)[1:10]
  for(k in 1:10){
    hefs_mn <- mean(hefs_forward[keysite,,ord_dist[k],ld])
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-rbeta(42,shape1 = par$alpha[srt_obs[i]],shape2 = par$beta[srt_obs[i]])*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    pred_mn <- mean(hefs_ens_pred)
    mn_scale <- hefs_mn/obs_forward_all_leads_hind[keysite,ord_dist[k],ld] * obs_forward_all_leads_hind[keysite,srt_obs[i],ld] / pred_mn
    hefs_ens_pred <- hefs_ens_pred * mn_scale
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  boxplot(samp_hefs_ecrps,ylab='eCRPS',border='orange',col='white',ylim=c(0,1))
  points(1,hefs_ecrps,col='blue',cex=3,pch=17)
}


#------------------------
#bias plots
ecrps_err<-c()

for(i in 1:100){
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  knn_dist <- sqrt((obs_forward_all_leads_hind[keysite,srt_obs[i],ld] - obs_forward_all_leads_hind[keysite,,ld])^2)
  ord_dist <- order(knn_dist)[1:10]
  samp_hefs_ecrps <- c()
  for(k in 1:10){
    hefs_mn <- mean(hefs_forward[keysite,,ord_dist[k],ld])
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-rbeta(42,shape1 = par$alpha[srt_obs[i]],shape2 = par$beta[srt_obs[i]])*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    pred_mn <- mean(hefs_ens_pred)
    mn_scale <- hefs_mn/obs_forward_all_leads_hind[keysite,ord_dist[k],ld] * obs_forward_all_leads_hind[keysite,srt_obs[i],ld] / pred_mn
    hefs_ens_pred <- hefs_ens_pred * mn_scale
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main='HEFS-frac;no-sort ens',ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)


#/////////////////////////////////////////////////////////////////////////////
#NN fractionation
#ensemble scatter plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  plot(0,obs_forward_all_leads_hind[keysite,srt_obs[i],ld],xlim=c(0,10),ylim=c(0,20),xlab='NN #',ylab='ens flow',cex=3,pch=18)
  points(rep(0,42),hefs_forward[keysite,,srt_obs[i],ld],col=alpha('blue',alpha=0.25))
  knn_dist <- sqrt((obs_forward_all_leads_hind[keysite,srt_obs[i],ld] - obs_forward_all_leads_hind[keysite,,ld])^2)
  ord_dist <- order(knn_dist)[1:10]
  knn_dist_frac <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist_frac <- order(knn_dist_frac)[1:10]
  for(k in 1:10){
    hefs_mn <- mean(hefs_forward[keysite,,ord_dist[k],ld])
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist_frac[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    pred_mn <- mean(hefs_ens_pred[,1])
    mn_scale <- hefs_mn/obs_forward_all_leads_hind[keysite,ord_dist[k],ld] * obs_forward_all_leads_hind[keysite,srt_obs[i],ld] / pred_mn
    #hefs_ens_pred <- hefs_ens_pred * mn_scale
    hefs_ens_pred <- hefs_ens_pred + (pred_mn * mn_scale - pred_mn);hefs_ens_pred[hefs_ens_pred<0]<-0
    points(rep(k,42),hefs_ens_pred[,ld],col=alpha('orange',alpha=0.25))
  }
}

#eCRPS plots
par(mfrow=c(4,3),mar=c(3,3,0,0))
for(i in 1:12){
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  samp_hefs_ecrps <- c()
  knn_dist <- sqrt((obs_forward_all_leads_hind[keysite,srt_obs[i],ld] - obs_forward_all_leads_hind[keysite,,ld])^2)
  ord_dist <- order(knn_dist)[1:10]
  knn_dist_frac <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  ord_dist_frac <- order(knn_dist_frac)[1:10]
  for(k in 1:10){
    hefs_mn <- mean(hefs_forward[keysite,,ord_dist[k],ld])
    hefs_cumul <- hefs_forward_cumul[keysite,,srt_obs[i]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist_frac[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    pred_mn <- mean(hefs_ens_pred[,1])
    mn_scale <- hefs_mn/obs_forward_all_leads_hind[keysite,ord_dist[k],ld] * obs_forward_all_leads_hind[keysite,srt_obs[i],ld] / pred_mn
    #hefs_ens_pred <- hefs_ens_pred * mn_scale
    hefs_ens_pred <- hefs_ens_pred + (pred_mn * mn_scale - pred_mn);hefs_ens_pred[hefs_ens_pred<0]<-0
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  boxplot(samp_hefs_ecrps,ylab='eCRPS',border='orange',col='white')
  points(1,hefs_ecrps,col='blue',cex=3,pch=17)
}

#//////////////////////////////////////////////////////////////////////////////
#compare methods with and without mean shift
my_power <- -1    #larger neg values = faster decay, more weight on early leads
decay <- (1:15)^my_power / sum(((1:15)^my_power)) #k dim vector
nn <- 20
#------------------------
#bias plots
par(mfrow=c(1,2))
ecrps_err<-c()

for(i in 1:100){
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  #knn_dist <- sqrt((obs_forward_all_leads_hind[keysite,srt_obs[i],ld] - obs_forward_all_leads_hind[keysite,,ld])^2)
  #ord_dist <- order(knn_dist)[1:30]
  #knn_dist <- sqrt(colSums(decay*(obs_frac[,srt_obs[i]] - obs_frac)^2))
  #ord_dist <- order(knn_dist)[1:30]
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:30]
  #knn_dist_frac <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  #ord_dist_frac <- order(knn_dist_frac)[1:10]
  samp_hefs_ecrps <- c()
  for(k in 1:nn){
    hefs_cumul <- hefs_forward_cumul[keysite,,ord_dist[k]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main=paste('No mean shift, NN=',nn),ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)


#bias plots
ecrps_err<-c()

for(i in 1:100){
  hefs_ecrps <- eCRPS(hefs_forward[keysite,,srt_obs[i],ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  #knn_dist <- sqrt((obs_forward_all_leads_hind[keysite,srt_obs[i],ld] - obs_forward_all_leads_hind[keysite,,ld])^2)
  #ord_dist <- order(knn_dist)[1:30]
  #knn_dist <- sqrt(colSums(decay*(obs_frac[,srt_obs[i]] - obs_frac)^2))
  #ord_dist <- order(knn_dist)[1:30]
  knn_dist <- sqrt(colSums(decay*(obs_forward_all_leads_hind[keysite,srt_obs[i],] - t(obs_forward_all_leads_hind[keysite,,]))^2))
  ord_dist <- order(knn_dist)[1:30]
  #knn_dist_frac <- sqrt(colSums(decay*(hefs_frac_mean[srt_obs[i],] - t(hefs_frac_mean))^2))
  #ord_dist_frac <- order(knn_dist_frac)[1:30]
  samp_hefs_ecrps <- c()
  for(k in 1:nn){
    hefs_mn <- mean(hefs_forward[keysite,,ord_dist[k],ld])
    hefs_cumul <- hefs_forward_cumul[keysite,,ord_dist[k]]
    hefs_ens_pred <-hefs_forward_frac[keysite,,ord_dist[k],]*matrix(rep(hefs_cumul,leads),ncol=leads,byrow=F)
    pred_mn <- mean(hefs_ens_pred[,1])
    mn_scale <- hefs_mn/obs_forward_all_leads_hind[keysite,ord_dist[k],ld] * obs_forward_all_leads_hind[keysite,srt_obs[i],ld] / pred_mn
    hefs_ens_pred <- hefs_ens_pred * mn_scale
    #hefs_ens_pred <- hefs_ens_pred + (pred_mn * mn_scale - pred_mn);hefs_ens_pred[hefs_ens_pred<0]<-0
    samp_hefs_ecrps[k] <- eCRPS(hefs_ens_pred[,ld],obs_forward_all_leads_hind[keysite,srt_obs[i],ld])
  }
  err <- hefs_ecrps - samp_hefs_ecrps
  ecrps_err <- c(ecrps_err,scale(err,center = F))
}

boxplot(ecrps_err,main=paste('mean shift, NN=',nn),ylab='eCRPS error',border='orange',col='white',ylim=c(-3,3))
abline(h=0,col='gray',lty=2,lwd=2)

##############################################END###############################
