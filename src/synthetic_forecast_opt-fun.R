
syn_opt <- function(obs_opt,obs_fwd_opt,hefs_fwd_opt,opt_strat,opt_samps,idx,loc,keysite_name,
                    cal_val_setup,pcnt,obj_pwr,scale_pwr,dif,lo,sig_a,sig_b,knn_val,knn_pwr){


#/////////////////////////////////////////
#Primary user defined settings
source('../Synthetic-Forecast_Verification/src/forecast_verification_functions.R')

kk <- as.integer(knn_val)

pcnt = pcnt

n_events <- round((1-pcnt) * length(obs_opt))

evt_loc <- order(obs_opt,decreasing=T)[1:n_events]
if(any(evt_loc<(leads+1))){evt_loc = evt_loc[-c(which(evt_loc<(leads+1)))]} #ensure no events too close to beginning of record

evt_yrs <- ixx_opt_wy$year[evt_loc]

if(cal_val_setup=='cal'){
  val_idx = evt_loc}

if(cal_val_setup=='5fold'){
  wy = 90:119
  wy_arr = array(NA,c(5,6))
  val_idx_lst = vector('list',5)
  set.seed(1)
  for(i in 1:5){
    samp = sample(wy,6,replace=F)
    wy_arr[i,] = samp + 1900
    wy = wy[!(wy%in%samp)]
    val_idx_lst[[i]] <- evt_loc[evt_yrs%in%(wy_arr[i,]-1900)]
  }
}

if(cal_val_setup=='5fold-test'){
  wy = 90:119
  wy_arr = array(NA,c(5,6))
  fit_arr = array(0,c(5,4))
  val_idx_lst = vector('list',5)
  set.seed(1)
  samp = sample(1:5,5,replace=F)
  while(any((1:5-samp)==0)==T){
    samp = sample(1:5,5,replace=F)
  }
  fit_arr[,1]<-samp
  set.seed(1)
  for(i in 1:5){
    samp = sample(wy,6,replace=F)
    wy_arr[i,] = samp + 1900
    wy = wy[!(wy%in%samp)]
    fit_arr[i,2:4]<-c(1:5)[-c(i,fit_arr[i,1])]
  }
  evt_loc_lst <- vector('list',5)
  for(i in 1:5){
    val_years = wy_arr[fit_arr[i,1],]
    val_idx = ixx_opt_wy$year%in%(val_years-1900)
    n_events <- round((1-pcnt) * length(which(val_idx==T)))
    evt_loc_sub <- order(obs_opt[val_idx],decreasing=T)[1:n_events]
    evt_loc <- (1:length(obs_opt))[obs_opt%in%obs_opt[val_idx][evt_loc_sub]]
    if(any(evt_loc<(leads+1))){evt_loc = evt_loc[-c(which(evt_loc<(leads+1)))]} #ensure no events too close to beginning of record
    evt_loc_lst[[i]] <- evt_loc
    evt_yrs <- ixx_opt_wy$year[evt_loc]
    val_idx_lst[[i]] <- evt_loc[evt_yrs%in%(val_years-1900)]}
}


#generate samples
if(cal_val_setup=='cal'){
  syn_hefs_val <- array(NA,c(opt_samps,dim(hefs_fwd_opt)[1],n_events,dim(hefs_fwd_opt)[3]))
  cal_idx = ixx_opt_wy$year%in%unique(ixx_opt_wy$year)
  for(s in 1:opt_samps){
    syn_hefs_out <- syn_gen(s,kk,keysite,knn_pwr,scale_pwr,dif,lo,sig_a,sig_b,cal_idx,val_idx,
                              obs_opt,obs_fwd_opt,hefs_fwd_opt,ixx_opt_wy)
    syn_hefs_val[s,,,] <- syn_hefs_out}
  rm(syn_hefs_out)
  gc()
}

if(cal_val_setup=='5fold'){
  syn_hefs_val <- array(NA,c(opt_samps,dim(hefs_fwd_opt)[1],n_events,dim(hefs_fwd_opt)[3]))
  for(s in 1:opt_samps){
    for(i in 1:5){
      val_years = wy_arr[i,]
      val_idx = val_idx_lst[[i]]
      cal_idx = !ixx_opt_wy$year%in%(val_years-1900)
      syn_hefs_out <- syn_gen(s,kk,keysite,knn_pwr,scale_pwr,dif,lo,sig_a,sig_b,cal_idx,val_idx,
                              obs_opt,obs_fwd_opt,hefs_fwd_opt,ixx_opt_wy)
      syn_hefs_val[s,,evt_loc%in%val_idx,] <- syn_hefs_out}
  }
  rm(syn_hefs_out)
  gc()
}

if(cal_val_setup=='5fold-test'){
  #evt_loc = evt_loc_lst[[idx]]
  val_idx = val_idx_lst[[idx]]
  cal_years = sort(wy_arr[fit_arr[idx,2:4],])
  cal_idx = ixx_opt_wy$year%in%(cal_years-1900)
  syn_hefs_val <- array(NA,c(opt_samps,dim(hefs_fwd_opt)[1],length(val_idx),dim(hefs_fwd_opt)[3]))
  for(s in 1:opt_samps){
    syn_hefs_out <- syn_gen(s,kk,keysite,knn_pwr,scale_pwr,dif,lo,sig_a,sig_b,cal_idx,val_idx,
                            obs_opt,obs_fwd_opt,hefs_fwd_opt,ixx_opt_wy)
    syn_hefs_val[s,,,] <- syn_hefs_out}
rm(syn_hefs_out)
gc()}

##obs_date_loc <- order(obs_val,decreasing=TRUE)[1:num_events_crps]  #index for maximum observation
#remove any obs dates that are too close to start of timeseries for lead 15 indexing
##if(any(obs_date_loc<=leads)){
  ##obs_date_loc <- obs_date_loc[-c(which(obs_date_loc<=leads))]}
#match obs dates to hefs/shefs indexing
##if(cal_val_setup!='cal' & cal_val_setup!='5fold'){
  ##obs_dates <- ixx_obs[ixx_obs%in%ixx_hefs][val_idx][obs_date_loc]}
##if(cal_val_setup=='cal' | cal_val_setup=='5fold'){
  ##obs_dates <- ixx_hefs[obs_date_loc]}
##obs_major_events <- obs_val[obs_date_loc]
##hefs_date_loc <- match(obs_dates,ixx_hefs)

#remove any VAL events that are not fully contained within the VAL index
##shefs_events <- syn_hefs_val[1,1,(hefs_date_loc-15),1]
##shefs_events <- syn_hefs_val[1,1,1,(hefs_date_loc-15),1]
##if(anyNA(shefs_events)==T){
  ##hefs_date_loc <- hefs_date_loc[!is.na(shefs_events)]
  ##obs_date_loc <- obs_date_loc[!is.na(shefs_events)]
##}

if(cal_val_setup!='5fold-test'){
hefs_ecrps_vec<-array(NA,c(length(evt_loc),leads))
shefs_ecrps_vec<-array(NA,c(opt_samps,length(evt_loc),leads))

for(s in 1:opt_samps){
  for(ld in 1:leads){
    for(i in 1:length(evt_loc)){
      hefs_idx <- evt_loc[i]-ld #need to back up by lds[ld] because forecasts are in 'forward' format
      HEFS <- hefs_fwd_opt[,hefs_idx,ld]
      hefs_ecrps_vec[i,ld] <- eCRPS(HEFS,obs_fwd_opt[hefs_idx,ld])
      SYN_HEFS <- syn_hefs_val[s,,i,ld]
      shefs_ecrps_vec[s,i,ld] <- eCRPS(SYN_HEFS,obs_fwd_opt[hefs_idx,ld])
    }
  }
}
}

if(cal_val_setup=='5fold-test'){
hefs_ecrps_vec<-array(NA,c(length(val_idx),leads))
shefs_ecrps_vec<-array(NA,c(opt_samps,length(val_idx),leads))

for(s in 1:opt_samps){
  for(ld in 1:leads){
    for(i in 1:length(val_idx)){
      hefs_idx <- val_idx[i]-ld #need to back up by lds[ld] because forecasts are in 'forward' format
      HEFS <- hefs_fwd_opt[,hefs_idx,ld]
      hefs_ecrps_vec[i,ld] <- eCRPS(HEFS,obs_fwd_opt[hefs_idx,ld])
      SYN_HEFS <- syn_hefs_val[s,,i,ld]
      shefs_ecrps_vec[s,i,ld] <- eCRPS(SYN_HEFS,obs_fwd_opt[hefs_idx,ld])
    }
  }
}
}
#MSE calculation
##MSE_HEFS <- c()
##MSE_syn_HEFS <- c()


##for (ld in 1:length(lds)) {
  ##obs_forward_idx <- obs_date_loc-ld
  ##hefs_idx <- hefs_date_loc-ld
  ##syn_idx <- hefs_date_loc-ld
  
  #ensemble mean
  ##HEFS_ens_mean <- apply(hefs_val[,hefs_idx,ld],FUN=mean,2)
  ##syn_HEFS_ens_mean <- apply(syn_hefs_val[1,,syn_idx,ld],2,FUN=mean)
  ##obs_events <- obs_fwd_val[obs_forward_idx,ld]
  
  ##MSE_HEFS[ld] <- mean(sqrt((HEFS_ens_mean - obs_events)^2))
  ##MSE_syn_HEFS[ld] <-mean(sqrt((syn_HEFS_ens_mean - obs_events)^2))   
##}

##rm(HEFS,SYN_HEFS,HEFS_ens_mean,syn_HEFS_ens_mean,syn_hefs_val,hefs_val,obs_val)
rm(HEFS,SYN_HEFS,syn_hefs_val)
gc()

my_power <- obj_pwr
w <- 1:leads
decay <- (w^my_power / sum(w^my_power))

#median diff approach
if(opt_strat=='ecrps-med'){
hefs_ecrps_med <- apply(hefs_ecrps_vec,2,median)
shefs_ecrps_med <- apply(shefs_ecrps_vec,3,median)
sse_ecrps_int <- (decay*(hefs_ecrps_med-shefs_ecrps_med)/hefs_ecrps_med)^2
sse_ecrps <- sum(sse_ecrps_int)}

#ks test approach
if(opt_strat=='ecrps-ks'){
ks_vec = sapply(1:leads, function(x){out = ks.test(hefs_ecrps_vec[,x],as.vector(shefs_ecrps_vec[,,x]),alternative='two.sided',exact=F);return(out$statistic)})
sse_ecrps <- 1e3*sum((decay*ks_vec)^2)}

if(opt_strat=='ecrps-dts'){
  dts_vec = sapply(1:leads, function(x){
    out = two_sample(hefs_ecrps_vec[,x]/sd(hefs_ecrps_vec[,x]),as.vector(shefs_ecrps_vec[,,x])/sd(hefs_ecrps_vec[,x]),keep.boots = T)
    return(out[1])})
sse_ecrps <- mean((decay*dts_vec)^2)}

##diff_vec_ecrps <- hefs_ecrps_vec - shefs_ecrps_vec
##sc_diff_vec_ecrps <- scale(diff_vec_ecrps,center = F)

##mse_ecrps_int <- apply(sc_diff_vec_ecrps,2,function(x){out<-mean(sqrt(x^2))})

##mse_ecrps <- mean(decay*mse_ecrps_int)

#mse_tot <- mse_ecrps + mean(decay*sqrt((MSE_HEFS - MSE_syn_HEFS)^2))

return(sse_ecrps)
}


###################################################END##################################################