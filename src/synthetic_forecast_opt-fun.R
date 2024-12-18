
syn_opt <- function(idx,loc,keysite_name,cal_val_setup,kk,knn_pwr,pcnt,obj_pwr,scale_pwr,diff,lo,sig_a,sig_b){

#/////////////////////////////////////////
#Primary user defined settings
source('./src/syn_gen_opt.R')
source('../Synthetic-Forecast_Verification/src/forecast_verification_functions.R')

wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

#primary hyperparameters
##loc = 'YRS'              #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM'             #number of samples to generate
##keysite_name = 'ORDC1'   #specific site ID for 'keysite' which conditions the kNN sampling
fit_gen_strategy = 'all'  #'all' fits to all available fit data and generate across all available observations
##cal_val_setup = 'val'
##kk <- 30           #sampling window
##knn_pwr <- 0      #larger negative values weights early lead times higher for sampling
##scale_pwr = -1  
##hi = 5  
##lo = 2      
##sig_a = 1  
##sig_b = 0 

#load in the prepared data
load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward_cumul,hefs_forward_frac,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid)
gc()

if(cal_val_setup=='cal'){
  leave_out_years = c()
}

if(cal_val_setup=='val'|cal_val_setup=='val-test'){
  lv_out_samps = 6          #how many validation years to leave out
  opt_leave_out = readRDS(paste('./data/',loc,'/opt_val_years_samp=',lv_out_samps,'.rds',sep=''))  #optimal validation subset
  leave_out_years = opt_leave_out #years from 'fit' period to leave out of fitting for model validation, use opt value or input vector of water years, e.g. c(1995,2000,2005,etc)
}

if(cal_val_setup=='wet'|cal_val_setup=='wet-test'){
  lv_out_samps = 6          #how many validation years to leave out
  opt_leave_out = readRDS(paste('./data/',loc,'/wet_val_years_samp=',lv_out_samps,'.rds',sep=''))  #optimal validation subset
  leave_out_years = opt_leave_out #years from 'fit' period to leave out of fitting for model validation, use opt value or input vector of water years, e.g. c(1995,2000,2005,etc)
}

if(cal_val_setup=='5fold'){
  wy = 90:119
  wy_arr = array(NA,c(5,6))
  set.seed(1)
  for(i in 1:5){
    samp = sample(wy,6,replace=F)
    wy_arr[i,] = samp + 1900
    wy = wy[!(wy%in%samp)]
  }
}

if(cal_val_setup=='5fold-test'){
  wy = 90:119
  wy_arr = array(NA,c(5,6))
  fit_arr = array(0,c(5,4))
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
}

#index for selected keysite
keysite <- which(site_names==keysite_name)

#date/time manipulation
if(fit_gen_strategy=='all'){
  ixx_hefs <- ixx_hefs[ixx_hefs%in%ixx_obs_forward]
  fit_start<-ixx_hefs[1]
  fit_end<-ixx_hefs[length(ixx_hefs)]
  ixx_hefs_wy<-wy_fun(ixx_hefs)
}


if(cal_val_setup=='cal'){
  val_years = leave_out_years
  cal_idx = !(unique(ixx_hefs_wy$year)%in%(val_years-1900))
  cal_years = unique(ixx_hefs_wy$year)[cal_idx]+1900
  syn_hefs_val <- syn_gen(1,kk,keysite,knn_pwr,scale_pwr,diff,lo,sig_a,sig_b,fit_start,fit_end,cal_years,val_years,
                 obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                 ixx_obs_forward)
  hefs_val <- hefs_forward[keysite,,,]
  obs_val <-obs[ixx_obs%in%ixx_hefs,keysite]
  obs_fwd_val <-obs_forward_all_leads[keysite,ixx_obs_forward%in%ixx_hefs,]
}

if(cal_val_setup=='val'|cal_val_setup=='wet'){
  val_years = leave_out_years
  cal_idx = !(unique(ixx_hefs_wy$year)%in%(val_years-1900))
  cal_years = unique(ixx_hefs_wy$year)[cal_idx]+1900
  syn_hefs_val <- array(NA,dim(hefs_forward[keysite,,,,drop=F]))
  syn_hefs_out <- syn_gen(1,kk,keysite,knn_pwr,scale_pwr,diff,lo,sig_a,sig_b,fit_start,fit_end,cal_years,val_years,
                 obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                 ixx_obs_forward)
  ixx_hefs_wy <- wy_fun(ixx_hefs)
  val_idx <- which(ixx_hefs_wy$year%in%(val_years-1900)==T)
  syn_hefs_val[,,val_idx,] <- syn_hefs_out
  hefs_val <- hefs_forward[keysite,,,]
  obs_val <-obs[ixx_obs%in%ixx_hefs,keysite][val_idx]
  obs_fwd_val <-obs_forward_all_leads[keysite,ixx_obs_forward%in%ixx_hefs,][val_idx,]
  rm(syn_hefs_out)
  gc()}

if(cal_val_setup=='val-test'|cal_val_setup=='wet-test'){
  tst_years = leave_out_years
  cal_idx = !(unique(ixx_hefs_wy$year)%in%(tst_years-1900))
  cal_years = unique(ixx_hefs_wy$year)[cal_idx]+1900
  val_years = sample(cal_years,6,replace=F)
  cal_years = cal_years[!(cal_years%in%val_years)]
  syn_hefs_val <- array(NA,dim(hefs_forward[keysite,,,,drop=F]))
  syn_hefs_out <- syn_gen(1,kk,keysite,knn_pwr,scale_pwr,diff,lo,sig_a,sig_b,fit_start,fit_end,cal_years,val_years,
                          obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                          ixx_obs_forward)
  ixx_hefs_wy <- wy_fun(ixx_hefs)
  val_idx <- which(ixx_hefs_wy$year%in%(val_years-1900)==T)
  syn_hefs_val[,,val_idx,] <- syn_hefs_out
  hefs_val <- hefs_forward[keysite,,,]
  obs_val <-obs[ixx_obs%in%ixx_hefs,keysite][val_idx]
  obs_fwd_val <-obs_forward_all_leads[keysite,ixx_obs_forward%in%ixx_hefs,][val_idx,]
  rm(syn_hefs_out)
  gc()}

if(cal_val_setup=='5fold-test'){
  val_years = wy_arr[fit_arr[idx,1],]
  cal_years = sort(wy_arr[fit_arr[idx,2:4],])
  syn_hefs_val <- array(NA,dim(hefs_forward[keysite,,,,drop=F]))
  syn_hefs_out <- syn_gen(1,kk,keysite,knn_pwr,scale_pwr,diff,lo,sig_a,sig_b,fit_start,fit_end,cal_years,val_years,
                          obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                          ixx_obs_forward)
  ixx_hefs_wy <- wy_fun(ixx_hefs)
  val_idx <- which(ixx_hefs_wy$year%in%(val_years-1900)==T)
  syn_hefs_val[,,val_idx,] <- syn_hefs_out
  hefs_val <- hefs_forward[keysite,,,]
  obs_val <-obs[ixx_obs%in%ixx_hefs,keysite][val_idx]
  obs_fwd_val <-obs_forward_all_leads[keysite,ixx_obs_forward%in%ixx_hefs,][val_idx,]
  rm(syn_hefs_out)
  gc()}

if(cal_val_setup=='5fold'){
  syn_hefs_val <- hefs_forward[keysite,,,,drop=FALSE]
  for(i in 1:5){
    leave_out_years = wy_arr[i,]
    val_years = leave_out_years
    cal_idx = !(unique(ixx_hefs_wy$year)%in%(val_years-1900))
    cal_years = unique(ixx_hefs_wy$year)[cal_idx]+1900
    syn_hefs_out <- syn_gen(1,kk,keysite,knn_pwr,scale_pwr,diff,lo,sig_a,sig_b,fit_start,fit_end,cal_years,val_years,
                 obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                 ixx_obs_forward)
    ixx_hefs_wy <- wy_fun(ixx_hefs)
    val_idx <- which(ixx_hefs_wy$year%in%(val_years-1900)==T)
    syn_hefs_val[,,val_idx,] <- syn_hefs_out}
  hefs_val <- hefs_forward[keysite,,,]
  obs_val <-obs[ixx_obs%in%ixx_hefs,keysite]
  obs_fwd_val <-obs_forward_all_leads[keysite,ixx_obs_forward%in%ixx_hefs,]
  rm(syn_hefs_out)
  gc()
}

rm(hefs_forward,obs)
gc()

pcnt = pcnt
lds = 1:leads

num_events_crps <- round((1-pcnt) * length(obs_val))

obs_date_loc <- order(obs_val,decreasing=TRUE)[1:num_events_crps]  #index for maximum observation
#remove any obs dates that are too close to start of timeseries for lead 15 indexing
if(any(obs_date_loc<=leads)){
  obs_date_loc <- obs_date_loc[-c(which(obs_date_loc<=leads))]}
#match obs dates to hefs/shefs indexing
if(cal_val_setup!='cal' & cal_val_setup!='5fold'){
  obs_dates <- ixx_obs[ixx_obs%in%ixx_hefs][val_idx][obs_date_loc]}
if(cal_val_setup=='cal' | cal_val_setup=='5fold'){
  obs_dates <- ixx_hefs[obs_date_loc]}
obs_major_events <- obs_val[obs_date_loc]
hefs_date_loc <- match(obs_dates,ixx_hefs)

#remove any VAL events that are not fully contained within the VAL index
shefs_events <- syn_hefs_val[1,1,(hefs_date_loc-15),1]
if(anyNA(shefs_events)==T){
  hefs_date_loc <- hefs_date_loc[!is.na(shefs_events)]
  obs_date_loc <- obs_date_loc[!is.na(shefs_events)]
}


hefs_ecrps_vec<-array(NA,c(length(obs_date_loc),length(lds)))
shefs_ecrps_vec<-array(NA,c(length(obs_date_loc),length(lds)))

for(ld in 1:length(lds)){
  for(i in 1:length(obs_date_loc)){
    hefs_idx <- hefs_date_loc[i]-ld #need to back up by lds[ld] because forecasts are in 'forward' format
    syn_idx <- hefs_date_loc[i]-ld #need to back up by lds[ld] because forecasts are in 'forward' format
    HEFS <- hefs_val[,hefs_idx,ld]
    hefs_ecrps_vec[i,ld] <- eCRPS(HEFS,obs_major_events[i])
    SYN_HEFS <- syn_hefs_val[1,,syn_idx,ld]
    shefs_ecrps_vec[i,ld] <- eCRPS(SYN_HEFS,obs_major_events[i])
  }
}

#MSE calculation
MSE_HEFS <- c()
MSE_syn_HEFS <- c()


for (ld in 1:length(lds)) {
  obs_forward_idx <- obs_date_loc-ld
  hefs_idx <- hefs_date_loc-ld
  syn_idx <- hefs_date_loc-ld
  
  #ensemble mean
  HEFS_ens_mean <- apply(hefs_val[,hefs_idx,ld],FUN=mean,2)
  syn_HEFS_ens_mean <- apply(syn_hefs_val[1,,syn_idx,ld],2,FUN=mean)
  obs_events <- obs_fwd_val[obs_forward_idx,ld]
  
  MSE_HEFS[ld] <- mean(sqrt((HEFS_ens_mean - obs_events)^2))
  MSE_syn_HEFS[ld] <-mean(sqrt((syn_HEFS_ens_mean - obs_events)^2))   
}

rm(HEFS,SYN_HEFS,HEFS_ens_mean,syn_HEFS_ens_mean,syn_hefs_val,hefs_val,obs_val)
gc()

#hefs_ecrps_med <- apply(hefs_ecrps_vec,2,median)
#shefs_ecrps_med <- apply(shefs_ecrps_vec,2,median)

#sse_ecrps <- mean(sqrt((hefs_ecrps_med-shefs_ecrps_med)^2))

my_power <- obj_pwr
w <- 1:leads
decay <- (w^my_power / sum(w^my_power))

diff_vec_ecrps <- hefs_ecrps_vec - shefs_ecrps_vec
sc_diff_vec_ecrps <- scale(diff_vec_ecrps,center = F)

mse_ecrps_int <- apply(sc_diff_vec_ecrps,2,function(x){out<-mean(sqrt(x^2))})

mse_ecrps <- mean(decay*mse_ecrps_int)

mse_tot <- mse_ecrps + mean(decay*sqrt((MSE_HEFS - MSE_syn_HEFS)^2))

return(mse_ecrps)
}


###################################################END##################################################