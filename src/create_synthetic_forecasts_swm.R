
args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

#load packages, functions, data
library(future)
library(future.apply)
#library(future.batchtools)
#library(fBasics)

source('./src/syn_gen_swm.R')

#/////////////////////////////////////////
#Primary user defined settings

loc = 'YRS'              #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM'
n_samp = 1               #number of samples to generate
swm_samps = 1
keysite_name = 'ORDC1'   #specific site ID for 'keysite' which conditions the kNN sampling
fit_gen_strategy = 'all'  #'all' fits to all available fit data and generate across all available observations
cal_val_setup = 'cal'
pcnt_opt = 0.99

#swm runs
use_bcf=F
k_swm=106
set.seed(1)

#load in the prepared data
load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward_cumul,hefs_forward_frac,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid)

#index for selected keysite
keysite <- which(site_names==keysite_name)

if(cal_val_setup!='5fold-test'){
  opt_params <- readRDS(paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt_opt,'_',keysite_name,'.rds',sep=''))}
if(cal_val_setup=='5fold-test'){
  opt_params <- readRDS(paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt_opt,'_',keysite_name,'-',idx,'.rds',sep=''))}

#orimary hyperparameters
kk <- 30           #sampling window
knn_pwr <- 0     #larger negative values weights early lead times higher for sampling

scale_pwr = opt_params[2] 
hi = opt_params[3]  
lo = opt_params[4]   
sig_a = opt_params[5] 
sig_b = opt_params[6] 

#parallel configuration
parllel=T          #use parallel processing
workrs=10        #number of workers (cores) to use

ix_swm<-seq(as.Date('2001-01-01'),as.Date('3008-01-08'),'day')
ix2_swm<-as.POSIXlt(ix_swm)

synq_sim_arr<-readRDS(paste('../syn-forecast_stochastic/out/synq_sim_arr_use-bcf=',use_bcf,'_k=',k_swm,'.rds',sep=''))
synq_sim_arr_cc<-readRDS(paste('../syn-forecast_stochastic/out/synq_sim_arr_cc_use-bcf=',use_bcf,'_k=',k_swm,'.rds',sep=''))

swm_samp<-synq_sim_arr[sample(1:dim(synq_sim_arr)[1],swm_samps),]
swm_samp_cc<-synq_sim_arr_cc[sample(1:dim(synq_sim_arr_cc)[1],swm_samps),]

rm(synq_sim_arr,synq_sim_arr_cc)

obs_forward_all_leads_hind <- obs_forward_all_leads_hind[keysite,,,drop=F]

hefs_forward <- hefs_forward[keysite,,,,drop=F]

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#If desiring to fit and generate for specific timeperiods:
#input fit_gen_strategy = 'specify' above and delineate specific dates in 'fit-start/end' and 'gen-start/end'

#date start and end for fitting; date doesn't matter if 'all' selected above
fit_start =   '1989-10-01'       
fit_end =     '2019-09-30'          
#date start for generation; date doesn't matter if 'all' selected above
gen_start =   '1989-10-01'      
gen_end =     '2021-10-01'      # note: generation only possible to end of 'ixx_obs_forward' index
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#////////////////////////////////////////

#------------Synthetic generation----------------------
print(paste('syngen start',Sys.time()))

#date/time manipulation
if(fit_gen_strategy=='all'){
  fit_start<-ixx_hefs[1]
  fit_end<-ixx_hefs[length(ixx_hefs)]
  gen_start<-ixx_obs_forward[1]
  gen_end<-ixx_obs_forward[length(ixx_obs_forward)]
}

ixx_gen <- as.POSIXlt(seq(as.Date(gen_start),as.Date(gen_end),by='day'),tz = "UTC") #ixx_obs_forward   or   ixx_hefs


#///////////////////////////////////////////////////////////////////////////////////////////////////
#synthetic generation for reference synthetic forecast based on actual obs
if(idx==1){
print(paste('swm-ref start',Sys.time()))
  
obs_forward_all_leads<-obs_forward_all_leads[keysite,,,drop=F]

keysite = 1
n_sites = 1

syngen_fun <- function(x){message("x : ", x,Sys.time()); out <- syn_gen(x,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,fit_start,fit_end,gen_start,gen_end,
                                                                        obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                                                                        ixx_obs_forward)}

n_samp = 100
syn_hefs_forward <- array(NA,c(n_samp,n_sites,n_ens,length(ixx_gen),leads))

if(parllel==TRUE){
  plan(multicore,workers=workrs)
  fut_vec <- future_lapply(as.list(1:n_samp),syngen_fun,future.seed=TRUE)
  
  #compile to multidimensional array and save
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-fut_vec[[m]]
  }
}

#sequential processing if no parallel flag
if(parllel==FALSE){
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-syngen_fun(m)
    print(paste('m=',m,Sys.time()))
  }
}


saveRDS(syn_hefs_forward,file=paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_swm-ref.rds',sep=''))

print(paste('swm-ref end',Sys.time()))
}
#///////////////////////////////////////////////////////////////////////////////////////////////////
#synthetic generation for SWM
#these are the dates that PRECEDE the 15 day cumulative totals, i.e., 
#on date t, here is the 15 day total over the NEXT 15 days
if(idx==2){
print(paste('swm start',Sys.time()))
  
ixx_obs <- ix2_swm
ixx_obs_forward <- ixx_obs[1:dim(obs_forward_all_leads)[2]]
n_obs_forward <- length(ixx_obs_forward)
n_obs <- length(ixx_obs)
n_sites <- 1
n_samp = 1
keysite = 1

#date/time manipulation
if(fit_gen_strategy=='all'){
  fit_start<-ixx_hefs[1]
  fit_end<-ixx_hefs[length(ixx_hefs)]
  gen_start<-ixx_obs_forward[1]
  gen_end<-ixx_obs_forward[length(ixx_obs_forward)]
}

ixx_gen <- as.POSIXlt(seq(as.Date(gen_start),as.Date(gen_end),by='day'),tz = "UTC") #ixx_obs_forward   or   ixx_hefs


obs_forward_all_leads <- array(NA,c(1,length(1:(length(ix2_swm)-leads)),leads))

for(i in 1:(length(ix2_swm)-leads)) {
  obs_forward_all_leads[1,i,] <- swm_samp[(i+1):(i+leads)]
}

syngen_fun <- function(x){message("x : ", x,Sys.time()); out <- syn_gen(x,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,fit_start,fit_end,gen_start,gen_end,
                                                                        obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                                                                        ixx_obs_forward)}

syn_hefs_forward <- array(NA,c(n_samp,n_sites,n_ens,length(ixx_gen),leads))

if(parllel==TRUE){
  plan(multicore,workers=workrs)
  fut_vec <- future_lapply(as.list(1:n_samp),syngen_fun,future.seed=TRUE)

#compile to multidimensional array and save
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-fut_vec[[m]]
  }
}

#sequential processing if no parallel flag
if(parllel==FALSE){
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-syngen_fun(m)
    print(paste('m=',m,Sys.time()))
  }
}

saveRDS(syn_hefs_forward,file=paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_swm.rds',sep=''))

print(paste('swm end',Sys.time()))
}

if(idx==3){
print(paste('swm-cc start',Sys.time()))
  
ixx_obs <- ix2_swm
ixx_obs_forward <- ixx_obs[1:dim(obs_forward_all_leads)[2]]
n_obs_forward <- length(ixx_obs_forward)
n_obs <- length(ixx_obs)
n_sites <- 1
n_samp = 1
keysite = 1
  
#date/time manipulation
if(fit_gen_strategy=='all'){
  fit_start<-ixx_hefs[1]
  fit_end<-ixx_hefs[length(ixx_hefs)]
  gen_start<-ixx_obs_forward[1]
  gen_end<-ixx_obs_forward[length(ixx_obs_forward)]
}
  
ixx_gen <- as.POSIXlt(seq(as.Date(gen_start),as.Date(gen_end),by='day'),tz = "UTC") #ixx_obs_forward   or   ixx_hefs
  
obs_forward_all_leads <- array(NA,c(1,length(1:(length(ix2_swm)-leads)),leads))

for(i in 1:(length(ix2_swm)-leads)) {
  obs_forward_all_leads[1,i,] <- swm_samp_cc[(i+1):(i+leads)]
}

syngen_fun <- function(x){message("x : ", x,Sys.time()); out <- syn_gen(x,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,fit_start,fit_end,gen_start,gen_end,
                                                                        obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                                                                        ixx_obs_forward)}

syn_hefs_forward <- array(NA,c(n_samp,n_sites,n_ens,length(ixx_gen),leads))

if(parllel==TRUE){
  plan(multicore,workers=workrs)
  fut_vec <- future_lapply(as.list(1:n_samp),syngen_fun,future.seed=TRUE)
  
  #compile to multidimensional array and save
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-fut_vec[[m]]
  }
}

#sequential processing if no parallel flag
if(parllel==FALSE){
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-syngen_fun(m)
    print(paste('m=',m,Sys.time()))
  }
}

saveRDS(syn_hefs_forward,file=paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_swm-cc.rds',sep=''))


saveRDS(ixx_gen,file=paste('out/',loc,'/ixx_gen_swm.rds',sep=''))
saveRDS(n_samp,file=paste('out/',loc,'/n_samp_swm.rds',sep=''))

print(paste('swm-cc end',Sys.time()))
}

print(paste('syngen end',Sys.time()))

rm(list = ls());gc()

###################################################END##################################################