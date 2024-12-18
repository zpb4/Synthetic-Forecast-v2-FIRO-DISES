
args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

#/////////////////////////////////////////
#Primary user defined settings

loc = 'YRS'              #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM'
n_samp = 100               #number of samples to generate
keysite_name = 'NBBC1'   #specific site ID for 'keysite' which conditions the kNN sampling
fit_gen_strategy = 'all'  #'all' fits to all available fit data and generate across all available observations
cal_val_setup = '5fold-test'
pcnt_opt = 0.99

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
workrs=5        #number of workers (cores) to use

if(cal_val_setup=='cal'){
  leave_out_years = c()
}

if(cal_val_setup=='val'){
  lv_out_samps = 6          #how many validation years to leave out
  opt_leave_out = readRDS(paste('./data/',loc,'/opt_val_years_samp=',lv_out_samps,'.rds',sep=''))  #optimal validation subset
  leave_out_years = opt_leave_out #years from 'fit' period to leave out of fitting for model validation, use opt value or input vector of water years, e.g. c(1995,2000,2005,etc)
}

if(cal_val_setup=='wet'){
  lv_out_samps = 6          #how many validation years to leave out
  opt_leave_out = readRDS(paste('./data/',loc,'/wet_val_years_samp=',lv_out_samps,'.rds',sep=''))  #optimal validation subset
  leave_out_years = opt_leave_out #years from 'fit' period to leave out of fitting for model validation, use opt value or input vector of water years, e.g. c(1995,2000,2005,etc)
}

if(cal_val_setup=='5fold' | cal_val_setup=='5fold-test'){
  wy = 90:119
  wy_arr = array(NA,c(5,6))
  set.seed(1)
  for(i in 1:5){
    samp = sample(wy,6,replace=F)
    wy_arr[i,] = samp + 1900
    wy = wy[!(wy%in%samp)]
  }
  leave_out_years = wy_arr[idx,]
}

print(leave_out_years)

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#If desiring to fit and generate for specific timeperiods:
#input fit_gen_strategy = 'specify' above and delineate specific dates in 'fit-start/end' and 'gen-start/end'

#date start and end for fitting; date doesn't matter if 'all' selected above
fit_start =   '1989-10-01'       
fit_end =     '2019-09-30'          
#date start for generation; date doesn't matter if 'all' selected above
gen_start =   '1979-10-02'       
gen_end =     '2021-10-01'      # note: generation only possible to end of 'ixx_obs_forward' index
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#////////////////////////////////////////

#------------Synthetic generation----------------------
print(paste('syngen start',Sys.time()))

#rm(list=ls())
#library(MTS)
varx_fun <- MTS::VARX
#library(doParallel)
#library(doMPI)
#library(Rmpi)
library(future)
library(future.apply)
library(future.batchtools)
library(fBasics)

source('./src/syn_gen.R')

#load in the prepared data
load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward_cumul,hefs_forward_frac,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid)

#index for selected keysite
keysite <- which(site_names==keysite_name)

#date/time manipulation
if(fit_gen_strategy=='all'){
  fit_start<-ixx_hefs[1]
  fit_end<-ixx_hefs[length(ixx_hefs)]
  gen_start<-ixx_obs_forward[1]
  gen_end<-ixx_obs_forward[length(ixx_obs_forward)]
}

ixx_gen <- as.POSIXlt(seq(as.Date(gen_start),as.Date(gen_end),by='day'),tz = "UTC") #ixx_obs_forward   or   ixx_hefs

#///////////////////////////////////////////////////////////////////////////////////////////////////
#synthetic generation

syngen_fun <- function(x){message("x : ", x,Sys.time()); out <- syn_gen(x,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,fit_start,fit_end,gen_start,gen_end,leave_out_years,
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

if(cal_val_setup != '5fold' & cal_val_setup!='5fold-test'){
  saveRDS(syn_hefs_forward,file=paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'.rds',sep=''))}
if(cal_val_setup == '5fold' | cal_val_setup =='5fold-test'){
  saveRDS(syn_hefs_forward,file=paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'-',idx,'.rds',sep=''))}

saveRDS(ixx_gen,file=paste('out/',loc,'/ixx_gen.rds',sep=''))
saveRDS(n_samp,file=paste('out/',loc,'/n_samp.rds',sep=''))

print(paste('syngen end',Sys.time()))

rm(list = ls());gc()

###################################################END##################################################