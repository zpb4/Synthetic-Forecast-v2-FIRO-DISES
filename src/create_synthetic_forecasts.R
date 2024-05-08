

#/////////////////////////////////////////
#Primary user defined settings

loc = 'LAM'              #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM'
n_samp = 10             #number of samples to generate
keysite_name = 'LAMC1'   #specific site ID for 'keysite' which conditions the kNN sampling
fit_gen_strategy = 'all'  #'all' fits to all available fit data and generate across all available observations

lv_out_samps = 6          #how many validation years to leave out
opt_leave_out = readRDS(paste('./data/',loc,'/opt_val_years_samp=',lv_out_samps,'.rds',sep=''))  #optimal validation subset
leave_out_years = opt_leave_out #years from 'fit' period to leave out of fitting for model validation, use opt value or input vector of water years, e.g. c(1995,2000,2005,etc)
#leave_out_years = c()

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

#/////////////////////////////////////////////////
#Other user defined settings that can be modified
correct_leads = F
lds_to_correct = 1:(leads-1)
fix_order = F
use_ar = F
refrac = F

#parameterization
parm <- 'a'

# 'a'  : kk=20; knn_pwr=-1; keysite='NHGC1'; obs0 + obs fwd sample, use_ar=F, fix_order = F, correct_leads = F 
# 'b'  : kk=10; knn_pwr=-1; keysite='NHGC1'; obs0 + obs fwd sample, use_ar=F, fix_order = F, correct_leads = F 
# 'c'  : kk=10; knn_pwr=-1; keysite='NHGC1'; obs0 + obs fwd sample, use_ar=F, fix_order = F, correct_leads = T, refrac=T; correct all leads, no decay 
# 'd'  : kk=10; knn_pwr=-1; keysite='NHGC1'; obs0 + obs fwd sample, use_ar=F, fix_order = F, correct_leads = T, refrac=F; correct all leads, no decay 
# 'e'  : kk=10; knn_pwr=-1; keysite='NHGC1'; obs0 + obs fwd sample, use_ar=F, fix_order = F, correct_leads = T, refrac=T; correct all leads, decay=-1.5 + 1.25 scale 
# 'f'  : kk=10; knn_pwr=-1; keysite='NHGC1'; obs0 + obs fwd sample, use_ar=F, fix_order = F, correct_leads = T, refrac=F; correct all leads, decay=-1.5 + 1.25 scale 

#k for KNN sampling
kk <- 20
knn_pwr <- -1  #larger negative values weights early lead times higher for sampling
parllel=T
workrs=10

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

#save.image("./out/model-fit_rdata.RData")
#-------------------------------synthetic forecast generation-------------------
#define function for parallel apply 
syngen_fun <- function(x){message("x : ", x,Sys.time()); out <- syn_gen(x,kk,keysite,knn_pwr,fit_start,fit_end,gen_start,gen_end,leave_out_years,mpi,obs_forward_all_leads_gen,obs_forward_all_leads_hind,hefs_forward_cumul,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,hefs_forward_frac,ixx_hefs,
                                                                        ixx_obs_forward,varx_fun,correct_leads,lds_to_correct,fix_order,use_ar,refrac); return(out)}

syn_hefs_forward <- array(NA,c(n_samp,n_sites,n_ens,length(ixx_gen),leads))

if(parllel==TRUE){
  plan(multicore,workers=workrs)
  fut_vec <- future_lapply(as.list(1:n_samp),syngen_fun,future.seed=TRUE)

#compile to multidimensional array and save
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-fut_vec[[m]]
  }
}


if(parllel==FALSE){
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-syngen_fun(m)
    print(paste('m=',m,Sys.time()))
  }
}

 
saveRDS(syn_hefs_forward,file=paste('out/',loc,'/syn_hefs_forward-',parm,'.rds',sep=''))
saveRDS(ixx_gen,file=paste('out/',loc,'/ixx_gen.rds',sep=''))
saveRDS(n_samp,file=paste('out/',loc,'/n_samp.rds',sep=''))

print(paste('syngen end',Sys.time()))

rm(list = ls());gc()

###################################################END##################################################