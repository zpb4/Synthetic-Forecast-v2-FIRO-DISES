
rm(list=ls())
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

print(paste('syngen start',Sys.time()))

#------------Synthetic generation----------------------

#load in the prepared data
load("out/data_prep_rdata.RData")

#parameterization
parm <- 'c'
# 'a'  : kk=20; knn_pwr=-0.5; keysite='NHGC1'; original method
# 'b'  : kk=20; knn_pwr=-0.5; keysite='NHGC1'; ss method 5 apr (empirical resid sampling)
# 'c'  : kk=20; knn_pwr=-0.5; keysite='NHGC1'; 20 preds (sq_rnk,sq_diff_mn,sq_diff_sd,sq_diff_kurt,sq_diff_skew,sq_obs_seq)
# 'd'  : kk=20; knn_pwr=-0.5; keysite='NHGC1'; lead 1-3 only + cumul sum


#number of synthetic ensembles
n_samp <- 10
nsamp_lst = as.list(1:n_samp)
#k for KNN sampling
kk <- 20
knn_pwr <- -0.5  #larger negative values weights early lead times higher for sampling
#the site to use as the key site for resampling; generally main reservoir inflow site
keysite <- which(site_names=="NHGC1")
parllel=F
workrs=35

#either input 'all' to fit to all available fit data and generate across all available observations
#or input 'specify' and delineate specific dates in 'fit-start/end' and 'gen-start/end'
fit_gen_strategy<- 'all'

#date start and end for fitting; date doesn't matter if 'all' selected above
fit_start<-   '1989-10-01'       #  '1989-10-01'    
fit_end<-     '2019-09-30'       #  '2019-09-30'    
#date start for generation; date doesn't matter if 'all' selected above
gen_start<-   '1979-10-02'      #  '1979-10-02'     
gen_end<-     '2021-10-01'      #  '2021-10-01' ; note: generation only possible to end of 'ixx_obs_forward' index

#years from 'fit' period to leave out of fitting for model validation
leave_out_years <- c(1995,2000,2005,2010,2015)

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
syngen_fun <- function(x){message("x : ", x,Sys.time()); out <- syn_gen(x,kk,keysite,knn_pwr,fit_start,fit_end,gen_start,gen_end,leave_out_years,mpi,obs_forward_all_leads_gen,obs_forward_all_leads_hind,hefs_forward_cumul,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,hefs_forward_frac,ixx_hefs,ixx_obs_forward,varx_fun); return(out)}

syn_hefs_forward <- array(NA,c(n_samp,n_sites,n_ens,length(ixx_gen),leads))

if(parllel==TRUE){
  plan(multicore,workers=workrs)
  fut_vec <- future_lapply(nsamp_lst,syngen_fun,future.seed=TRUE)

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

 
saveRDS(syn_hefs_forward,file=paste('out/syn_hefs_forward-',parm,'.rds',sep=''))
saveRDS(ixx_gen,file='out/ixx_gen.rds')
saveRDS(n_samp,file='out/n_samp.rds')

print(paste('syngen end',Sys.time()))

rm(list = ls());gc()

###################################################END##################################################