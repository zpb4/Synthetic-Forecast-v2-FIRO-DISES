
rm(list=ls())
library(MTS)
library(doParallel)

source('./src/syn_gen.R')

print(paste('syngen start',Sys.time()))
#parallelization code
n.cores <- parallel::detectCores()
my.cluster<-parallel::makeCluster(n.cores,type = 'PSOCK')
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

#------------Synthetic generation----------------------

#load in the prepared data
load("out/data_prep_rdata.RData")

#number of synthetic ensembles
n_samp <- 100
#k for KNN sampling
kk <- 20
#the site to use as the key site for resampling; generally main reservoir inflow site
keysite <- which(site_names=="NHGC1") 
parllel=F

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

save.image("./out/model-fit_rdata.RData")
#-------------------------------synthetic forecast generation-------------------
#array to store 
syn_hefs_forward <- array(NA,c(n_samp,n_sites,n_ens,length(ixx_gen),leads))

if(parllel==TRUE){
  #parallel 'foreach' implementation of synthetic forecast generation
  syn_hefs_out<-foreach(m = 1:n_samp,.inorder=F)%dopar%{
    syn_hefs<- syn_gen(1,kk,keysite,fit_start,fit_end,gen_start,gen_end,leave_out_years)
    return(syn_hefs)
  }

  #compile to multidimensional array and save
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-syn_hefs_out[[m]]
  }
}

if(parllel==FALSE){
  for(m in 1:n_samp){
    syn_hefs_forward[m,,,,]<-syn_gen(1,kk,keysite,fit_start,fit_end,gen_start,gen_end,leave_out_years)
  }
}

saveRDS(syn_hefs_forward,file='out/syn_hefs_forward.rds')
saveRDS(ixx_gen,file='out/ixx_gen.rds')
saveRDS(n_samp,file='out/n_samp.rds')

rm(list = ls());gc()

print(paste('syngen end',Sys.time()))

stopCluster()