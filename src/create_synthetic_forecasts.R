
rm(list=ls())
library(MTS)
library(doParallel)


source('./src/syn_gen.R')

#parallelization code
parallel::detectCores()
n.cores <- parallel::detectCores()
my.cluster<-parallel::makeCluster(n.cores,type = 'PSOCK')
print(my.cluster)
doParallel::registerDoParallel(cl = my.cluster)
foreach::getDoParRegistered()

#------------Synthetic generation----------------------

#load in the prepared data
load("out/data_prep_rdata.RData")

#number of synthetic ensembles
n_samp <- 10  #some of the diagnostic plots require at least 2 samples
#k for KNN sampling
kk <- 20
#the site to use as the key site for resampling
keysite <- which(site_names=="NHGC1") 
#the period over which to generate synthetic forecasts. good options include the entire obs record or the hindcast record
ixx_sim <- as.POSIXlt(seq(as.Date('1989-10-01'),as.Date('2019-09-15'),by='day'),tz = "UTC") #ixx_obs_forward   or   ixx_hefs
n_sim <- length(ixx_sim)

syn_hefs_forward <- array(NA,c(n_samp,n_sites,n_ens,n_sim,leads))

#parallel 'foreach' implementation
syn_hefs_out<-foreach(m = 1:n_samp,.combine='c',.inorder=F)%dopar% {
  print(m)
  syn_hefs<- syn_gen(n_samp,kk,keysite,ixx_sim)
  return(syn_hefs)
}

#compile to multidimensional array and save
for(m in 1:n_samp){
  syn_hefs_forward[m,,,,]<-syn_hefs_out[[m]]
}

saveRDS(syn_hefs_forward,file='out/syn_hefs_forward.rds')
saveRDS(ixx_sim,file='out/ixx_sim.rds')
saveRDS(n_samp,file='out/n_samp.rds')

