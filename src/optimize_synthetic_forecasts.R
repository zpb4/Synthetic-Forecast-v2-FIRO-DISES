args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

print(paste('opt begin',Sys.time()))

#library(optimParallel)
library(DEoptim)
library(twosamples)
library(lubridate)

#/////////////////////////////////////////
#Primary user defined settings

loc = 'SOD'             #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM'             #number of samples to generate
keysite_name = 'SRWC1'   #specific site ID for 'keysite' which conditions the kNN sampling
cal_val_setup = '5fold'

opt_strat = 'ecrps-dts'
opt_samps = 1
pcnt = 0.9901     #percentile to target for eCRPS loss function; e.g. 0.9 = 90-100% of obs inflows
obj_pwr = 0       #decay function value to weight leads differently, neg values to weight shorter leads, pos to weight longer leads

#optimized parameters (reasonable values shown for testing the function, actual optimization per)
kk = 30            #sampling window
knn_pwr = 0        #sampling decay
scale_pwr = 0      #decay scaling between hi and lo
dif = 3            #amount to add to lo to get hi end (lead 1) threshold scaling envelope, must be positive
lo = 2              #lo end (lead K) threshold scaling envelope
sig_a = 1           #slope of sigmoid obs dependence function
sig_b = 0           #location of sigmoid obs dependence function

source('./src/synthetic_forecast_opt-fun.R')
source('./src/syn_gen_opt.R')
source('../Synthetic-Forecast_Verification/src/forecast_verification_functions.R')
load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))
rm(hefs_forward_cumul,hefs_forward_frac,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid)
gc()

#index for selected keysite
keysite <- which(site_names==keysite_name)

#date/time manipulation
ixx_opt <- as.POSIXlt(seq(as.Date('1990-10-01'),as.Date('2019-09-15'),by='day'))
if(loc=='LAM'){
ixx_opt <- as.POSIXlt(seq(as.Date('1990-10-01'),as.Date('2019-08-15'),by='day'))}
ixx_opt_wy <- wy_fun(ixx_opt)

obs_opt <- obs[ixx_obs%in%ixx_opt,keysite]
obs_fwd_opt <- obs_forward_all_leads[keysite,ixx_obs_forward%in%ixx_opt,]
hefs_fwd_opt <- hefs_forward[keysite,,ixx_hefs%in%ixx_opt,]

rm(obs,obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward)
gc()

synopt_fun <- function(pars){
  out <- syn_opt(obs_opt,obs_fwd_opt,hefs_fwd_opt,opt_strat,opt_samps,idx,loc,keysite_name,cal_val_setup,pcnt,obj_pwr,pars[1],pars[2],pars[3],pars[4],pars[5],pars[6],pars[7]);return(out)}


#test the function to be optimized
system.time(synopt_fun(c(scale_pwr,dif,lo,sig_a,sig_b,kk,knn_pwr)))

#initial population of random uniform perturbed values around reasonable values specified above
#st = c(scale_pwr,dif,lo,sig_a,sig_b)
#st = c(scale_pwr,dif,lo,sig_a,sig_b,kk,knn_pwr)
#init_pop <- matrix(rep(st,10*length(st)),ncol=length(st),byrow=T)+matrix((runif(10*length(st)^2)-0.5),ncol=length(st))

#constraints for opt params
##upr = c(1,30,10,5,2,100,0)
##lwr = c(-1,0.1,1,0.1,-2,10,-2)

#constraints for fixed kk and knn_pwr
upr = c(1,30,10,5,2,30,0)
lwr = c(-1,0.1,1,0.1,-2,30,0)

#init_pop <- matrix(rep(st,10*length(st)),ncol=length(st),byrow=T)+matrix((runif(10*length(st)^2)-0.5),ncol=length(st))

if(tolower(.Platform$OS.type)=='windows'){
  cl <- makeCluster(detectCores(),type='PSOCK')
  clusterExport(cl=cl,varlist=list('synopt_fun','syn_opt','loc','keysite_name','cal_val_setup','kk','knn_pwr','ixx_obs','leads',
                                   'n_ens'))}

if(tolower(.Platform$OS.type)=='unix'){
  cl <- makeCluster(detectCores(),type='FORK')
  #cl <- makeCluster(detectCores(),type='PSOCK',setup_timeout=300)
  }

setDefaultCluster(cl=cl)
print(cl)

#out <- optimParallel(par=st,fn=synopt_fun,lower = lwr,upper = upr,control=list(maxit=10))
#out_DE <- DEoptim(fn=synopt_fun,lower=lwr,upper=upr,DEoptim.control(VTR=0.3,cluster=cl,itermax=10,parallelType='parallel',trace=TRUE))
#out_DE <- DEoptim(fn=synopt_fun,lower=lwr,upper=upr,DEoptim.control(cluster=cl,itermax=100,initialpop=init_pop,parallelType='parallel',trace=TRUE))
out_DE <- DEoptim(fn=synopt_fun,lower=lwr,upper=upr,DEoptim.control(cluster=cl,itermax=100,parallelType='parallel',trace=T))


ecrps_sse <- out_DE$optim$bestval
best_par <- out_DE$optim$bestmem
hi <- best_par[3]+best_par[2]
best_par[2] <- hi

opt_out <- c(ecrps_sse,best_par)
names(opt_out) <- c('eCRPS_sse','scale_pwr','hi','lo','sig_a','sig_b','kk','knn_pwr')

if(cal_val_setup!='5fold-test'){
  saveRDS(out_DE,paste('./out/',loc,'/DEopt_',cal_val_setup,'_pcnt=',pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'.rds',sep=''))
  saveRDS(opt_out,paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'.rds',sep=''))
  write.csv(t(opt_out),paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'.csv',sep=''),row.names = F)}

if(cal_val_setup=='5fold-test'){
  saveRDS(out_DE,paste('./out/',loc,'/DEopt_',cal_val_setup,'_pcnt=',pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'-',idx,'.rds',sep=''))
  saveRDS(opt_out,paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'-',idx,'.rds',sep=''))
  write.csv(t(opt_out),paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'-',idx,'.csv',sep=''),row.names = F)}

print(paste('opt end',Sys.time()))

#stopCluster()

rm(list = ls());gc()

###################################################END##################################################