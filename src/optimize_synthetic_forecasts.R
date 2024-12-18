args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

print(paste('opt begin',Sys.time()))

#library(optimParallel)
library(DEoptim)

#/////////////////////////////////////////
#Primary user defined settings

loc = 'YRS'             #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM'             #number of samples to generate
keysite_name = 'NBBC1'   #specific site ID for 'keysite' which conditions the kNN sampling
cal_val_setup = '5fold-test'
kk <- 30            #sampling window
knn_pwr = 0        #sampling decay
pcnt = 0.99         #percentile to target for eCRPS loss function; e.g. 0.9 = 90-100% of obs inflows
obj_pwr = 0       #decay function value to weight leads differently, neg values to weight shorter leads, pos to weight longer leads

#optimized parameters (reasonable values shown for testing the function, actual optimization per)
scale_pwr = -1      #decay scaling between hi and lo
diff = 3            #amount to add to lo to get hi end (lead 1) threshold scaling envelope, must be positive
lo = 2              #lo end (lead K) threshold scaling envelope
sig_a = 1           #slope of sigmoid obs dependence function
sig_b = 0           #location of sigmoid obs dependence function

source('./src/synthetic_forecast_opt-fun.R')
load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))


synopt_fun <- function(pars){
  out <- syn_opt(idx=idx,loc=loc,keysite_name=keysite_name,cal_val_setup=cal_val_setup,kk=kk,knn_pwr=knn_pwr,pcnt=pcnt,obj_pwr=obj_pwr,pars[1],pars[2],pars[3],pars[4],pars[5]);
  return(out)}

#test the function to be optimized
system.time(synopt_fun(c(scale_pwr,diff,lo,sig_a,sig_b)))

#initial population of random uniform perturbed values around reasonable values specified above
st = c(scale_pwr,diff,lo,sig_a,sig_b)
init_pop <- matrix(rep(st,10*length(st)),ncol=length(st),byrow=T)+matrix((runif(10*length(st)^2)-0.5),ncol=length(st))

#constraints for opt params
upr = c(-0.1,30,10,5,2)
lwr = c(-5,0.1,1,0.1,-2)

init_pop <- matrix(rep(st,10*length(st)),ncol=length(st),byrow=T)+matrix((runif(10*length(st)^2)-0.5),ncol=length(st))

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
out_DE <- DEoptim(fn=synopt_fun,lower=lwr,upper=upr,DEoptim.control(cluster=cl,itermax=100,initialpop=init_pop,parallelType='parallel',trace=TRUE))


ecrps_sse <- out_DE$optim$bestval
best_par <- out_DE$optim$bestmem
hi <- best_par[3]+best_par[2]
best_par[2] <- hi

opt_out <- c(ecrps_sse,best_par)
names(opt_out) <- c('eCRPS_sse','scale_pwr','hi','lo','sig_a','sig_b')

if(cal_val_setup!='5fold-test'){
  saveRDS(out_DE,paste('./out/',loc,'/DEopt_',cal_val_setup,'_pcnt=',pcnt,'_',keysite_name,'.rds',sep=''))
  saveRDS(opt_out,paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt,'_',keysite_name,'.rds',sep=''))
  write.csv(t(opt_out),paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt,'_',keysite_name,'.csv',sep=''),row.names = F)}

if(cal_val_setup=='5fold-test'){
  saveRDS(out_DE,paste('./out/',loc,'/DEopt_',cal_val_setup,'_pcnt=',pcnt,'_',keysite_name,'-',idx,'.rds',sep=''))
  saveRDS(opt_out,paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt,'_',keysite_name,'-',idx,'.rds',sep=''))
  write.csv(t(opt_out),paste('./out/',loc,'/DE-opt-params_',cal_val_setup,'_pcnt=',pcnt,'_',keysite_name,'-',idx,'.csv',sep=''),row.names = F)}

print(paste('opt end',Sys.time()))

#stopCluster()

rm(list = ls());gc()

###################################################END##################################################