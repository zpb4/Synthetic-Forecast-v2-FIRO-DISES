
rm(list=ls())
library(lubridate)
library(ncdf4)


#----------------------------------------
loc = 'YRS'
keysite_name = 'ORDC1'
pcnt_opt = 0.99
cal_val_setup = '5fold-test'

#scaling specifications
scale_site = 'ORDC1'
evt_to_scale = 1  #which event, by rank, to scale
scale_rng = 15    # +/- days to implement linear scaling framework to data
return_period = 500

load(paste('out/',loc,'/data_prep_rdata_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'.RData',sep=''))

#remove previous files if existing (netcdf will not overwrite files)
unlink(paste('./out/',loc,'/Qf-hefs_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'.nc',sep=''),recursive=TRUE)
unlink(paste('out/',loc,'/Qf-syn_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'.nc',sep=''),recursive=TRUE)

syn_hefs_forward <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'.rds',sep=''))
ixx_gen <- readRDS(paste('./out/',loc,'/ixx_gen.rds',sep='')) 
n_samp <- readRDS(paste('./out/',loc,'/n_samp.rds',sep='')) 
#syn_hefs_forward <- syn_hefs_forward[1:10,,,,]

#add a single entry dimension to match synthetic forecasts
hefs_fwd<-array(NA,c(1,dim(hefs_forward)))
hefs_fwd[1,,,,]<-hefs_forward

#rearrange dimensions to match NHG model in Python
hefs_out<-aperm(hefs_fwd,c(1,2,5,3,4))

#---------------------------------------------------
#hefs
#define the dimensions
ens_dim<-ncdim_def('ensemble','',0:(dim(hefs_fwd)[1]-1))
site_dim<-ncdim_def('site','',0:(dim(hefs_fwd)[2]-1))
trace_dim<-ncdim_def('trace','',0:(dim(hefs_fwd)[3]-1))
date_dim<-ncdim_def('date',paste('days since',as.character(ixx_hefs[1]),'00:00:00'),0:(dim(hefs_fwd)[4]-1),calendar = 'proleptic_gregorian')
ld_dim<-ncdim_def('lead','',0:(dim(hefs_fwd)[5]-1))

#write the variable to the netcdf file and save
hefs_var<-ncvar_def('hefs','kcfs',dim=list(ens_dim,site_dim,ld_dim,trace_dim,date_dim))
hefs_nc<-nc_create(paste('./out/',loc,'/Qf-hefs_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'.nc',sep=''),hefs_var,force_v4 = F)
ncvar_put(hefs_nc,hefs_var,hefs_out)
nc_close(hefs_nc)

#-------------------------------------------------------------
#syn-hefs
shefs_out<-aperm(syn_hefs_forward,c(1,2,5,3,4))
#define the dimensions
ens_dim<-ncdim_def('ensemble','',0:(dim(syn_hefs_forward)[1]-1))
site_dim<-ncdim_def('site','',0:(dim(syn_hefs_forward)[2]-1))
trace_dim<-ncdim_def('trace','',0:(dim(syn_hefs_forward)[3]-1))
date_dim<-ncdim_def('date',paste('days since',as.character(ixx_gen[1]),'00:00:00'),0:(dim(syn_hefs_forward)[4]-1),calendar = 'proleptic_gregorian')
ld_dim<-ncdim_def('lead','',0:(dim(syn_hefs_forward)[5]-1))

#write the variable to the netcdf file and save
shefs_var<-ncvar_def('syn','kcfs',dim=list(ens_dim,site_dim,ld_dim,trace_dim,date_dim))
shefs_nc<-nc_create(paste('out/',loc,'/Qf-syn_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'.nc',sep=''),shefs_var,force_v4 = F)
ncvar_put(shefs_nc,shefs_var,shefs_out)
nc_close(shefs_nc)

rm(list=ls());gc()

#################################END#############################################

