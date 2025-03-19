#script to slice out smaller subset of output array for plotting on local machine
loc = 'YRS'
keysite_name = 'ORDC1'
pcnt_opt = 0.99
cal_val_setup = '5fold-test'

#scaling specifications
scale_site = 'ORDC1'
evt_to_scale = 1  #which event, by rank, to scale
scale_rng = 15    # +/- days to implement linear scaling framework to data
return_period = 500

syn_hefs_fwd <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'.rds',sep=''))
syn_hefs_forward <- syn_hefs_fwd[1:10,,,,,drop=FALSE]

saveRDS(syn_hefs_forward,paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'_plot-ens.rds',sep=''))

rm(list=ls());gc()


#########################################END##########################