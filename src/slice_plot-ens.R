#script to slice out smaller subset of output array for plotting on local machine
loc = 'LAM'
keysite_name = 'LAMC1'
pcnt_opt = 0.99
cal_val_setup = '5fold'

syn_hefs_fwd <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'.rds',sep=''))
syn_hefs_forward <- syn_hefs_fwd[1:10,,,,,drop=FALSE]

saveRDS(syn_hefs_forward,paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_plot-ens.rds',sep=''))

rm(list=ls());gc()


#########################################END##########################