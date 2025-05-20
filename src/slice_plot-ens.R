#script to slice out smaller subset of output array for plotting on local machine
loc = 'YRS'
keysite_name = 'ORDC1'
pcnt_opt = 0.99
cal_val_setup = '5fold-test'
obj_pwr = 0
opt_strat = 'ecrps-dts'

syn_hefs_fwd <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'.rds',sep=''))
syn_hefs_forward <- syn_hefs_fwd[1:10,,,,,drop=FALSE]

saveRDS(syn_hefs_forward,paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'_plot-ens.rds',sep=''))

rm(list=ls());gc()


#########################################END##########################