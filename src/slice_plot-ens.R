#script to slice out smaller subset of output array for plotting on local machine
loc = 'LAM'
parm = 'r'

syn_hefs_fwd <- readRDS(paste('out/',loc,'/syn_hefs_forward-',parm,'.rds',sep=''))
syn_hefs_forward <- syn_hefs_fwd[1:10,,,,,drop=FALSE]

saveRDS(syn_hefs_forward,paste('out/',loc,'/syn_hefs_forward-',parm,'_plot-ens.rds',sep=''))

rm(list=ls());gc()


#########################################END##########################