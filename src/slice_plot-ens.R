#script to slice out smaller subset of output array for plotting on local machine
parm = 'c'

syn_hefs_fwd <- readRDS(paste('out/syn_hefs_forward-',parm,'.rds',sep=''))
syn_hefs_forward <- syn_hefs_fwd[1:10,,,,]

saveRDS(syn_hefs_forward,paste('out/syn_hefs_forward-',parm,'_plot-ens.rds',sep=''))

rm(list=ls());gc()


#########################################END##########################