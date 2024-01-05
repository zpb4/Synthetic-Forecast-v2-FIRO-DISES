#script to slice out smaller subset of output array for plotting on local machine

syn_hefs_fwd <- readRDS('out/syn_hefs_forward.rds')
syn_hefs_forward <- syn_hefs_fwd[1:10,,,,]

saveRDS(syn_hefs_forward,'out/syn_hefs_forward_plot-ens.rds')

rm(list=ls());gc()


#########################################END##########################