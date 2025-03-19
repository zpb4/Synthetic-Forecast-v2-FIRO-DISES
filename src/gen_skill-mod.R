
args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

library(ncdf4)


#----------------------------------------
loc = 'YRS'
keysite_name = 'ORDC1'
pcnt_opt = 0.99
cal_val_setup = '5fold'

#Rouge et al. 2023, 'Forecast families..'
#S = 1 - k
#F_skill-mod = (1-k)O_t + kF_orig
skill_mods = c(0,0.2,0.4,0.6,0.8,1)
skill_mod = skill_mods[idx]
skill_dcy = 0.1
skill_tail = 0.5

load(paste('./out/',loc,'/data_prep_rdata.RData',sep=''))

#remove previous files if existing (netcdf will not overwrite files)
unlink(paste('./out/',loc,'/Qf-hefs_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_dcy,'.nc',sep=''),recursive=TRUE)
unlink(paste('out/',loc,'/Qf-syn_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_dcy,'_v3.nc',sep=''),recursive=TRUE)

syn_hefs_forward <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_v3.rds',sep=''))
print(dim(syn_hefs_forward))
#ixx_gen <- readRDS(paste('./out/',loc,'/ixx_gen.rds',sep='')) 
ixx_gen <- ixx_obs_forward
#n_samp <- readRDS(paste('./out/',loc,'/n_samp.rds',sep='')) 
n_samp <- dim(syn_hefs_forward)[1]

#add a single entry dimension to match synthetic forecasts
hefs_fwd<-array(NA,c(1,dim(hefs_forward)))
hefs_fwd[1,,,,]<-hefs_forward


#apply skill scaling
w <- leads:1
#decay <- w^skill_dcy / sum(w^skill_dcy)
#dcy_out <- (decay - min(decay))/(max(decay)-min(decay)) * (1-skill_tail) + skill_tail

dcy = (exp(skill_dcy*w)-exp(skill_dcy)) / (exp(skill_dcy*w[length(w)-1])-exp(skill_dcy))
dcy_out = dcy/max(dcy) * (skill_mod - skill_mod*skill_tail) + (skill_mod*skill_tail)

##plot(1:15,dcy_out,type='l',ylim=c(0,1),main=paste('Decay function, dcy=',skill_dcy,' tail=',skill_tail,sep=''),xlab='lead',ylab='eCRPSS scaling factor')

#S = skill_mod * dcy_out
S = dcy_out
k = 1 - S

#par(mfrow=c(1,2))
#plot(1:15,S,type='l',ylim=c(0,1),main=paste('eCRPSS, mod=',skill_mod,' dcy=',skill_dcy,' tail=',skill_tail,sep=''),xlab='lead',ylab='eCRPSS')
#plot(1:15,k,type='l',ylim=c(0,1),main=paste('k factor, mod=',skill_mod,' dcy=',skill_dcy,' tail=',skill_tail,sep=''),xlab='lead',ylab='k factor')

#recalculate HEFS with updated skill
hefs_mod <- array(NA,dim(hefs_fwd))
#HEFS w/modified skill
for(s in 1:dim(hefs_forward)[1]){
  for(e in 1:dim(hefs_forward)[2]){
    hefs_mod[1,s,e,,] <- matrix(rep((1-k),dim(hefs_forward)[3]),nrow=dim(hefs_forward)[3],byrow = T) * obs_forward_all_leads_hind[s,,] + 
      matrix(rep(k,dim(hefs_forward)[3]),nrow=dim(hefs_forward)[3],byrow = T) * hefs_forward[s,e,,]
  }
}

saveRDS(hefs_mod, paste('out/',loc,'/hefs_forward_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'.rds',sep=''))

#recalculate syn-HEFS with updated skill
shefs_mod <- array(NA,dim(syn_hefs_forward))
#HEFS w/modified skill
for(n in 1:n_samp){
  for(s in 1:dim(hefs_forward)[1]){
    for(e in 1:dim(hefs_forward)[2]){
      shefs_mod[n,s,e,,] <- matrix(rep((1-k),dim(syn_hefs_forward)[4]),nrow=dim(syn_hefs_forward)[4],byrow = T) * obs_forward_all_leads[s,,] + 
        matrix(rep(k,dim(syn_hefs_forward)[4]),nrow=dim(syn_hefs_forward)[4],byrow = T) * syn_hefs_forward[n,s,e,,]
    }
  }
}

saveRDS(shefs_mod, paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_v3.rds',sep=''))
saveRDS(shefs_mod[1:10,,,,], paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_plot-ens_v3.rds',sep=''))

#rearrange dimensions to match NHG model in Python
hefs_out<-aperm(hefs_mod,c(1,2,5,3,4))

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
hefs_nc<-nc_create(paste('./out/',loc,'/Qf-hefs_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'.nc',sep=''),hefs_var,force_v4 = F)
ncvar_put(hefs_nc,hefs_var,hefs_out)
nc_close(hefs_nc)

#-------------------------------------------------------------
#syn-hefs
shefs_out<-aperm(shefs_mod,c(1,2,5,3,4))
#define the dimensions
ens_dim<-ncdim_def('ensemble','',0:(dim(syn_hefs_forward)[1]-1))
site_dim<-ncdim_def('site','',0:(dim(syn_hefs_forward)[2]-1))
trace_dim<-ncdim_def('trace','',0:(dim(syn_hefs_forward)[3]-1))
date_dim<-ncdim_def('date',paste('days since',as.character(ixx_gen[1]),'00:00:00'),0:(dim(syn_hefs_forward)[4]-1),calendar = 'proleptic_gregorian')
ld_dim<-ncdim_def('lead','',0:(dim(syn_hefs_forward)[5]-1))

#write the variable to the netcdf file and save
shefs_var<-ncvar_def('syn','kcfs',dim=list(ens_dim,site_dim,ld_dim,trace_dim,date_dim))
shefs_nc<-nc_create(paste('out/',loc,'/Qf-syn_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_v3.nc',sep=''),shefs_var,force_v4 = F)
ncvar_put(shefs_nc,shefs_var,shefs_out)
nc_close(shefs_nc)

rm(list=ls());gc()

#################################END#############################################

