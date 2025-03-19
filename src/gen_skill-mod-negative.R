
args = commandArgs(trailingOnly=TRUE)
print(paste('task #',args[1]))
idx = as.numeric(args[1])

#idx=2

library(ncdf4)


#----------------------------------------
loc = 'YRS'
keysite_name = 'ORDC1'
pcnt_opt = 0.99
cal_val_setup = 'cal'

#Rouge et al. 2023, 'Forecast families..'
#S = 1 - k
#F_skill-mod = (1-k)O_t + kF_orig
skill_mods = c(-0.2,-0.4,-0.6,-0.8,-1)
skill_mod = skill_mods[idx]
skill_dcy = 0.1
skill_tail = 0.5

load(paste('./out/',loc,'/data_prep_rdata.RData',sep=''))

#remove previous files if existing (netcdf will not overwrite files)
unlink(paste('./out/',loc,'/Qf-hefs_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_dcy,'.nc',sep=''),recursive=TRUE)
unlink(paste('out/',loc,'/Qf-syn_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_dcy,'.nc',sep=''),recursive=TRUE)

syn_hefs_forward <- readRDS(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'.rds',sep=''))
##print(dim(syn_hefs_forward))
#ixx_gen <- readRDS(paste('./out/',loc,'/ixx_gen.rds',sep='')) 
ixx_gen <- ixx_obs_forward
#n_samp <- readRDS(paste('./out/',loc,'/n_samp.rds',sep='')) 
n_samp <- dim(syn_hefs_forward)[1]

#add a single entry dimension to match synthetic forecasts
hefs_fwd<-array(NA,c(1,dim(hefs_forward)))
hefs_fwd[1,,,,]<-hefs_forward

lobs_fnf <- read.csv('./data/YRS/ORDC1_fnf_daily.csv')
lobs_date <- lobs_fnf$Date..ending.12.UTC.
lobs_dtg = seq(as.Date(lobs_date[1]),as.Date(lobs_date[length(lobs_date)]),'day')
lobs_dtg_ix = as.POSIXlt(lobs_dtg)
lobs = c(0,lobs_fnf$Flow..KCFS.[-c(length(lobs_dtg))])
lobs_hefs = lobs[lobs_date%in%as.character(ixx_hefs)]

dowy_fun <- function(date_vec){
  ref_yr = seq(as.Date('2023-10-01'),as.Date('2024-09-30'),'day')
  ref_yr_idx = as.POSIXlt(ref_yr)
  dowy_idx = rep(NA,length(ref_yr_idx))
  #set any feb 29 (leap year) to same dowy as feb 28 for convenience
  dowy_idx[c(-which(ref_yr_idx$mo==1 & ref_yr_idx$mday==29))] <- 1:(length(ref_yr_idx)-1)
  dowy_idx_na <- dowy_idx
  dowy_idx[which(ref_yr_idx$mo==1 & ref_yr_idx$mday==29)] <- dowy_idx[which(ref_yr_idx$mo==1 & ref_yr_idx$mday==28)]
  match_vec = rbind(ref_yr_idx$mo,ref_yr_idx$mday,dowy_idx)
  match_vec_na = rbind(ref_yr_idx$mo,ref_yr_idx$mday,dowy_idx_na)
  dowy_out = rep(NA,length(date_vec))
  dowy_out_na = rep(NA,length(date_vec))
  for(i in 1:length(ref_yr_idx)){
    dowy_out[which(date_vec$mo==match_vec[1,i] & date_vec$mday==match_vec[2,i])] = match_vec[3,i]
    dowy_out_na[which(date_vec$mo==match_vec_na[1,i] & date_vec$mday==match_vec_na[2,i])] = match_vec_na[3,i]
  }
  return(list(dowy_out,dowy_out_na))
}

climo_idx = which(lobs_dtg=='1982-10-01'):which(lobs_dtg=='2023-09-30') #subtract a day for synch
climo_obs = lobs[climo_idx]
climo_dtg = dowy_fun(lobs_dtg_ix[climo_idx])[[2]]
climo_forecast_mat = sapply(1:365,function(x){out = climo_obs[which(climo_dtg==x)]})

climo_hefs_dowy = dowy_fun(ixx_hefs)[[1]]
climo_hefs_dowy = c(climo_hefs_dowy,1:leads)
climo_hefs = array(1,dim(hefs_forward))
for(i in 1:length(ixx_hefs)){
  climo_hefs[3,,i,] <- climo_forecast_mat[,climo_hefs_dowy[(i+1):(i+15)]]
}
climo_hefs[climo_hefs<0] <- 0

saveRDS(climo_hefs, paste('out/',loc,'/hefs_climo.rds',sep=''))

climo_shefs_dowy = dowy_fun(ixx_gen)[[1]]
climo_shefs_dowy = c(climo_shefs_dowy,1:leads)
climo_shefs = array(1,dim(syn_hefs_forward)[2:5])
for(i in 1:length(ixx_gen)){
  climo_shefs[3,,i,] <- climo_forecast_mat[,climo_shefs_dowy[(i+1):(i+15)]]
}
climo_shefs[climo_shefs<0] <- 0

saveRDS(climo_shefs, paste('out/',loc,'/synhefs_climo.rds',sep=''))

#which(ixx_hefs$year==97&ixx_hefs$mo==0&ixx_hefs$mday==1)

#hefs_vec = hefs_forward[3,,2283,]
#climo_vec = climo_hefs[3,,2283,]

#plot(0:15,lobs_hefs[2283:2298],type='l',ylim=c(0,350),lwd=2)
#for(i in 1:n_ens){
  #lines(0:15,c(lobs_hefs[2283],hefs_vec[i,]),col='gray')
#}

#plot(0:15,lobs_hefs[2283:2298],type='l',ylim=c(0,350),lwd=2)
#for(i in 1:n_ens){
  #lines(0:15,c(lobs_hefs[2283],climo_vec[i,]),col='gray')
#}

roll_sum <- function(x){
  out <- c()
  for(i in 1:length(x)){
    out[i] <- sum(x[1:i])
  }
  return(out)
}

skill_collapse <- function(tgt_forc,base_forc,mod){
  base_forc_skillmod<-rep(NA,length(base_forc))
  vec_base <- sort(base_forc,index.return=T)
  vec_tgt <- sort(tgt_forc,index.return=T)
  base_forc_skillmod[vec_base$ix] <- vec_base$x + (vec_tgt$x - vec_base$x) * mod
  return(base_forc_skillmod)
}

skill_collapse_mat <- function(tgt_forc,base_forc,mod_vec){
  sapply(1:dim(tgt_forc)[2],function(x){out = skill_collapse(tgt_forc[,x],base_forc[,x],mod_vec[x])})
}

skill_collapse_cumul <- function(tgt_forc,base_forc,mod_vec){
  lds = dim(base_forc)[2]
  cumul_base <- t(apply(base_forc,1,roll_sum))
  cumul_tgt <- t(apply(tgt_forc,1,roll_sum))
  vec_base_c <- sort(cumul_base[,lds],index.return=T)
  vec_tgt_c <- sort(cumul_tgt[,lds],index.return=T)
  cumul_out <- array(NA,dim(cumul_base))
  cumul_out[vec_base_c$ix,] <- cumul_base[vec_base_c$ix,] + (cumul_tgt[vec_tgt_c$ix,] - cumul_base[vec_base_c$ix,]) * matrix(rep(mod_vec,n_ens),ncol=lds,byrow=T)
  base_mod <- cbind(cumul_out[,1],t(apply(cumul_out,1,diff)))
  return(base_mod)
}

#v = skill_collapse_cumul(climo_vec,hefs_vec,rep(.5,leads))

#plot(0:15,lobs_hefs[2283:2298],type='l',ylim=c(0,350),lwd=2)
#for(i in 1:n_ens){
  #lines(0:15,c(lobs_hefs[2283],hefs_vec[i,]),col='gray')
#}

#plot(0:15,lobs_hefs[2283:2298],type='l',ylim=c(0,350),lwd=2)
#for(i in 1:n_ens){
  #lines(0:15,c(lobs_hefs[2283],v[i,]),col='gray')
#}

#plot(0:15,lobs_hefs[2283:2298],type='l',ylim=c(0,350),lwd=2)
#for(i in 1:n_ens){
  #lines(0:15,c(lobs_hefs[2283],climo_vec[i,]),col='gray')
#}

#apply skill scaling
w <- leads:1
dcy = (exp(skill_dcy*w)-exp(skill_dcy)) / (exp(skill_dcy*w[length(w)-1])-exp(skill_dcy))
dcy_out = dcy/max(dcy) * ((-skill_mod) - (-skill_mod)*skill_tail) + ((-skill_mod)*skill_tail)

#recalculate HEFS with updated skill
hefs_mod <- array(NA,dim(hefs_fwd))
#HEFS w/modified skill
for(s in 1:dim(hefs_forward)[1]){
  for(i in 1:dim(hefs_forward)[3]){
    hefs_mod[1,s,,i,] <- skill_collapse_cumul(climo_hefs[s,,i,],hefs_forward[s,,i,],dcy_out)
  }
}

saveRDS(hefs_mod, paste('out/',loc,'/hefs_forward_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'.rds',sep=''))

#recalculate syn-HEFS with updated skill
shefs_mod <- array(NA,dim(syn_hefs_forward))
#HEFS w/modified skill
for(n in 1:n_samp){
  for(s in 1:dim(hefs_forward)[1]){
    for(i in 1:dim(syn_hefs_forward)[4]){
      shefs_mod[n,s,,i,] <- skill_collapse_cumul(climo_shefs[s,,i,],syn_hefs_forward[n,s,,i,],dcy_out)
    }
  }
}

saveRDS(shefs_mod, paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'.rds',sep=''))
saveRDS(shefs_mod[1:10,,,,], paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'_plot-ens.rds',sep=''))

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
shefs_nc<-nc_create(paste('out/',loc,'/Qf-syn_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_mod=',skill_mod,'_dcy=',skill_dcy,'_tail=',skill_tail,'.nc',sep=''),shefs_var,force_v4 = F)
ncvar_put(shefs_nc,shefs_var,shefs_out)
nc_close(shefs_nc)

rm(list=ls());gc()

#################################END#############################################

