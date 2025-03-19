#/////////////////////////////////////////
#Primary user defined settings

loc = 'YRS'              #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM' 'ADO'

keysite_name = 'NBBC1'   #specific site ID for 'keysite' which conditions the kNN sampling

#////////////////////////////////////////
load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))
library(ks)

ixx_hefs <- ixx_hefs[ixx_hefs%in%ixx_obs]

wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

ixx_hefs_wy <- wy_fun(ixx_hefs)
val_size <- 6


obs_trn <- obs[ixx_obs%in%ixx_hefs,keysite_name]

ann_tot_trn <- sapply(unique(ixx_hefs_wy$year),function(x){ob<-obs_trn[ixx_hefs_wy$year%in%x];return(sum(ob))})
wet_idx <- order(ann_tot_trn,decreasing = T)[1:val_size]
wet_years <- unique(ixx_hefs_wy$year)[wet_idx]

val_vec <- ixx_hefs_wy$year%in%wet_years
cal_vec <- !ixx_hefs_wy$year%in%wet_years

ann_tot_val <- sapply(wet_years,function(x){ob<-obs_trn[ixx_hefs_wy$year%in%x];return(sum(ob))})
ann_tot_cal <- sapply(unique(ixx_hefs_wy$year[cal_vec]),function(x){ob<-obs_trn[ixx_hefs_wy$year%in%x];return(sum(ob))})

hist(ann_tot_cal,xlim=c(0,max(ann_tot_val)))
hist(ann_tot_val,col='orange',add=T)

cal_ks <- kde(ann_tot_cal,eval.points = seq(min(ann_tot_trn),max(ann_tot_trn),(max(ann_tot_trn)-min(ann_tot_trn))/100))
val_ks <- kde(ann_tot_val,eval.points = seq(min(ann_tot_trn),max(ann_tot_trn),(max(ann_tot_trn)-min(ann_tot_trn))/100))

png(filename=paste('./data/',loc,'/wet-years_kde-plot.png',sep=''))
plot(cal_ks$eval.points,cal_ks$estimate,type='l',lwd=3,xlim=c(0,max(ann_tot_val)),xlab='Annual flow (kcfs)',ylab='Density')
lines(val_ks$eval.points,val_ks$estimate,col='orange',lwd=2)
legend('topright',c('Cal','Val'),col=c('black','orange'))
dev.off()

saveRDS(sort(wet_years+1900),paste('./data/',loc,'/wet_val_years_samp=',val_size,'.rds',sep=''))
write.csv(sort(wet_years+1900),paste('./data/',loc,'/wet_val_years_samp=',val_size,'.csv',sep=''))

#######################################END##################################################