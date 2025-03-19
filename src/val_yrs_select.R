#/////////////////////////////////////////
#Primary user defined settings

loc = 'NHG'              #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM' 'ADO'

keysite_name = 'NHGC1'   #specific site ID for 'keysite' which conditions the kNN sampling

#////////////////////////////////////////
library(ks)

load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))

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
samps <- 100
val_yr_vec <- array(NA,c(samps,val_size))
diff <- c()
set.seed(1)

for(i in 1:100){

  val_samp <- sample(unique(ixx_hefs_wy$year),val_size,replace = F)
  val_yr_vec[i,]<-val_samp

  val_vec <- ixx_hefs_wy$year%in%val_samp
  cal_vec <- !ixx_hefs_wy$year%in%val_samp

  obs_trn <- obs[ixx_obs%in%ixx_hefs,keysite_name]

  ann_tot_trn <- sapply(unique(ixx_hefs_wy$year),function(x){ob<-obs_trn[ixx_hefs_wy$year%in%x];return(sum(ob))})
  ann_tot_val <- sapply(val_samp,function(x){ob<-obs_trn[ixx_hefs_wy$year%in%x];return(sum(ob))})
  ann_tot_cal <- sapply(unique(ixx_hefs_wy$year[cal_vec]),function(x){ob<-obs_trn[ixx_hefs_wy$year%in%x];return(sum(ob))})

  cal_ks <- kde(ann_tot_cal,eval.points = seq(min(ann_tot_trn),max(ann_tot_trn),(max(ann_tot_trn)-min(ann_tot_trn))/100))
  val_ks <- kde(ann_tot_val,eval.points = seq(min(ann_tot_trn),max(ann_tot_trn),(max(ann_tot_trn)-min(ann_tot_trn))/100))

  diff[i] <- sqrt(sum((cal_ks$estimate - val_ks$estimate)^2))}

#plot comparison
#val_years <- readRDS(paste('./data/',loc,'/opt_val_years_samp=',val_size,'.rds',sep=''))-1900

val_years <- val_yr_vec[which.min(diff),]

val_vec <- ixx_hefs_wy$year%in%val_years
cal_vec <- !ixx_hefs_wy$year%in%val_years

ann_tot_val <- sapply(val_years,function(x){ob<-obs_trn[ixx_hefs_wy$year%in%x];return(sum(ob))})
ann_tot_cal <- sapply(unique(ixx_hefs_wy$year[cal_vec]),function(x){ob<-obs_trn[ixx_hefs_wy$year%in%x];return(sum(ob))})

hist(ann_tot_cal)
hist(ann_tot_val,col='yellow',add=T)

cal_ks <- kde(ann_tot_cal,eval.points = seq(min(ann_tot_trn),max(ann_tot_trn),(max(ann_tot_trn)-min(ann_tot_trn))/100))
val_ks <- kde(ann_tot_val,eval.points = seq(min(ann_tot_trn),max(ann_tot_trn),(max(ann_tot_trn)-min(ann_tot_trn))/100))

png(filename=paste('./data/',loc,'/val-years_kde-plot.png',sep=''))
plot(cal_ks$eval.points,cal_ks$estimate,type='l',lwd=3,xlab='Annual flow (kcfs)',ylab='Density')
lines(val_ks$eval.points,val_ks$estimate,col='yellow',lwd=2)
legend('topright',c('Cal','Val'),c('black','yellow'))
dev.off()

saveRDS(sort(val_years+1900),paste('./data/',loc,'/opt_val_years_samp=',val_size,'.rds',sep=''))
write.csv(sort(val_years+1900),paste('./data/',loc,'/opt_val_years_samp=',val_size,'.csv',sep=''))

#######################################END##################################################