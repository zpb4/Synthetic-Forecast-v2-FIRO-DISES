#Script to process csv HEFS data from RFC

library(lubridate)
library(stringr)
library(extRemes)

print(paste('datapro start',Sys.time()))

#input main location ID
loc = 'LAM'  #current options: 'NHG' 'YRS' 'LAM'
load(paste('out/',loc,'/data_prep_rdata.RData',sep=''))

#scaling specifications
scale_site = 'LAMC1'
evt_to_scale = 1  #which event, by rank, to scale
scale_rng = 15    # +/- days to implement linear scaling framework to data
return_period = 100  #return period value or 'largest' to scale another event to largest on record

if(return_period=='largest' & evt_to_scale==1){
  stop('cannot scale largest event to itself')
}

#obs
#obs <- read.csv(paste('./data/',loc,'/observed_flows.csv',sep=''),header=TRUE)
#obs <- data.frame(obs)
#ixx_obs <- as.POSIXlt(ymd(obs[,1]),tz = "UTC")  
#site_names <- sort(colnames(obs)[-1])
#n_sites <- length(site_names)
#n_obs <- nrow(obs)
#col_names <- colnames(obs)
#col_names <- col_names[col_names!='Date']
#col_names <- sort(col_names) #rearrange alphabetically to jive w/'folder_list' below
#obs <- subset(obs,select=col_names)

#HEFS
#folder_list <- list.files(paste('./data/',loc,'/HEFS',sep=''))
#all_files <- list.files(paste('./data/',loc,'/HEFS/',folder_list[1],'/',sep=""))

#HEFS_dates <- sort(as.numeric(gsub("\\D", "", all_files)))
#start_date <- ymd(substr(HEFS_dates[1],1,8))
#end_date <- ymd(substr(HEFS_dates[length(HEFS_dates)],1,8))
#ixx_hefs <- as.POSIXlt(seq(as.Date(start_date),as.Date(end_date),'day'))
#ixx_hefs_out <- as.POSIXlt(seq(as.Date(start_date),as.Date(end_date+15),'day'))
#n_hefs <- length(ixx_hefs)

#hefs_forward <- readRDS(paste('./out/',loc,'/hefs_forward.rds',sep=''))

#leads <- min(dim(hefs_forward)[4],15) #only want max of 15d dynamical forecasts (some have 30 leads)
#n_ens <- dim(hefs_forward)[2]

#data process
obs_hefs <- obs[ixx_obs%in%ixx_hefs,]

wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

ixx_hefs_wy <- wy_fun(ixx_hefs)
ixx_obs_wy <- wy_fun(ixx_obs)

yrs <- unique(ixx_obs_wy$year)

ann_max <- array(NA,c(length(yrs),n_sites))

for(i in 1:length(yrs)){
  ann_max[i,] <- apply(obs[which(ixx_obs_wy$year==yrs[i]),],2,max)
}

fit<-fevd(ann_max[,which(col_names==scale_site)],type='GEV')

if(return_period != 'largest'){
  dsn<-qevd((1-1/return_period),loc=fit$results$par[1],scale=fit$results$par[2],shape=fit$results$par[3],type='GEV')}
if(return_period == 'largest'){
  lrg_idx <- order(obs_hefs[,which(col_names==scale_site)],decreasing = T)[1]
  dsn <- obs_hefs[lrg_idx,which(col_names==scale_site)]
}

#hist(ann_max[,which(col_names==scale_site)],xlim=c(0,1000),freq = F)
#lines(seq(1,1000,1),devd(seq(1,1000,1),loc=fit$results$par[1],scale=fit$results$par[2],shape=fit$results$par[3],type='GEV'))
#abline(v=dsn,col='red')

evt_idx <- order(obs_hefs[,which(col_names==scale_site)],decreasing = T)[evt_to_scale]
evt_idx_obs <- which(obs[,which(col_names==scale_site)]==obs_hefs[evt_idx,which(col_names==scale_site)])
evt_date <- ixx_hefs[evt_idx]
evt_to_scale_ratio <- dsn / obs_hefs[evt_idx,which(col_names==scale_site)]

scale_window <- c(seq(1,evt_to_scale_ratio,(evt_to_scale_ratio-1)/scale_rng)[-(scale_rng+1)],evt_to_scale_ratio,rev(seq(1,evt_to_scale_ratio,(evt_to_scale_ratio-1)/scale_rng)[-(scale_rng+1)]))


#rearrange fwd looking forecast to obs-synch
fwd_forecast_rearrange<-function(forecast){
  forecast_out<-array(0,dim(forecast))
  for(i in 1:dim(forecast)[3]){
    forecast_out[,(i+1):dim(forecast)[2],i]<-forecast[,1:(dim(forecast)[2]-i),i]
  }
  return(forecast_out)
}

rearrange_to_fwd_forecast<-function(forecast){
  forecast_out<-array(0,dim(forecast))
  for(i in 1:dim(forecast)[3]){
    forecast_out[,1:(dim(forecast)[2]-i),i]<-forecast[,(i+1):dim(forecast)[2],i]
  }
  return(forecast_out)
}

hefs<-array(NA,dim(hefs_forward))

for(i in 1:n_sites){
  hefs[i,,,]<-fwd_forecast_rearrange(hefs_forward[i,,,])
}

#scale obs
obs_scale <- obs
obs_scale[(evt_idx_obs-scale_rng):(evt_idx_obs+scale_rng),] <- obs[(evt_idx_obs-scale_rng):(evt_idx_obs+scale_rng),] * matrix(rep(scale_window,n_sites),ncol = n_sites,byrow = F)

#scale hefs
for(i in 1:n_sites){
  for(k in 1:dim(hefs_forward)[2]){
    hefs[i,k,(evt_idx-scale_rng):(evt_idx+scale_rng),] <- hefs[i,k,(evt_idx-scale_rng):(evt_idx+scale_rng),] * matrix(rep(scale_window,dim(hefs_forward)[4]),ncol = dim(hefs_forward)[4],byrow = F) 
  }
}

ld = 5

png(paste('./',loc,'-',scale_site,'_evt-rnk=',evt_to_scale,'_scale-event-comp.png',sep=''),height = 1500,width = 2500,res=300)

par(cex.axis=2,cex.lab=2.5,mar=c(4.5,5,0.5,0.5))
plot(-15:15,obs_scale[(evt_idx_obs-scale_rng):(evt_idx_obs+scale_rng),which(col_names==scale_site)],type='l',col='darkorange3',lwd=3,xlab='days from event',ylab='flow (kcfs)',ylim=c(0,1.05*max(obs_scale[(evt_idx_obs-ld):(evt_idx_obs+15-ld),which(col_names==scale_site)])))
     #main=c('scaled event comparison',paste(evt_date,'-',return_period,'yr event')),ylim=c(0,15))
lines(-15:15,obs[(evt_idx_obs-scale_rng):(evt_idx_obs+scale_rng),which(col_names==scale_site)],type='l',col='gray20',lwd=3)
legend('topright',c('unscaled','scaled'),lwd=c(3,3),col=c('gray20','darkorange3'),bty ='n',cex=2)
dev.off()

evt_idx_hist <- order(obs[,which(col_names==scale_site)],decreasing = T)[1]
evt_date_hist <- ixx_obs[evt_idx_hist]

png(paste('./',loc,'-',scale_site,'_evt-rnk=',evt_to_scale,'_hist-event-comp.png',sep=''),height = 1500,width = 2500,res=300)
par(cex.axis=2,cex.lab=2.5,mar=c(4.5,5,0.5,0.5))
plot(-15:15,obs[(evt_idx_hist-scale_rng):(evt_idx_hist+scale_rng),which(col_names==scale_site)],type='l',col='darkorange3',lwd=3,xlab='days from event',ylab='flow (kcfs)',ylim = c(0,1.05*max(obs_scale[(evt_idx_obs-ld):(evt_idx_obs+15-ld),which(col_names==scale_site)])))
     #main=c('hist event comparison',paste(evt_date,'-',evt_date_hist)),ylim=c(0,15))
lines(-15:15,obs[(evt_idx_obs-scale_rng):(evt_idx_obs+scale_rng),which(col_names==scale_site)],type='l',col='gray20',lwd=3)
legend('topright',c(as.character(evt_date),as.character(evt_date_hist)),lwd=c(3,3),col=c('gray20','darkorange3'),cex=2,bty='n')
dev.off()

hefs_fwd<-array(NA,dim(hefs_forward))

for(i in 1:n_sites){
  hefs_fwd[i,,,]<-rearrange_to_fwd_forecast(hefs[i,,,])
}

#correct output array for data lost at end of array due to back and forth swap
for(i in 1:dim(hefs_forward)[4]){
  hefs_fwd[,,(dim(hefs_forward)[3]-i):dim(hefs_forward)[3],]<-hefs_forward[,,(dim(hefs_forward)[3]-i):dim(hefs_forward)[3],]
}

##Plotting

png(paste('./',loc,'-',scale_site,'_evt-rnk=',evt_to_scale,'_unscale-hefs.png',sep=''),height = 1500,width = 2500,res=300)
par(cex.axis=2,cex.lab=2.5,mar=c(4.5,5,0.5,0.5))
plot(0:leads,obs[(evt_idx_obs-ld):(evt_idx_obs+leads-ld),which(col_names==scale_site)],type='l',col='gray20',xlab='days',ylab='flow (kcfs)',ylim=c(0,1.05*max(obs_scale[(evt_idx_obs-ld):(evt_idx_obs+15-ld),which(col_names==scale_site)])),lwd=3)
     #xlab='days',ylab='flow (kcfs)',main=c('unscaled',paste(evt_date,'-',ld,'d lead')))
for(i in 1:n_ens){
  lines(0:leads,c(obs[(evt_idx_obs-ld),which(col_names==scale_site)],hefs_forward[which(col_names==scale_site),i,evt_idx-ld,]),col='gray80')
}
lines(0:leads,obs[(evt_idx_obs-ld):(evt_idx_obs+leads-ld),which(col_names==scale_site)],col='gray20',lwd=3)
legend('topright',c('obs','HEFS'),lwd=c(3,1),col=c('gray20','gray80'),cex=2,bty='n')
abline(v=ld,lty=2,col='darkorange3')
dev.off()

png(paste('./',loc,'-',scale_site,'_evt-rnk=',evt_to_scale,'_scale-hefs.png',sep=''),height = 1500,width = 2500,res=300)
par(cex.axis=2,cex.lab=2.5,mar=c(4.5,5,0.5,0.5))
plot(0:leads,obs_scale[(evt_idx_obs-ld):(evt_idx_obs+leads-ld),which(col_names==scale_site)],type='l',col='gray20',ylim=c(0,1.05*max(obs_scale[(evt_idx_obs-ld):(evt_idx_obs+15-ld),which(col_names==scale_site)])),lwd=3,xlab='days',ylab='flow (kcfs)')
     #xlab='days',ylab='flow',main=c(paste('scaled -',return_period,'yr event'),paste(evt_date,'-',ld,'d lead')))
for(i in 1:n_ens){
  lines(0:leads,c(obs_scale[(evt_idx_obs-ld),which(col_names==scale_site)],hefs_fwd[which(col_names==scale_site),i,evt_idx-ld,]),col='gray80')
}
lines(0:leads,obs_scale[(evt_idx_obs-ld):(evt_idx_obs+leads-ld),which(col_names==scale_site)],lwd=3,col='gray20')
legend('topright',c('obs','HEFS'),lwd=c(3,1),col=c('gray20','gray80'),cex=2,bty='n')
abline(v=ld,lty=2,col='darkorange3')
dev.off()

png(paste('./',loc,'-',scale_site,'_evt-rnk=',evt_to_scale,'_hist-hefs.png',sep=''),height = 1500,width = 2500,res=300)
par(cex.axis=2,cex.lab=2.5,mar=c(4.5,5,0.5,0.5))
plot(0:leads,obs[(evt_idx_hist-ld):(evt_idx_hist+leads-ld),which(col_names==scale_site)],type='l',col='gray20',ylim=c(0,1.05*max(obs_scale[(evt_idx_obs-ld):(evt_idx_obs+15-ld),which(col_names==scale_site)])),lwd=3,xlab='days',ylab='flow')
     #xlab='days',ylab='flow',main=c('pre-hindcast extreme',paste(evt_date_hist,'-',ld,'d lead')))
legend('topright',c('obs','HEFS??'),lwd=c(3,1),col=c('gray20','gray80'),cex=2,bty='n')
abline(v=ld,lty=2,col='darkorange3')
dev.off()

hefs_forward <- hefs_fwd
obs <- obs_scale

csv_obs <- cbind(ixx_obs,obs)
colnames(csv_obs)[1]<-'Date'
write.csv(csv_obs,paste('data/',loc,'/observed_flows_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'.csv',sep=''),row.names = F)

rm(hefs,hefs_fwd,obs_scale,csv_obs)

#----------------------Cumulative stats for obs----------------------------------------

#####prepare observed data######
#calculate cumulative observed flow totals
#these are forward looking total, not including the current day
#therefore, we cannot have a value for the first time entry. 
#we only retain these values for dates where there are a full 15 days that follow
obs_forward_all_leads <- array(NA,c(n_sites,length(1:(n_obs-leads)),leads))

for (j in 1:n_sites) {
  for(i in 1:(n_obs-leads)) {
    obs_forward_all_leads[j,i,] <- obs[(i+1):(i+leads),j]
  }
}

#these are the dates that PRECEDE the 15 day cumulative totals, i.e., 
#on date t, here is the 15 day total over the NEXT 15 days
ixx_obs_forward <- ixx_obs[1:(n_obs-leads)]
n_obs_forward <- length(ixx_obs_forward)

#subset cumul obs from to hindcast period
obs_forward_all_leads_hind <- obs_forward_all_leads[,ixx_obs_forward%in%ixx_hefs,,drop=FALSE]


save.image(paste('out/',loc,'/data_prep_rdata_scaled_',scale_site,'_evt=',evt_to_scale,'_rtn=',return_period,'.RData',sep=''))

print(paste('datapro end',Sys.time()))

rm(list=ls());gc()

###############################################END################################################