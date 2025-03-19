#/////////////////////////////////////////
#Primary user defined settings

loc = 'LAM'              #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM' 'ADO'
keysite_name = 'LAMC1'
pcnt_opt = 0.99
cal_val_setup = 'cal'

wy = 90:119
wy_arr = array(NA,c(5,6))
set.seed(1)
for(i in 1:5){
  samp = sample(wy,6,replace=F)
  wy_arr[i,] = samp + 1900
  wy = wy[!(wy%in%samp)]
}


ixx_gen<-readRDS(file=paste('out/',loc,'/ixx_gen_hourly.rds',sep=''))
n_samp<-readRDS(file=paste('out/',loc,'/n_samp_hourly.rds',sep=''))

#function to convert an input date vector to water year reference vector
wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

ixx_gen <- wy_fun(ixx_gen)

syn_hefs_forward<-readRDS(file=paste('out/',loc,'/syn_hefs_forward_hourly_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'.rds',sep=''))

print(dim(syn_hefs_forward))

syn_hefs_forward_plot <- syn_hefs_forward[1:10,,,,1:72]

saveRDS(syn_hefs_forward_plot,file=paste('out/',loc,'/syn_hefs_forward_hourly_pcnt=',pcnt_opt,'_',keysite_name,'_',cal_val_setup,'_plot-ens.rds',sep=''))


rm(list=ls());gc()

#########################################END#####################################################