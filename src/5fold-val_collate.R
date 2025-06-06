#/////////////////////////////////////////
#Primary user defined settings

loc = 'YRS'              #main hindcast location ID, current options: 'NHG' 'YRS' 'LAM' 'ADO'
keysite_name = 'ORDC1'
pcnt_opt = 0.99
cal_val_setup = '5fold-test'
obj_pwr = 0
opt_strat = 'ecrps-dts'

wy = 90:119
wy_arr = array(NA,c(5,6))
set.seed(1)
for(i in 1:5){
  samp = sample(wy,6,replace=F)
  wy_arr[i,] = samp + 1900
   wy = wy[!(wy%in%samp)]
}


ixx_gen<-readRDS(file=paste('out/',loc,'/ixx_gen.rds',sep=''))
n_samp<-readRDS(file=paste('out/',loc,'/n_samp.rds',sep=''))

#function to convert an input date vector to water year reference vector
wy_fun<-function(date_vec){
  wy_vec <- date_vec$year
  wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
  date_vec_wy <- date_vec
  date_vec_wy$year <- wy_vec
  return(date_vec_wy)
}

ixx_gen <- wy_fun(ixx_gen)

syn_hefs_forward<-readRDS(file=paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'-',1,'.rds',sep=''))

print(dim(syn_hefs_forward))

for(i in 2:5){
  leave_out_years = wy_arr[i,]
  val_idx <- ixx_gen$year%in%(leave_out_years-1900)
  syn_hefs_fwd<-readRDS(file=paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'-',i,'.rds',sep=''))
  syn_hefs_forward[,,,val_idx,]<-syn_hefs_fwd[,,,val_idx,]
}

print(dim(syn_hefs_forward))

saveRDS(syn_hefs_forward,file=paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'.rds',sep=''))

#remove/delete 5 separate split files
unlink(paste('out/',loc,'/syn_hefs_forward_pcnt=',pcnt_opt,'_objpwr=',obj_pwr,'_optstrat=',opt_strat,'_',keysite_name,'_',cal_val_setup,'-',1:5,'.rds',sep=''),recursive=TRUE)

#rm(list=ls());gc()

#########################################END#####################################################