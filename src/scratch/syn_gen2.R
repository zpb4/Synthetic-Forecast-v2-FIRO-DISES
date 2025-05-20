
##syn_gen <- function (seed,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,fit_start,fit_end,gen_start,gen_end,leave_out_years,
                     ##obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                     ##ixx_obs_forward) {

##syn_gen <- function (seed,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,t1,t2,t3,fit_start,fit_end,gen_start,gen_end,leave_out_years,
                       ##obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                       ##ixx_obs_forward) {
  
##syn_gen <- function (seed,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,t1,t2,t3,a1,a2,a3,b1,b2,b3,fit_start,fit_end,gen_start,gen_end,leave_out_years,
                       ##obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                       ##ixx_obs_forward) {
  
syn_gen <- function (seed,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,bias_hi,bias_pwr,var_bias_hi,var_bias_pwrs,l1_mn,l1_var,fit_start,fit_end,gen_start,gen_end,leave_out_years,
                       obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                       ixx_obs_forward) {
    
  set.seed(seed)
  
  ixx_gen <- as.POSIXlt(seq(as.Date(gen_start),as.Date(gen_end),by='day'),tz = "UTC") 
  ixx_obs_forward <- as.POSIXlt(seq(as.Date(ixx_obs_forward[1]),as.Date(ixx_obs_forward[length(ixx_obs_forward)]),by='day'),tz = "UTC")
  ixx_hefs <- as.POSIXlt(seq(as.Date(ixx_hefs[1]),as.Date(ixx_hefs[length(ixx_hefs)]),by='day'),tz = "UTC")
  ixx_obs <- as.POSIXlt(seq(as.Date(ixx_obs[1]),as.Date(ixx_obs[length(ixx_obs)]),by='day'),tz = "UTC")
  n_gen <- length(ixx_gen)
  ixx_fit_all <- as.POSIXlt(seq(as.Date(fit_start),as.Date(fit_end),by='day'),tz = "UTC")
  
  #function to convert an input date vector to water year reference vector
  wy_fun<-function(date_vec){
    wy_vec <- date_vec$year
    wy_vec[date_vec$mo%in%c(9,10,11)] <- wy_vec[date_vec$mo%in%c(9,10,11)]+1
    date_vec_wy <- date_vec
    date_vec_wy$year <- wy_vec
    return(date_vec_wy)
  }
  
  ixx_fit_all_wy <- wy_fun(ixx_fit_all)
  
  trn_idx <- !(ixx_fit_all_wy$year%in%(leave_out_years-1900))
  ixx_fit <- ixx_fit_all[trn_idx] #fit years excluding leave out years
  ixx_fit <- ixx_fit[ixx_fit%in%ixx_obs_forward] #ensure all fit years are in the 'ixx_obs_forward' index, remove if not
  ixx_hefs_obs_fwd <- ixx_hefs[ixx_hefs%in%ixx_obs_forward]
  
  #some error catching
  if (n_samp < 1 | !is.numeric(n_samp) | n_samp%%1!=0) {stop("n_samp is not a valid entry")}
  if (kk < 1 | !is.numeric(kk) | kk%%1!=0) {stop("k is not a valid entry")}
  if (ixx_gen[1] < ixx_obs_forward[1] | ixx_gen[length(ixx_gen)] > ixx_obs_forward[length(ixx_obs_forward)]) {
    stop("simulation period outside available observational period")
  }
  
  
  ################################set of the KNN###################################
  
  #knn parameters to use throughout
  tot <- sum(rep(1,kk) / 1:kk)
  wts <- rep(1,kk) / 1:kk / rep(tot,kk)

  #for knn, weigh earlier leads more than later leads
  w <- 1:leads
  ##w <- 1:leads
  decay <- w^knn_pwr / sum(w^knn_pwr)
  decay <- c(decay[1],decay)
  
  #set the key observation-based covariate values for the simulation
  obs_forward_all_leads_gen <- obs_forward_all_leads[,ixx_obs_forward%in%ixx_gen,,drop=FALSE]
  n_gen <- length(ixx_gen)
  obs_forward_all_leads_fit <- obs_forward_all_leads_hind[,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE]
  
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  ##gen_knn_data <- t(obs_forward_all_leads_gen[keysite,,])
  ##fit_knn_data <- t(obs_forward_all_leads_fit[keysite,,])
  gen_knn_data <- t(cbind(obs[ixx_obs%in%ixx_gen,keysite],obs_forward_all_leads_gen[keysite,,]))
  fit_knn_data <- t(cbind(obs[ixx_obs%in%ixx_fit,keysite],obs_forward_all_leads_fit[keysite,,]))

  #calculate distance
  knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    sqrt(colSums(decay*(gen_knn_data[,x] - fit_knn_data)^2))
  })

  #the resampled locations viaa KNN to use
  #note: vectors that have zero distance will not be resampled; prevents resampling of same event where fit and gen intersect
  knn_lst <- apply(knn_dist,2,function(x) {x[x==0]<-NA;sample(order(x)[1:kk],size=1,prob=wts)})

  hefs_resamp_vec <- ixx_fit[knn_lst]
  
  rm(gen_knn_data,fit_knn_data)
  gc()
  #######################################################################################################

    
  #######################################################################################################
  
  #the fractions to sample from
  hefs_forward_resamp <- hefs_forward[,,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE]
  hefs_mn_resamp <- apply(hefs_forward[,,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE],c(1,3,4),mean)
  
  #final array for the ensemble forecasts at all lead times
  all_leads_gen <- array(NA,c(n_sites,n_ens,n_gen,leads))
  
  #w <- 1:leads
  #decay <- (w^pwr / sum(w^pwr)) 
  #hi = hi
  #lo = lo
  #dcy <- (decay-min(decay)) * (hi-lo)/(max(decay)-min(decay)) + lo
  
  scale_decay_fun <- function(hi,lo,pwr,lds){
    w = 1:lds
    if(pwr!=0){
      win = rev(w)
      dcy = (exp(pwr*win)-exp(pwr)) / (exp(pwr*win[length(win)-1])-exp(pwr))
      dcy_out = dcy/max(dcy) * (hi-lo) + lo}
    if(pwr==0){
      dcy = (hi - lo)/(length(w)-1)
      dcy_out = hi - 0:(length(w)-1)*dcy}
    return(dcy_out)
  }
  
  dcy <- scale_decay_fun(hi,lo,scale_pwr,leads)
  dcy_bias <- scale_decay_fun(bias_hi,0,bias_pwr,leads)
  dcy_bias[1] <- l1_mn
  dcy_bias_var <- scale_decay_fun(var_bias_hi,0,var_bias_pwr,leads)
  dcy_bias_var[1] <- l1_var
  ##dcy_in <- scale_decay_fun(hi,lo,scale_pwr,(leads-3))
  ##dcy <- c(t1,t2,t3,dcy_in)

  sigmoid_fun <- function(x,a,b){out <- 1/(1+exp(-x*a+b));return(out)}
  
  ##sigas <- c(a1,a2,a3,rep(sig_a,(leads-3)))
  ##sigbs <- c(b1,b2,b3,rep(sig_b,(leads-3)))
  
  #plot(seq(-2,2,0.1),sigmoid_fun(seq(-2,2,0.1),1,0),type='l')
  
  HEFS_scale_out <- array(NA,c(n_sites,n_gen,leads))
  
  for (j in 1:n_sites) {
    #what to add to adjust the resampled HEFS ensemble
    HEFS_scale <-  obs_forward_all_leads_gen[j,,]/obs_forward_all_leads_fit[j,knn_lst,] 
    HEFS_mn_shift <- array(NA,dim(HEFS_scale))
    #HEFS_scale[which(is.na(HEFS_scale) | HEFS_scale==Inf | HEFS_scale > ratio_threshold)] <- 1
    for(k in 1:leads){
      HEFS_scale[which(is.na(HEFS_scale[,k]) | HEFS_scale[,k]==Inf),k] <- 1
      HEFS_sc_in <- HEFS_scale[,k]
      scale_up_vec <- rep(0,length(HEFS_sc_in))
      scale_up_vec[HEFS_sc_in>1] <- 1
      ##HEFS_sc <- HEFS_scale[,k]
      
      obs_sc <- obs_forward_all_leads_gen[j,,k]
      obs_sc[obs_sc<=0]<-min(obs_sc[obs_sc>0])
      obs_sc<-log(obs_sc)
      obs_scale <- scale(obs_sc)
      
      ##bias_val<-(obs_forward_all_leads_gen[j,,k] - min(obs_forward_all_leads_gen[j,,k])) / (max(obs_forward_all_leads_gen[j,,k]) - min(obs_forward_all_leads_gen[j,,k]))
      bias_val<-(obs_sc - min(obs_sc)) / (max(obs_sc) - min(obs_sc)) 
      ##bias_scale<- 1 - (dcy_bias[k] * bias_val)
      ##bias_scale<- 1 - (dcy_bias[k] * bias_val * scale_up_vec)
      ##bias_scale<- 1 - (rnorm(length(bias_Val),mean=dcy_bias[k],sd=dcy_bias[k]) * bias_val * scale_up_vec);bias_scale[bias_scale<=0]<-1
      HEFS_mn_shift[,k] <- -(dcy_bias[k] * bias_val * scale_up_vec) * hefs_mn_resamp[j,knn_lst,k]
      var_bias_scale<- 1 + (dcy_bias_var[k] * bias_val * scale_up_vec)
      
      ##HEFS_sc <- HEFS_scale[,k] * bias_scale
      HEFS_sc <- HEFS_scale[,k] * var_bias_scale
      
      ratio_threshold <- sigmoid_fun(obs_scale,sig_a,sig_b) * (dcy[k]-1) + 1
      ##ratio_threshold <- sigmoid_fun(obs_scale,sigas[k],sigbs[k]) * (dcy[k]-1) + 1
      #ratio_threshold <- dcy[k]
      
      #allow unlimited scaling for largest 10 events
      ##abv_thresh <- order(obs_forward_all_leads_fit[j,knn_lst,k],decreasing = T)[10]
      ##ext_scale <- obs_forward_all_leads_gen[j,,k] / obs_forward_all_leads_fit[j,knn_lst,k][abv_thresh]
      ##ext_scale[ext_scale < 1] <- 1
      ##ratio_threshold <- ratio_threshold * ext_scale
      #HEFS_sc[which(obs_sc > obs_sc[abv_thresh])] <- HEFS_scale[which(obs_sc > obs_sc[abv_thresh]),k]
      
      #accept/reject format for ratios exceeding threshold
      HEFS_sc[which(HEFS_sc > ratio_threshold)]<-1
      ##HEFS_sc[which(HEFS_sc > ratio_threshold)]<-ratio_threshold[which(HEFS_sc > ratio_threshold)]
      ##HEFS_sc[which(HEFS_sc > ratio_threshold)]<-runif(length(which(HEFS_sc > ratio_threshold))) * (ratio_threshold[which(HEFS_sc > ratio_threshold)]-1) + 1
      
      HEFS_scale[,k]<-HEFS_sc
    }
    OBS_scale <- obs[ixx_obs%in%ixx_gen,j] / obs[ixx_obs%in%ixx_fit,j][knn_lst]
    ##HEFS_scale_sm <- t(apply(HEFS_scale,1,function(x){out=ksmooth(1:leads,x,bandwidth=2,x.points=1:leads);return(out$y)}))
    HEFS_scale_sm <- t(apply(cbind(OBS_scale,HEFS_scale,HEFS_scale[,leads]),1,function(x){out=ksmooth(1:(leads+2),x,bandwidth=1,x.points=1:(leads+2));return(out$y)}))
    #make adjustments for each ensemble member
    for(e in 1:n_ens) {
      ##all_leads_gen[j,e,,] <- (hefs_forward_resamp[j,e,knn_lst,])*HEFS_scale
      ##all_leads_gen[j,e,,] <- (hefs_forward_resamp[j,e,knn_lst,])*HEFS_scale_sm
      all_leads_gen[j,e,,] <- (hefs_forward_resamp[j,e,knn_lst,]) * HEFS_scale_sm[,2:(leads+1)] + HEFS_mn_shift
    }
    ##HEFS_scale_out[j,,] <- HEFS_scale
    ##HEFS_scale_out[j,,] <- HEFS_scale_sm
    HEFS_scale_out[j,,] <- HEFS_scale_sm[,2:(leads+1)]
  }
  rm(hefs_forward_resamp,hefs_forward_resamp_sub,hefs_forward,obs_forward_all_leads_fit,
     obs_forward_all_leads_gen,obs_forward_all_leads_hind,obs_forward_all_leads)
  gc()  
  
    
  return(list(all_leads_gen,hefs_resamp_vec,HEFS_scale_out))  
    
}


#######################################################################################################

