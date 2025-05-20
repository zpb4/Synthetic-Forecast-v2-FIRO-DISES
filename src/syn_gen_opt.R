
syn_gen <- function (seed,kk,keysite,knn_pwr,scale_pwr,dif,lo,sig_a,sig_b,cal_idx,val_idx,
                     obs_opt,obs_fwd_opt,hefs_fwd_opt,ixx_opt_wy) {

  set.seed(seed)
  n_gen <- length(val_idx)
  
  ################################set of the KNN###################################
  
  #knn parameters to use throughout
  tot <- sum(rep(1,kk) / 1:kk)
  wts <- rep(1,kk) / 1:kk / rep(tot,kk)

  #for knn, weigh earlier leads more than later leads
  my_power <- knn_pwr
  w <- 1:leads
  decay <- w^my_power / sum(w^my_power) 
  decay <- c(rep(decay[1],1),decay)
  
  val_idx_lds <- rep(-1:(-leads),each=length(val_idx)) + rep(val_idx,leads)
  #set the key observation-based covariate values for the simulation
  obs_forward_all_leads_gen <- obs_fwd_opt[val_idx_lds,]
  #obs_forward_all_leads_fit <- obs_forward_all_leads_hind[,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE]
  obs_forward_all_leads_fit <- obs_fwd_opt[cal_idx,]
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  ##gen_knn_data <- t(obs_forward_all_leads_gen)
  ##fit_knn_data <- t(obs_forward_all_leads_fit)
  
  gen_knn_data <- t(cbind(obs_opt[val_idx_lds],obs_forward_all_leads_gen))
  fit_knn_data <- t(cbind(obs_opt[cal_idx],obs_forward_all_leads_fit))
  
  ##gen_knn <- t(cbind(obs_opt[val_idx_lds],obs_forward_all_leads_gen))
  ##fit_knn <- t(cbind(obs_opt[cal_idx],obs_forward_all_leads_fit))
  
  ##gen_knn_data <- apply(gen_knn,2,function(x){x / sum(x)})
  ##fit_knn_data <- apply(fit_knn,2,function(x){x / sum(x)})
  
  ##gen_knn_data <- rbind((obs_opt[val_idx_lds]/max(obs_opt)),gen_knn_data)
  ##fit_knn_data <- rbind((obs_opt[cal_idx]/max(obs_opt)),fit_knn_data)
  
  ##gen_knn_data <- rbind((apply(gen_knn,2,max)/max(obs_opt)),gen_knn_data)
  ##fit_knn_data <- rbind((apply(fit_knn,2,max)/max(obs_opt)),fit_knn_data)
  
  ##gen_knn_data <- rbind(gen_knn_data,(apply(gen_knn,2,max)/max(obs_opt)))
  ##fit_knn_data <- rbind(fit_knn_data,(apply(fit_knn,2,max)/max(obs_opt)))
  

  #calculate distance
  knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    sqrt(colSums(decay*(gen_knn_data[,x] - fit_knn_data)^2))
  })
  #diag(knn_dist) <- 10*max(knn_dist)   #so that we don't sample "in sample" events

  #the resampled locations viaa KNN to use
  #note: vectors that have zero distance will not be resampled; prevents resampling of same event where fit and gen intersect
  knn_lst <- apply(knn_dist,2,function(x) {x[x==0]<-NA;sample(order(x)[1:kk],size=1,prob=wts)})
  rm(gen_knn_data,fit_knn_data)
  gc()
  #######################################################################################################

    
  #######################################################################################################
  
  #the fractions to sample from
  hefs_forward_resamp <- hefs_fwd_opt[,cal_idx,]
  
  #final array for the ensemble forecasts at all lead times
  all_leads_gen <- array(NA,c(1,n_ens,n_gen,leads))
  
  hi = lo + dif
  
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

  sigmoid_fun <- function(x,a,b){out <- 1/(1+exp(-x*a+b));return(out)}

  #what to add to adjust the resampled HEFS ensemble
  gen_scale <- obs_forward_all_leads_gen
  fit_scale <- obs_forward_all_leads_fit[knn_lst,]
  gen_scale[gen_scale==0] <- min(gen_scale[gen_scale>0])
  fit_scale[fit_scale==0] <- min(gen_scale[gen_scale>0])
  HEFS_scale_mat <-  gen_scale/fit_scale
  #HEFS_scale[which(is.na(HEFS_scale) | HEFS_scale==Inf | HEFS_scale > ratio_threshold)] <- 1
  for(k in 1:leads){
    mat_idx <- (1+((k-1)*n_gen)):(k*n_gen)
    HEFS_scale <- HEFS_scale_mat[mat_idx,]
    HEFS_scale[which(is.na(HEFS_scale[,k]) | HEFS_scale[,k]==Inf | HEFS_scale[,k]==0),k] <- 1
    HEFS_sc <- HEFS_scale[,k]
      
    obs_sc <- obs_forward_all_leads_gen[mat_idx,k]
    obs_sc[obs_sc<=0]<-min(obs_sc[obs_sc>0]) # set all 0 to min flow, log can't take 0s
    obs_sc<-log(obs_sc) #log of obs to make approx normal for scaling
    obs_scale <- scale(obs_sc)  #scale obs
    
    #ratio threshold is value between hi and lo threshold dependent on obs size via the sigmoid function
    ratio_threshold <- sigmoid_fun(obs_scale,sig_a,sig_b) * (dcy[k]-1) + 1
    
    #HEFS_sc[which(HEFS_sc > ratio_threshold)]<-1
    HEFS_sc[which(HEFS_sc > ratio_threshold)]<-ratio_threshold[which(HEFS_sc > ratio_threshold)]
    #HEFS_sc[which(HEFS_sc > ratio_threshold)]<-runif(length(which(HEFS_sc > ratio_threshold))) * (ratio_threshold[which(HEFS_sc > ratio_threshold)]-1) + 1
    
    #alternate approach to set ratios that exceed the threshold to the threshold value
    #HEFS_sc[which(obs_forward_all_leads_gen[j,,k]>pcntile_val)]<-HEFS_scale[which(obs_forward_all_leads_gen[j,,k]>pcntile_val),k]  
    HEFS_scale[,k]<-HEFS_sc
    
    for(e in 1:n_ens) {
      all_leads_gen[1,e,,k] <- (hefs_forward_resamp[e,knn_lst[mat_idx],k])*HEFS_scale[,k]
    }
  }
  #make adjustments for each ensemble member
  
    
  rm(hefs_forward_resamp,obs_forward_all_leads_fit,obs_forward_all_leads_gen,obs_fwd_opt,hefs_fwd_opt,obs_opt)
  gc()  
  
    
  return(all_leads_gen)  
    
}


#######################################################################################################

