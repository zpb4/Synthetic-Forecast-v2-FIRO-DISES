
syn_gen <- function (seed,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,fit_start,fit_end,gen_start,gen_end,
                     obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                     ixx_obs_forward) {
  
  set.seed(seed)
  
  ixx_gen <- as.POSIXlt(seq(as.Date(gen_start),as.Date(gen_end),by='day'),tz = "UTC") 
  ixx_obs_forward <- as.POSIXlt(seq(as.Date(ixx_obs_forward[1]),as.Date(ixx_obs_forward[length(ixx_obs_forward)]),by='day'),tz = "UTC")
  ixx_hefs <- as.POSIXlt(seq(as.Date(ixx_hefs[1]),as.Date(ixx_hefs[length(ixx_hefs)]),by='day'),tz = "UTC")
  ixx_obs <- as.POSIXlt(seq(as.Date(ixx_obs[1]),as.Date(ixx_obs[length(ixx_obs)]),by='day'),tz = "UTC")
  
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
  my_power <- knn_pwr
  ##w <- 1:(leads+1)
  w <- 1:leads
  decay <- w^my_power / sum(w^my_power)
  
  #set the key observation-based covariate values for the simulation
  obs_forward_all_leads_gen <- obs_forward_all_leads[,,,drop=FALSE]
  n_gen <- length(ixx_gen)
  obs_forward_all_leads_fit <- obs_forward_all_leads_hind[,,,drop=FALSE]
  
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  gen_knn_data <- t(obs_forward_all_leads_gen[keysite,,])
  fit_knn_data <- t(obs_forward_all_leads_fit[keysite,,])

  #calculate distance
  knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    sqrt(colSums(decay*(gen_knn_data[,x] - fit_knn_data)^2))
  })

  #the resampled locations viaa KNN to use
  #note: vectors that have zero distance will not be resampled; prevents resampling of same event where fit and gen intersect
  knn_lst <- apply(knn_dist,2,function(x) {x[x==0]<-NA;sample(order(x)[1:kk],size=1,prob=wts)})

  rm(gen_knn_data,fit_knn_data)
  gc()
  #######################################################################################################

  #final array for the ensemble forecasts at all lead times
  all_leads_gen <- array(NA,c(n_sites,n_ens,n_gen,leads))
  
  my_power <- scale_pwr
  w <- 1:leads
  decay <- (w^my_power / sum(w^my_power)) 
  hi = hi
  lo = lo
  dcy <- (decay-min(decay)) * (hi-lo)/(max(decay)-min(decay)) + lo

  sigmoid_fun <- function(x,a,b){out <- 1/(1+exp(-x*a+b));return(out)}
  
  #plot(seq(-2,2,0.1),sigmoid_fun(seq(-2,2,0.1),1,0),type='l')
  
  for (j in 1:n_sites) {
    #calculate scaling value
    HEFS_scale <-  obs_forward_all_leads_gen[j,,]/obs_forward_all_leads_fit[j,knn_lst,] 
    #HEFS_scale[which(is.na(HEFS_scale) | HEFS_scale==Inf | HEFS_scale > ratio_threshold)] <- 1
    for(k in 1:leads){
      HEFS_scale[which(is.na(HEFS_scale[,k]) | HEFS_scale[,k]==Inf),k] <- 1
      HEFS_sc <- HEFS_scale[,k]
      
      obs_sc <- obs_forward_all_leads_gen[j,,k]
      obs_sc[obs_sc<=0]<-min(obs_sc[obs_sc>0])
      obs_sc<-log(obs_sc)
      obs_scale <- scale(obs_sc)
      
      ratio_threshold <- sigmoid_fun(obs_scale,sig_a,sig_b) * (dcy[k]-1) + 1
      
      #allow unlimited scaling for largest 10 events
      ##abv_thresh <- order(obs_forward_all_leads_fit[j,knn_lst,k],decreasing = T)[10]
      ##ext_scale <- obs_forward_all_leads_gen[j,,k] / obs_forward_all_leads_fit[j,knn_lst,k][abv_thresh]
      ##ext_scale[ext_scale < 1] <- 1
      ##ratio_threshold <- ratio_threshold * ext_scale
      #HEFS_sc[which(obs_sc > obs_sc[abv_thresh])] <- HEFS_scale[which(obs_sc > obs_sc[abv_thresh]),k]
      
      #accept/reject format for ratios exceeding threshold
      HEFS_sc[which(HEFS_sc > ratio_threshold)]<-1
      
      HEFS_scale[,k]<-HEFS_sc
    }
    #make adjustments for each ensemble member
    for(e in 1:n_ens) {
      all_leads_gen[j,e,,] <- (hefs_forward[j,e,knn_lst,])*HEFS_scale
    }
    
  }
  rm(hefs_forward,obs_forward_all_leads_fit,
     obs_forward_all_leads_gen,obs_forward_all_leads_hind,obs_forward_all_leads)
  gc()  
  
    
  return(all_leads_gen)  
    
}


#######################################################################################################

