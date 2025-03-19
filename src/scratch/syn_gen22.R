
syn_gen <- function (seed,kk,keysite,knn_pwr,ratio_thresh,fit_start,fit_end,gen_start,gen_end,leave_out_years,
                     obs_forward_all_leads,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
                     ixx_obs_forward) {
  
  #n_samp = the number of synthetic ensembles to create
  #kk = the k used in KNN resampling
  #keysite = an index for which site should be used to drive resampling  
  #ixx_gen = date object that spans the length of the desired synthetic forecasts
  #ixx_fit = date object from hindcast period to fit model, excluding and 'leave out years' specified in model setup
  
  #dimensional refs:
  #k = # of leads
  #ne = # of ensembles
  #nsite = # of sites
  #n = time dimension
  #library(MTS)
  #library(future)
  #library(parallel)
  #library(doParallel)
  #library(doMPI)
  #load('./out/data_prep_rdata.RData')
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
  my_power <- knn_pwr
  ##w <- 1:(leads+1)
  w <- 1:leads
  decay <- w^my_power / sum(w^my_power)

  #the threshold for the ratio of resampled and current obs flows beyond which no mean adjusts are made to HEFS resampled forecasts
  ratio_threshold <- ratio_thresh
  
  #set the key observation-based covariate values for the simulation
  obs_forward_all_leads_gen <- obs_forward_all_leads[,ixx_obs_forward%in%ixx_gen,,drop=FALSE]
  n_gen <- length(ixx_gen)
  obs_forward_all_leads_fit <- obs_forward_all_leads_hind[,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE]
  
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  gen_knn_data <- t(obs_forward_all_leads_gen[keysite,,])
  fit_knn_data <- t(obs_forward_all_leads_fit[keysite,,])
  #add on the previous observed day, since this might also help get the rising limb of a hydrograph
  ##gen_knn_data <- rbind(
    ##gen_knn_data[1,c(1,1:(ncol(gen_knn_data)-1)),drop=FALSE],
    ##gen_knn_data)
  ##fit_knn_data <- rbind(
    ##fit_knn_data[1,c(1,1:(ncol(fit_knn_data)-1)),drop=FALSE],
    ##fit_knn_data)
  #calculate distance
  knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    sqrt(colSums(decay*(gen_knn_data[,x] - fit_knn_data)^2))
  })
  #diag(knn_dist) <- 10*max(knn_dist)   #so that we don't sample "in sample" events

  #the resampled locations viaa KNN to use
  #note: vectors that have zero distance will not be resampled; prevents resampling of same event where fit and gen intersect
  knn_lst <- apply(knn_dist,2,function(x) {x[x==0]<-NA;sample(order(x)[1:kk],size=1,prob=wts)})
  #set.seed(seed+1)
  #knn_lst2 <- apply(knn_dist,2,function(x) {x[x==0]<-NA;sample(order(x)[1:kk],size=1,prob=wts)})
  rm(gen_knn_data,fit_knn_data)
  gc()
  #######################################################################################################

    
  #######################################################################################################
  
  #the fractions to sample from
  hefs_forward_resamp_sub <- hefs_forward[,,ixx_hefs%in%ixx_obs_forward,,drop=FALSE]
  hefs_forward_resamp <- hefs_forward[,,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE]
  
  #final array for the ensemble forecasts at all lead times
  all_leads_gen <- array(NA,c(n_sites,n_ens,n_gen,leads))
  
  my_power <- -0.5
  w <- 1:leads
  decay <- (w^my_power / sum(w^my_power)) + 1

  for (j in 1:n_sites) {
    HEFS_ens_mean_resamp_old <- apply(hefs_forward_resamp[j,,knn_lst,],FUN=mean,c(2,3))
    #what to add to adjust the resampled HEFS ensemble
    HEFS_scale <-  obs_forward_all_leads_gen[j,,]/obs_forward_all_leads_fit[j,knn_lst,] 
    #HEFS_scale[which(is.na(HEFS_scale) | HEFS_scale==Inf | HEFS_scale > ratio_threshold)] <- 1
    for(k in 1:leads){
      #HEFS_scale[which(is.na(HEFS_scale[,k]) | HEFS_scale[,k]==Inf | HEFS_scale[,k] > ratio_threshold[k]),k] <- 1
      HEFS_scale[which(is.na(HEFS_scale[,k]) | HEFS_scale[,k]==Inf),k] <- 1
      HEFS_sc <- HEFS_scale[,k]
      #ratio_threshold <- 1+scale(obs_forward_all_leads_gen[j,,k],center=F)/max(scale(obs_forward_all_leads_gen[j,,k],center=F))
      pcntile_idx <- round(0.99999*length(obs_forward_all_leads_gen[j,,k]))
      sel_idx <- order(obs_forward_all_leads_gen[j,,k])[pcntile_idx]
      pcntile_val <- obs_forward_all_leads_gen[j,sel_idx,k]
      #ratio_threshold <- 1+(obs_forward_all_leads_gen[j,,k]/(pcntile_val))
      ratio_threshold <- (1+(obs_forward_all_leads_gen[j,,k]/(pcntile_val)))*decay
      HEFS_sc[which(HEFS_sc > ratio_threshold)]<-ratio_threshold[which(HEFS_sc > ratio_threshold)]
      HEFS_scale[,k]<-HEFS_sc
      #HEFS_scale[which(HEFS_scale[,k] > ratio_threshold[k]),k]<-1
      #rand_ratio <- rnorm(length(HEFS_scale[,k]),mean=ratio_threshold[k],sd=(0.1*ratio_threshold[k]))
      #rand_ratio[rand_ratio<=0]<-1
      #HEFS_scale[which(HEFS_scale[,k]>rand_ratio),k]<-rand_ratio[HEFS_scale[,k]>rand_ratio]
      #HEFS_scale[which(HEFS_scale[,k]>rand_ratio),k]<-1
    }
    #make adjustments for each ensemble member
    for(e in 1:n_ens) {
      all_leads_gen[j,e,,] <- (hefs_forward_resamp[j,e,knn_lst,])*HEFS_scale
    }
    
  }
  rm(hefs_forward_resamp,hefs_forward_resamp_sub,hefs_forward,obs_forward_all_leads_fit,
     obs_forward_all_leads_gen,obs_forward_all_leads_hind,obs_forward_all_leads)
  gc()  
  
    
  return(all_leads_gen)  
    
}


#######################################################################################################

