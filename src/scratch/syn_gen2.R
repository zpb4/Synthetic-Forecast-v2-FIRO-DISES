
#syn_gen <- function (n_samp=1,kk=20,keysite=1,ixx_sim) {
  
syn_gen <- function (seed,kk,keysite,knn_pwr,fit_start,fit_end,gen_start,gen_end,leave_out_years,
                       obs_forward_all_leads_gen,obs_forward_all_leads_hind,hefs_forward,ixx_hefs,
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
  my_power <- -1
  w <- 1:(leads+1)
  decay <- w^my_power / sum(w^my_power)

  #the threshold for the ratio of resampled and current obs flows beyond which no mean adjusts are made to HEFS resampled forecasts
  ratio_threshold <- 2
  
  #set the key observation-based covariate values for the simulation
  obs_forward_all_leads_gen <- obs_forward_all_leads[,ixx_obs_forward%in%ixx_gen,,drop=FALSE]
  n_gen <- length(ixx_gen)
  
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  sim_knn_data <- t(obs_forward_all_leads_gen[keysite,,])
  hind_knn_data <- t(obs_forward_all_leads_hind[keysite,ixx_hefs%in%ixx_fit,])
  #add on the previous observed day, since this might also help get the rising limb of a hydrograph
  sim_knn_data <- rbind(
    sim_knn_data[1,c(1,1:(ncol(sim_knn_data)-1)),drop=FALSE],
    sim_knn_data)
  hind_knn_data <- rbind(
    hind_knn_data[1,c(1,1:(ncol(hind_knn_data)-1)),drop=FALSE],
    hind_knn_data)
  #calculate distance
  knn_dist <- sapply(1:ncol(sim_knn_data),function(x){
    sqrt(colSums(decay*(sim_knn_data[,x] - hind_knn_data)^2))
  })
  diag(knn_dist) <- 10*max(knn_dist)   #so that we don't sample "in sample" events

  #the resampled locations viaa KNN to use
  knn_lst <- apply(knn_dist,2,function(x) {sample(c(0,order(x)[1:kk]),size=1,prob=c(0,wts))})
  #######################################################################################################

    
  #######################################################################################################
  
  #the fractions to sample from
  hefs_forward_resamp <- hefs_forward[,,ixx_hefs%in%ixx_obs_forward,,drop=FALSE]
  
  #final array for the ensemble forecasts at all lead times
  all_leads_sim <- array(NA,c(n_sites,n_ens,n_sim,leads))

  for (j in 1:n_sites) {
    HEFS_ens_mean_resamp_old <- apply(hefs_forward_resamp[j,,knn_lst,],FUN=mean,c(2,3))
    #what to add to adjust mean
    resamp_ratio <-  obs_forward_all_leads_sim[j,,]/obs_forward_all_leads_hind[j,knn_lst,] 
    resamp_ratio[which(is.na(resamp_ratio) | resamp_ratio==Inf | resamp_ratio > ratio_threshold)] <- 1
    HEFS_ens_mean_resamp_new <- HEFS_ens_mean_resamp_old * (1 +  (resamp_ratio - 1)) 
    #how to scale to adjust STDEV
    HEFS_scale <- HEFS_ens_mean_resamp_new/HEFS_ens_mean_resamp_old
    HEFS_scale[is.na(HEFS_scale)==TRUE | HEFS_scale==Inf] <- 1
    #make adjustments for each ensemble member
    for(e in 1:n_ens) {
      all_leads_sim[j,e,,] <- (hefs_forward_resamp[j,e,knn_lst,])*HEFS_scale
    }
    
  }
  gc()  
  
    
  return(all_leads_sim)  
    
}


#######################################################################################################

