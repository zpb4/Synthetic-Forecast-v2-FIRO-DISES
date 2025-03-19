
syn_gen <- function (seed,kk,keysite,knn_pwr,scale_pwr,hi,lo,sig_a,sig_b,fit_start,fit_end,gen_start,gen_end,leave_out_years,
                     obs_forward_all_leads,obs_forward_all_leads_hourly,obs_forward_all_leads_hind,obs_forward_all_leads_hind_hourly,hefs_forward,hefs_forward_hourly,ixx_hefs,
                     ixx_obs_forward,ixx_obs_forward_hourly,ixx_obs_hourly) {
  
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
  ixx_obs_forward_hourly <- as.POSIXlt(seq(as.Date(ixx_obs_forward_hourly[1]),as.Date(ixx_obs_forward_hourly[length(ixx_obs_forward_hourly)]),by='day'),tz = "UTC")
  ixx_hefs <- as.POSIXlt(seq(as.Date(ixx_hefs[1]),as.Date(ixx_hefs[length(ixx_hefs)]),by='day'),tz = "UTC")
  ixx_obs <- as.POSIXlt(seq(as.Date(ixx_obs[1]),as.Date(ixx_obs[length(ixx_obs)]),by='day'),tz = "UTC")
  ixx_obs_hourly <- as.POSIXlt(seq(as.Date(ixx_obs_hourly[1]),as.Date(ixx_obs_hourly[length(ixx_obs_hourly)]),by='day'),tz = "UTC")
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

  #adjust weights of sampling; neg values prioritize earlier leads, 0 is equal weighting
  lds_hourly <- leads*4
  my_power <- knn_pwr
  ##w <- 1:(leads+1)
  w <- 1:lds_hourly
  decay <- w^my_power / sum(w^my_power)
  #obs and 6 hr forecast have equal weights to capture rising or falling limbs
  decay <- c(decay[1],decay)
  
  #set the key observation-based covariate values for the simulation
  obs_forward_all_leads_gen <- obs_forward_all_leads[,ixx_obs_forward%in%ixx_gen,,drop=FALSE]
  obs_forward_all_leads_gen_hourly <- obs_forward_all_leads_hourly[,ixx_obs_forward_hourly%in%ixx_gen,,drop=FALSE]
  obs_forward_all_leads_gen_hourly_m1 <- obs_forward_all_leads_gen_hourly[,,24,drop=F]
  obs_forward_all_leads_gen_hourly_m1[,2:dim(obs_forward_all_leads_gen_hourly)[2],1]<-obs_forward_all_leads_gen_hourly[,1:(dim(obs_forward_all_leads_gen_hourly)[2]-1),24,drop=FALSE]

  n_gen <- length(ixx_gen)
  
  #obs_forward_all_leads_hind_hourly<-abind(obs_forward_all_leads_hind_hourly,obs_forward_all_leads_hind_hourly[,11195,,drop=F],along=2)
  
  obs_forward_all_leads_fit <- obs_forward_all_leads_hind[,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE]
  obs_forward_all_leads_fit_hourly <- obs_forward_all_leads_hind_hourly[,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE]
  obs_forward_all_leads_hind_hourly_m1 <- obs_forward_all_leads_hind_hourly[,,24,drop=F]
  obs_forward_all_leads_hind_hourly_m1[,2:dim(obs_forward_all_leads_hind_hourly)[2],1]<-obs_forward_all_leads_hind_hourly[,1:(dim(obs_forward_all_leads_hind_hourly)[2]-1),24,drop=FALSE]
  obs_forward_all_leads_fit_hourly_m1 <- obs_forward_all_leads_hind_hourly_m1[,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE]
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  #gen_knn_data <- t(obs_forward_all_leads_gen[keysite,,])
  #fit_knn_data <- t(obs_forward_all_leads_fit[keysite,,])

  hr6_index <- seq(6,leads*24,6)
  
  #gen_knn_data <- t(obs_forward_all_leads_gen_hourly[keysite,,hr6_index])
  #fit_knn_data <- t(obs_forward_all_leads_fit_hourly[keysite,,hr6_index])
  
  #bind the 6-hourly forecasts with observation (24 hr value from previous day) to include obs in the sampling routine
  gen_knn_data <- rbind(obs_forward_all_leads_gen_hourly_m1[keysite,,1],t(obs_forward_all_leads_gen_hourly[keysite,,hr6_index]))
  fit_knn_data <- rbind(obs_forward_all_leads_fit_hourly_m1[keysite,,1],t(obs_forward_all_leads_fit_hourly[keysite,,hr6_index]))
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
  hefs_forward_resamp_sub_hourly <- hefs_forward_hourly[,,ixx_hefs%in%ixx_obs_forward,,drop=FALSE]
  hefs_forward_resamp_hourly <- hefs_forward_resamp_sub_hourly[,,ixx_hefs_obs_fwd%in%ixx_fit,,drop=FALSE]
  rm(hefs_forward_resamp_sub_hourly)
  gc()
  
  #final array for the ensemble forecasts at all lead times
  all_leads_gen <- array(NA,c(n_sites,n_ens,n_gen,leads*24))
  
  my_power <- scale_pwr
  w <- 1:leads
  decay <- (w^my_power / sum(w^my_power)) 
  hi = hi
  lo = lo
  dcy <- (decay-min(decay)) * (hi-lo)/(max(decay)-min(decay)) + lo

  sigmoid_fun <- function(x,a,b){out <- 1/(1+exp(-x*a+b));return(out)}
  
  hr6_inc <- seq(6,24*leads,6)
  hr6_mat_idx <- matrix(1:(leads*4),ncol=4,byrow = T)
  
  for (j in 1:n_sites) {
    HEFS_scale_hr6 <-  obs_forward_all_leads_gen_hourly[j,,hr6_inc]/obs_forward_all_leads_fit_hourly[j,knn_lst,hr6_inc] 
    HEFS_scale_dly <-  obs_forward_all_leads_gen[j,,]/obs_forward_all_leads_fit[j,knn_lst,]
    HEFS_scale_hourly <- array(1,c(dim(HEFS_scale_dly)[1],(leads*24)))
    #HEFS_scale[which(is.na(HEFS_scale) | HEFS_scale==Inf | HEFS_scale > ratio_threshold)] <- 1
    #for(k in 1:(leads*24)){
    #for(k in 1:leads){
    for(k in 1:leads){
      HEFS_scale <- HEFS_scale_dly
      HEFS_scale[which(is.na(HEFS_scale[,k]) | HEFS_scale[,k]==Inf),k] <- 1
      HEFS_sc <- HEFS_scale[,k]
      
      #obs_sc <- obs_forward_all_leads_gen_hourly[j,,hr6_inc[k]]
      obs_sc <- obs_forward_all_leads_gen[j,,k]
      obs_sc[obs_sc<=0]<-min(obs_sc[obs_sc>0])
      obs_sc<-log(obs_sc)
      obs_scale <- scale(obs_sc)
      
      ratio_threshold <- sigmoid_fun(obs_scale,sig_a,sig_b) * (dcy[k]-1) + 1
      
      #accept/reject format for ratios exceeding threshold
      HEFS_sc[which(HEFS_sc > ratio_threshold)]<-1
      
      #allow unlimited scaling for largest 10 events
      #abv_thresh <- order(obs_forward_all_leads_fit_hourly[j,knn_lst,hr6_inc[k]],decreasing = T)[10]
      #ext_scale <- obs_forward_all_leads_gen_hourly[j,,hr6_inc[k]] / obs_forward_all_leads_fit_hourly[j,knn_lst,hr6_inc[k]][abv_thresh]
      ##abv_thresh <- order(obs_forward_all_leads_gen[j,,k],decreasing = T)[1:10]
      ##ext_scale <- HEFS_sc[abv_thresh]
      ##HEFS_sc[abv_thresh]<-ext_scale
      
      #daily scaling matrix at lead k set to updated scaling (based on thresholding)
      HEFS_scale_dly[,k]<-HEFS_sc
      #disaggregate to 6-hourly scaling
      HEFS_scale_hr6[,hr6_mat_idx[k,]] <- matrix(rep(HEFS_sc,4),ncol=4,byrow = F) 
    }

    #linearly interpolate scaling between hour 18 of forecast t-1 to hour 6 of forecast to smooth transition
    for(k in 1:(leads-1)){
      lin_interp <- ((HEFS_scale_hr6[,(k*4+2)] - HEFS_scale_hr6[,(k*4-1)]) / 3)
      hr6_interp <- matrix(rep(1:2),length(lin_interp),ncol=2,byrow = T) * lin_interp + HEFS_scale_hr6[,(k*4-2)]
      HEFS_scale_hr6[,(k*4):(k*4+1)] <- hr6_interp}
    
    #linear interpolation disaggregation of 6 hour values to hourly forecasts
    for(e in 1:n_ens) {
      #6-hourly forecast values
      sset <- cbind(c(obs_forward_all_leads_gen_hourly[j,1,1],obs_forward_all_leads_gen_hourly[j,-c(length(knn_lst)),24]),hefs_forward_resamp_hourly[j,e,knn_lst,hr6_inc]*HEFS_scale_hr6)
      #linear interpolation to hourly
      disagg_scale <- t(apply(sset,1,function(x){approx(c(0,hr6_inc),x,xout=(0:(leads*24)))$y}))
      #remove observed value (24 hr obs from previous day) when output to final array
      all_leads_gen[j,e,,1:(leads*24)] <- disagg_scale[,-c(1)]
      #if desiring to back-propagate 6 hourly value across leads 1-6 hr, apply this:
      #all_leads_gen[j,e,,1:6] <- matrix(rep(disagg_scale[,7],6),ncol=6,byrow=F)
    }
    
  }
  rm(hefs_forward_resamp_hourly,hefs_forward_hourly,obs_forward_all_leads_fit,
     obs_forward_all_leads_gen,obs_forward_all_leads_hind,obs_forward_all_leads,obs_forward_all_leads_fit_hourly,
     obs_forward_all_leads_gen_hourly,obs_forward_all_leads_hind_hourly,obs_forward_all_leads_hourly)
  gc()  
  
    
  return(list(all_leads_gen,hefs_resamp_vec))  
    
}


#######################################################################################################

