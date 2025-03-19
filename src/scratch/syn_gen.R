
syn_gen <- function (seed,kk,keysite,knn_pwr,fit_start,fit_end,gen_start,gen_end,leave_out_years,use_mpi,
obs_forward_all_leads_gen,obs_forward_all_leads_hind,hefs_forward_cumul,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,hefs_forward_frac,ixx_hefs,
ixx_obs_forward,varx_fun,correct_leads,lds_to_correct,fix_order,use_ar,refrac) {
  
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
  
  ################################Step 0: set of the KNN###################################
  
  #knn parameters to use throughout
  tot <- sum(rep(1,kk) / 1:kk)
  wts <- rep(1,kk) / 1:kk / rep(tot,kk)
  
  #for knn, weigh earlier leads more than later leads
  my_power <- knn_pwr    #larger neg values = faster decay, more weight on early leads
  decay <- (1:(leads+1))^my_power / sum(((1:(leads+1))^my_power)) #k dim vector
  
  #set the key observation-based covariate values for the simulation
  obs_forward_all_leads_gen <- obs_forward_all_leads[,ixx_obs_forward%in%ixx_gen,,drop=FALSE]
  n_gen <- length(ixx_gen)
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  gen_knn_data <- t(obs_forward_all_leads_gen[keysite,,])
  fit_knn_data <- t(obs_forward_all_leads_hind[keysite,ixx_hefs_obs_fwd%in%ixx_fit,]) #k x n matrix

  gen_knn_data <- rbind(
    gen_knn_data[1,c(1,1:(ncol(gen_knn_data)-1)),drop=FALSE],
    gen_knn_data)
  fit_knn_data <- rbind(
    fit_knn_data[1,c(1,1:(ncol(fit_knn_data)-1)),drop=FALSE],
    fit_knn_data)
  
  #apply function over all indices in the generation obs
  knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    sq_diff <- (gen_knn_data[,x] - fit_knn_data)^2
    decay_diff <- sq_diff *decay
    #take square root of difference summed over leads to get Euclidean distance modified by decay function
    sqrt(colSums(decay_diff))
  })
  #knn_dist is a n_fit x n_gen matrix of distances
  
  #the resampled locations via KNN to use
  knn_lst <- apply(knn_dist,2,function(x){
    #order the distances along the fit dimension yielding a value for each of the n_gen dimensions
    ord_dist <- order(x)[1:kk] #order function sorts in ascending order and returns indices for the top 1:kk distances
    #output value is a distance weighted sampling of the 1:kk ordered indices
    sample(c(0,ord_dist),size=1,prob=c(0,wts))})
  
  #knn_lst is a n_gen length vector of indices to sample from the fitted period in following routines
  
  rm(gen_knn_data,fit_knn_data)
  gc()
  #####################Step 1: fit and simulate the ensembe average model##########################
  
  #fit var model set to only use obs flows from keysite, and assumes a lag-1 structure
  varx_x <- obs_forward_all_leads_hind[keysite,ixx_hefs_obs_fwd%in%ixx_fit,]
  #var_model <- MTS::VARX(zt=t(hefs_forward_cumul_ens_avg[,ixx_hefs%in%ixx_fit,drop=FALSE]), p=1, xt = varx_x)
  var_model <- varx_fun(zt=t(hefs_forward_cumul_ens_avg[,ixx_hefs%in%ixx_fit,drop=FALSE]), p=1, xt = varx_x,include.mean=FALSE,output=FALSE)
  #zt = a T-by-k data matrix of a k-dim time series
  #p VAR order
  #xt = a T-by-kx matrix of kx exogenous variables
  
  #simulate from var model
  varx_x <- obs_forward_all_leads_gen[keysite,,]
  varx_resid <- rbind(rep(0,n_sites),var_model$residuals) #assume zero residuals for first time step
  exg_pred <- varx_x%*%t(var_model$beta)
  innovations <- varx_resid[knn_lst,,drop=FALSE]
  
  #mv nomr innov model
  #mv_fit = MGMM::FitGMM(var_model$residuals)
  #innovations <- mvtnorm::rmvnorm(length(knn_lst),mean=mv_fit@Mean,sigma=mv_fit@Covariance)
  
  #simulations of ensemble average of cumulative flows
  cumul_ens_avg_gen <- array(0,c(n_gen,n_sites))
  for (i in 2:n_gen) {
    cumul_ens_avg_gen[i,] <- exg_pred[i,] + var_model$Phi%*%cumul_ens_avg_gen[i-1,] + innovations[i,]
    cumul_ens_avg_gen[i,]  <- sapply(cumul_ens_avg_gen[i,],function(x){max(x,0)})
  }
  rm(varx_x,var_model,varx_resid,exg_pred,innovations)
  gc()

  ######Step 2: Model the cumulative ensemble spread around the cumulative ensemble mean###################
  
  cumul_ens_resid_gen <- array(0,c(n_sites,n_ens,n_gen))
  cumul_gen <- array(NA,c(n_sites,n_ens,n_gen))
  
  #Individual AR version
  ar_ens_resid_models <- list(n_sites)
  for (j in 1:n_sites) {
    ar_ens_resid_models[[j]] <- list(n_ens)
    for(e in 1:n_ens) {
      #resample residuals of ensemble members around ensemble mean
      hefs_forward_cumul_ens_resid_fit <- hefs_forward_cumul_ens_resid[,,ixx_hefs%in%ixx_fit,drop=FALSE]
      if(use_ar==FALSE){
        cumul_ens_resid_gen[j,e,] <- hefs_forward_cumul_ens_resid_fit[j,e,knn_lst]}
      #fit AR model to ensemble residual, then simulate from that model with resampled innovations
      if(use_ar==TRUE){
        ar_ens_resid_models[[j]][[e]] <- arima(hefs_forward_cumul_ens_resid_fit[j,e,],order=c(1,0,0),include.mean=FALSE)
        innovations <- ar_ens_resid_models[[j]][[e]]$residuals[knn_lst]
        cumul_ens_resid_gen[j,e,] <- arima.sim(n_gen,model=list('ar'=ar_ens_resid_models[[j]][[e]]$coef),innov=innovations)}
      #cumul_ens_resid_gen[j,e,] <- arima.sim(n_gen,model=list('ar'=ar_ens_resid_models[[j]][[e]]$coef))
      #use these resampled residuals to estimate individual cumulative ensemble members
      cumul_gen[j,e,] <- cumul_ens_resid_gen[j,e,] + cumul_ens_avg_gen[,j]
      #address any negative values
      cumul_gen[j,e,] <- sapply(cumul_gen[j,e,],function(x) {max(0,x)})
      
    }
  }
  if(use_ar==T){rm(ar_ens_resid_models,innovations,cumul_ens_resid_gen,hefs_forward_cumul_ens_resid_fit,hefs_forward_cumul_ens_resid)}
  if(use_ar==F){rm(cumul_ens_resid_gen,hefs_forward_cumul_ens_resid_fit,hefs_forward_cumul_ens_resid)}
  gc()
  
  ######Step 3: Model the fractionation of each ensemble across lead times###################
  #kk=10
  tot <- sum(rep(1,kk) / 1:kk)
  wts <- rep(1,kk) / 1:kk / rep(tot,kk)
  
  #decay vector for correction by leads
  #define with decay function
  #my_power2 <- -1    #larger neg values = faster decay, more weight on early leads
  #decay2 <- (1:(leads))^my_power2 / sum(((1:leads)^my_power2)) * c(rep #k dim vector
  #manually tune
  decay2 <- c(0.7,0.55,0.45,0.43,0.42,0.2,0.2,0.1,0,0,0,0,0,0,0)
  
  estBetaParams <- function(mu, var) {
    alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
    beta <- alpha * (1 / mu - 1)
    return(params = list(alpha = alpha, beta = beta))
  }
  
  #the resampled locations via KNN to use
  set.seed(seed+1)
  knn_lst <- apply(knn_dist,2,function(x) {sample(c(0,order(x)[1:kk]),size=1,prob=c(0,wts))})

  #final array for the ensemble forecasts at all lead times
  all_leads_gen_frac <- array(NA,c(dim(cumul_gen),leads))
  all_leads_gen <- array(NA,c(dim(cumul_gen),leads))
  hefs_fwd_cumul_fit <- hefs_forward_cumul[,,ixx_hefs%in%ixx_fit,drop=FALSE]
  hefs_fwd_cumul <- hefs_fwd_cumul_fit[,,knn_lst,drop=FALSE]
  
  for (j in 1:n_sites) {
    #resorts randomly generated cumulative ensemble to match original ensemble sorting
    if(fix_order==TRUE){
      for(k in 1:length(ixx_gen)){
        rnk<-rank(hefs_fwd_cumul[j,,k])
        cumul_gen[j,,k]<-sort(cumul_gen[j,,k])[rnk]
      }
    }
    for(e in 1:n_ens) {
      #knn to resample fractions
      hefs_fit_forward_frac <- hefs_forward_frac[,,ixx_hefs%in%ixx_fit,,drop=FALSE]
      all_leads_gen_frac[j,e,,] <- hefs_fit_forward_frac[j,e,knn_lst,]
      #set undefined fractions (due to zero totals) to all zeros
      all_leads_gen_frac[j,e,,][which(is.na(all_leads_gen_frac[j,e,,]))] <- 0
      #calculate final 15-day ahead forecasts
      all_leads_gen[j,e,,] <- sapply(1:leads,function(x) {all_leads_gen_frac[j,e,,x]*cumul_gen[j,e,]})
    }
    #arrays for mean and stdev shift implementation
    hefs_fwd_fit <- hefs_forward[j,,ixx_hefs%in%ixx_fit,]
    hefs_fwd_cumul_fit <- hefs_forward_cumul[j,,ixx_hefs%in%ixx_fit]
    obs_fwd_fit <- obs_forward_all_leads_hind[j,ixx_hefs_obs_fwd%in%ixx_fit,]
    cumul_gen_step <- cumul_gen[j,,]
    
    #code to correct leads by mean/sdev shift
    if(correct_leads==TRUE){
      for(k in lds_to_correct){
        #calculate syn-hefs mean scale param
        hefs_mn <- apply(hefs_fwd_fit[,knn_lst,k],2,mean)
        pred_mn <- apply(all_leads_gen[j,,,k],2,mean)
        mn_scale <- hefs_mn / obs_fwd_fit[knn_lst,k] * obs_forward_all_leads_gen[j,,k] / pred_mn  
        mn_scale[is.na(mn_scale)==T|mn_scale==0|mn_scale==Inf|mn_scale==-Inf]<-1
        #calculate the amount to add to each member of the ensemble for the mean shift and apply it
        #beta distribution model for variability in skill correction; not yet implemented or tested
        ##decay2 <- rbeta(length(pred_mn),shape1 = estBetaParams(decay2[k],.01)$alpha,shape2 = estBetaParams(decay2[k],.01)$beta)
        
        #additive value to apply across all ensemble members
        mn_add <- (pred_mn * mn_scale - pred_mn) * decay2[k]
        all_leads_gen[j,,,k] <- all_leads_gen[j,,,k] + matrix(rep(mn_add,each=n_ens),ncol=length(ixx_gen),byrow=F)
        
        #calculate standard dev scaling
        sd_scale <- (apply(hefs_fwd_fit[,knn_lst,k],2,function(x){sd(x,na.rm=TRUE)})/hefs_mn) * (pred_mn+mn_add) / apply(all_leads_gen[j,,,k],2,function(x){sd(x,na.rm=TRUE)}) 
        sd_scale[is.na(sd_scale)==T|sd_scale==0|sd_scale==Inf|sd_scale==-Inf]<-1
        #order of magnitude envelope for stability
        env = 10
        sd_scale[abs(sd_scale)>env]<-1
        #scale ensemble spread by the sd_scale param
        #Note: this line subtracts the ensemble mean from the ensemble at lead k, multiplies those values by the sd_scale, and then adds them back to the ensemble mean 
        #to result in the scaled ensemble
        all_leads_gen[j,,,k] <- (all_leads_gen[j,,,k] - matrix(rep((pred_mn + mn_add),each=n_ens),ncol=length(ixx_gen),byrow=F)) * matrix(rep(sd_scale,each=n_ens),ncol=length(ixx_gen),byrow=F) + matrix(rep((pred_mn + mn_add),each=n_ens),ncol=length(ixx_gen),byrow=F)
        
        #iteratively refractionate ensemble to maintain generated cumulative ensemble
        #no refrac necessary for k = leads
        if(refrac==T & k<leads){
          #refractionate remaining cumulative ensemble
          #subtract cumulative ensemble up to lead k from total cumul ensemble 
          cumul_gen_step[,] = cumul_gen[j,,] - apply(all_leads_gen[j,,,1:k],c(1,2),sum)
          cumul_gen_step[cumul_gen_step<0]<-0
          #remaining cumul array to multiply against fractionation array
          cumul_gen_mat <- array(NA,c(dim(cumul_gen_step),(leads-k)))
          for(i in 1:(leads-k)){
            cumul_gen_mat[,,i] <- cumul_gen_step
          }
          #recalculate fractionation on remaining leads and multiply by remaining cumul ensemble
          if(k<(leads-1)){
            hefs_frac_step <- aperm(apply(hefs_fwd_fit[,knn_lst,(k+1):leads],c(1,2),function(x){x/sum(x)}),c(2,3,1))}
          #apply function drops a dimension with size 1, so need to modify at (leads-1)
          if(k==(leads-1)){
            hefs_fr <- apply(hefs_fwd_fit[,knn_lst,(k+1):leads],c(1,2),function(x){x/sum(x)})
            hefs_frac_step <- array(NA,c(dim(hefs_fr),1))
            hefs_frac_step[,,1] <- hefs_fr}
          hefs_frac_step[is.na(hefs_frac_step)|hefs_frac_step==Inf|hefs_frac_step==-Inf|hefs_frac_step<0]<-0
          all_leads_gen[j,,,(k+1):leads] <- hefs_frac_step[,,] * cumul_gen_mat[,,]}
      }
    }
  }
  rm(all_leads_gen_frac,cumul_gen,hefs_fit_forward_frac,hefs_forward_frac,knn_lst,hefs_fwd_fit,hefs_fwd_cumul_fit,obs_fwd_fit,cumul_gen_step)
  gc()

  all_leads_gen[is.na(all_leads_gen)|all_leads_gen==Inf|all_leads_gen==-Inf|all_leads_gen<0]<-0
    
  return(all_leads_gen)  
  
  rm(list = ls());gc()
}

