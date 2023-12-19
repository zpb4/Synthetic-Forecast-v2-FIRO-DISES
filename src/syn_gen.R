
syn_gen <- function (n_samp=1,kk=20,keysite=1,fit_start,fit_end,gen_start,gen_end,leave_out_years) {
  
  #n_samp = the number of synthetic ensembles to create
  #kk = the k used in KNN resampling
  #keysite = an index for which site should be used to drive resampling  
  #ixx_gen = date object that spans the length of the desired synthetic forecasts
  #ixx_fit = date object from hindcast period to fit model, excluding and 'leave out years' specified in model setup
  
  ixx_gen <- as.POSIXlt(seq(as.Date(gen_start),as.Date(gen_end),by='day'),tz = "UTC") 
  n_gen <- length(ixx_gen)
  ixx_fit_all <- as.POSIXlt(seq(as.Date(fit_start),as.Date(fit_end),by='day'),tz = "UTC")
  trn_idx <- !(ixx_fit_all$year%in%(leave_out_years-1900))
  ixx_fit <- ixx_fit_all[trn_idx] #fit years excluding leave out years
  ixx_fit <- ixx_fit[ixx_fit%in%ixx_obs_forward] #ensure all fit years are in the 'ixx_obs_forward' index, remove if not
  
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
  my_power <- -.5
  decay <- (1:leads)^my_power / sum(((1:leads)^my_power))
  
  #set the key observation-based covariate values for the simulation
  obs_forward_all_leads_gen <- obs_forward_all_leads[,ixx_obs_forward%in%ixx_gen,,drop=FALSE]
  n_gen <- length(ixx_gen)
  
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  gen_knn_data <- t(obs_forward_all_leads_gen[keysite,,])
  fit_knn_data <- t(obs_forward_all_leads_hind[keysite,ixx_hefs%in%ixx_fit,])
  knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    sqrt(colSums(decay*(gen_knn_data[,x] - fit_knn_data)^2))
  })
  
  #the resampled locations viaa KNN to use
  knn_lst <- apply(knn_dist,2,function(x) {sample(c(0,order(x)[1:kk]),size=1,prob=c(0,wts))})
  #######################################################################################################

    
  
  #####################Step 1: fit and simulate the ensembe average model##########################
  
  #fit var model set to only use obs flows from keysite, and assumes a lag-1 structure
  varx_x <- obs_forward_all_leads_hind[keysite,ixx_hefs%in%ixx_fit,]
  var_model <- MTS::VARX(zt=t(hefs_forward_cumul_ens_avg[,ixx_hefs%in%ixx_fit,drop=FALSE]), p=1, xt = varx_x)
  #zt = a T-by-k data matrix of a k-dim time series
  #p VAR order
  #xt = a T-by-kx matrix of kx exogenous variables
  
  #simulate from var model
  varx_x <- obs_forward_all_leads_gen[keysite,,]
  varx_resid <- rbind(rep(0,n_sites),var_model$residuals) #assume zero residuals for first time step
  exg_pred <- varx_x%*%t(var_model$beta)
  innovations <- varx_resid[knn_lst,,drop=FALSE]
  
  #simulations of ensemble average of cumulative flows
  cumul_ens_avg_gen <- array(0,c(n_gen,n_sites))
  for (i in 2:n_gen) {
    cumul_ens_avg_gen[i,] <- exg_pred[i,] + var_model$Phi%*%cumul_ens_avg_gen[i-1,] + innovations[i,]
    cumul_ens_avg_gen[i,]  <- sapply(cumul_ens_avg_gen[i,],function(x){max(x,0)})
  }
  gc()
  #######################################################################################################
  
  
  ######Step 2: Model the cumulative ensemble spread around the cumulative ensemble mean###################
  
  cumul_ens_resid_gen <- array(0,c(n_sites,n_ens,n_gen))
  cumul_gen <- array(NA,c(n_sites,n_ens,n_gen))
  
  #Individual AR version
  ar_ens_resid_models <- list(n_sites)
  for (j in 1:n_sites) {
    ar_ens_resid_models[[j]] <- list(n_ens)
    for(e in 1:n_ens) {
      #fit AR model to ensemble residual, then simulate from that model with resampled innovations
      ar_ens_resid_models[[j]][[e]] <- arima(hefs_forward_cumul_ens_resid[j,e,ixx_hefs%in%ixx_fit],order=c(1,0,0),include.mean=FALSE)
      innovations <- ar_ens_resid_models[[j]][[e]]$residuals[knn_lst]
      cumul_ens_resid_gen[j,e,] <- arima.sim(n_gen,model=list('ar'=ar_ens_resid_models[[j]][[e]]$coef),innov=innovations)
      
      #use these resampled residuals to estimate individual cumulative ensemble members
      cumul_gen[j,e,] <- cumul_ens_resid_gen[j,e,] + cumul_ens_avg_gen[,j]
      #address any negative values
      cumul_gen[j,e,] <- sapply(cumul_gen[j,e,],function(x) {max(0,x)})
      
    }
  }
  rm(ar_ens_resid_models)
  gc()
  
  ##simple residual resample version
  #for (j in 1:n_sites) {
  #  for (e in 1:n_ens) {
  #    #resample residuals of ensemble members around ensemble mean
  #    #cumul_ens_resid_sim[j,e,] <- hefs_forward_cumul_ens_resid[j,e,knn_lst]
  #    #use these resampled residuals to estimate individual cumulative ensemble members
  #    cumul_sim[j,e,] <- cumul_ens_resid_sim[j,e,] + cumul_ens_avg_sim[,j]
  #    #address any negative values
  #    cumul_sim[j,e,] <- sapply(cumul_sim[j,e,],function(x) {max(0,x)})
      
  #  }
  #}     
      
      
  # #VAR version
  # for (j in 1:n_sites) {
  #   
  #   var_model <- MTS::VAR(x=t(hefs_forward_cumul_ens_resid[j,1:(n_ens-1),]), p=1)
  #   var_resid <- rbind(rep(0,n_ens-1),var_model$residuals) #assume zero residuals for first time step
  #   innovations <- var_resid[knn_lst,]
  #   for (i in 2:n_sim) {
  #     cumul_ens_resid_sim[j,1:(n_ens-1),i] <- var_model$Phi%*%cumul_ens_resid_sim[j,1:(n_ens-1),i-1] + innovations[i,]
  #   }
  #   #the ensemble member residuals are centered around zero by definition
  #   #(they are the difference between the ensemble members and the mean of the ensemble). 
  #   #so to get the last ensemble member, subtract off the sum of the rest of the ensemble member residuals from  0 
  #   cumul_ens_resid_sim[j,n_ens,] <- 0 - apply(cumul_ens_resid_sim[j,1:(n_ens-1),],FUN=sum,2)
  #   
  #   for (e in 1:n_ens) {
  #     #use these resampled residuals to estimate individual cumulative ensemble members
  #     cumul_sim[j,e,] <- cumul_ens_resid_sim[j,e,] + cumul_ens_avg_sim[,j]
  #     #address any negative values
  #     cumul_sim[j,e,] <- sapply(cumul_sim[j,e,],function(x) {max(0,x)})
  #   }
  # }
  
    
  #######################################################################################################
    
    
   
  ######Step 3: Model the fractionation of each ensemble across lead times###################
  
  #knn again to add in some variaiability
  knn_lst <- apply(knn_dist,2,function(x) {sample(c(0,order(x)[1:kk]),size=1,prob=c(0,wts))})
    
  #final array for the ensemble forecasts at all lead times
  all_leads_gen_frac <- array(NA,c(dim(cumul_gen),leads))
  all_leads_gen <- array(NA,c(dim(cumul_gen),leads))
  for (j in 1:n_sites) {
    for(e in 1:n_ens) {
      #knn to resample fractions
      hefs_fit_forward_frac <- hefs_forward_frac[,,ixx_hefs%in%ixx_fit,]
      all_leads_gen_frac[j,e,,] <- hefs_fit_forward_frac[j,e,knn_lst,]
      #set undefined fractions (due to zero totals) to all zeros
      all_leads_gen_frac[j,e,,][which(is.na(all_leads_gen_frac[j,e,,]))] <- 0
      #calculate final 15-day ahead forecasts
      all_leads_gen[j,e,,] <- sapply(1:leads,function(x) {all_leads_gen_frac[j,e,,x]*cumul_gen[j,e,]})
    }
  }
  rm(all_leads_gen_frac)
  gc()
  #######################################################################################################
    
  return(all_leads_gen)  
    
}


#######################################################################################################


# 
# if (verbose) {
#   par(mfrow=c(3,2),mar=c(2,2,1,1))
#   plot(lm_cumul_ens_avg$fitted.values,hefs_forward_cumul_ens_avg)
#   abline(0,1)
#   plot(lm_cumul_ens_avg_final_pred,hefs_forward_cumul_ens_avg)
#   abline(0,1)
#   plot(lm_cumul_ens_avg$fitted.values,lm_cumul_ens_avg_resid)
#   plot(hefs_forward_cumul_ens_avg,lm_cumul_ens_avg_resid)
#   #remaining heterosckadasticity in innovations to be resampled based on obs cumul flows
#   plot(lm_cumul_ens_avg$fitted.values,lm_cumul_ens_avg_resid_model_resid)
#   plot(hefs_forward_cumul_ens_avg,lm_cumul_ens_avg_resid_model_resid)
# }
# 
# if (verbose) {
#   par(mfrow=c(2,2),mar=c(2,2,1,1))
#   keep <- as.POSIXlt(seq(as.Date('1985-10-01'),as.Date('1998-09-30'),'day'))
#   keep_hind <- ixx_hefs_forward%in%keep
#   keep_sim <- ixx_sim%in%keep
#   ymax <- max(cumul_ens_avg_sim,cumul_ens_avg_sim_final,hefs_forward_cumul_ens_avg)
#   plot(ixx_sim[keep_sim],cumul_ens_avg_sim[keep_sim],col="black",lty=1,type="l",ylim=c(0,ymax),lwd=2)
#   lines(ixx_hefs_forward[keep_hind],hefs_forward_cumul_ens_avg[keep_hind],col="red",lty=2,lwd=1)
#   lines(ixx_sim[keep_sim],cumul_ens_avg_sim_final[keep_sim],col='blue',lty=3,lwd=1)
#   
#   keep.sim <- ixx_sim%in%keep
#   keep.obs <- ixx_obs_forward%in%keep
#   plot(ixx_sim[keep.sim],cumul_ens_avg_sim[keep.sim],type="l")
#   lines(ixx_obs_forward[keep.obs],obs_keysite_forward_cumul[keep.obs],col="red",lty=2)
#   lines(ixx_sim[keep.obs],cumul_ens_avg_sim_final[keep.sim],col="blue",lty=2)
#   plot(ixx_sim[keep.sim],simulated_innovations[keep.sim],type="l")
#   plot(ixx_sim[keep.sim],cumul_ens_avg_resid_sim[keep.sim],type="l")
#   
# }
# 
# if (verbose) {
#   par(mfrow=c(2,1),mar=c(2,2,1,1))
#   keep <- as.POSIXlt(seq(as.Date('1994-10-01'),as.Date('1998-09-30'),'day'))
#   keep_hind <- ixx_hefs_forward%in%keep
#   keep_sim <- ixx_sim%in%keep
#   
#   ymax <- max(hefs_forward_cumul_ens_avg,cumul_ens_avg_sim_final,hefs_forward_cumul,hefs_forward_cumul_sim)
#   #original hindcast
#   plot(ixx_hefs_forward[keep_hind],hefs_forward_cumul_ens_avg[keep_hind],type="l",ylim=c(0,ymax))
#   for(j in 1:n_ens) {
#     lines(ixx_hefs_forward[keep_hind],hefs_forward_cumul[j,keep_hind],col="red")
#   }
#   lines(ixx_hefs_forward[keep_hind],hefs_forward_cumul_ens_avg[keep_hind],lwd=2)
#   
#   #synthetic hindcast
#   plot(ixx_sim[keep_sim],cumul_ens_avg_sim_final[keep_sim],type="l",ylim=c(0,ymax))
#   for(j in 1:n_ens) {
#     lines(ixx_sim[keep_sim],hefs_forward_cumul_sim[j,keep_sim],col="grey")
#   }
#   lines(ixx_sim[keep_sim],cumul_ens_avg_sim_final[keep_sim],lwd=2)
# }
# 







# 
# kk <- 3  #k for kNN 
# tot <- sum(rep(1,kk) / 1:kk)
# wts <- rep(1,kk) / 1:kk / rep(tot,kk)
# 
# hefs_forward_cumul_sim_resid <- hefs_forward_cumul
# hefs_forward_cumul_sim <- hefs_forward_cumul
# 
# for (e in 1:n_ens) {
#   
#   cur_sim <- cumul_ens_avg_sim_final  #cumul_ens_avg_sim_final or hefs_forward_cumul_ens_avg
#   cur_hind_knn <- cbind(hefs_forward_cumul_ens_avg[2:n_hind_forward],obs_keysite_forward_cumul_hind[2:n_hind_forward],hefs_forward_cumul_ens_resid[e,1:(n_hind_forward-1)])
#   
#   #start for time = 1
#   cur_sim_knn <- c(cur_sim[1],obs_keysite_forward_cumul_hind[1],0)
#   cur_dist <- sqrt(rowSums((cur_hind_knn - cur_sim_knn)^2))
#   knn_lst <- sample(c(0,order(cur_dist)[1:kk]),size=1,prob=c(0,wts))  
#   hefs_forward_cumul_sim_resid[e,1] <- hefs_forward_cumul_ens_resid[e,knn_lst]
#   temp_hefs_forward <- hefs_forward_cumul_sim_resid[e,1] + cur_sim[1]
#   if (temp_hefs_forward<0) {hefs_forward_cumul_sim_resid[e,1] <- hefs_forward_cumul_sim_resid[e,1] - temp_hefs_forward} 
#   hefs_forward_cumul_sim[e,i] <- hefs_forward_cumul_sim_resid[e,1] + cur_sim[1]
#   lag_sampled_resid <- hefs_forward_cumul_sim_resid[e,1]
#   #loop through rest of the time steps
#   for (i in 2:n_hind_forward) {
#     cur_sim_knn <- c(cur_sim[i],obs_keysite_forward_cumul_hind[i],lag_sampled_resid)
#     cur_dist <- sqrt(rowSums((cur_hind_knn - cur_sim_knn)^2))
#     knn_lst <- sample(c(0,order(cur_dist)[1:kk]),size=1,prob=c(0,wts))  
#     hefs_forward_cumul_sim_resid[e,i] <- hefs_forward_cumul_ens_resid[e,knn_lst]
#     temp_hefs_forward <- hefs_forward_cumul_sim_resid[e,i] + cur_sim[i]
#     if (temp_hefs_forward<0) {hefs_forward_cumul_sim_resid[e,i] <- hefs_forward_cumul_sim_resid[e,i] - temp_hefs_forward} 
#     hefs_forward_cumul_sim[e,i] <- hefs_forward_cumul_sim_resid[e,i] + cur_sim[i]
#     lag_sampled_resid <- hefs_forward_cumul_sim_resid[e,i]
#   }
#   
# }
