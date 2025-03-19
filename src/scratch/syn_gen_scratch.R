
syn_gen <- function (seed,kk,keysite,knn_pwr,fit_start,fit_end,gen_start,gen_end,leave_out_years,use_mpi,
obs_forward_all_leads_gen,obs_forward_all_leads_hind,hefs_forward_cumul,hefs_forward_cumul_ens_avg,hefs_forward_cumul_ens_resid,hefs_forward_frac,ixx_hefs,
ixx_obs_forward,varx_fun,correct_leads,lds_to_correct,fix_order,use_ar) {
  
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
  decay <- (1:(leads+1))^my_power / sum(((1:leads)^my_power)) #k dim vector
  
  #plot(1:leads,decay,type='l',ylim=c(0,1),lwd=2,xlab='Leads',ylab='decay function value')
  #idx = c(-seq(0.5,5,1))
  #for (i in 1:length(idx)){
    #my_power <- idx[i]
    #decay <- (1:leads)^my_power / sum(((1:leads)^my_power))
    #lines(1:leads,decay,ylim=c(0,0.5),col=i,lwd=2)
  #}
  #legend('topright',paste('D =',idx),col=1:length(idx),lwd=2)
  
  #set the key observation-based covariate values for the simulation
  obs_forward_all_leads_gen <- obs_forward_all_leads[,ixx_obs_forward%in%ixx_gen,,drop=FALSE]
  n_gen <- length(ixx_gen)
  
  #syn_hefs_out<-foreach(m = 1:n_samp,.packages=c('MTS'),.inorder=T)%dopar%{
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  gen_knn_data <- t(obs_forward_all_leads_gen[keysite,,])
  fit_knn_data <- t(obs_forward_all_leads_hind[keysite,ixx_hefs_obs_fwd%in%ixx_fit,]) #k x n matrix
  ##gen_knn_data <- apply(obs_forward_all_leads_gen[keysite,,],1,cumsum) #k x n matrix
  ##fit_knn_data <- apply(obs_forward_all_leads_hind[keysite,ixx_hefs_obs_fwd%in%ixx_fit,],1,cumsum) #k x n matrix
  ##gen_knn_data <- rbind(obs$NHGC1[ixx_obs%in%ixx_gen],apply(obs_forward_all_leads_gen[keysite,,],1,cumsum)) #k x n matrix
  ##fit_knn_data <- rbind(obs$NHGC1[ixx_obs%in%ixx_fit],apply(obs_forward_all_leads_hind[keysite,ixx_hefs_obs_fwd%in%ixx_fit,],1,cumsum)) #k x n matrix
  ##gen_knn_data <- rbind(obs$NHGC1[ixx_obs%in%ixx_fit],c(0,obs$NHGC1[ixx_obs%in%ixx_fit][-length(ixx_fit)]),apply(obs_forward_all_leads_gen[keysite,,],1,cumsum)) #k x n matrix
  ##fit_knn_data <- rbind(obs$NHGC1[ixx_obs%in%ixx_fit],c(0,obs$NHGC1[ixx_obs%in%ixx_fit][-length(ixx_fit)]),apply(obs_forward_all_leads_hind[keysite,ixx_hefs_obs_fwd%in%ixx_fit,],1,cumsum)) #k x n matrix
  
  gen_knn_data <- rbind(
    gen_knn_data[1,c(1,1:(ncol(gen_knn_data)-1)),drop=FALSE],
    gen_knn_data)
  fit_knn_data <- rbind(
    fit_knn_data[1,c(1,1:(ncol(fit_knn_data)-1)),drop=FALSE],
    fit_knn_data)
  

  ##frac_fun <- function(x){out<-x / sum(x);return(out)}
  #apply function over all indices in the generation obs
  knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    #squared differences between the k x 1 vector of gen obs data and the k x n vector of fit data
    ##frac_fit <- apply(fit_knn_data,2,frac_fun)
    ##gen_knn_data[,x] <- frac_fun(gen_knn_data[,x])
    ##frac_diff <- (gen_knn_data[,x] - frac_fit)^2
    sq_diff <- (gen_knn_data[,x] - fit_knn_data)^2
    ##sq_diff_scale <- t(apply(sq_diff,1,function(x){x_=scale(x,center=F);return(x_)}))  
    ##sq_cumul_diff <- t(scale((sum(gen_knn_data[,x]) - apply(fit_knn_data,2,sum)^2),center=F))
    #multiply squared differences by decay function to accentuate differences at certain leads (typically emphasize early leads)
    ##comb_diff <- rbind(frac_diff,sq_cumul_diff)
    #decay_diff <- comb_diff*decay[1:4]
    decay_diff <- sq_diff *decay
    #take square root of difference summed over leads to get Euclidean distance modified by decay function
    sqrt(colSums(decay_diff))
    #colSums(decay_diff)
    #sqrt(colSums(decay*(gen_knn_data[,x] - fit_knn_data)^2))
  })
  #knn_dist is a n_fit x n_gen matrix of distances
  
  #apply function over all indices in the generation obs
  ##knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    #squared differences between the k x 1 vector of gen obs data and the k x n vector of fit data
    ##abs_diff <- abs(sum(gen_knn_data[,x]) - apply(fit_knn_data,2,sum)) 
  ##})
  #knn_dist is a n_fit x n_gen matrix of distances
  
  #the resampled locations via KNN to use
  knn_lst <- apply(knn_dist,2,function(x){
    #order the distances along the fit dimension yielding a value for each of the n_gen dimensions
    ord_dist <- order(x)[1:kk] #order function sorts in ascending order and returns indices for the top 1:kk distances
    #output value is a distance weighted sampling of the 1:kk ordered indices
    sample(c(0,ord_dist),size=1,prob=c(0,wts))})
    ##order(x)[1]})
  
  #knn_lst is a n_gen length vector of indices to sample from the fitted period in following routines
  
  #rm(gen_knn_data,fit_knn_data)
  #gc()
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
  
  #a new decay function, this time just focusing on the first 2 days of observed flow, and a heavy weight on the HEFS cumulative ensemble mean forecast
  ##my_power <- knn_pwr    #larger neg values = faster decay, more weight on early leads
  ##decay <- (1:leads)^my_power / sum(((1:leads)^my_power)) #k dim vector
  ##decay <- c(1,1,rep(0,13),2)
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #we also include the cumulative ensemble forecast total as another variable on which to resample
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  ##gen_knn_data <- rbind(t(obs_forward_all_leads_gen[keysite,,]),cumul_ens_avg_gen[,keysite])
  ##hind_knn_data <- rbind(t(obs_forward_all_leads_hind[keysite,ixx_hefs%in%ixx_fit,]),t(hefs_forward_cumul_ens_avg[,ixx_hefs%in%ixx_fit,drop=FALSE])[,keysite])
  
  ##gen_knn_data <- t(apply(gen_knn_data,1,function(x) {scale(x)[,1]}))
  ##hind_knn_data <- t(apply(hind_knn_data,1,function(x) {scale(x)[,1]}))
  
  ##knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    ##sqrt(colSums(decay*(gen_knn_data[,x] - hind_knn_data)^2))
  ##})
  
  #the resampled locations via KNN to use
  ##knn_lst <- apply(knn_dist,2,function(x) {sample(c(0,order(x)[1:kk]),size=1,prob=c(0,wts))})
  
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
  #a new decay function, this time just focusing on the first 2 days of observed flow, and a heavy weight on the HEFS cumulative ensemble mean forecast
  #decay <- c(1,1,rep(0,13),2)
  #my_power <- knn_pwr    #larger neg values = faster decay, more weight on early leads
  #y = 20
  #decay <- (1:y)^my_power / sum(((1:y)^my_power)) #k dim vector
  
  #this calculates the distance between each 1:leads set of obs in the simulation period and the hindcast period
  #we also include the cumulative ensemble forecast total as another variable on which to resample
  #knn_dist is a matrix of distances, with dimensions (n_hind_forward,n_sim)
  #note, this distance is only calculated for the keysite
  ##gen_knn_data <- rbind(t(obs_forward_all_leads_gen[keysite,,])) #,cumul_ens_avg_gen[,keysite])
  ##hind_knn_data <- rbind(t(obs_forward_all_leads_hind[keysite,ixx_hefs_obs_fwd%in%ixx_fit,])) #,t(hefs_forward_cumul_ens_avg[,ixx_hefs%in%ixx_fit,drop=FALSE])[,keysite])
  
  ##gen_knn_data <- t(apply(gen_knn_data,1,function(x) {scale(x)[,1]}))
  ##hind_knn_data <- t(apply(hind_knn_data,1,function(x) {scale(x)[,1]}))
  
  ##knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
    ##sqrt(colSums(decay*(gen_knn_data[,x] - hind_knn_data)^2))
  ##})
  
  ###hefs_forward_frac[is.na(hefs_forward_frac)==T]<-0
  
  ###frac_fun <- function(x){
    ###out<-x / sum(x)
    ###out[is.na(out)==T]<-0
    ###return(out)}
  
  ###obs_frac_gen <- apply(obs_forward_all_leads_gen[keysite,,],1,frac_fun)
  ###hefs_frac_mean_fit <- apply(hefs_forward_frac[keysite,,ixx_hefs%in%ixx_fit,],c(2,3),function(x){mean(x,na.rm=T)})
  ###gen_cumul_sd <- apply(cumul_gen[keysite,,],2,sd)
  ###fit_cumul_sd <- apply(hefs_forward_cumul[keysite,,ixx_hefs%in%ixx_fit],2,sd)
  ###fit_cumul_obs_fwd <- apply(obs_forward_all_leads_hind[keysite,ixx_hefs%in%ixx_fit,],1,sum)
  ###gen_cumul_obs_fwd <- apply(obs_forward_all_leads_gen[keysite,,],1,sum)
  
  ###gen_knn_data <- rbind(gen_cumul_obs_fwd,cumul_ens_avg_gen[,keysite],gen_cumul_sd,obs_frac_gen)
  ###hind_knn_data <- rbind(fit_cumul_obs_fwd,t(hefs_forward_cumul_ens_avg[,ixx_hefs%in%ixx_fit,drop=FALSE])[,keysite],fit_cumul_sd,t(hefs_frac_mean_fit))
  
  ###gen_knn_data <- t(apply(gen_knn_data,1,function(x) {scale(x)[,1]}))
  ###hind_knn_data <- t(apply(hind_knn_data,1,function(x) {scale(x)[,1]}))
  
  ###knn_dist <- sapply(1:ncol(gen_knn_data),function(x){
  ###sqrt(colSums((gen_knn_data[,x] - hind_knn_data)^2))
  ###})
  
  #the resampled locations via KNN to use
  ##knn_lst <- apply(knn_dist,2,function(x) {sample(c(0,order(x)[1:kk]),size=1,prob=c(0,wts))})
  
  
  ##fit_knn_hefs_data <- hefs_forward_cumul[keysite,,ixx_hefs%in%ixx_fit]
  ##gen_knn_shefs_data <- cumul_gen[keysite,,]
  ##fit_cumul_obs_fwd <- apply(obs_forward_all_leads_hind[keysite,ixx_hefs_obs_fwd%in%ixx_fit,],1,sum)
  ##gen_cumul_obs_fwd <- apply(obs_forward_all_leads_gen[keysite,,],1,sum)
  ##fit_cumul_hefs <- hefs_forward_cumul_ens_avg[keysite,ixx_hefs%in%ixx_fit]
  
  ##rnk_fun <- function(x){
    ##rnk_cdf <- rank(x,ties.method = 'random') / (length(x)+1)
    ##return(rnk_cdf[1])
  ##}
  
  ##plan(multicore)
  
  ##knn_dist <- future_sapply(1:ncol(gen_knn_data),function(x){
    ##sq_diff_skew <- (t(scale(skewness(gen_knn_shefs_data[,x]) - apply(fit_knn_hefs_data,2,skewness))))^2
    ##sq_diff_kurt <- (t(scale(kurtosis(gen_knn_shefs_data[,x]) - apply(fit_knn_hefs_data,2,kurtosis))))^2
    ##sq_diff_sd <- (t(scale(sd(gen_knn_shefs_data[,x]) - apply(fit_knn_hefs_data,2,sd))))^2
    ##sq_diff_mn <- (t(scale(mean(gen_knn_shefs_data[,x]) - apply(fit_knn_hefs_data,2,mean))))^2
    ##ob_hefs_gen <- c(gen_cumul_obs_fwd[x],gen_knn_shefs_data[,x])
    ##rnk_cdf_gen <- rnk_fun(ob_hefs_gen)
    ##rnk_cdf_fit <- apply(rbind(fit_cumul_obs_fwd,fit_knn_hefs_data),2,rnk_fun)
    ##sq_rnk <- (rnk_cdf_gen - rnk_cdf_fit)^2
    ##sq_obs_seq <- (gen_knn_data[,x] - hind_knn_data)^2
    
    ##sq_diff <- rbind(sq_rnk,sq_diff_mn,sq_diff_sd,sq_diff_kurt,sq_diff_skew,sq_obs_seq)
    ##sq_diff[is.na(sq_diff)==T] <- 0

    ##decay_diff <- sq_diff *decay
    ##sqrt(colSums(decay_diff))
  ##})

  
  #the resampled locations via KNN to use
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
    #1d lead mean shift
    hefs_fwd_fit <- hefs_forward[j,,ixx_hefs%in%ixx_fit,]
    hefs_fwd_cumul_fit <- hefs_forward_cumul[j,,ixx_hefs%in%ixx_fit]
    obs_fwd_fit <- obs_forward_all_leads_hind[j,ixx_hefs_obs_fwd%in%ixx_fit,]
    cumul_gen_step <- cumul_gen[j,,]
    
    if(correct_leads==TRUE){
      for(k in lds_to_correct){
        hefs_mn <- apply(hefs_fwd_fit[,knn_lst,k],2,mean)
        pred_mn <- apply(all_leads_gen[j,,,k],2,mean)
        mn_scale <- hefs_mn / obs_fwd_fit[knn_lst,k] * obs_forward_all_leads_gen[j,,k] / pred_mn
        mn_scale[is.na(mn_scale)==T|mn_scale==0|mn_scale==Inf|mn_scale==-Inf]<-1
        mn_add <- pred_mn * mn_scale - pred_mn
        ##all_leads_gen[j,,,k] <- all_leads_gen[j,,,k] * matrix(rep(mn_scale,each=n_ens),ncol=length(ixx_gen),byrow=F)
        all_leads_gen[j,,,k] <- all_leads_gen[j,,,k] + matrix(rep(mn_add,each=n_ens),ncol=length(ixx_gen),byrow=F)
        ##obs_scale <- obs_forward_all_leads_gen[j,,k] / obs_fwd_fit[knn_lst,k]
        ##obs_scale[is.na(obs_scale)==T|obs_scale==0|obs_scale==Inf|obs_scale==-Inf]<-1
        ##cumul_obs_scale <- cumul_ens_avg_gen[,j] / apply(hefs_fwd_cumul_fit[,knn_lst],2,mean)
        ##cumul_obs_scale[is.na(cumul_obs_scale)==T|cumul_obs_scale==0|cumul_obs_scale==Inf|cumul_obs_scale==-Inf]<-1
        ##rng_scale <- (apply(hefs_fwd_fit[,knn_lst,k],2,function(x){out = diff(range(x))})/hefs_mn) * (pred_mn+mn_add) / apply(all_leads_gen[j,,,k],2,function(x){out = diff(range(x))})
        ##rng_scale[is.na(rng_scale)==T|rng_scale==0|rng_scale==Inf|rng_scale==-Inf]<-1
        sd_scale <- (apply(hefs_fwd_fit[,knn_lst,k],2,function(x){sd(x,na.rm=TRUE)})/hefs_mn) * (pred_mn+mn_add) / apply(all_leads_gen[j,,,k],2,function(x){sd(x,na.rm=TRUE)})
        sd_scale[is.na(sd_scale)==T|sd_scale==0|sd_scale==Inf|sd_scale==-Inf]<-1
    
        env = 10
        #rng_scale[abs(rng_scale)>env]<-1
        sd_scale[abs(sd_scale)>env]<-1
        #scale variance by different scaling techniques
        #all_leads_gen[j,,,k] <- (all_leads_gen[j,,,k] - matrix(rep((pred_mn + mn_add),each=n_ens),ncol=length(ixx_gen),byrow=F)) * matrix(rep(rng_scale,each=n_ens),ncol=length(ixx_gen),byrow=F) + matrix(rep((pred_mn + mn_add),each=n_ens),ncol=length(ixx_gen),byrow=F)
        all_leads_gen[j,,,k] <- (all_leads_gen[j,,,k] - matrix(rep((pred_mn + mn_add),each=n_ens),ncol=length(ixx_gen),byrow=F)) * matrix(rep(sd_scale,each=n_ens),ncol=length(ixx_gen),byrow=F) + matrix(rep((pred_mn + mn_add),each=n_ens),ncol=length(ixx_gen),byrow=F)
        
        #refractionate remaining cumulative ensemble
        cumul_gen_step[,] = cumul_gen[j,,] - apply(all_leads_gen[j,,,1:k],c(1,2),sum)
        cumul_gen_step[cumul_gen_step<0]<-0
        cumul_gen_mat <- array(NA,c(dim(cumul_gen_step),(leads-k)))
        for(i in 1:(leads-k)){
          cumul_gen_mat[,,i] <- cumul_gen_step
        }
        hefs_frac_step <- aperm(apply(hefs_fwd_fit[,knn_lst,(k+1):leads],c(1,2),function(x){x/sum(x)}),c(2,3,1))
        hefs_frac_step[is.na(hefs_frac_step)|hefs_frac_step==Inf|hefs_frac_step==-Inf|hefs_frac_step<0]<-0
        all_leads_gen[j,,,(k+1):leads] <- hefs_frac_step[,,] * cumul_gen_mat[,,]}
    }
  }
  rm(all_leads_gen_frac,cumul_gen,hefs_fit_forward_frac,hefs_forward_frac,knn_lst,hefs_fwd_fit,hefs_fwd_cumul_fit,obs_fwd_fit,cumul_gen_step)
  gc()

  all_leads_gen[is.na(all_leads_gen)|all_leads_gen==Inf|all_leads_gen==-Inf|all_leads_gen<0]<-0
    
  return(all_leads_gen)  
  
  rm(list = ls());gc()
}

