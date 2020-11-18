library(apastats)
load.libs(c('circular','doSNOW','foreach','doParallel','mgcv','ggplot2','data.table','DEoptim','bbmle','egg','stringr','plyr','dplyr','scales','Hmisc', 'ggthemes'))
source('get_preds_funcs_mle2.R')
mode_from_dens<-function(x){
  dens<-density(x)
  dens$x[dens$y==max(dens$y)]
}
get_template_mu<-function(exp, peaks){
  if (length(peaks)>1) {
    peaks = peaks[1]
    
  }
  if (length(exp)>1) {
    exp = exp[1]
  }
  if (peaks=='avgd'){
    template_mu<-0
  } else if (peaks=='single'){
    template_mu<-ifelse(exp=='e1', 30, 25)
  } else{
    template_mu<-case_when(exp=='e1'~ c(-30,30), T~c(-25,25))
  }
  template_mu
}
cbind.fill <- function(...){
  nm <- list(...) 
  dfdetect <- grepl("data.frame|matrix", unlist(lapply(nm, function(cl) paste(class(cl), collapse = " ") )))
  # first cbind vectors together 
  vec <- data.frame(nm[!dfdetect])
  n <- max(sapply(nm[dfdetect], nrow)) 
  vec <- data.frame(lapply(vec, function(x) rep(x, n)))
  if (nrow(vec) > 0) nm <- c(nm[dfdetect], list(vec))
  nm <- lapply(nm, as.data.frame)
  
  do.call(cbind, lapply(nm, function (df1) 
    rbind(df1, as.data.frame(matrix(NA, ncol = ncol(df1), nrow = n-nrow(df1), dimnames = list(NULL, names(df1))))) )) 
}


st_data_e2<-single_target_data[expName=='Main Experiment'&!is.na(lrt)]
st_data_e1<-single_target_data[expName=='Suppl. Experiment'&!is.na(lrt)]

boot_mle2<-function(exp, peaks, resample = T){
  if (resample == T){
      if (exp=='e1'){
        data_s <- st_data_e1[,.SD[sample(.N, .N, replace = T)], by = .(subjectId,stimTypef)]
        distractors <- c(runif(18, -40, -20), runif(18, 20, 40))
        
      } else if (exp=='e2') {
        data_s <- st_data_e2[,.SD[sample(.N, .N, replace = T)], by = .(subjectId,stimTypef)]
        distractors <- c(runif(18, -30, -20), runif(18, 20, 30))
      }
  } else {
    print('No resampling')
    if (exp=='e1'){
      data_s <- st_data_e1
      distractors <- c(runif(18, -40, -20), runif(18, 20, 40))
      
    } else if (exp=='e2') {
      data_s <- st_data_e2
      distractors <- c(runif(18, -30, -20), runif(18, 20, 30))
    }
  }
  
  if (peaks == 'single')
    LL_fun<-LL_singlepeak2_mle
  else if (peaks == 'two')
    LL_fun<-LL_twopeaks_mle
  else if (peaks=='avgd')
    LL_fun<-LL_averaged_mle
  
  res<-mle2(LL_fun,
            data = list(ctpd=data_s$ctpd, lrtc=data_s$lrtc, exp=exp), 
            method = 'L-BFGS-B', 
            start = list(template_sd = sd(distractors)*0.5, scaling = 10+rnorm(1,10,5), sigma = mysd(data_s$lrtc), mu=0),  
            lower = c(template_sd = 1e-5, scaling = -1000, sigma=1e-5, mu=-10), 
            upper = c(template_sd = 180, scaling = 1000, sigma=10, mu=10), 
            vecpar = F, trace=T)
  res_dt<-as.data.table(as.list(coef(res)))
  res_dt$ll<-as.numeric(logLik(res))
  res_dt$bic<-as.numeric(BIC(res))
  res_dt$exp<-exp
  res_dt$peaks<-peaks
  res_dt
}

nrep = 500

par_df <- expand.grid(n=1:nrep, exp = c('e1','e2'), peaks=c('single', 'two', 'avgd')) 

cl<-makeCluster(8)
registerDoSNOW(cl)

pb <- txtProgressBar(max = nrow(par_df), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)

res_mle_three_models<-foreach(n=par_df$n, exp=par_df$exp, peaks=par_df$peaks, .combine = rbind, .packages = c('data.table','mgcv','apastats','bbmle'), .options.snow = list(progress=progress)) %dopar% boot_mle2(exp, peaks)

### get the p_first predictions 

for_scaling<-datam2_fp[,.(p_first=mymean( p_first),.N),by=.(ctpd, ctpd_other,targets_fp, stimType, subjectId)]

get_dens_ratio<-function(template_mu, template_sd, ctpd, ctpd_other){
  rowMeans(sapply(template_mu, function(mu) dnorm(ctpd_other, mu, template_sd)))/rowMeans(sapply(template_mu, function(mu) dnorm(ctpd, mu, template_sd)))
}
get_p_first_probs<-function(template_mu, template_sd, scaling, ctpd, ctpd_other){
  probs<-0.5+scaling*log(get_dens_ratio(template_mu, template_sd, ctpd, ctpd_other))
  probs
}
p_first_LL<-function(template_mu, template_sd, scaling){
  probs<-get_p_first_probs(template_mu, template_sd, scaling, ctpd, ctpd_other)  
  -sum(dbinom(round(p_first*N), N, probs, log=T))
}

template_sds<-res_mle_three_models[,.(template_sd=mode_from_dens(template_sd), scaling=mean(scaling), mu=mean(mu), sigma=mean(sigma)) ,by=.(exp,peaks)]

p_first_scaling<-ddply(template_sds[exp=='e2'], .variables = .(peaks, exp), .fun = function (df){
  template_mu<-get_template_mu(df$exp, df$peaks)
  
  data<-list(ctpd=for_scaling$ctpd, ctpd_other=for_scaling$ctpd_other, template_sd=template_sds[exp==df$exp&peaks==df$peaks, template_sd], p_first=for_scaling$p_first, template_mu=template_mu, N=for_scaling$N)
  res<-as.vector(coef(m<-mle2(p_first_LL, start = list(scaling = 0),  data = data, vecpar = F)))
  true_bic <- log(nrow(for_scaling))*2+deviance(m)
  
  
  cbind(for_scaling, 
        true_bic = true_bic,
        dens_ratio = get_dens_ratio(template_mu=template_mu, template_sd=data$template_sd, ctpd=data$ctpd, ctpd_other=data$ctpd_other),
        log_dens_ratio = log(get_dens_ratio(template_mu=template_mu, template_sd=data$template_sd, ctpd=data$ctpd, ctpd_other=data$ctpd_other)),
        pf_scaling = rep((res), nrow(for_scaling)),
        prob = get_p_first_probs(template_mu=template_mu, template_sd=data$template_sd, scaling=mean(res), ctpd=data$ctpd, ctpd_other=data$ctpd_other))
})

setDT(p_first_scaling)
p_first_scaling[,typefp:=interaction(substr(stimType,1,1),targets_fp, sep='|')]

### get single target RT predictions

for_st_rt_scaling<-single_target_data[,.(mean_rt=(lrtc), grandmean), keyby=.(short_exp, subjectId, stimTypef, mean_ctpd)]

rt_scaling<-ddply(template_sds, .variables = .(peaks, exp), .fun = myfun<-function (df){
  template_mu <- get_template_mu(df$exp, df$peaks)
  
  data<-list(ctpd=for_st_rt_scaling[short_exp==df$exp,mean_ctpd], 
             template_sd=df$template_sd, 
             lrtc=for_st_rt_scaling[short_exp==df$exp,mean_rt], 
             template_mu=template_mu, 
             exp=df$exp)
  
  if (df$peaks == 'single'){
    LL_fun <- LL_singlepeak2_mle
    pred_fun <- get_sp_preds_mle
  } else if (df$peaks == 'two'){
    LL_fun <- LL_twopeaks_mle
    pred_fun <- get_twopeaks_preds_mle
  } else if (df$peaks=='avgd'){
    LL_fun <- LL_averaged_mle
    pred_fun <- get_averaged_preds_mle
  }
  
  res<-coef(m<-mle2(LL_fun,
                 data = data,
                 method = 'L-BFGS-B', 
                 start = list(scaling = 0, sigma = 1e3, mu=0),  
                 lower = c(scaling = -1000, sigma=1e-5, mu=-10), 
                 upper = c(scaling = 1000, sigma=10, mu=10), 
                 vecpar = F, trace=F))
  
  mod_dev <- deviance(m)
  
  # 4 comes from the number of parameters (templpate_sd, scaling, sigma, mu)
  true_bic <- log(length(data$lrtc))*4+mod_dev
  
  cbind(for_st_rt_scaling[short_exp==df$exp], 
        true_bic = true_bic,
        t(res),
        pred_fun(template_sd=data$template_sd, scaling=res['scaling'], sigma=res['sigma'], mu=res['mu'], ctpd=data$ctpd, exp=data$exp))
})

setDT(rt_scaling)

### get two-targets RT predictions

for_two_t_rt_scaling<-data_r_e2[totCorrect==2&trial==0,.(mean_rt=(lrt_2c), grandmean), keyby=.(subjectId, short_exp, targets_fp, ctpd_1=pmin(ctpd_1, ctpd_2), ctpd_2=pmax(ctpd_1, ctpd_2))]

get_dens_sum<-function(template_mu, template_sd, ctpd, ctpd_other){
  0.5*(rowMeans(sapply(template_mu, function(mu) dnorm(ctpd_other, mu, template_sd))) + 
       rowMeans(sapply(template_mu, function(mu) dnorm(ctpd, mu, template_sd))))
}
get_two_rt_preds<-function(template_mu, template_sd, scaling, mu, ctpd, ctpd_other){
  pred_rts<-mu+scaling*get_dens_sum(template_mu, template_sd, ctpd, ctpd_other)
  pred_rts
}
two_rt_LL<-function(template_mu, template_sd, scaling, mu, sigma, exp){
  pred_rts<-get_two_rt_preds(template_mu, template_sd, scaling, mu, ctpd_1, ctpd_2)  
  -sum(dnorm(lrtc-pred_rts, 0, sigma, log=T))
}


rt_two_t_scaling <- ddply(template_sds[exp=='e2'], .variables = .(peaks, exp), .fun = myfun<-function (df){
  template_mu <- get_template_mu(df$exp, df$peaks)
  dens_sum<-get_dens_sum(template_mu = template_mu, template_sd = df$template_sd, ctpd = for_two_t_rt_scaling$ctpd_1, ctpd_other = for_two_t_rt_scaling$ctpd_2)
  
  data<-list(ctpd_1 = for_two_t_rt_scaling$ctpd_1, 
             ctpd_2 = for_two_t_rt_scaling$ctpd_2, 
             template_sd = df$template_sd, 
             lrtc = for_two_t_rt_scaling[short_exp==df$exp,mean_rt], 
             template_mu = template_mu, 
             dens_sum = dens_sum,
             exp = df$exp)
  
  res<-coef(m<-mle2(two_rt_LL,
                 data = data,
                 method = 'L-BFGS-B', 
                 start = list(scaling = 0, sigma = 10, mu=0),  
                 lower = c(scaling = -1000, sigma=1e-5, mu=-30), 
                 upper = c(scaling = 1000, sigma=20, mu=20), 
                 vecpar = F, trace=T, control=list(maxit = 8000)))
  bic<-BIC(m)
  mod_ll <- logLik(m)
  mod_dev <- deviance(m)
  
  # 4 comes from the number of parameters (templpate_sd, scaling, sigma, mu)
  true_bic <- log(length(data$lrtc))*4+mod_dev
  
  cbind(for_two_t_rt_scaling[short_exp==df$exp], 
        t(res),
        exp_probs = dens_sum,
        true_bic = true_bic,
        pred_rt = get_two_rt_preds(template_mu = template_mu, 
                                   template_sd = data$template_sd, 
                                   scaling = res['scaling'], 
                                   mu = res['mu'], 
                                   ctpd = data$ctpd_1, ctpd_other = data$ctpd_2))
})
setDT(rt_two_t_scaling)
rt_two_t_scaling[,targets_fp:=factor(targets_fp, levels= c('O+B','P+O','P+B','P+P'))]
rt_two_t_scaling[order(peaks, targets_fp)]

bic_tbl <- merge(merge(rt_scaling[,.(bic_rt_st = true_bic[1]),keyby=.(exp, peaks)],  p_first_scaling[,.(bic_pf=true_bic[1]),keyby=.(exp, peaks)], all = T), rt_two_t_scaling[,.(bic_rt_2t = true_bic[1]),keyby=.(exp, peaks)], all=T)

save(res_mle_three_models, rt_scaling, rt_two_t_scaling, p_first_scaling, template_sds, bic_tbl, file = 'simulations_results.RData')