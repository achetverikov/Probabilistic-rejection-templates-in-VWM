LL_singlepeak2_mle <- function( template_sd, scaling, sigma, mu, perc_noise_sd=0, exp) {
  n_ctpd <- length(ctpd)
  dmin = 20
  dmax = ifelse(exp=='e2', 30, 40)
  
  dmean <- (dmin+dmax)/2#ifelse(runif(n_ctpd)>=0.5, (dmin+dmax)/2, -(dmin+dmax)/2)
  probs <- dnorm(ctpd+rnorm(n_ctpd, 0, perc_noise_sd), dmean, template_sd)
  pred_rts <- mu+scaling*(probs)
  R = dnorm(lrtc-pred_rts, 0, sigma, log=T)
  if (any(is.infinite(R))){
    print(list(template_sd, mu, sigma, scaling, perc_noise_sd))
  }
  -sum(R)
}

get_sp_preds_mle<-function( template_sd, scaling, sigma, mu, perc_noise_sd=0, ctpd, exp='e2'){
  n_ctpd <- length(ctpd)
  dmin = 20
  dmax = ifelse(exp=='e2', 30, 40)
  dmean <- (dmin+dmax)/2#ifelse(runif(n_ctpd)>=0.5, (dmin+dmax)/2, -(dmin+dmax)/2)
  probs <- dnorm(ctpd+rnorm(n_ctpd, 0, perc_noise_sd), dmean, template_sd)
  pred_rts <- mu+scaling*probs
  data.frame(y = pred_rts, ymin = pred_rts-1.96*sigma,  ymax = pred_rts+1.96*sigma)
}
LL_averaged_mle <- function(template_sd, mu, sigma, scaling, perc_noise_sd=0) {
  n_ctpd <- length(ctpd)
  probs <- dnorm(ctpd+rnorm(n_ctpd, 0, perc_noise_sd), 0, template_sd)
  pred_rts <- mu+scaling*probs
  R = dnorm(lrtc-pred_rts, 0, sigma, log=T)
  if (any(is.infinite(R))){
    print(list(template_sd, mu, sigma, scaling, perc_noise_sd))
  }
  -sum(R)
}
get_averaged_preds_mle<-function( template_sd, scaling, sigma, mu, perc_noise_sd=0, ctpd, exp='e2'){
  n_ctpd <- length(ctpd)
  probs <- dnorm(ctpd+rnorm(n_ctpd, 0, perc_noise_sd), 0, template_sd)
  pred_rts <- mu+scaling*probs
  
  data.frame(y = pred_rts, ymin = pred_rts-1.96*sigma,  ymax = pred_rts+1.96*sigma)
}
LL_twopeaks_mle <- function( template_sd, scaling, sigma, mu, perc_noise_sd=0, exp) {
  n_ctpd <- length(ctpd)
  dmin = 20
  dmax = ifelse(exp=='e2', 30, 40)
  
  dmean <- (dmin+dmax)/2
  probs <- (dnorm(ctpd+rnorm(n_ctpd, 0, perc_noise_sd), dmean, template_sd)+dnorm(ctpd+rnorm(n_ctpd, 0, perc_noise_sd), -dmean, template_sd))/2
  pred_rts <- mu+scaling*(probs)
  R = dnorm(lrtc-pred_rts, 0, sigma, log=T)
  -sum(R)
}

get_twopeaks_preds_mle<-function( template_sd, scaling, sigma, mu, perc_noise_sd=0, ctpd, exp='e2'){
  n_ctpd <- length(ctpd)
  dmin = 20
  dmax = ifelse(exp=='e2', 30, 40)
  
  dmean <- (dmin+dmax)/2
  probs <- (dnorm(ctpd+rnorm(n_ctpd, 0, perc_noise_sd), dmean, template_sd)+dnorm(ctpd+rnorm(n_ctpd, 0, perc_noise_sd), -dmean, template_sd))/2
  pred_rts <- mu+scaling*(probs)
  
  data.frame(y = pred_rts, ymin = pred_rts-1.96*sigma,  ymax = pred_rts+1.96*sigma)
}