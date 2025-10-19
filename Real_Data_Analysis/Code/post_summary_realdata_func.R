# post_summary_realdata_func.R
post_summary_RD_fn = function(mu_samples = NULL){
  pmean = mean(mu_samples)
  pvar = var(mu_samples)
  pquan = quantile(mu_samples, probs = c(0.025, 0.975))
  # pcov = unname((theta0 > pquan[1] & theta0 < pquan[2])*1)
  plen = unname(pquan[2] - pquan[1])
  psummary = c(pmean, pvar, pquan, plen)
  
  # remove unnecessary stuff
  remove(pmean, pvar, pquan, plen)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  return(list(psummary = psummary))
}