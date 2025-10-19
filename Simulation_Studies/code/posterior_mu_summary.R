#posterior summary of induced posterior of mu obtained using regression
# posterior_mu_summary.R

post_mu_summary = function(mu.samples = NULL, mu_true = NULL){
  p.mean = mean(mu.samples)
  p.var = var(mu.samples)
  p.quan = quantile(mu.samples, probs = c(0.025, 0.975))
  p.cover = mu_true > p.quan[1] & mu_true < p.quan[2]
  
  p.summary = c(p.mean, p.var, p.quan, p.cover)
  return(p.summary = p.summary)
}