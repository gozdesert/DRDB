
# oracle_estimator_mu.R 02/26 only linear case
oracle_est_mu_func = function(Y = NULL, X = NULL, trt = NULL, alpha1 = NULL, alpha2 = NULL, 
                                           beta1 = NULL,
                                           beta2 = NULL, beta12 =NULL, beta3 = NULL,
                                           model_type = NULL, coeff = NULL){
  if(model_type == "L-L"){
    a = as.vector(X%*%beta3) - 0.08
    piX = exp(a)/(1 + exp(a))
    # trt arm part
    m1Xstar = 2*X%*%beta1 + alpha1 
    ipw1 = trt*(Y - m1Xstar)/piX
    fstar1 = m1Xstar + ipw1
    oracle_mu1 = mean(fstar1)
    # control arm part
    m0Xstar = alpha2 + X%*%beta2
    ipw0 = (1-trt)*(Y - m0Xstar)/(1-piX)
    fstar0 = m0Xstar + ipw0
    oracle_mu0 = mean(fstar0)
    fstar = fstar1 - fstar0
    #final step of IF
    oracle_mu = oracle_mu1 - oracle_mu0
    sd_oracle = sqrt(mean((fstar - oracle_mu)^2))/sqrt(length(Y))
  }else if(model_type == "Q-L"){
    a = as.vector(X%*%beta3) - 0.08
    piX = exp(a)/(1 + exp(a))
    m1Xstar = alpha1 + 2*X%*%beta1 + coeff*(X^2)%*%beta12 
    ipw1 = trt*(Y - m1Xstar)/piX
    fstar1 = m1Xstar + ipw1
    oracle_mu1 = mean(fstar1)
    # control arm part
    m0Xstar = alpha2 + X%*%beta2
    ipw0 = (1-trt)*(Y - m0Xstar)/(1-piX)
    fstar0 = m0Xstar + ipw0
    oracle_mu0 = mean(fstar0)
    fstar = fstar1 - fstar0
    #final step of IF
    oracle_mu = oracle_mu1 - oracle_mu0
    sd_oracle = sqrt(mean((fstar - oracle_mu)^2))/sqrt(length(Y))
  }
  return(list(oracle_mu = oracle_mu, sd_oracle = sd_oracle))
}

