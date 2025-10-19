
# post_m_ps_nlp_samples.R

# updated the modelSelection to remove the message.
library(mombf)

post_nlp_samples = function(Ytrain = NULL, Xtrain = NULL, trt_tr = NULL, 
                            ind1_tr = NULL, ind0_tr = NULL, M = NULL){
  # PS and phat estimation 
  fit_pi = modelSelection(y = trt_tr, x = Xtrain, center = TRUE, scale = TRUE, family = "binomial", verbose = FALSE, priorCoef=momprior(tau=1/3), priorVar=igprior(.01,.01))
  ps_nlp_samples = as.matrix(rnlp(msfit = fit_pi, niter = M))
  phat = sum(trt_tr)/length(trt_tr)
  
  # m1 estimation
  Y1 = Ytrain[ind1_tr]; X1 = as.matrix(Xtrain[ind1_tr, ])
  fit_m1 = modelSelection(y = Y1, x = X1, center = TRUE, scale = TRUE, family = "normal", verbose = FALSE, priorCoef=momprior(tau=0.348), priorVar=igprior(.01,.01))
  out_m1 = as.matrix(rnlp(msfit = fit_m1, niter = M))
  m1_nlp_samples = out_m1[, -ncol(out_m1)]  # Mx(p+1) matrix
  remove(Y1, fit_m1, out_m1, fit_pi)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  # m0 estimation
  Y0 = Ytrain[ind0_tr]; X0 = as.matrix(Xtrain[ind0_tr, ])
  fit_m0 = modelSelection(y = Y0, x = X0, center = TRUE, scale = TRUE, family = "normal", verbose = FALSE, priorCoef=momprior(tau=0.348), priorVar=igprior(.01,.01))
  out_m0 = as.matrix(rnlp(msfit = fit_m0, niter = M))
  m0_nlp_samples = out_m0[, -ncol(out_m0)]  # Mx(p+1) matrix
  remove(Y0, fit_m0, out_m0)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  return(list(ps_nlp_samples = ps_nlp_samples, phat = phat, m1_nlp_samples = m1_nlp_samples,
              m0_nlp_samples = m0_nlp_samples)) 
}