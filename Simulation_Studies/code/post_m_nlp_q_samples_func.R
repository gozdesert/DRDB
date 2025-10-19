# post_m_nlp_q_samples_func.R

# updated the modelSelection to remove the message.
library(mombf)


post_m1_nlp_q_samples = function(Ytrain = NULL, Xtrain = NULL, 
                                    ind1_tr = NULL, M = NULL){
  # m1 estimation
  Y1 = Ytrain[ind1_tr]; X1 = as.matrix(Xtrain[ind1_tr, ])
  W = cbind(X1, X1^2)
  fit_m1 = modelSelection(y = Y1, x = W, center = TRUE, scale = TRUE, family = "normal", verbose = FALSE, priorCoef=momprior(tau=0.348), priorVar=igprior(.01,.01))
  out_m1 = as.matrix(rnlp(msfit = fit_m1, niter = M))
  m1_nlp_samples = out_m1[, -ncol(out_m1)]  # Mx(p+1) matrix
  remove(W, X1, Y1, fit_m1, out_m1)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  return(list(m1_nlp_samples = m1_nlp_samples)) 
}