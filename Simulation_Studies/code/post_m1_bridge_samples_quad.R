#obtaining samples of beta using Bridge
# post_m1_bridge_samples_quad.R
# lambda.min is extracted to calculate bridge posterior for m1
library(glmnet)


post_m1_bridge_samples_quad_func = function(Ytrain = NULL, Xtrain = NULL,
                                      ind1_tr = NULL, M = NULL){
  
  ###########################################################################
  # posterior for m1 based on bridge
  ###########################################################################
  ## obtaining the treatment group
  Y1 = Ytrain[ind1_tr]; W = as.matrix(Xtrain[ind1_tr, ])
  n1 = nrow(W) 
  X1 = cbind(W, W^2)
  p = ncol(Xtrain) # note 2p many coefficients should be in X1
  # standardization of X1 
  X_centered = apply(X1, 2, function(x) x - mean(x)) # centered X
  X1s = apply(X_centered, 2, function(x) x / sqrt(sum(x^2) / n1)) # scaled X: Xs
  n1 = nrow(X1s)
  A = cbind(0, diag(2*p))
  Ahat = t(A)%*%A
  X1cv = X1s 
  # obtain lambda.ridge
  cv_model.ridge = cv.glmnet(x = X1cv, y = Y1, alpha = 0, standardize = TRUE, thresh=1e-20)
  #find optimal lambda value that minimizes test MSE
  lambda.ridge = cv_model.ridge$lambda.min
  
  Xtilde = cbind(1, X1s)
  Xhat = t(Xtilde)%*%Xtilde # (2p+1)x(2p+1) matrix
  # scale of Y
  Y1_sd = sqrt(((n1 - 1)/n1))*sd(Y1) # we need this edit to use lambda.ridge
  #obtain lambda.hat from glmnet
  lambda.hat = (n1*lambda.ridge)/Y1_sd #lambda is for original scale Xn1!
  #### Ridge regression sampling
  #marginal posterior of sigma 
  #obtain K samples for sig2 from InvGam dist 
  a.ridge = (n1 - 1)/2 #scale parameter
  Q.ridge = Xhat + lambda.hat*Ahat  #(2p+1)x(2p+1) matrix
  Q.inv = solve(Q.ridge) #(2p+1)x(2p+1) matrix
  PX.ridge = Xtilde%*%Q.inv%*%t(Xtilde)
  r.ridge = (t(Y1)%*%(diag(n1) - PX.ridge)%*%Y1)/2
  sig2.ridge = 1/rgamma(n = M, shape = a.ridge, scale = 1/r.ridge)
  #conditional posterior of gamma give sigma2
  L.ridge = t(chol(Q.ridge)) #(2p+1)x(2p+1) matrix
  yr.ridge = solve(t(L.ridge),t(mvrnorm(n = M, mu = rep(0, 2*p+1), Sigma = diag(2*p+1))))
  B.ridge = t(sqrt(1/sig2.ridge)%*% t(as.vector(t(Xtilde)%*%Y1)))
  Vr.ridge = solve(L.ridge, B.ridge)            # (2p+1)xM 
  Thetar.ridge = solve(t(L.ridge), Vr.ridge)    # (2p+1)xM 
  #obtain one sample of gamma by choosing a random column from Pm1.samples.bols since its colums are samples we have M many of them
  Pm1.samples.ridge = as.matrix(t(t(Thetar.ridge + yr.ridge)*sqrt(sig2.ridge))) # (2p+1)xM
  X_mean <- colMeans(X1)
  X_sd <- apply(X_centered, 2, function(x) sqrt(sum(x^2) / n1))
  mu_sd = X_mean/X_sd
  
  Pm1.samples.ridge_betas = Pm1.samples.ridge[-1 ,]/X_sd
  Pm1.samples.ridge_beta0 = Pm1.samples.ridge[1 ,] - colMeans(Pm1.samples.ridge[-1 ,]*mu_sd)*(2*p) # it should be sum not colmeans!
  Pm1.samples.ridge_original = matrix(rbind(Pm1.samples.ridge_beta0, Pm1.samples.ridge_betas), ncol = M)
  
  rm(Y1, W, X1, X_centered, X1s, n1, X1cv, cv_model.ridge, lambda.ridge, Xtilde, Xhat, Y1_sd, 
     a.ridge, Q.ridge, Q.inv, PX.ridge, r.ridge, sig2.ridge, L.ridge, yr.ridge, B.ridge, 
     Vr.ridge, Thetar.ridge, Pm1.samples.ridge, X_mean, X_sd, mu_sd, Pm1.samples.ridge_betas,
     Pm1.samples.ridge_beta0, lambda.hat)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  return(list(m1.samples.ridge_quad = Pm1.samples.ridge_original))
}