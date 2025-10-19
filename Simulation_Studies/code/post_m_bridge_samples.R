#obtaining samples of beta using Bridge
# post_m_bridge_samples.R
# lambda.min is extracted to calculate bridge posterior for m1
library(glmnet)


post_m_bridge_samples_func = function(Ytrain = NULL, Xtrain = NULL,
                                    ind1_tr = NULL, ind0_tr = NULL, M = NULL){
  p = ncol(Xtrain)
  ###########################################################################
  # posterior for m1 based on bridge
  ###########################################################################
  ## obtaining the treatment group
  Y1 = Ytrain[ind1_tr]; X1 = as.matrix(Xtrain[ind1_tr, ])
  n1 = nrow(X1) 
  # standardization of X1 
  X_centered = apply(X1, 2, function(x) x - mean(x)) # centered X
  X1s = apply(X_centered, 2, function(x) x / sqrt(sum(x^2) / n1)) # scaled X: Xs
  n1 = nrow(X1s)
  A = cbind(0, diag(p))
  Ahat = t(A)%*%A
  X1cv = X1s 
  # obtain lambda.ridge
  cv_model.ridge = cv.glmnet(x = X1cv, y = Y1, alpha = 0, standardize = TRUE, thresh=1e-20)
  #find optimal lambda value that minimizes test MSE
  lambda.ridge = cv_model.ridge$lambda.min
      
  Xtilde = cbind(1, X1s)
  Xhat = t(Xtilde)%*%Xtilde # (p+1)x(p+1) matrix
  # scale of Y
  Y1_sd = sqrt(((n1 - 1)/n1))*sd(Y1) # we need this edit to use lambda.ridge
  #obtain lambda.hat from glmnet
  lambda.hat = (n1*lambda.ridge)/Y1_sd 
  #### Ridge regression sampling
  #marginal posterior of sigma 
  #obtain K samples for sig2 from InvGam dist 
  a.ridge = (n1 - 1)/2 #scale parameter
  Q.ridge = Xhat + lambda.hat*Ahat  #(p+1)x(p+1) matrix
  Q.inv = solve(Q.ridge) #(p+1)x(p+1) matrix
  PX.ridge = Xtilde%*%Q.inv%*%t(Xtilde)
  r.ridge = (t(Y1)%*%(diag(n1) - PX.ridge)%*%Y1)/2
  sig2.ridge = 1/rgamma(n = M, shape = a.ridge, scale = 1/r.ridge)
  #conditional posterior of gamma given sigma2
  L.ridge = t(chol(Q.ridge)) #(p+1)x(p+1) matrix
  yr.ridge = solve(t(L.ridge),t(mvrnorm(n = M, mu = rep(0, p+1), Sigma = diag(p+1))))
  B.ridge = t(sqrt(1/sig2.ridge)%*% t(as.vector(t(Xtilde)%*%Y1)))
  Vr.ridge = solve(L.ridge, B.ridge)            # (p+1)xM 
  Thetar.ridge = solve(t(L.ridge), Vr.ridge)    # (p+1)xM 
  #obtain one sample of gamma by choosing a random column from Pm1.samples.bols since its columns are samples we have M many of them
  Pm1.samples.ridge = as.matrix(t(t(Thetar.ridge + yr.ridge)*sqrt(sig2.ridge))) # (p+1)xM
  X_mean <- colMeans(X1)
  X_sd <- apply(X_centered, 2, function(x) sqrt(sum(x^2) / n1))
  mu_sd = X_mean/X_sd
  
  Pm1.samples.ridge_betas = Pm1.samples.ridge[-1 ,]/X_sd
  Pm1.samples.ridge_beta0 = Pm1.samples.ridge[1 ,] - colMeans(Pm1.samples.ridge[-1 ,]*mu_sd)*p # it should be sum not colmeans!
  Pm1.samples.ridge_original = matrix(rbind(Pm1.samples.ridge_beta0, Pm1.samples.ridge_betas), ncol = M)
  
  rm(Y1, X1, X_centered, X1s, n1, X1cv, cv_model.ridge, lambda.ridge, Xtilde, Xhat, Y1_sd, 
     a.ridge, Q.ridge, Q.inv, PX.ridge, r.ridge, sig2.ridge, L.ridge, yr.ridge, B.ridge, 
     Vr.ridge, Thetar.ridge, Pm1.samples.ridge, X_mean, X_sd, mu_sd, Pm1.samples.ridge_betas,
     Pm1.samples.ridge_beta0, lambda.hat)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  
  
  ###########################################################################
  # posterior for m0 based on bridge
  ###########################################################################
  ## obtaining the control group
  Y0 = Ytrain[ind0_tr]; X0 = as.matrix(Xtrain[ind0_tr, ])
  n0 = nrow(X0)
  # standardization of X1 
  X0_centered = apply(X0, 2, function(x) x - mean(x)) # centered X
  X0s = apply(X0_centered, 2, function(x) x / sqrt(sum(x^2) / n0)) # scaled X: Xs
  n0 = nrow(X0s)
  # obtain lambda.ridge
  cv_model.ridge0 = cv.glmnet(x = X0s, y = Y0, alpha = 0, standardize = TRUE, thresh=1e-20)
  #find optimal lambda value that minimizes test MSE
  lambda.ridge0 = cv_model.ridge0$lambda.min
  
  X0tilde = cbind(1, X0s)
  X0hat = t(X0tilde)%*%X0tilde # (p+1)x(p+1) matrix
  # scale of Y
  Y0_sd = sqrt(((n0 - 1)/n0))*sd(Y0) # we need this edit to use lambda.ridge
  #obtain lambda.hat from glmnet
  lambda.hat0 = (n0*lambda.ridge0)/Y0_sd #lambda is for original scale Xn1!
  #### Ridge regression sampling
  a.ridge0 = (n0 - 1)/2 #scale parameter
  Q.ridge0 = X0hat + lambda.hat0*Ahat  #(p+1)x(p+1) matrix
  Q.inv0 = solve(Q.ridge0) #(p+1)x(p+1) matrix
  PX.ridge0 = X0tilde%*%Q.inv0%*%t(X0tilde)
  r.ridge0 = (t(Y0)%*%(diag(n0) - PX.ridge0)%*%Y0)/2
  sig2.ridge0 = 1/rgamma(n = M, shape = a.ridge0, scale = 1/r.ridge0)
  #conditional posterior of gamma given sigma2
  L.ridge0 = t(chol(Q.ridge0)) #(p+1)x(p+1) matrix
  yr.ridge0 = solve(t(L.ridge0),t(mvrnorm(n = M, mu = rep(0, p+1), Sigma = diag(p+1))))
  B.ridge0 = t(sqrt(1/sig2.ridge0)%*% t(as.vector(t(X0tilde)%*%Y0)))
  Vr.ridge0 = solve(L.ridge0, B.ridge0)            # (p+1)xK 
  Thetar.ridge0 = solve(t(L.ridge0), Vr.ridge0)    # (p+1)xK 
  #obtain one sample of gamma by choosing a random column from Pm0.samples.bols
  Pm0.samples.ridge = as.matrix(t(t(Thetar.ridge0 + yr.ridge0)*sqrt(sig2.ridge0))) # (p+1)xM
  X0_mean <- colMeans(X0)
  X0_sd <- apply(X0_centered, 2, function(x) sqrt(sum(x^2) / n0))
  mu_sd0 = X0_mean/X0_sd
  
  Pm0.samples.ridge_betas = Pm0.samples.ridge[-1 ,]/X0_sd
  Pm0.samples.ridge_beta0 = Pm0.samples.ridge[1 ,] - colMeans(Pm0.samples.ridge[-1 ,]*mu_sd0)*p
  
  Pm0.samples.ridge_original = matrix(rbind(Pm0.samples.ridge_beta0, Pm0.samples.ridge_betas), ncol = M)
  
  rm(Y0, X0, X0_centered, X0s, n0, cv_model.ridge0, lambda.ridge0, X0tilde, X0hat, Y0_sd, 
     a.ridge0, Q.ridge0, Q.inv0, PX.ridge0, r.ridge0, sig2.ridge0, L.ridge0, yr.ridge0, 
     B.ridge0, Vr.ridge0, Thetar.ridge0, Pm0.samples.ridge, X0_mean, X0_sd, mu_sd0, 
     Pm0.samples.ridge_betas, Pm0.samples.ridge_beta0, lambda.hat0)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  return(list(m1.samples.ridge = Pm1.samples.ridge_original,
              m0.samples.ridge = Pm0.samples.ridge_original))
}

