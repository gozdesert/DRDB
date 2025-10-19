# post_m_barts_samples_func.R
# lambda.min is extracted to calculate bridge posterior for m1
library(BART)

post_m_barts_samples = function(Ytrain = NULL, Xtrain = NULL, Xtest = NULL,
                               ind1_tr = NULL, ind0_tr = NULL, M = NULL){
  # obtain m1X and m0X based on bart
  X1t = Xtrain[ind1_tr, ]
  Y1t = Ytrain[ind1_tr]
  m1_fit = wbart(x.train = X1t, y.train = Y1t, x.test = Xtest, nskip=1000, ndpost=M, ntree = 200, printevery=1000, sparse = TRUE)
  m1X_samples = t(m1_fit$yhat.test) # nkxM
  remove(X1t, Y1t, m1_fit)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  X0t = Xtrain[ind0_tr, ]
  Y0t = Ytrain[ind0_tr]
  m0_fit = wbart(x.train = X0t, y.train = Y0t, x.test = Xtest, nskip=1000, ndpost=M, ntree = 200, printevery=1000, sparse = TRUE)
  m0X_samples = t(m0_fit$yhat.test) # nkxM
  remove(X0t, Y0t, m0_fit, Xtest)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  return(list(m1X_samples = m1X_samples, m0X_samples = m0X_samples))
}