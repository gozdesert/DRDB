#posterior calculation for mu1: 
# post_mu1_samples.R
post_mu_samples = function(m1samples = NULL, m0samples = NULL, ps_samples = NULL, phat = NULL, Xtest = NULL, Ytest = NULL, trt_test = NULL, ind1_test = NULL, ind0_test = NULL, M = NULL){
  
  Xtest = as.matrix(Xtest)
  n1_test = length(ind1_test)
  n0_test = length(ind0_test)
  n = n0_test + n1_test
  # treated data
  Ytest1 = Ytest[ind1_test]
  X1 = as.matrix(Xtest[ind1_test, ]) # test data on T = 1
  X1 = cbind(1, X1) # test data on T = 1
  # control data
  Ytest0 = Ytest[ind0_test]
  X0 = as.matrix(Xtest[ind0_test, ]) # test data on T = 0
  X0 = cbind(1, X0) # test data on T = 0
  # whole data for mu
  Xtilde = cbind(1, Xtest) # all X values on the test data
  
  ##########################################################################
  #.............posterior for bias1 part: T=1 data.........................#
  ##########################################################################
  m1_X_t1 = X1 %*% m1samples  # n1xM matrix each column is m_1(X)
  #obtain lambda tilde
  X1beta = X1 %*% ps_samples # n1xM matrix each column is sample
  eX1beta = exp(X1beta) # n1xM matrix each column is sample
  piX1 = eX1beta/(1 + eX1beta) # n1xM matrix each column is sample
  lambdatilde_X_t1 = phat/piX1  # n1xM matrix each column is sample
  #obtaining samples for the first bias term
  df_b1 = n1_test - 1  # degrees of freedom
  f_X_Y_1 = lambdatilde_X_t1*(Ytest1 - m1_X_t1)  # n1xM matrix each column is sample
  eta_b1 = colMeans(f_X_Y_1) # M many center 
  s2_b1 = colMeans((f_X_Y_1 - eta_b1)^2)/(n1_test - 1) # M many scales 
  # t-distribution for bias1 
  b1_samples = rt(M, df = df_b1)*sqrt(s2_b1) + eta_b1
  remove(m1_X_t1,X1beta,eX1beta,piX1,lambdatilde_X_t1,df_b1,f_X_Y_1,eta_b1,s2_b1)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  ##########################################################################
  #.............posterior for bias2 part: T=0 data.........................#
  ##########################################################################
  m0_X_t0 = X0 %*% m0samples  # n0xM matrix each column is m_0(X)
  #obtain lambda tilde
  X0beta = X0 %*% ps_samples # n0xM matrix each column is sample
  eX0beta = exp(X0beta) # n0xM matrix each column is sample
  piX0 = eX0beta/(1 + eX0beta) # n0xM matrix each column is sample
  lambdatilde_X_t0 = (1-phat)/(1-piX0)  # n0xM matrix each column is sample
  #obtaining samples for the first bias term
  df_b0 = n0_test - 1  # degrees of freedom
  f_X_Y_0 = lambdatilde_X_t0*(Ytest0 - m0_X_t0)  # n0xM matrix each column is sample
  eta_b0 = colMeans(f_X_Y_0) # M many center 
  s2_b0 = colMeans((f_X_Y_0 - eta_b0)^2)/(n0_test - 1) # M many scales 
  # t-distribution for bias0 
  b0_samples = rt(M, df = df_b0)*sqrt(s2_b0) + eta_b0
  remove(m0_X_t0,X0beta,eX0beta,piX0,lambdatilde_X_t0,df_b0,f_X_Y_0,eta_b0,s2_b0)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  # obtain the final samples for the bias: b1-b0
  b_samples = b1_samples - b0_samples
  
  #### obtaining samples for mu1tilde
  m_samples = m1samples - m0samples
  mX = Xtilde %*% m_samples # n many m1(X)
  df_m = n - 1 
  mX_bar = colMeans(mX)
  eta_m = mX_bar  + b_samples
  s2_m = colMeans((mX - mX_bar)^2)/(n - 1)
  
  # t-distribution for mu1tilde
  mutilde_samples = rt(M, df = df_m)*sqrt(s2_m) + eta_m
  
  # mu1 samples M many and its mean and variance
  mu_samples = mutilde_samples
  
  return(list(b1_samples = b1_samples, b0_samples = b0_samples, mu_samples = mu_samples))
}
