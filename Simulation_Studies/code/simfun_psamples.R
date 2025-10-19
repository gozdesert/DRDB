# obtaining posterior samples for imp and mu_DR
# bart, barts and freq version of DR and imp are added
#simfunc_psamples.R
remove(list = ls())
source("posterior_for_mu/functions/gen_data_different_models.R")
source("posterior_for_mu/functions/true_mu.R")
source("posterior_for_mu/functions/oracle_estimator_mu.R")
source("posterior_for_mu/functions/post_m_bridge_samples.R")
source("posterior_for_mu/functions/post_m_ps_nlp_samples_func.R")
source("posterior_for_mu/functions/posterior_mu_summary.R")
source("posterior_for_mu/functions/post_mu_samples.R")
source("posterior_for_mu/functions/post_m1_bridge_samples_quad.R")
source("posterior_for_mu/functions/post_mu_samples_m1_q.R")
source("posterior_for_mu/functions/post_m_nlp_q_samples_func.R")
source("posterior_for_mu/functions/post_mu_samples_bart_func.R")
source("posterior_for_mu/functions/post_m_bart_samples_func.R")
source("posterior_for_mu/functions/post_m_barts_samples_func.R")
source("posterior_for_mu/functions/post_mu_oracle_samples_func.R")


# data generating mechanism for mu1   ####           
n = 1000; p = 200; s = 14; d = 2; #p2 for modeling propensity score
beta1 = c(rep(1,(s)/2), rep(0.5,(s)/2), rep(0, p-s));
beta2 = c(rep(1,(s)/2), rep(0.5,(s)/2), rep(0, p-s))
alpha1 = 5; alpha2 = 3
beta12 = c(0.5, 1, 1, rep(0,p-3))  #coeff for "Q-L", changed 060925
mu = rep(0, p)
Sig = diag(c(1, 1, rep(1, (p-2))))
# rho = 0.7
# Sig = matrix(0, nrow = p, ncol = p)
# for (i in 1:p) {
#   for (j in 1:p) {
#     if (i == j) { Sig[i, j] <- 1} else {Sig[i, j] <- rho^(abs(i-j))} #for correlated X
#   }}
#Omega = 2*Sig^2 # covariance  of W=(X1^2,X_2^2,...)
#diag(Omega) = 2
#mu = rep(0, p); Sig = diag(c(1, 1, rep(1, (p-2))))
Omega = matrix(rep(1, p^2), nrow = p)
diag(Omega) = 3
ratio = 3
# 4 comes from L part of m1X and 5 is the variance ratio of L and Q parts in m1X
coeff = as.numeric(sqrt(4*t(beta1)%*%Sig%*%beta1)/sqrt(ratio*t(beta12)%*%Omega%*%beta12)); 
beta3 = c(rep(0.35, d), rep(0, p-d)) 
model_type = "L-L" # "L-L"; "Q-L"
Kfoldset = c(5)
M = 1000;  snr = 5; nsample = 1000; niter = 20
K = length(Kfoldset)
Kf = as.numeric(paste0(Kfoldset, collapse = ""))
input = c(p, s, d, n, model_type, M, niter, Kf, ratio, coeff)

#calculate the true ATE mu #####
mu_true = true_mu_func(n = 500000, p = p, s = s, alpha1 = alpha1, alpha2 = alpha2,
                       beta1 = beta1, beta12 = beta12, beta2 = beta2, beta3 = beta3, 
                       model_type = model_type, coeff = coeff,
                       mu = mu, Sig = Sig, snr = snr)$mu

print("Finished step 1")


simfunc = function(i, mu_true){
  # generate the data
  DG = GD_different_models_func(n = n, p = p, s = s, alpha1 = alpha1, alpha2 = alpha2,
                                beta1 = beta1, beta12 = beta12, beta2 = beta2, beta3 = beta3, 
                                model_type = model_type, coeff = coeff,
                                mu = mu, Sig = Sig, snr = snr)
  Y = DG$Y; X = DG$X; Y1 = DG$Y1; Y0 = DG$Y0
  trt = DG$trt; ind1 = DG$indx1; ind0 = DG$indx0
  
  #############################################################################
  #.....................    DR Oracle Approach    ...........................#
  #############################################################################  
  out_mu_DR_oracle = post_mu_oracle_samples(Y = Y, X = X, trt = trt,
                                            ind1= ind1, ind0 = ind0,
                                            alpha1 = alpha1, alpha2 = alpha2, beta1 = beta1,
                                                       beta2 = beta2, beta3 = beta3,
                                                       model_type = NULL, coeff = NULL, M = M)
  mu_DR_oracle_psamples = out_mu_DR_oracle$mu_oracle
  
  #############################################################################
  #.....................Imputation Approach for NLP...........................#
  #############################################################################
  out_imp_nlp = post_nlp_samples(Ytrain = Y, Xtrain = X, trt_tr = trt, 
                                 ind1_tr = ind1, ind0_tr = ind0, M = M)
  m_nlp_samples = out_imp_nlp$m1_nlp_samples - out_imp_nlp$m0_nlp_samples # Mx(p+1) matrix
  mX_imp_nlp_samples = cbind(1, X)%*%t(m_nlp_samples) #nxM matrix
  mu_imp_nlp_est = colMeans(mX_imp_nlp_samples)
  c2_imp_nlp = colMeans((mX_imp_nlp_samples - mu_imp_nlp_est)^2)/(n-1)
  mu_imp_nlp_psamples = rt(n = M, df = n-1)*sqrt(c2_imp_nlp) + mu_imp_nlp_est
  #mu_imp_nlp = post_mu_summary(mu.samples = mu_imp_nlp_psamples, mu_true = mu_true)
  remove(m_nlp_samples, mX_imp_nlp_samples, mu_imp_nlp_est,
         c2_imp_nlp)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  #############################################################################
  #...............Imputation Approach for NLP with quad.......................#
  #############################################################################
  out_imp_nlp_quad = post_m1_nlp_q_samples(Ytrain = Y, Xtrain = X, ind1_tr = ind1, M = M)
  mX_imp_nlp_samples = cbind(1, X, X^2)%*%t(out_imp_nlp_quad$m1_nlp_samples) - cbind(1, X)%*%t(out_imp_nlp$m0_nlp_samples)
  mu_imp_nlp_est = colMeans(mX_imp_nlp_samples)
  c2_imp_nlp = colMeans((mX_imp_nlp_samples - mu_imp_nlp_est)^2)/(n-1)
  mu_imp_nlp_psamples = rt(n = M, df = n-1)*sqrt(c2_imp_nlp) + mu_imp_nlp_est
  #mu_imp_nlp_quad = post_mu_summary(mu.samples = mu_imp_nlp_psamples, mu_true = mu_true)
  remove(out_imp_nlp_quad, mX_imp_nlp_samples, mu_imp_nlp_est,
         c2_imp_nlp)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  #############################################################################
  #...................Imputation Approach for Bridge .........................#
  #############################################################################
  out_imp_br = post_m_bridge_samples_func(Ytrain = Y, Xtrain = X, 
                                          ind1_tr = ind1, ind0_tr = ind0, M = M)
  m_brig_samples = out_imp_br$m1.samples.ridge - out_imp_br$m0.samples.ridge # Mx(p+1) matrix
  mX_imp_brig_samples = cbind(1, X)%*%m_brig_samples #nxM matrix
  mu_imp_brig_est = colMeans(mX_imp_brig_samples)
  c2_imp_brig = colMeans((mX_imp_brig_samples - mu_imp_brig_est)^2)/(n-1)
  mu_imp_brig_psamples = rt(n = M, df = n-1)*sqrt(c2_imp_brig) + mu_imp_brig_est
  #mu_imp_brid = post_mu_summary(mu.samples = mu_imp_brig_psamples, mu_true = mu_true)
  remove(m_brig_samples, mX_imp_brig_samples, mu_imp_brig_est,
         c2_imp_brig)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  #############################################################################
  #............. Imputation Approach for Bridge Quad .........................#
  #############################################################################
  out_imp_br_q = post_m1_bridge_samples_quad_func(Ytrain = Y, Xtrain = X, ind1_tr = ind1, M = M)
  mX_imp_brig_samples = cbind(1, X, X^2)%*%out_imp_br_q$m1.samples.ridge_quad - cbind(1, X)%*%out_imp_br$m0.samples.ridge 
  mu_imp_brig_est = colMeans(mX_imp_brig_samples)
  c2_imp_brig = colMeans((mX_imp_brig_samples - mu_imp_brig_est)^2)/(n-1)
  mu_imp_brig_quad_psamples = rt(n = M, df = n-1)*sqrt(c2_imp_brig) + mu_imp_brig_est
  #mu_imp_brid_quad = post_mu_summary(mu.samples = mu_imp_brig_quad_psamples, mu_true = mu_true)
  remove(out_imp_br, out_imp_br_q, mX_imp_brig_samples, mu_imp_brig_est,
         c2_imp_brig)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  
  #############################################################################
  #................... Imputation Approach for BART ..........................#
  #############################################################################
  # out_imp_bart = post_m_bart_samples(Ytrain = Y, Xtrain = X, Xtest = X, 
  #                                    ind1_tr = ind1,ind0_tr = ind0, M = M)
  # mX_bart_samples = out_imp_bart$m1X_samples - out_imp_bart$m0X_samples
  # mu_imp_bart_est = colMeans(mX_bart_samples)
  # c2_imp_bart = colMeans((mX_bart_samples - mu_imp_bart_est)^2)/(n-1)
  # mu_imp_bart_psamples = rt(n = M, df = n-1)*sqrt(c2_imp_bart) + mu_imp_bart_est
  # #mu_imp_bart = post_mu_summary(mu.samples = mu_imp_bart_psamples, mu_true = mu_true)
  # remove(out_imp_bart, mX_bart_samples, mu_imp_bart_est, c2_imp_bart)
  # gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  # 
  # 
  #############################################################################
  #............. Imputation Approach for BART with sparse ....................#
  #############################################################################
  out_imp_barts = post_m_barts_samples(Ytrain = Y, Xtrain = X, Xtest = X, 
                                       ind1_tr = ind1,ind0_tr = ind0, M = M)
  mX_barts_samples = out_imp_barts$m1X_samples - out_imp_barts$m0X_samples
  mu_imp_barts_est = colMeans(mX_barts_samples)
  c2_imp_barts = colMeans((mX_barts_samples - mu_imp_barts_est)^2)/(n-1)
  mu_imp_barts_psamples = rt(n = M, df = n-1)*sqrt(c2_imp_barts) + mu_imp_barts_est
  #mu_imp_barts = post_mu_summary(mu.samples = mu_imp_barts_psamples, mu_true = mu_true)
  remove(out_imp_barts, mX_barts_samples, mu_imp_barts_est, c2_imp_barts)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  #############################################################################
  #                                                                           #
  # ...................Our Method with Cross-Fitting ........................ #
  #                                                                           #
  #############################################################################
  
 
    Kfold = as.numeric(Kfoldset)
    ind_fold = cut(1:n, Kfold, labels = FALSE) # K different folds
    ind = sample(ind_fold) # randomize data splitting
    
    mu_DR_nlp_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
    mu_DR_nlp_q_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
    mu_DR_bridge_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
    mu_DR_bridge_q_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
    mu_DR_bart_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
    mu_DR_barts_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
    
    for (k in 1:Kfold) {
      Xtest = as.matrix(X[ind == k, ]) # (n/Kfold)xp matrix
      Ytest = Y[ind == k] # n/Kfold vector
      trt_test = trt[ind == k]; ind1_test=which(trt_test == 1); ind0_test=which(trt_test == 0)
      #n_test = length(Ytest)
      # create training data
      Xtrain = as.matrix(X[ind != k, ]); Ytrain = Y[ind != k]; trt_tr = trt[ind != k]
      ind1_tr = which(trt_tr == 1); ind0_tr = which(trt_tr == 0)
      #n1_tr = length(ind1_tr); n0_tr = length(ind0_tr)
      
      ##########################################################################
      #..................post for mu_DR for NLP & NLP-Q .......................#
      ##########################################################################
      out_DR_nlp = post_nlp_samples(Ytrain = Ytrain, Xtrain = Xtrain, trt_tr = trt_tr,
                                    ind1_tr = ind1_tr, ind0_tr = ind0_tr, M = M)
      ## posterior samples for nuisance parameters based on nlp
      ps_nlp_samples = t(out_DR_nlp$ps_nlp_samples) #(p+1)xM matrix
      ps_nlp_sample = ps_nlp_samples[, sample(1:M, 1)] #(p+1)xM matrix
      phat = out_DR_nlp$phat
      m1_nlp_samples = t(out_DR_nlp$m1_nlp_samples) #(p+1)xM matrix
      m1_nlp_sample =  m1_nlp_samples[, sample(1:M, 1)]
      m0_nlp_samples = t(out_DR_nlp$m0_nlp_samples) #(p+1)xM matrix
      m0_nlp_sample =  m0_nlp_samples[, sample(1:M, 1)]
      out_DR_nlp_q = post_m1_nlp_q_samples(Ytrain = Ytrain, Xtrain = Xtrain,
                                           ind1_tr = ind1_tr, M = M)
      m1_nlp_q_samples = t(out_DR_nlp_q$m1_nlp_samples)
      m1_nlp_q_sample = m1_nlp_q_samples[, sample(1:M, 1)]
      
      ##########################################################################      
      ##-- posterior for mu DR--one sample
      out_mu_DR_nlp = post_mu_samples(m1samples = m1_nlp_sample, m0samples = m0_nlp_sample,
                                      ps_samples = ps_nlp_sample, phat = phat, Xtest = Xtest, 
                                      Ytest = Ytest, trt_test = trt_test, ind1_test = ind1_test,
                                      ind0_test = ind0_test, M = M)
      mu_DR_nlp_psamples_foldK[k, ] = out_mu_DR_nlp$mu_samples
      remove(out_mu_DR_nlp, m0_nlp_samples, m1_nlp_samples)
      gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
      ##########################################################################
      ##-- posterior for mu DR with nlp with quadratic m1
      out_mu_DR_nlp_q = post_mu_samples_m1_q_func(m1samples = m1_nlp_q_sample, m0samples = m0_nlp_sample, ps_samples = ps_nlp_sample, phat = phat, Xtest = Xtest, Ytest = Ytest, trt_test = trt_test, ind1_test = ind1_test, ind0_test = ind0_test, M = M)
      mu_DR_nlp_q_psamples_foldK[k, ] = out_mu_DR_nlp_q$mu_samples
      remove(out_mu_DR_nlp_q, m0_nlp_sample, m1_nlp_sample, m1_nlp_q_samples, m1_nlp_q_sample)
      gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
      
      
      ##########################################################################
      #.................posterior for mu_DR for Bridge.........................#
      ##########################################################################
      out_DR_br = post_m_bridge_samples_func(Ytrain = Ytrain, Xtrain = Xtrain, 
                                             ind1_tr = ind1_tr, ind0_tr = ind0_tr, M = M)
      m1_brig_samples = out_DR_br$m1.samples.ridge #(p+1)xM
      m1_brig_sample = m1_brig_samples[, sample(1:M, 1)] #(p+1) vector
      m0_brig_samples = out_DR_br$m0.samples.ridge #(p+1)xM
      m0_brig_sample = m0_brig_samples[, sample(1:M, 1)] #(p+1) vector
      remove(out_DR_br)
      gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
      
      out_DR_br_q = post_m1_bridge_samples_quad_func(Ytrain = Ytrain, Xtrain = Xtrain, 
                                                     ind1_tr = ind1_tr, M = M)
      m1_br_q_samples = out_DR_br_q$m1.samples.ridge_quad #(2p+1)xM
      m1_br_q_sample = m1_br_q_samples[, sample(1:M, 1)] #(2p+1) vector
      remove(out_DR_br_q)
      gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
      ##########################################################################      
      ##-- posterior for mu DR-one sample
      out_mu_DR_br = post_mu_samples(m1samples = m1_brig_sample, m0samples = m0_brig_sample,
                                     ps_samples = ps_nlp_sample, phat = phat, Xtest = Xtest, 
                                     Ytest = Ytest, trt_test = trt_test, ind1_test = ind1_test,
                                     ind0_test = ind0_test, M = M)
      mu_DR_bridge_psamples_foldK[k, ] = out_mu_DR_br$mu_samples
      remove(out_mu_DR_br, m0_brig_samples, m1_brig_samples)
      gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
      ##########################################################################
      ##-- posterior for mu DR-one sample with quadratic
      out_mu_DR_br_q = post_mu_samples_m1_q_func(m1samples = m1_br_q_sample, m0samples = m0_brig_sample,
                                                 ps_samples = ps_nlp_sample, phat = phat, Xtest = Xtest, 
                                                 Ytest = Ytest, trt_test = trt_test, ind1_test = ind1_test,
                                                 ind0_test = ind0_test, M = M)
      mu_DR_bridge_q_psamples_foldK[k, ] = out_mu_DR_br_q$mu_samples
      remove(out_mu_DR_br_q, m0_brig_sample, m1_brig_sample)
      gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
      
      ##########################################################################
      #.................posterior for mu_DR for BART.........................#
      ##########################################################################
      # out_DR_bart = post_m_bart_samples(Ytrain = Ytrain, Xtrain = Xtrain, 
      #                                   Xtest = Xtest, ind1_tr = ind1_tr,
      #                                   ind0_tr = ind0_tr, M = M)
      # m1X_samples_DR = out_DR_bart$m1X_samples #nkxM
      # m1X_sample_DR = m1X_samples_DR[, sample(1:M, 1)]
      # m0X_samples_DR = out_DR_bart$m0X_samples #nkxM
      # m0X_sample_DR = m0X_samples_DR[, sample(1:M, 1)]
      # 
      # out_mu_DR_bart = post_mu_samples_bart(m1X_samples = m1X_sample_DR, 
      #                                       m0X_samples = m0X_sample_DR, 
      #                                       Xtest = Xtest, Ytest = Ytest,
      #                                       ind1_test = ind1_test, 
      #                                       ind0_test = ind0_test, 
      #                                       ps_samples = ps_nlp_sample, 
      #                                       phat = phat, M = M)
      # mu_DR_bart_psamples_foldK[k, ] = out_mu_DR_bart$mu_samples
      # remove(out_DR_bart, out_mu_DR_bart, m1X_samples_DR, m1X_sample_DR, 
      #        m0X_samples_DR, m0X_sample_DR)
      # gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
      
      
      ##########################################################################
      #...........posterior for mu_DR for Bart with sparsity..................#
      ##########################################################################
      out_DR_barts = post_m_barts_samples(Ytrain = Ytrain, Xtrain = Xtrain, 
                                          Xtest = Xtest, ind1_tr = ind1_tr,
                                          ind0_tr = ind0_tr, M = M)
      m1Xs_samples_DR = out_DR_barts$m1X_samples #nkxM
      m1Xs_sample_DR = m1Xs_samples_DR[, sample(1:M, 1)]
      m0Xs_samples_DR = out_DR_barts$m0X_samples #nkxM
      m0Xs_sample_DR = m0Xs_samples_DR[, sample(1:M, 1)]
      
      out_mu_DR_barts = post_mu_samples_bart(m1X_samples = m1Xs_sample_DR, 
                                             m0X_samples = m0Xs_sample_DR, 
                                             Xtest = Xtest, Ytest = Ytest,
                                             ind1_test = ind1_test, 
                                             ind0_test = ind0_test, 
                                             ps_samples = ps_nlp_sample, 
                                             phat = phat, M = M)
      mu_DR_barts_psamples_foldK[k, ] = out_mu_DR_barts$mu_samples
      remove(out_DR_barts, out_mu_DR_barts, m1Xs_samples_DR, m1Xs_sample_DR, 
             m0Xs_samples_DR, m0Xs_sample_DR)
      gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
      
      
      
      
    }
    ###########################################################################
    #...........Aggregate posteriors obtained from each data fold.............#
    ###########################################################################
    mu_DR_nlp_psamples = colMeans(mu_DR_nlp_psamples_foldK)
    mu_DR_nlp_q_psamples = colMeans(mu_DR_nlp_q_psamples_foldK)
    mu_DR_brid_psamples = colMeans(mu_DR_bridge_psamples_foldK)
    mu_DR_brid_q_psamples = colMeans(mu_DR_bridge_q_psamples_foldK)
    mu_DR_bart_psamples = colMeans(mu_DR_bart_psamples_foldK)
    mu_DR_barts_psamples = colMeans(mu_DR_barts_psamples_foldK)
    
    remove(ind_fold, ind, Kfold, mu_DR_nlp_psamples_foldK, mu_DR_nlp_q_psamples_foldK,
           mu_DR_bridge_psamples_foldK, mu_DR_bridge_q_psamples_foldK, mu_DR_bart_psamples_foldK,
           mu_DR_barts_psamples_foldK)
    
    gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  print("Finished step 2")
  
  remove(Y, X, trt, Y1, Y0, ind1, ind0)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  print(paste0("Print i: ",i)) 
  
  list(mu_DR_oracle_psamples = mu_DR_oracle_psamples,
       mu_imp_nlp_psamples = mu_imp_nlp_psamples,
       mu_imp_nlp_psamples = mu_imp_nlp_psamples,
       mu_imp_brig_psamples = mu_imp_brig_psamples,
       mu_imp_brig_quad_psamples = mu_imp_brig_quad_psamples,
       #mu_imp_bart_psamples = mu_imp_bart_psamples,
       mu_imp_barts_psamples = mu_imp_barts_psamples,
       mu_DR_nlp_psamples = mu_DR_nlp_psamples,
       mu_DR_nlp_q_psamples = mu_DR_nlp_q_psamples,
       mu_DR_brid_psamples = mu_DR_brid_psamples,
       mu_DR_brid_q_psamples = mu_DR_brid_q_psamples,
       #mu_DR_bart_psamples = mu_DR_bart_psamples,
       mu_DR_barts_psamples =mu_DR_barts_psamples,
       mu_true = mu_true, input = input)
}


#############################################################################
#......................Run the code in parallel.............................#
#############################################################################

library(parallel)
RNGkind("L'Ecuyer-CMRG")
ncore = 1

timestamp()
out <- mclapply(1:niter, simfunc, mu_true = mu_true, mc.cores = ncore)
timestamp()

filename = paste0("posterior_for_mu/SimPlots/ATE_psamples_n",n, "p",p, "s",s, "d",d, "Kfold",Kf, "model_type_",model_type, "niter", niter, "alpha1_", alpha1, "alpha2_", alpha2, "coeff_", coeff,
                  format(Sys.time(), "_%Y-%m_%d_%I-%p"), ".RDS")
saveRDS(out, file = filename)

beepr::beep()