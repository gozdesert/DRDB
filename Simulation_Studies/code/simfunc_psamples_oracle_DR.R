# obtaining posterior samples only for oracle and mu_DR
# simfunc_psamples_oracle_DR
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
n = 1000; p = 50; s = 13; d = 2; #p2 for modeling propensity score
# beta1 = c(rep(1,(s)/2), rep(0.5,(s)/2), rep(0, p-s));
beta1 = c(rep(1,(s+1)/2), rep(0.5,(s-1)/2), rep(0, p-s))   
beta2 = c(rep(1,(s+1)/2), rep(0.5,(s-1)/2), rep(0, p-s))  # s = odd
alpha1 = 5; alpha2 = 3
beta12 = c(0.5, 1, 1, rep(0,p-3))  #coeff for "Q-L", changed 060925
mu = rep(0, p)
Sig = diag(c(1, 1, rep(1, (p-2))))
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
  out_orc = oracle_est_mu_func(Y = Y, X = X, trt = trt, alpha1 = alpha1, alpha2 = alpha2, 
                               beta1 = beta1, beta2 = beta2, beta12 =beta12, beta3 = beta3,
                               model_type = model_type, coeff = coeff)
  pmean_oracle = out_orc$oracle_mu
  psd_oracle = out_orc$sd_oracle
  mu_DR_oracle_psamples = rnorm(n = M, mean = pmean_oracle, sd = psd_oracle)
  remove(DG, out_orc)
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
  mu_DR_bridge_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
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
  mu_DR_brid_psamples = colMeans(mu_DR_bridge_psamples_foldK)
  mu_DR_barts_psamples = colMeans(mu_DR_barts_psamples_foldK)
  
  print("Finished step 2")
  
  remove(Y, X, trt, Y1, Y0, ind1, ind0)
  gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  print(paste0("Print i: ",i)) 
  
  list(mu_DR_oracle_psamples = mu_DR_oracle_psamples,
       mu_DR_nlp_psamples = mu_DR_nlp_psamples,
       mu_DR_brid_psamples = mu_DR_brid_psamples,
       mu_DR_barts_psamples =mu_DR_barts_psamples,
       mu_true = mu_true, input = input)
}


#############################################################################
#......................Run the code in parallel.............................#
#############################################################################

library(parallel)
RNGkind("L'Ecuyer-CMRG")
ncore = 1

t1 = timestamp()
out <- mclapply(1:niter, simfunc, mu_true = mu_true, mc.cores = ncore)
t2 = timestamp()

filename = paste0("posterior_for_mu/SimPlots/ATE_psamples_oracle_normal_n",n, "p",p, "s",s, "d",d, "Kfold",Kf, "model_type_",model_type, "niter", niter, "alpha1_", alpha1, "alpha2_", alpha2, "coeff_", coeff,
                  format(Sys.time(), "_%Y-%m_%d_%I-%p"), ".RDS")
saveRDS(out, file = filename)

beepr::beep()