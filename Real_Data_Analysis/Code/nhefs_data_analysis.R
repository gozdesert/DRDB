remove(list = ls())
# Load required functions
source("posterior_for_mu/Real_Data/NHEFS/nhefs_functions.R")
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
source("posterior_for_mu/functions/post_m_barts_samples_func.R")
source("posterior_for_mu/functions/post_mu_samples_bart_func.R")
source("posterior_for_mu/functions/mu_freq_lasso_ridge_res.R")
source("posterior_for_mu/Real_Data/post_summary_realdata_func.R")

# Load the dataset
nhefs <- read.csv("posterior_for_mu/Real_Data/NHEFS/nhefs.csv")
#ind_miss=which(is.na(nhefs$wt82_71)|(nhefs$alcoholpy==2)) #added this line to remove missing rows
# Apply proper factor conversion first
nhefs_clean <- prepare_nhefs_data(nhefs)
# Check the structure
str(nhefs_clean[c("sex", "race", "education", "asthma", "exercise")])

# Create different analysis datasets
# Y = wt82 ; trt = qsmk
cov_set = c("sex", "race", "age", "education", "smokeintensity", "smokeyrs", "exercise", "active", "wt71", 
            "sbp")
  
# c("sex", "race", "age", "education","smokeintensity", "smokeyrs", 
#             "exercise", "active", "wt71", "asthma", "cholesterol",
#             "dbp","ht","price71_82","sbp","tax71_82", "allergies")
  
#c("sex", "race", "age", "education","smokeintensity","smokeyrs", "exercise", "active", "wt71")

# Check if this set of covariates satisfies positivity condition 
issues <- quick_positivity_check(nhefs_clean[c(cov_set, "qsmk")])

# Create model matrices for given set of covariates
#create_model_matrices <- function(data, covariate_set, treatment_var, outcome_var)
gen_data <- create_model_matrices(data = nhefs_clean, covariate_set = cov_set, 
                                  treatment_var = "qsmk", outcome_var = "wt82_71")
X = gen_data$X_no_int
dim(X) 
Y = gen_data$Y
trt = gen_data$A
ind1 = which(trt == 1)
ind0 = which(trt == 0)
M = 1000 # number of posterior samples
n = length(Y)
Kfoldset = c(5); K = length(Kfoldset)
mu_true = mean(Y)

#############################################################################
#........... ATE estimation based on naive estimator  ...............#
#############################################################################
# Calculate means and sample sizes
mean_treated <- mean(Y[trt == 1])
mean_control <- mean(Y[trt == 0])
n_treated <- sum(trt == 1)
n_control <- sum(trt == 0)
# Calculate variances
var_treated <- var(Y[trt == 1])
var_control <- var(Y[trt == 0])
# Empirical difference
empdiff <- mean_treated - mean_control
# Standard error
se_empdiff <- sqrt(var_treated / n_treated + var_control / n_control)
# 95% confidence interval
alpha <- 0.05
z <- qnorm(1 - alpha/2)
ci_lower <- empdiff - z * se_empdiff
ci_upper <- empdiff + z * se_empdiff
mu_naive = cbind(empdiff, se_empdiff^2, ci_lower, ci_upper, ci_upper-ci_lower)

#############################################################################
#........... mu freq ncf estimators based on lasso and ridge ...............#
#############################################################################
out_freq = mu_freq_lasso_ridge_func(Y=Y, X=X, trt=trt, ind1=ind1, ind0=ind0)
mu_DR_ncf_lasso = out_freq$mu_DR_lasso
mu_DR_ncf_lasso[5] = mu_DR_ncf_lasso[4] - mu_DR_ncf_lasso[3] 
mu_DR_ncf_ridge = out_freq$mu_DR_ridge
mu_DR_ncf_ridge[5] = mu_DR_ncf_ridge[4] - mu_DR_ncf_ridge[3]
mu_imp_lasso = out_freq$mu_imp_lasso
mu_imp_lasso[5] = mu_imp_lasso[4] - mu_imp_lasso[3]
mu_imp_ridge = out_freq$mu_imp_ridge
mu_imp_ridge[5] = mu_imp_ridge[4] - mu_imp_ridge[3]
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
mu_imp_nlp = post_summary_RD_fn(mu_samples = mu_imp_nlp_psamples)$psummary
remove(m_nlp_samples, mX_imp_nlp_samples, mu_imp_nlp_est,
       c2_imp_nlp, mu_imp_nlp_psamples)

#############################################################################
#...............Imputation Approach for NLP with quad.......................#
#############################################################################
out_imp_nlp_quad = post_m1_nlp_q_samples(Ytrain = Y, Xtrain = X, ind1_tr = ind1, M = M)
mX_imp_nlp_samples = cbind(1, X, X^2)%*%t(out_imp_nlp_quad$m1_nlp_samples) - cbind(1, X)%*%t(out_imp_nlp$m0_nlp_samples)
mu_imp_nlp_est = colMeans(mX_imp_nlp_samples)
c2_imp_nlp = colMeans((mX_imp_nlp_samples - mu_imp_nlp_est)^2)/(n-1)
mu_imp_nlp_psamples = rt(n = M, df = n-1)*sqrt(c2_imp_nlp) + mu_imp_nlp_est
mu_imp_nlp_q = post_summary_RD_fn(mu_samples = mu_imp_nlp_psamples)$psummary
remove(out_imp_nlp_quad, mX_imp_nlp_samples, mu_imp_nlp_est,
       c2_imp_nlp, mu_imp_nlp_psamples)

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
mu_imp_brid = post_summary_RD_fn(mu_samples = mu_imp_brig_psamples)$psummary
remove(m_brig_samples, mX_imp_brig_samples, mu_imp_brig_est,
       c2_imp_brig, mu_imp_brig_psamples)
#gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)

#############################################################################
#............. Imputation Approach for Bridge Quad .........................#
#############################################################################
out_imp_br_q = post_m1_bridge_samples_quad_func(Ytrain = Y, Xtrain = X, ind1_tr = ind1, M = M)
mX_imp_brig_samples = cbind(1, X, X^2)%*%out_imp_br_q$m1.samples.ridge_quad - cbind(1, X)%*%out_imp_br$m0.samples.ridge 
mu_imp_brig_est = colMeans(mX_imp_brig_samples)
c2_imp_brig = colMeans((mX_imp_brig_samples - mu_imp_brig_est)^2)/(n-1)
mu_imp_brig_psamples = rt(n = M, df = n-1)*sqrt(c2_imp_brig) + mu_imp_brig_est
mu_imp_brid_q = post_summary_RD_fn(mu_samples = mu_imp_brig_psamples)$psummary
remove(out_imp_br, out_imp_br_q, mX_imp_brig_samples, mu_imp_brig_est,
       c2_imp_brig, mu_imp_brig_psamples)

#############################################################################
#............. Imputation Approach for BART with sparse ....................#
#############################################################################
out_imp_barts = post_m_barts_samples(Ytrain = Y, Xtrain = X, Xtest = X, 
                                     ind1_tr = ind1,ind0_tr = ind0, M = M)
mX_barts_samples = out_imp_barts$m1X_samples - out_imp_barts$m0X_samples
mu_imp_barts_est = colMeans(mX_barts_samples)
c2_imp_barts = colMeans((mX_barts_samples - mu_imp_barts_est)^2)/(n-1)
mu_imp_barts_psamples = rt(n = M, df = n-1)*sqrt(c2_imp_barts) + mu_imp_barts_est
mu_imp_bart = post_summary_RD_fn(mu_samples = mu_imp_barts_psamples)$psummary
remove(out_imp_barts, mX_barts_samples, mu_imp_barts_est, c2_imp_barts, mu_imp_barts_psamples) 

#############################################################################
#                                                                           #
# ...................Our Method with Cross-Fitting ........................ #
#                                                                           #
#############################################################################

Kfold = Kfoldset[1]
ind_fold = cut(1:n, Kfold, labels = FALSE) # K different folds
ind = sample(ind_fold) # randomize data splitting

mu_DR_nlp_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
mu_DR_nlp_q_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
mu_DR_bridge_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
mu_DR_bridge_q_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)
mu_DR_bart_psamples_foldK = matrix(rep(0, M*Kfold), nrow = Kfold)

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
  #gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  ##########################################################################
  ##-- posterior for mu DR with nlp with quadratic m1
  out_mu_DR_nlp_q = post_mu_samples_m1_q_func(m1samples = m1_nlp_q_sample, m0samples = m0_nlp_sample, ps_samples = ps_nlp_sample, phat = phat, Xtest = Xtest, Ytest = Ytest, trt_test = trt_test, ind1_test = ind1_test, ind0_test = ind0_test, M = M)
  mu_DR_nlp_q_psamples_foldK[k, ] = out_mu_DR_nlp_q$mu_samples
  remove(out_mu_DR_nlp_q, m0_nlp_sample, m1_nlp_sample, m1_nlp_q_samples, m1_nlp_q_sample)
  #gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  
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
  #gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
  out_DR_br_q = post_m1_bridge_samples_quad_func(Ytrain = Ytrain, Xtrain = Xtrain, 
                                                 ind1_tr = ind1_tr, M = M)
  m1_br_q_samples = out_DR_br_q$m1.samples.ridge_quad #(2p+1)xM
  m1_br_q_sample = m1_br_q_samples[, sample(1:M, 1)] #(2p+1) vector
  remove(out_DR_br_q)
  #gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  ##########################################################################      
  ##-- posterior for mu DR-one sample
  out_mu_DR_br = post_mu_samples(m1samples = m1_brig_sample, m0samples = m0_brig_sample,
                                 ps_samples = ps_nlp_sample, phat = phat, Xtest = Xtest, 
                                 Ytest = Ytest, trt_test = trt_test, ind1_test = ind1_test,
                                 ind0_test = ind0_test, M = M)
  mu_DR_bridge_psamples_foldK[k, ] = out_mu_DR_br$mu_samples
  remove(out_mu_DR_br, m0_brig_samples, m1_brig_samples)
  #gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  ##########################################################################
  ##-- posterior for mu DR-one sample with quadratic
  out_mu_DR_br_q = post_mu_samples_m1_q_func(m1samples = m1_br_q_sample, m0samples = m0_brig_sample,
                                             ps_samples = ps_nlp_sample, phat = phat, Xtest = Xtest, 
                                             Ytest = Ytest, trt_test = trt_test, ind1_test = ind1_test,
                                             ind0_test = ind0_test, M = M)
  mu_DR_bridge_q_psamples_foldK[k, ] = out_mu_DR_br_q$mu_samples
  remove(out_mu_DR_br_q, m0_brig_sample, m1_brig_sample)
  #gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
  
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
  mu_DR_bart_psamples_foldK[k, ] = out_mu_DR_barts$mu_samples
  remove(out_DR_barts, out_mu_DR_barts, m1Xs_samples_DR, m1Xs_sample_DR, 
         m0Xs_samples_DR, m0Xs_sample_DR)
  
}


###########################################################################
#...........Aggregate posteriors obtained from each data fold.............#
###########################################################################
mu_DR_nlp_psamples_agg = colMeans(mu_DR_nlp_psamples_foldK)
mu_DR_nlp_q_psamples_agg = colMeans(mu_DR_nlp_q_psamples_foldK)
mu_DR_brid_psamples_agg = colMeans(mu_DR_bridge_psamples_foldK)
mu_DR_brid_q_psamples_agg = colMeans(mu_DR_bridge_q_psamples_foldK)
mu_DR_bart_psamples_agg = colMeans(mu_DR_bart_psamples_foldK)

mu_DR_nlp = post_summary_RD_fn(mu_samples = mu_DR_nlp_psamples_agg)$psummary
mu_DR_nlp_q = post_summary_RD_fn(mu_samples = mu_DR_nlp_q_psamples_agg)$psummary
mu_DR_brid = post_summary_RD_fn(mu_samples = mu_DR_brid_psamples_agg)$psummary  
mu_DR_brid_q = post_summary_RD_fn(mu_samples = mu_DR_brid_q_psamples_agg)$psummary 
mu_DR_bart = post_summary_RD_fn(mu_samples = mu_DR_bart_psamples_agg)$psummary  

remove(ind_fold, ind, Kfold, mu_DR_nlp_psamples_agg, mu_DR_nlp_q_psamples_agg,
       mu_DR_brid_psamples_agg, mu_DR_brid_q_psamples_agg)
#gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)
remove(Y, trt, ind1, ind0)
#gc(verbose = getOption("verbose"), reset = FALSE, full = TRUE)

#####################
summary_data = rbind(mu_naive,
                     mu_imp_nlp,
                     mu_imp_nlp_q,
                     mu_imp_brid,
                     mu_imp_brid_q,
                     mu_imp_lasso,
                     mu_imp_ridge,
                     mu_imp_bart,
                     mu_DR_ncf_lasso,
                     mu_DR_ncf_ridge,
                     mu_DR_nlp, # last component is for Kfold
                     mu_DR_nlp_q,
                     mu_DR_brid,
                     mu_DR_brid_q,
                     mu_DR_bart)
colnames(summary_data) = c("pmean", "pvar", "2.5%", "97.5%", "CILeng")
rownames(summary_data)[1] = c("mu_naive")
beepr::beep()
dim(X); cov_set
cap = c('ATE for NHEFS with cov:', cov_set, "dim of X:", dim(X))

knitr::kable(summary_data, format = "pipe",
             caption = cat(cap), align = "c",
             booktabs = TRUE, valign = 't', linesep = c("", "", "", "", "", "", "\\addlinespace"))
