# calculate the true mu1, mu0

source("posterior_for_mu/functions/gen_data_different_models.R") # generate data with different models

true_mu_func = function(n = NULL, p = NULL, s = NULL, alpha1 = NULL, alpha2 = NULL,
                       beta1 = NULL, beta12 = NULL,beta2 = NULL, beta3 = NULL, 
                       model_type = NULL, 
                                    coeff = NULL,
                                    mu = NULL, Sig = NULL, snr = NULL){
  
 DG = GD_different_models_func(n = n, p = p, s = s, alpha1 = alpha1, alpha2 = alpha2, 
                          beta1 = beta1,
                          beta12 = beta12, beta2 = beta2, beta3 = beta3, 
                          model_type = model_type, 
                          coeff = coeff, mu = mu, Sig = Sig, snr = snr) 
 mu0 = mean(DG$Y0)
 mu1 = mean(DG$Y1)
 mu_true = mu1 - mu0
 
 sd0 = sd(DG$Y0)
 sd1 = sd(DG$Y1)
 sdY = sd(DG$Y)
 
 DF <- data.frame(Y = DG$Y1, X = DG$X)
 
 # ---- Fitted linear regression (with intercept) ----
 model_with <- lm(Y ~ ., data=DF)
 betastar <- unname(model_with$coefficients)
 
 # ---- Fitted linear regression (WITHOUT intercept) ----
 model_without <- lm(Y ~ . - 1, data=DF)
 betastar_no_int <- unname(model_without$coefficients) 
 
 DG_additional = GD_different_models_func(n = n, p = p, s = s, alpha1 = alpha1, alpha2 = alpha2, 
                                          beta1 = beta1,
                                          beta12 = beta12, beta2 = beta2, beta3 = beta3, 
                                          model_type = model_type, 
                                          coeff = coeff, mu = mu, Sig = Sig, snr = snr) 
 X = DG_additional$X
 Xtilde = cbind(1, X)
 
 # calculate the targetted regresion function 
 mstar_X = Xtilde%*%betastar
 mstar_no_int_X = X%*%betastar_no_int
 
 # Compare means
 mean_mstar <- mean(mstar_X)
 mean_mstar_no_int <- mean(mstar_no_int_X)
 
 return(list(mu0 = mu0, mu1 = mu1, mu = mu_true,
             sd0 = sd0, sd1 = sd1, sdY = sdY,
             mean_mstar = mean_mstar, mean_mstar_no_int = mean_mstar_no_int))
 
}

