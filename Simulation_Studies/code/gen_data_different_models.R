

library(MASS)
# GD_different_models.R
# alpha1 added to Q-L model 02/26/2025
GD_different_models_func = function(n = NULL, p = NULL, s = NULL, alpha1 = NULL, alpha2 = NULL,
                                    beta1 = NULL, beta12 = NULL,
                                    beta2 = NULL, beta3 = NULL, model_type = NULL, coeff = NULL,
                                    mu = NULL, Sig = NULL, snr = NULL){
  X = mvrnorm(n = n, mu = mu, Sigma = Sig)  # nxp-matrix
  if(model_type == "L-L"){
    norm2.1 =  as.numeric(t(beta1)%*%Sig%*%beta1) # norm^2 to calculate sd1
    sd1 = sqrt(norm2.1/snr)
    norm2.2 =  as.numeric(t(beta2)%*%Sig%*%beta2) # norm^2 to calculate sd2
    sd2 = sqrt(norm2.2/snr)
    m1X = X%*%beta1
    Y1 = alpha1 + 2*m1X + rnorm(n, sd = 2*sd1) # I just added 5: 06/13
    m0X = alpha2 + X%*%beta2
    Y0 = m0X + rnorm(n, sd = sd2)
    #generate treatment indicator T|X #calculate pi(X)
    a = as.vector(X%*%beta3) - 0.08 # I just added -0.08 for s = 7, d= 2: 06/15
    # I just added 0.08 for s = 15, d= 5: 06/15
    px = exp(a)/(1 + exp(a))
    trt = rbinom(n, size = 1, prob = px)
    indx1 = which(trt == 1)
    indx0 = which(trt == 0)
    Y = trt*Y1 + (1 - trt)*Y0
  }else if(model_type == "L-C"){ # quadratic m(X) and logistic p(X)
    norm2.1 =  as.numeric(t(beta1)%*%Sig%*%beta1) # norm^2 to calculate sd1
    sd1 = sqrt(norm2.1/snr)
    norm2.2 =  as.numeric(t(beta2)%*%Sig%*%beta2) # norm^2 to calculate sd2
    sd2 = sqrt(norm2.2/snr)
    m1X = X%*%beta1
    Y1 = alpha1 + 2*m1X + rnorm(n, sd = 2*sd1) # I just added 5: 06/13
    m0X = alpha2 + X%*%beta2
    Y0 = m0X + rnorm(n, sd = sd2)
    #generate treatment indicator T|X #calculate pi(X)
    a = sin(as.vector(X%*%beta3) - 0.08)
    px = exp(a)/(1 + exp(a))
    trt = rbinom(n, size = 1, prob = px)
    indx1 = which(trt == 1)
    indx0 = which(trt == 0)
    Y = trt*Y1 + (1 - trt)*Y0  
  }else if(model_type == "L-ONL"){ # quadratic m(X) and osculated p(X)
    norm2.1 =  as.numeric(t(beta1)%*%Sig%*%beta1) # norm^2 to calculate sd1
    sd1 = sqrt(norm2.1/snr)
    norm2.2 =  as.numeric(t(beta2)%*%Sig%*%beta2) # norm^2 to calculate sd2
    sd2 = sqrt(norm2.2/snr)
    m1X = X%*%beta1
    Y1 = 2*m1X + rnorm(n, sd = 2*sd1)
    m0X = X%*%beta2
    Y0 = m0X + rnorm(n, sd = sd2)
    #generate treatment indicator T|X
    #calculate pi(X)
    ex = exp(sin(X))
    px = ex/(1 + ex)
    trt = rbinom(n, size = 1, prob = px)
    indx1 = which(trt == 1)
    indx0 = which(trt == 0)
    #observed Y values
    Y = trt*Y1 + (1 - trt)*Y0
  }else if(model_type == "Q-L"){ # quadratic m(X) and logistic p(X)
    #generate Y(1)|X and Y(0)|X
    m1X = alpha1 + 2*X%*%beta1 + coeff*(X^2)%*%beta12
    norm2.1 =  var(m1X) # norm^2 to calculate sd1
    sd1 = sqrt(norm2.1/snr)
    Y1 = m1X + rnorm(n, sd = sd1) #coeff is calculated before
    # mean((X%*%beta1)^2)/mean((coeff^2)*(X%*%beta12)^4)
    m0X = alpha2 + X%*%beta2
    norm2.2 =  var(m0X) # norm^2 to calculate sd2
    sd2 = sqrt(norm2.2/snr)
    Y0 = m0X + rnorm(n, sd = sd2)
    #generate treatment indicator T|X #calculate pi(X)
    a = as.vector(X%*%beta3) - 0.08
    px = exp(a)/(1 + exp(a))
    
    trt = rbinom(n, size = 1, prob = px)
    indx1 = which(trt == 1)
    indx0 = which(trt == 0)
    Y = trt*Y1 + (1 - trt)*Y0
    }else if(model_type == "Q-C"){ # quadratic m(X) and logistic p(X)
      #generate Y(1)|X and Y(0)|X
      m1X = alpha1 + 2*X%*%beta1 + coeff*(X^2)%*%beta12
      norm2.1 =  var(m1X) # norm^2 to calculate sd1
      sd1 = sqrt(norm2.1/snr)
      Y1 = m1X + rnorm(n, sd = sd1) #coeff is calculated before
      # mean((X%*%beta1)^2)/mean((coeff^2)*(X%*%beta12)^4)
      m0X = alpha2 + X%*%beta2
      norm2.2 =  var(m0X) # norm^2 to calculate sd2
      sd2 = sqrt(norm2.2/snr)
      Y0 = m0X + rnorm(n, sd = sd2)
      #generate treatment indicator T|X #calculate pi(X)
      a = sin(as.vector(X%*%beta3) - 0.08)
      px = exp(a)/(1 + exp(a))
      trt = rbinom(n, size = 1, prob = px)
      indx1 = which(trt == 1)
      indx0 = which(trt == 0)
      Y = trt*Y1 + (1 - trt)*Y0  
  }else if(model_type == "Q-ONL"){ # quadratic m(X) and osculated p(X)
    #generate Y(1)|X and Y(0)|X
    norm2.1 =  as.numeric(t(beta1)%*%Sig%*%beta1) # norm^2 to calculate sd1
    sd1 = sqrt(norm2.1/snr)
    Y1 = X%*%beta1 + coeff*(X%*%beta12)^2 + rnorm(n, sd = sd1) #coeff is calculated before
    norm2.2 =  as.numeric(t(beta2)%*%Sig%*%beta2) # norm^2 to calculate sd2
    sd2 = sqrt(norm2.2/snr)
    Y0 = X%*%beta2 + rnorm(n, sd = sd2)
    #generate treatment indicator T|X
    #calculate pi(X)
    ex = exp(sin(X))
    px = ex/(1 + ex)
    trt = rbinom(n, size = 1, prob = px)
    indx1 = which(trt == 1)
    indx0 = which(trt == 0)
    #observed Y values
    Y = trt*Y1 + (1 - trt)*Y0
  }else if(model_type == "SNL-L"){
    m1X = as.vector(X%*%beta1)*abs(as.vector(X%*%beta1)) # f(x) = x|x|
    norm2.1 = var(m1X)
    sd1 = sqrt(norm2.1/snr)
    Y1 = m1X + rnorm(n, sd = sd1)
    norm2.2 =  as.numeric(t(beta2)%*%Sig%*%beta2) # norm^2 to calculate sd2
    sd2 = sqrt(norm2.2/snr)
    m0X = X%*%beta2
    Y0 = m0X + rnorm(n, sd = sd2)
    #generate treatment indicator T|X #calculate pi(X)
    a = as.vector(X%*%beta3) + 0.1
    px = exp(a)/(1 + exp(a))
    
    trt = rbinom(n, size = 1, prob = px)
    indx1 = which(trt == 1)
    indx0 = which(trt == 0)
    Y = trt*Y1 + (1 - trt)*Y0
  }else if(model_type == "NSNL-L"){
    m1X = as.vector(1*(abs(X%*%beta1)<=1)) + as.vector(X%*%beta1) # f(x) = I(|X| \leq 1)
    norm2.1 = var(m1X)
    sd1 = sqrt(norm2.1/snr)
    Y1 = m1X + rnorm(n, sd = sd1)
    norm2.2 =  as.numeric(t(beta2)%*%Sig%*%beta2) # norm^2 to calculate sd2
    sd2 = sqrt(norm2.2/snr)
    m0X = X%*%beta2
    Y0 = m0X + rnorm(n, sd = sd2)
    #generate treatment indicator T|X #calculate pi(X)
    a = as.vector(X%*%beta3) + 0.1
    px = exp(a)/(1 + exp(a))
    
    trt = rbinom(n, size = 1, prob = px)
    indx1 = which(trt == 1)
    indx0 = which(trt == 0)
    Y = trt*Y1 + (1 - trt)*Y0
  }else{
    m1X = sin(2*pi*(X%*%beta1))
    norm2.1 = var(m1X)
    sd1 = sqrt(norm2.1/snr)
    Y1 = m1X + rnorm(n, sd = sd1)
    norm2.2 =  as.numeric(t(beta2)%*%Sig%*%beta2) # norm^2 to calculate sd2
    sd2 = sqrt(norm2.2/snr)
    m0X = X%*%beta2
    Y0 = m0X + rnorm(n, sd = sd2)
    #generate treatment indicator T|X #calculate pi(X)
    a = as.vector(X%*%beta3) + 0.1
    px = exp(a)/(1 + exp(a))
    
    trt = rbinom(n, size = 1, prob = px)
    indx1 = which(trt == 1)
    indx0 = which(trt == 0)
    Y = trt*Y1 + (1 - trt)*Y0
  }
  
  return( list(Y = as.vector(Y), X = as.matrix(X), Y1 = as.vector(Y1), Y0 = as.vector(Y0), trt = as.vector(trt), indx1 = indx1, indx0 = indx0) )
}
