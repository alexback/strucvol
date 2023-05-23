#'@title simstrucsystem
#'@description Simulate a bivariate structural volatility model.
#'@importFrom MASS mvrnorm 
#'@export

simstrucsystem <- function(len = 2000, pars = c(-10, 0.95, 0.3, -12, 0.9, 0.2, 0.9),
                           Ain = 100, Ein = 80, K = 20, r = 0.001, uv = 0.0005, ttv = 5)
  {
   
    mu1 <- pars[1]
    beta1 <- pars[2]
    sigma1 <- pars[3]
    
    mu2 <- pars[4]
    beta2 <- pars[5]
    sigma2 <- pars[6]
  
    rho <- pars[7]

    len <- as.integer(len)
    
    h_1 <- rep_len(0, length.out = len)
    h_2 <- rep_len(0, length.out = len)
    
    h0_1 <- rnorm(1, mean = mu1, sd = sigma1 / sqrt(1 - beta1^2))
    h0_2 <- rnorm(1, mean = mu2, sd = sigma2 / sqrt(1 - beta2^2))
    
    
    eta1 <- rnorm(len)
    eta2 <- rnorm(len)
    
    
    sigmaeps <- matrix(c(1, rho, rho, 1), byrow = T, ncol = 2, nrow = 2) 
    eps <- mvrnorm(n = len, mu = c(0,0), Sigma = sigmaeps)
    
    eps1 <- eps[,1]
    eps2 <- eps[,2]
    
    
    h_1[1] <- mu1 + beta1 * (h0_1 - mu1) + sigma1 * eta1[1]  # Initialize
    h_2[1] <- mu2 + beta2 * (h0_2 - mu2) + sigma2 * eta2[1]  # Initialize
    
    for (i in 2:length(h_1)){
    h_1[i] <- mu1 + beta1 * (h_1[i-1] - mu1) + sigma1 * eta1[i]
    h_2[i] <- mu2 + beta2 * (h_2[i-1] - mu2) + sigma2 * eta2[i]
    
    }
    
    ra <- exp(h_1 / 2) * eps1
    rm <- exp(h_2 / 2) * eps2
    
    
    A <- cumprod(c(Ain, exp(ra)))

    AK <- A / K

    d1 <- (log(AK) + (r + 0.5 * sqrt(uv * 252)^2 ) * ttv) / (sqrt(uv * 252 * ttv))
    
    
    E <- rep(0,len)
    E[1] <- Ein
    lm <- rep(0, len-1)
    
    re <- rep(0, len)
    re[1] <- ra[1]
    
    for (j in 2:len){
    lm[j -1] <- pnorm(d1[j -1]) * A[j - 1] / E[j - 1]
    re[j] <- lm[j -1] * ra[j]
    
    E[j] <- E[j - 1] * exp(re[j])
    }
    
    EK = E / K
    return(list(retseries = cbind(ra, re, rm), lm = lm, E = E, A = A, EK = EK))
    
  
}

