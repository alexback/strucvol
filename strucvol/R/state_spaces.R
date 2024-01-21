# Standard ARSV model state space maker. 1: MCL, 2: testing.
spmaker <- function(param, yt, Ht){
  
  ya <- log(yt^2)
  l <- length(ya)
  
  Tt <- diag(2) * c(param[1], 1) # OK
  
  Qt <- diag(2) * c(param[2]^2, 0)
  
  Zt <- array(dim = c(1,2,l))
  Zt[1,1,] <- 1
  Zt[1,2,] <- 1
  
  a0 <- c(0, 0)
  a <- array(dim = c(2,1,l))
  a[,,1] <- a0
  
  
  p <- array(dim = c(2,2,l))
  P0 <- matrix(diag(c(param[2]^2 / (1 - param[1]^2),100000)), nrow = 2, ncol = 2)
  
  
  
  return(list(Tt = Tt, Qt = Qt, Zt = Zt, a0 = a0, P0 = P0, yt = ya))
}


spmaker2 <- function(param, yt, Ht){
  
  ya <- log(yt^2) + 1.27
  l <- length(ya)
  
  Tt <- diag(2) * c(param[1], 1) # OK
  
  Qt <- diag(2) * c(param[2]^2, 0)
  
  Zt <- array(dim = c(1,2,l))
  Zt[1,1,] <- 1
  Zt[1,2,] <- 1
  
  a0 <- c(0, 0)
  a <- array(dim = c(2,1,l))
  a[,,1] <- a0
  
  
  p <- array(dim = c(2,2,l))
  P0 <- matrix(diag(c(param[2]^2 / (1 - param[1]^2),100000)), nrow = 2, ncol = 2)
  
  
  
  return(list(Tt = Tt, Qt = Qt, Zt = Zt, a0 = a0, P0 = P0, yt = ya))
}

# Structural stochastic volatility state space maker.
sspmaker <- function(param, yt, Ht, xt){
  
  ya <- log(yt^2)
  l <- length(ya)
  
  Tt <- diag(3) * c(param[1], 1, 1) # OK
  
  Qt <- diag(3) * c(param[2]^2, 0, 0)
  
  Zt <- array(dim = c(1,3,l))
  Zt[1,1,] <- 1
  Zt[1,2,] <- 1
  Zt[1,3,] <- xt 
  
  a0 <- c(0, 0, 0)
  a <- array(dim = c(3,1,l))
  a[,,1] <- a0
  
  
  p <- array(dim = c(3,3,l))
  P0 <- matrix(diag(c(param[2]^2 / (1 - param[1]^2),100000, 100000)), nrow = 3, ncol = 3)
  
  
  
  return(list(Tt = Tt, Qt = Qt, Zt = Zt, a0 = a0, P0 = P0, yt = ya))
}

# Return approximation based on Stirlings formula for large j to avoid "inf".
# Needed in the multivariate state space formulation.
gammaratio <- function(x){
  if (x <= 171){
    
    return(gamma(x) / gamma(0.5 + x))
    
  }
  
  else{
    
    return(x^(-0.5))
  }
  
}

# Multivariate state space maker for structural stochastic volatility.
msspmaker <- function(param, yt, xt){
  ya <- log(yt^2) + 1.27
  l <- length(ya[,1])
  
  ind <- ifelse((sum(yt[,1] * yt[,2] >= 0) / l) > 0.5, 1,-1)
  
  beta1 <- param[1]
  beta2 <- param[2]
  
  sigma1 <- param[3]
  sigma2 <- param[4]
  
  hrho <- param[5]
  
  Tt <- diag(5) * c(beta1, beta2, 1, 1, 1) # OK
  
  Qt <- diag(5) * c(sigma1^2, sigma2^2, 0, 0, 0)
  
  Zt <- array(dim = c(2,5,l))
  Zt[1,1,] <- 1
  Zt[1,2,] <- 0
  Zt[1,3,] <- 1
  Zt[1,4,] <- 0
  Zt[1,5,] <- xt
  
  Zt[2,1,] <- 0
  Zt[2,2,] <- 1
  Zt[2,3,] <- 0
  Zt[2,4,] <- 1
  Zt[2,5,] <- 0
  
  a0 <- c(0, 0, 0, 0, 0)
  
  P0 <- diag(5) * c(sigma1^2 / (1 - beta1^2), sigma2^2 / (1 - beta2^2), 
                    10^5, 10^5, 10^5)
  
  
  j <- 1:10000
  
  hrhorec <- sum(gamma(0.5) * (sapply(j, gammaratio) / j) * hrho^(2 * j))
  
  Ht <- matrix(c(pi^2 / 2, hrhorec, hrhorec, pi^2 / 2), nrow = 2, ncol = 2, byrow = T)
  
  return(list(Tt = Tt, Qt = Qt, Zt = Zt, Ht = Ht, a0 = a0, P0 = P0, yt = ya, ind = ind))
}

#################### Asymmetric model state space maker ########################

asyspmaker <- function(param, yt){
  
  l <- length(yt)
  Tt <- diag(2) * c(param[1], 1) # OK
  A <- 0.7979 * param[3] * param[2]
  B <- 1.1061 * param[3] * param[2]
  
  st <- sign(yt)
  
  Qt <- diag(2) * c(param[2]^2 - A^2, 0)
  
  Zt <- array(dim = c(1,2,l))
  Zt[1,1,] <- 1
  Zt[1,2,] <- 1
  
  
  Gt <- array(dim = c(2,1,l))
  Gt[1,1,] <- B * st
  Gt[2,1,] <- 0
  
  
  ct <- array(dim = c(2,1,l))
  ct[1,1,] <- A * st
  ct[2,1,] <- 0
  
  
  a0 <- c(0, 0)
  a <- array(dim = c(2,1,l))
  a[,,1] <- a0
  
  Ht <- rep(pi^2/2, l)
  
  P0 <- matrix(diag(c(param[2]^2 / (1 - param[1]^2), 10000000)),
               nrow = 2, ncol = 2)
  yt <- log(yt^2) + 1.27
  return(list(Tt = Tt, Qt = Qt, Zt = Zt, Gt = Gt, Ht = Ht, ct = ct, yt = yt, a0 = a0,P0 = P0))
}

################################################################################

############ Asymmetric model with covariates state space maker ################

asyspmakerx <- function(param, yt, x){
  
  l <- length(yt)
  Tt <- diag(3) * c(param[1], 1, 1) # OK
  A <- 0.7979 * param[3] * param[2]
  B <- 1.1061 * param[3] * param[2]
  
  st <- sign(yt)
  
  Qt <- diag(3) * c(param[2]^2 - A^2, 0, 0)
  
  Zt <- array(dim = c(1,3,l))
  Zt[1,1,] <- 1
  Zt[1,2,] <- 1
  Zt[1,3,] <- x
  
  Gt <- array(dim = c(3,1,l))
  Gt[1,1,] <- B * st
  Gt[2,1,] <- 0
  Gt[3,1,] <- 0
  
  ct <- array(dim = c(3,1,l))
  ct[1,1,] <- A * st
  ct[2,1,] <- 0
  ct[3,1,] <- 0
  
  a0 <- c(0, 0, 0)
  a <- array(dim = c(3,1,l))
  a[,,1] <- a0
  
  Ht <- rep(pi^2/2, l)
  
  P0 <- matrix(diag(c(param[2]^2 / (1 - param[1]^2), 10000000, 10000000)),
               nrow = 3, ncol = 3)
  yt <- log(yt^2) + 1.27
  return(list(Tt = Tt, Qt = Qt, Zt = Zt, Gt = Gt, Ht = Ht, ct = ct, yt = yt, a0 = a0,P0 = P0))
}

################################################################################