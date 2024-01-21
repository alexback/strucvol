###################### Standard ARSV model estimation ##########################


######################## Log likelihood function ###############################

svloglik2 <- function(param, yt, Ht, fit = F){
  
  
  sp <- spmaker2(param = param, yt = yt, Ht = Ht)
  
  kf <- svkf(Tt = sp$Tt, Qt = sp$Qt, Zt = sp$Zt, Ht = Ht, a0 = sp$a0,
             P0 = sp$P0, yt = sp$yt)
  
  loglik <- -0.5 * (log(2 * pi) * length(yt[-c(1:3)]) +
                      sum(log(abs(kf$Ft[-c(1:3)])) + kf$vt[-c(1:3)] * (1 / kf$Ft[-c(1:3)]) * kf$vt[-c(1:3)]))
  if (fit == T){
    return(-loglik)
  }
  
  else{
    
    return(list(loglik = loglik, kf = kf))
  }
  
  
  
}



################################################################################

######################## Log likelihood function ###############################

svloglik <- function(param, yt, Ht, fit = F){
  
  
  sp <- spmaker(param = param, yt = yt, Ht = Ht)
  
  kf <- svkf(Tt = sp$Tt, Qt = sp$Qt, Zt = sp$Zt, Ht = Ht, a0 = sp$a0,
             P0 = sp$P0, yt = sp$yt)
  
  loglik <- -0.5 * (log(2 * pi) * length(yt[-c(1:3)]) +
                      sum(log(abs(kf$Ft[-c(1:3)])) + kf$vt[-c(1:3)] * (1 / kf$Ft[-c(1:3)]) * kf$vt[-c(1:3)]))
  if (fit == T){
    return(-loglik)
  }
  
  else{
    
    return(list(loglik = loglik, kf = kf))
  }
  
  
  
}



################################################################################

################################ Ht maker ######################################

Htm <- function(param, yt, tol = 0.001){
  l <- length(yt)
  crit <- 1000
  Ht <- rep(pi^2 / 2, l)
  
  i <- 0
  
  while(crit > tol & i<100){
    i <- i + 1
    
    fsp <- spmaker(param = param, yt = yt, Ht = Ht)
    # Filter and smooth to obtain et 
    kf <- svkf(Tt = fsp$Tt, Qt = fsp$Qt, Zt = fsp$Zt, Ht = Ht,
               a0 = fsp$a0, P0 = fsp$P0, yt = fsp$yt)
    ks <- svks(vt = kf$vt, Ft = kf$Ft, Kt = kf$Kt, Tt = fsp$Tt,
               Qt = fsp$Qt, Zt = fsp$Zt,
               Ht = Ht, a0 = fsp$a0, P0 = fsp$P0, yt = fsp$yt)
    
    et <- ks$epsilont
    
    
    tempHt <- (2 * et) / (exp(et) - 1)  
    
    crit <- mean(abs(tempHt - Ht))
    
    Ht <- tempHt
    
  }
  
  kf <- svkf(Tt = fsp$Tt, Qt = fsp$Qt, Zt = fsp$Zt, Ht = Ht,
             a0 = fsp$a0, P0 = fsp$P0, yt = fsp$yt)
  ks <- svks(vt = kf$vt, Ft = kf$Ft, Kt = kf$Kt, Tt = fsp$Tt,
             Qt = fsp$Qt, Zt = fsp$Zt,
             Ht = Ht, a0 = fsp$a0, P0 = fsp$P0, yt = fsp$yt)

  
  return(list(Ht = Ht, fsp = fsp, fkf = kf, epsilont = et))
  
}

################################################################################ 
#'@title Calculate the importance weights
#'@export
wcalc <- function(tempeps, epsilont, Ht){
  
  antitempeps <- 2 * epsilont - tempeps
  
  wtempeps <- 0.5 * (log(Ht) + tempeps - exp(tempeps) + tempeps^2 / Ht)
  
  wantitempeps <- 0.5 * (log(Ht) + antitempeps - exp(antitempeps) + antitempeps^2 / Ht)
  
  w <- exp(sum((wtempeps + wantitempeps) / 2))
  
  
  return(w)  
  
}



#################### Sandmann and Koopman algorithm ############################
#'@title Calculate the Sandmann and Koopman log likelihood
#'@import Brobdingnag
#'@export
skloglik <- function(param, yt, N = 5, fit = T){
  
  l <- length(yt)
  stageone <- Htm(param = param, yt = yt)
  
  
  Ht <- stageone$Ht
  
  fsp <- stageone$fsp
  
  fkf <- stageone$fkf
  
  epsilont <- stageone$epsilont
  
  
  loglik <- svloglik(param = param, yt = yt, Ht = Ht, fit = T)
  
  w <- numeric(N)
  
  set.seed(123)
  tempepsdf <- replicate(N, svdsim(vt = fkf$vt, Ft = fkf$Ft, Kt = fkf$Kt, Tt = fsp$Tt, Qt = fsp$Qt, Zt = fsp$Zt, Ht = Ht, fsp$a0, fsp$P0,
                                   yt = fsp$yt))
  
  w <- apply(tempepsdf, 2, wcalc, epsilont = epsilont, Ht = Ht)
  
  #for(i in 1:N){
  
  
  #tempeps <- svdsim(vt = fkf$vt, Ft = fkf$Ft, Kt = fkf$Kt, Tt = fsp$Tt, Qt = fsp$Qt, Zt = fsp$Zt, Ht = Ht, fsp$a0, fsp$P0,
  # yt = fsp$yt)
  
  
  #antitempeps <- 2 * epsilont - tempeps
  
  #wtempeps <- 0.5 * (log(Ht) + tempeps - exp(tempeps) + tempeps^2 / Ht)
  
  #wantitempeps <- 0.5 * (log(Ht) + antitempeps - exp(antitempeps) + antitempeps^2 / Ht)
  
  #w[i] <- exp(sum((wtempeps + wantitempeps) / 2))
  
  
  #}
  
  w <- as.brob(w)
  wbar <- sum(w) / N
  swsquared <- (sum((w - wbar)^2)) / (N - 1) 
  print(swsquared)
  num <- swsquared
  denom <- 2 * N * wbar^2
  
  
  remainder <- log(wbar) + num / denom
  remainder <- as.numeric(remainder)
  
  if(fit == T){

    return(loglik - remainder)
  }
  else {
    
    return(list(w = w, Ht = Ht, kf = stageone$fkf, epsilont = epsilont))
  }
  
}




######################### Structural model estimation ##########################
######################### Log likelihood #######################################
# The Kalman filter (quasi-) log likelihood.

ssvloglik <- function(param, yt, Ht, xt, fit = T){
  
  
  sp <- sspmaker(param = param, yt = yt, Ht = Ht, xt = xt)
  
  kf <- ssvkf(Tt = sp$Tt, Qt = sp$Qt, Zt = sp$Zt, Ht = Ht, a0 = sp$a0,
             P0 = sp$P0, yt = sp$yt)
  
  loglik <- -0.5 * (log(2 * pi) * length(yt[-c(1:3)]) +
                      sum(log(abs(kf$Ft[-c(1:3)])) + kf$vt[-c(1:3)] * (1 / kf$Ft[-c(1:3)]) * kf$vt[-c(1:3)]))
  if (fit == T){
    return(-loglik)
  }
  
  else{
    
    return(list(loglik = loglik, kf = kf))
  }
}


################################ Ht maker ######################################

# Calculate values for Ht to be used in the importance sampling weights.

sHtm <- function(param, yt, xt, tol = 0.001){
  l <- length(yt)
  crit <- 1000
  Ht <- rep(pi^2 / 2, l)
  
  i <- 0
  
  while(crit > tol & i<100){
    i <- i + 1
    
    fsp <- sspmaker(param = param, yt = yt, Ht = Ht, xt = xt)
    # Filter and smooth to obtain et 
    kf <- ssvkf(Tt = fsp$Tt, Qt = fsp$Qt, Zt = fsp$Zt, Ht = Ht,
               a0 = fsp$a0, P0 = fsp$P0, yt = fsp$yt)
    ks <- ssvks(vt = kf$vt, Ft = kf$Ft, Kt = kf$Kt, Tt = fsp$Tt,
               Qt = fsp$Qt, Zt = fsp$Zt,
               Ht = Ht, a0 = fsp$a0, P0 = fsp$P0, yt = fsp$yt)
    
    et <- ks$epsilont
    
    
    tempHt <- (2 * et) / (exp(et) - 1)  
    
    crit <- mean(abs(tempHt - Ht))

    Ht <- tempHt
    
  }
  
  fsp <- sspmaker(param = param, yt = yt, Ht = Ht, xt = xt)
  
  kf <- ssvkf(Tt = fsp$Tt, Qt = fsp$Qt, Zt = fsp$Zt, Ht = Ht,
             a0 = fsp$a0, P0 = fsp$P0, yt = fsp$yt)
  ks <- ssvks(vt = kf$vt, Ft = kf$Ft, Kt = kf$Kt, Tt = fsp$Tt,
             Qt = fsp$Qt, Zt = fsp$Zt,
             Ht = Ht, a0 = fsp$a0, P0 = fsp$P0, yt = fsp$yt)
  
  return(list(Ht = Ht, fsp = fsp, fkf = kf, epsilont = et))
  
}

################################################################################ 

##################### Importance weight calculator #############################

# Helper function, calculates the importance sampling weights. 
swcalc <- function(tempeps, epsilont, Ht){
  
  antitempeps <- 2 * epsilont - tempeps
  
  wtempeps <- 0.5 * (log(Ht) + tempeps - exp(tempeps) + tempeps^2 / Ht)
  
  wantitempeps <- 0.5 * (log(Ht) + antitempeps - exp(antitempeps) + antitempeps^2 / Ht)
  
  w <- exp(sum((wtempeps + wantitempeps) / 2))
  
  
  return(w)  
  
}

################################################################################


#################### Sandmann and Koopman algorithm ############################

# Sandmann and Koopmann algorithm for importance sampling log likelihood.
# Commented out for loop illustrates what the function wcalc() does.
#'@title Calculate the Sandmann and Koopman log likelihood
#'@import Brobdingnag
#'@export
sskloglik <- function(param, yt, xt, N = 5, fit = T){
  
  l <- length(yt)
  stageone <- sHtm(param = param, yt = yt, xt = xt)
  
  
  Ht <- stageone$Ht
  
  fsp <- stageone$fsp
  
  fkf <- stageone$fkf
  
  epsilont <- stageone$epsilont
  
  
  loglik <- ssvloglik(param, yt, Ht = Ht, xt = xt)
  
  w <- numeric(N)
  
  set.seed(123)
  tempepsdf <- replicate(N, ssvdsim(vt = fkf$vt, Ft = fkf$Ft, Kt = fkf$Kt,
                                   Tt = fsp$Tt, Qt = fsp$Qt, Zt = fsp$Zt, 
                                   Ht = Ht, fsp$a0, fsp$P0,
                                   yt = fsp$yt))
  
  w <- apply(tempepsdf, 2, swcalc, epsilont = epsilont, Ht = Ht)
  w <- as.brob(w)
  
  wbar <- sum(w) / N
  swsquared <- (sum((w - wbar)^2)) / (N - 1) 
  print(swsquared)
  num <- swsquared
  denom <- 2 * N * wbar^2
  
  remainder <- log(wbar) + (num / denom)
  remainder <- as.numeric(remainder)
  print(remainder)
  
  

  
  if(fit == T){
    
    return(loglik - remainder)
  }
  else {
    
    return(list(w = w, Ht = Ht, kf = stageone$fkf, ks = stageone$fks))
  }
}

################################################################################


############## Multivariate stochastic volatility  Log likelihood ##############
mssvloglik <- function(param, yt, xt, fit = T){
  
  
  sp <- msspmaker(param = param, yt = yt, xt = xt)
  
  # kf <- svkf(Tt = sp$Tt, Qt = sp$Qt, Zt = sp$Zt, Ht = sp$Ht, a0 = sp$a0,
  #  P0 = sp$P0, yt = sp$yt)
  
  kf <- fkf(a0 = sp$a0, P0 = sp$P0, dt = c(0,0,0,0,0), ct = c(0,0),
            Tt = sp$Tt, Zt = sp$Zt, HHt = sp$Qt, GGt = sp$Ht, 
            yt = t(sp$yt))
  
  if (fit == T){
    return(-sum(kf$logLik))
  }
  
  else{
    
    return(kf)
  }
}
################################################################################

###################### Asymmetric model log likelihoods ########################
# Remove initial observations from ll due to diffuse init. Remove d obs.,
# where d = dim(state vector).
asysvloglik <- function(param, yt){
  
  
  sp <- asyspmaker(param = param, yt = yt)
  
  kf <- asysvkf(sp$Tt, sp$Qt, sp$Zt, Ht = sp$Ht, sp$Gt, sp$ct, sp$a0, sp$P0, sp$yt)
  
  loglik <- -sum(log(abs(kf$Ft[2:length(kf$Ft)])) +
                   kf$vt[2:length(kf$vt)] *
                   (1 / kf$Ft[2:length(kf$Ft)]) * kf$vt[2:length(kf$vt)])
  
  return(-loglik)
  
  
}
asysvloglikt <- function(param, yt){
  
  
  sp <- asyspmaker(param = param, yt = yt)
  
  kf <- asysvkf(sp$Tt, sp$Qt, sp$Zt, sp$Ht, sp$Gt, sp$ct, sp$a0, sp$P0, sp$yt)
  
  loglik <- -(log(abs(kf$Ft[2:length(kf$Ft)])) +
                kf$vt[2:length(kf$vt)] *
                (1 / kf$Ft[2:length(kf$Ft)]) * kf$vt[2:length(kf$vt)])
  
  return(-loglik)
  
  
}

################################################################################

############# Asymmetric model with covariates log likelihoods #################
# Remove initial observations from ll due to diffuse init. Remove d obs.,
# where d = dim(state vector).
asysvloglikx <- function(param, yt, x){
  
  
  sp <- asyspmakerx(param = param, yt = yt, x = x)
  
  kf <- asysvkfx(sp$Tt, sp$Qt, sp$Zt, Ht = sp$Ht, sp$Gt, sp$ct, sp$a0, sp$P0, sp$yt)
  
  loglik <- -sum(log(abs(kf$Ft[3:length(kf$Ft)])) +
                   kf$vt[3:length(kf$vt)] *
                   (1 / kf$Ft[3:length(kf$Ft)]) * kf$vt[3:length(kf$vt)])
  
  return(-loglik)
  
  
}

asysvlogliktx <- function(param, yt, x){
  
  
  sp <- asyspmakerx(param = param, yt = yt, x = x)
  
  kf <- asysvkfx(sp$Tt, sp$Qt, sp$Zt, sp$Ht, sp$Gt, sp$ct, sp$a0, sp$P0, sp$yt)
  
  loglik <- -(log(abs(kf$Ft[3:length(kf$Ft)])) +
                kf$vt[3:length(kf$vt)] *
                (1 / kf$Ft[3:length(kf$Ft)]) * kf$vt[3:length(kf$vt)])
  
  return(-loglik)
  
  
}
################################################################################