###################### Standard ARSV model estimation ##########################


######################## Log likelihood function ###############################

svloglik <- function(param, yt, Ht, fit = F){
  
  
  sp <- spmaker(param = param, yt = yt, Ht = Ht)
  
  kf <- svkf(Tt = sp$Tt, Qt = sp$Qt, Zt = sp$Zt, Ht = Ht, a0 = sp$a0,
             P0 = sp$P0, yt = sp$yt)
  
  loglik <- -0.5 * (log(2 * pi) * length(yt) +
                      sum(log(abs(kf$Ft)) + kf$vt * (1 / kf$Ft) * kf$vt))
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
wcalc <- function(tempeps, epsilont, Ht){
  
  antitempeps <- 2 * epsilont - tempeps
  
  wtempeps <- 0.5 * (log(Ht) + tempeps - exp(tempeps) + tempeps^2 / Ht)
  
  wantitempeps <- 0.5 * (log(Ht) + antitempeps - exp(antitempeps) + antitempeps^2 / Ht)
  
  w <- exp(sum((wtempeps + wantitempeps) / 2))
  
  
  return(w)  
  
}

#################### Sandmann and Koopman algorithm ############################

skloglik <- function(param, yt, N = 5, fit = T){
  
  l <- length(yt)
  stageone <- Htm(param = param, yt = yt)
  
  
  Ht <- stageone$Ht
  
  fsp <- stageone$fsp
  
  fkf <- stageone$fkf
  
  epsilont <- stageone$epsilont
  
  
  loglik <- svloglik(param, yt, Ht = Ht, fit = T)
  
  w <- numeric(N)
  
  set.seed(8493)
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
  
  wbar <- mean(w)
  swsquared <- var(w)
  num <- swsquared
  denom <- 2 * N * wbar^2
  
  
  remainder <- log(wbar) + num / denom
  
  if(fit == T){
    
    return(loglik - remainder)
  }
  else {
    
    return(list(w = w, Ht = Ht, kf = stageone$fkf, ks = stageone$fks))
  }
  
}





######################### Structural model estimation ##########################
######################### Log likelihood #######################################
# The Kalman filter (quasi-) log likelihood.

ssvloglik <- function(param, yt, Ht, xt, fit = T){
  
  
  sp <- sspmaker(param = param, yt = yt, Ht = Ht, xt = xt)
  
  kf <- ssvkf(Tt = sp$Tt, Qt = sp$Qt, Zt = sp$Zt, Ht = Ht, a0 = sp$a0,
             P0 = sp$P0, yt = sp$yt)
  
  loglik <- -0.5 * (log(2 * pi) * length(yt) +
                      sum(log(abs(kf$Ft)) + kf$vt * (1 / kf$Ft) * kf$vt))
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

sskloglik <- function(param, yt, xt, N = 5, fit = T){
  
  l <- length(yt)
  stageone <- sHtm(param = param, yt = yt, xt = xt)
  
  
  Ht <- stageone$Ht
  
  fsp <- stageone$fsp
  
  fkf <- stageone$fkf
  
  epsilont <- stageone$epsilont
  
  
  loglik <- ssvloglik(param, yt, Ht = Ht, xt = xt)
  
  w <- numeric(N)
  
  set.seed(8493)
  tempepsdf <- replicate(N, ssvdsim(vt = fkf$vt, Ft = fkf$Ft, Kt = fkf$Kt,
                                   Tt = fsp$Tt, Qt = fsp$Qt, Zt = fsp$Zt, 
                                   Ht = Ht, fsp$a0, fsp$P0,
                                   yt = fsp$yt))
  
  w <- apply(tempepsdf, 2, swcalc, epsilont = epsilont, Ht = Ht)
  
  wbar <- mean(w)
 
  swsquared <- var(w)

  num <- swsquared
  denom <- 2 * N * wbar^2
  
  
  remainder <- log(wbar) + num / denom
  
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