
######################### Standard ARSV model ##################################
########################## The Kalman filter ###################################


svkf <- function(Tt, Qt, Zt, Ht, Gt, a0, P0, yt){
  l <- length(yt)
  
  vt <- numeric(l)
  Kt <- array(dim = c(2,1,l))
  at <- array(dim = c(2,1,l + 1))
  at[,,1] <- a0
  Ft <- numeric(l)
  Lt <- array(dim = c(2,2,l))
  Pt <- array(dim = c(2,2,l + 1))
  Pt[,,1] <- P0
  
  
  
  for (i in 1:l){
    
    
    Qti <- Qt
    
    Zti <- t(matrix(Zt[,,i]))
    
    
    Pti <- matrix(Pt[,,i], nrow = 2, ncol = 2)
    
    vt[i] <- yt[i] - Zti%*%at[,,i]
    
    Ft[i] <- Zti%*%Pti%*%t(Zti) + Ht[i]
    
    Kti <- (Tt%*%Pti%*%t(Zti)) / Ft[i]
    
    Kt[,,i] <- Kti
    
    Lti <- Tt - Kti%*%Zti
    
    Lt[,,i] <- Lti
    
    at[,,i + 1] <- Tt%*%at[,,i] + Kti%*%vt[i]
    
    Pt[,,i + 1] <- Tt%*%Pti%*%t(Lti) + Qti 
    
  }
  
  return(list(vt = vt, Kt = Kt, at = at[,,1:l], Ft = Ft, Lt = Lt, Pt = Pt[,,1:l]))  
}



################################################################################



########################### The Kalman smoother ################################

svks <- function(vt, Ft, Kt, Tt, Qt, Zt, Ht, a0, P0, yt){
  
  l <- length(yt)
  et <- numeric(l)
  epsilont <- numeric(l)
  rt <- array(dim = c(2,1,l))
  rt[1,1,l] <- 0
  rt[2,1,l] <- 0
  
  
  Dt <- numeric(l)
  Nt <- array(dim = c(2,2,l))
  Nt[,,l] <- 0
  
  Ct <- numeric(l)
  
  
  for (i in l:1){
    
    
    Zti <- t(matrix(Zt[,,i]))
    Kti <- matrix(Kt[,,i])
    Lti <- Tt - Kti%*%Zti
    Nti <- matrix(Nt[,,i], ncol = 2, nrow = 2, byrow = F)
    rti <- matrix(rt[,,i])
    
    et[i] <- (1 / Ft[i]) * vt[i] - t(Kti)%*%rti
    
    rt[,,i - 1] <- t(Zti) %*% (1 / Ft[i]) %*% vt[i] + t(Lti)%*%rti
    
    Dt[i] <- (1 / Ft[i]) + t(Kti)%*%Nti%*%Kti
    
    Nt[,,i - 1] <- t(Zti) %*% (1 / Ft[i])%*%Zti + t(Lti)%*%Nti%*%Lti
    
    epsilont[i] <- Ht[i] * et[i]
    
    Ct[i] <- Ht[i] - Ht[i] * Dt[i] * Ht[i]
    # Use nearest symmetric positive-definite matrix! Function in pracma.   
  }
  
  return(list(rt = rt, Nt = Nt, Dt = Dt, epsilont = epsilont,  Ct = Ct))
  
  
}

################################################################################

################ Disturbance simulation smoother ###############################
svdsim <- function(vt, Ft, Kt, Tt, Qt, Zt, Ht, a0, P0, yt){
  
  l <- length(yt)
  
  et <- numeric(l)
  epsilont <- numeric(l)
  rt <- array(dim = c(2,1,l))
  rt[1,1,l] <- 0
  rt[2,1,l] <- 0
  
  
  Dt <- numeric(l)
  Nt <- array(dim = c(2,2,l))
  Nt[,,l] <- 0
  
  Mt <- array(dim = c(1,2,l))
  Mt[,,l] <- 0
  
  
  Ct <- numeric(l)
  
  
  for (i in l:1){
    
    Zti <- t(matrix(Zt[,,i]))
    Kti <- matrix(Kt[,,i])
    Lti <- Tt - Kti%*%Zti
    
    
    Nti <- matrix(Nt[,,i], ncol = 2, nrow = 2)
    rti <- matrix(rt[,,i])
    
    
    et[i] <- (1 / Ft[i]) %*% vt[i] - t(Kti)%*%rti
    
    
    Dt[i] <- (1 / Ft[i]) + t(Kti)%*%Nti%*%Kti
    
    
    Mti <- Ht[i]%*%Dt[i]%*%Zti - Ht[i]%*%t(Kti)%*%Nti%*%Tt
    
    
    Ct[i] <- Ht[i] - Ht[i] %*% Dt[i] %*% Ht[i]
    
    
    uti <- rnorm(1, mean = 0, sd = sqrt(Ct[i]))
    
    rt[,,i - 1] <- t(Zti) %*% (1 / Ft[i]) %*% vt[i] - t(Mti)%*%(1 / Ct[i])%*%uti + t(Lti)%*%rti
    
    Nt[,,i - 1] <- t(Zti) %*% (1 / Ft[i])%*%Zti + t(Mti)%*%(1 / Ct[i])%*%Mti + t(Lti)%*%Nti%*%Lti
    
    epsilont[i] <-  Ht[i] %*% et[i] +  uti 
  }
  
  
  
  return(epsilont)
  
}


################################################################################



######################### Structural model  ####################################
######################### The Kalman filter ####################################


ssvkf <- function(Tt, Qt, Zt, Ht, Gt, a0, P0, yt){
  l <- length(yt)
  
  vt <- numeric(l)
  Kt <- array(dim = c(3,1,l))
  at <- array(dim = c(3,1,l + 1))
  at[,,1] <- a0
  Ft <- numeric(l)
  Lt <- array(dim = c(3,3,l))
  Pt <- array(dim = c(3,3,l + 1))
  Pt[,,1] <- P0
  
  
  
  for (i in 1:l){
    
    
    Qti <- Qt
    
    Zti <- t(matrix(Zt[,,i]))
    
    
    Pti <- matrix(Pt[,,i], nrow = 3, ncol = 3)
    
    vt[i] <- yt[i] - Zti%*%at[,,i]
    
    Ft[i] <- Zti%*%Pti%*%t(Zti) + Ht[i]
    
    Kti <- (Tt%*%Pti%*%t(Zti)) / Ft[i]
    
    Kt[,,i] <- Kti
    
    Lti <- Tt - Kti%*%Zti
    
    Lt[,,i] <- Lti
    
    at[,,i + 1] <- Tt%*%at[,,i] + Kti%*%vt[i]
    
    Pt[,,i + 1] <- Tt%*%Pti%*%t(Lti) + Qti 
    
  }
  
  return(list(vt = vt, Kt = Kt, at = at[,,1:l], Ft = Ft, Lt = Lt, Pt = Pt[,,1:l]))  
}

################################################################################

######################### The Kalman smoother ##################################

ssvks <- function(vt, Ft, Kt, Tt, Qt, Zt, Ht, a0, P0, yt){
  
  l <- length(yt)
  et <- numeric(l)
  epsilont <- numeric(l)
  rt <- array(dim = c(3,1,l))
  rt[1,1,l] <- 0
  rt[2,1,l] <- 0
  rt[3,1,l] <- 0
  
  
  Dt <- numeric(l)
  Nt <- array(dim = c(3,3,l))
  Nt[,,l] <- 0
  
  Ct <- numeric(l)
  
  
  for (i in l:1){
    
    
    Zti <- t(matrix(Zt[,,i]))
    Kti <- matrix(Kt[,,i])
    Lti <- Tt - Kti%*%Zti
    Nti <- matrix(Nt[,,i], ncol = 3, nrow = 3, byrow = F)
    rti <- matrix(rt[,,i])
    
    et[i] <- (1 / Ft[i]) * vt[i] - t(Kti)%*%rti
    
    rt[,,i - 1] <- t(Zti) %*% (1 / Ft[i]) %*% vt[i] + t(Lti)%*%rti
    
    Dt[i] <- (1 / Ft[i]) + t(Kti)%*%Nti%*%Kti
    
    Nt[,,i - 1] <- t(Zti) %*% (1 / Ft[i])%*%Zti + t(Lti)%*%Nti%*%Lti
    
    epsilont[i] <- Ht[i] * et[i]
    
    Ct[i] <- Ht[i] - Ht[i] * Dt[i] * Ht[i]
    
  }
  
  return(list(rt = rt, Nt = Nt, Dt = Dt, epsilont = epsilont,  Ct = Ct))
  
  
}

################################################################################

################## Disturbance simulation smoother #############################

ssvdsim <- function(vt, Ft, Kt, Tt, Qt, Zt, Ht, a0, P0, yt){
  
  l <- length(yt)
  
  et <- numeric(l)
  epsilont <- numeric(l)
  rt <- array(dim = c(3,1,l))
  rt[1,1,l] <- 0
  rt[2,1,l] <- 0
  rt[3,1,l] <- 0
  
  
  Dt <- numeric(l)
  Nt <- array(dim = c(3,3,l))
  Nt[,,l] <- 0
  
  Mt <- array(dim = c(1,3,l))
  Mt[,,l] <- 0
  
  
  Ct <- numeric(l)
  
  
  for (i in l:1){
    
    Zti <- t(matrix(Zt[,,i]))
    Kti <- matrix(Kt[,,i])
    Lti <- Tt - Kti%*%Zti
    
    
    Nti <- matrix(Nt[,,i], ncol = 3, nrow = 3)
    rti <- matrix(rt[,,i])
    
    
    et[i] <- (1 / Ft[i]) %*% vt[i] - t(Kti)%*%rti
    
    
    Dt[i] <- (1 / Ft[i]) + t(Kti)%*%Nti%*%Kti
    
    
    Mti <- Ht[i]%*%Dt[i]%*%Zti - Ht[i]%*%t(Kti)%*%Nti%*%Tt
    
    
    Ct[i] <- Ht[i] - Ht[i] %*% Dt[i] %*% Ht[i]
    
    
    uti <- rnorm(1, mean = 0, sd = sqrt(Ct[i]))
    
    rt[,,i - 1] <- t(Zti) %*% (1 / Ft[i]) %*% vt[i] - t(Mti)%*%(1 / Ct[i])%*%uti + t(Lti)%*%rti
    
    Nt[,,i - 1] <- t(Zti) %*% (1 / Ft[i])%*%Zti + t(Mti)%*%(1 / Ct[i])%*%Mti + t(Lti)%*%Nti%*%Lti
    
    epsilont[i] <-  Ht[i] %*% et[i] +  uti 
  }
  
  
  
  return(epsilont)
  
}

################################################################################

######################## Asymmetric model without covariates ###################

########################## The Kalman filter ###################################

asysvkf <- function(Tt, Qt, Zt, Ht, Gt, ct, a0, P0, yt){
  l <- length(yt)
  
  vt <- numeric(l)
  Kt <- array(dim = c(2,1,l))
  at <- array(dim = c(2,1,l + 1))
  at[,,1] <- a0
  Ft <- numeric(l)
  Lt <- array(dim = c(2,2,l))
  Pt <- array(dim = c(2,2,l + 1))
  Pt[,,1] <- P0
  
  
  
  for (i in 1:l){
    
    
    Gti <- matrix(Gt[,,i])
    Qti <- Qt
    
    Zti <- t(matrix(Zt[,,i]))
    
    cti <- matrix(ct[,,i])
    
    Pti <- matrix(Pt[,,i], nrow = 2, ncol = 2)
    
    vt[i] <- yt[i] - Zti%*%at[,,i]
    
    Ft[i] <- Zti%*%Pti%*%t(Zti) + Ht[i]
    
    Kti <- (Tt%*%Pti%*%t(Zti) + Gti) / Ft[i]
    
    Kt[,,i] <- Kti
    
    Lti <- Tt - Kti%*%Zti
    
    Lt[,,i] <- Lti
    
    at[,,i + 1] <- Tt%*%at[,,i] + cti + Kti%*%vt[i]
    
    Pt[,,i + 1] <- Tt%*%Pti%*%t(Lti) + Qti - Gti%*%t(Kti)
    
  }
  
  return(list(vt = vt, Kt = Kt, at = at[,,1:l], Ft = Ft, Lt = Lt, 
              Pt = Pt[,,1:l]))  
}

################################################################################

######################### Asymmetric model with covariates #####################

########################## The Kalman filter ###################################

asysvkfx <- function(Tt, Qt, Zt, Ht, Gt, ct, a0, P0, yt){
  l <- length(yt)
  
  vt <- numeric(l)
  Kt <- array(dim = c(3,1,l))
  at <- array(dim = c(3,1,l + 1))
  at[,,1] <- a0
  Ft <- numeric(l)
  Lt <- array(dim = c(3,3,l))
  Pt <- array(dim = c(3,3,l + 1))
  Pt[,,1] <- P0
  
  
  
  for (i in 1:l){
    
    
    Gti <- matrix(Gt[,,i])
    Qti <- Qt
    
    Zti <- t(matrix(Zt[,,i]))
    
    cti <- matrix(ct[,,i])
    
    Pti <- matrix(Pt[,,i], nrow = 3, ncol = 3)
    
    vt[i] <- yt[i] - Zti%*%at[,,i]
    
    Ft[i] <- Zti%*%Pti%*%t(Zti) + Ht[i]
    
    Kti <- (Tt%*%Pti%*%t(Zti) + Gti) / Ft[i]
    
    Kt[,,i] <- Kti
    
    Lti <- Tt - Kti%*%Zti
    
    Lt[,,i] <- Lti
    
    at[,,i + 1] <- Tt%*%at[,,i] + cti + Kti%*%vt[i]
    
    Pt[,,i + 1] <- Tt%*%Pti%*%t(Lti) + Qti - Gti%*%t(Kti)
    
  }
  
  return(list(vt = vt, Kt = Kt, at = at[,,1:l], Ft = Ft, Lt = Lt, 
              Pt = Pt[,,1:l]))  
}

################################################################################
