############################################################
# REVISED TABLE 7
# UNIVERSAL INFERENCE E-VALUES
# Section 2.4 Implementation
#
# Option 3:
# Evaluate entire D2 under both alternative fits
#
# Outputs:
# Mean(EUI)
# Median(EUI)
# P(EUI >= 20)
# P(EUI >= 100)
############################################################

rm(list=ls())

set.seed(2025)

############################################################
# 1. CDF
############################################################

CDF_y <- function(y, alpha, beta, lambda1, lambda2){
  
  A <- 1 - lambda1*(lambda2 - 1)
  H <- lambda1*(2*lambda2 - 1)
  K <- lambda1*lambda2
  
  t <- pmax(1 - y^(-alpha),0)
  
  A*t^beta + H*t^(2*beta) - K*t^(3*beta)
}

############################################################
# 2. IIT SAMPLER
############################################################

IIT_fast <- function(n, alpha, beta, lambda1, lambda2){
  
  ys <- c(
    seq(1.0001,10,length=400),
    seq(10.01,100,length=400),
    exp(seq(log(100.1),log(1e4),length=400))
  )
  
  Fs <- CDF_y(
    ys,
    alpha,beta,
    lambda1,lambda2
  )
  
  Fs <- cummax(Fs)
  Fs <- pmin(pmax(Fs,0),1)
  
  u <- runif(n)
  
  approx(
    Fs,
    ys,
    xout=u,
    ties="ordered"
  )$y
}

###########
# 3. LOGLIK
###########

loglik_segment <- function(y, p){
  
  alpha <- p[1]
  beta  <- p[2]
  l1    <- p[3]
  l2    <- p[4]
  
  if(alpha<=0 || beta<=0) return(-Inf)
  if(any(y<=0)) return(-Inf)
  
  t <- 1 - y^(-alpha)
  
  if(any(t<=0)) return(-Inf)
  
  A <- 1 - l1*(l2 - 1)
  H <- l1*(2*l2 - 1)
  K <- l1*l2
  
  term <- A*t^(beta-1) +
    2*H*t^(2*beta-1) -
    3*K*t^(3*beta-1)
  
  if(any(term<=0)) return(-Inf)
  
  sum(
    log(alpha*beta) -
      (alpha+1)*log(y) +
      log(term)
  )
}
##########################
# 4. MLE
# MULTI-START OPTIMIZATION
##########################

mle_segment <- function(y){
  
  ##########################
  # Multiple starting values
  ##########################
  
  starts <- list(
    
    c(3,2,-0.3,1.1),
    
    c(2,1,0.2,0.8),
    
    c(5,3,0.5,1.5),
    
    c(1,1,-0.5,0.5),
    
    c(8,5,-0.2,1.8)
  )
  
  ########################
  # Negative log-likelihood
  #########################
  
  negll <- function(p){
    
    ll <- loglik_segment(y,p)
    
    if(is.finite(ll))
      return(-ll)
    
    1e12
  }
  
  ################
  # Store best fit
  ################
  
  best_fit <- NULL
  
  best_ll <- -Inf
  
  ##########################
  # Multi-start optimization
  ##########################
  
  for(st in starts){
    
    fit <- tryCatch(
      
      optim(
        par = st,
        fn = negll,
        method = "L-BFGS-B",
        
        lower = c(
          1e-6,
          1e-6,
          -1,
          0
        ),
        
        upper = c(
          100,
          100,
          1,
          2
        )
      ),
      
      error=function(e) NULL
    )
    
    #################
    # Update best fit
    #################
    
    if(!is.null(fit)){
      
      ll <- -fit$value
      
      if(is.finite(ll) &&
         ll > best_ll){
        
        best_ll <- ll
        
        best_fit <- fit
      }
    }
  }
  
  #############
  # No valid fit
  ##############
  
  if(is.null(best_fit)){
    
    return(
      list(
        par = rep(NA,4),
        ll = -Inf
      )
    )
  }
  
  #################
  # Return best fit
  #################
  
  list(
    par = best_fit$par,
    ll  = best_ll
  )
}

##########################
# 5. D1 CHANGEPOINT SEARCH
##########################

find_cp_D1 <- function(y, lmin=max(5,floor(0.1*length(y)))){
  
  n <- length(y)
  
  ll0 <- mle_segment(y)$ll
  
  bestT <- -Inf
  bestk <- NULL
  
  for(k in lmin:(n-lmin)){
    
    ll1 <- mle_segment(y[1:k])$ll +
      mle_segment(y[(k+1):n])$ll
    
    T <- 2*(ll1-ll0)
    
    if(is.finite(T) && T > bestT){
      
      bestT <- T
      bestk <- k
    }
  }
  
  bestk
}
#######################################
# 6. STRICT UNIVERSAL INFERENCE E-VALUE
#######################################

UI_evalue <- function(y){
  
  n <- length(y)
  
  m <- floor(n/2)
  
  ############
  # Split data
  ############
  
  D1 <- y[1:m]
  
  D2 <- y[(m+1):n]
  
  ################################
  # STEP 1:
  # Find changepoint using D1 only
  ################################
  
  khat <- find_cp_D1(D1)
  
  if(is.null(khat))
    return(NA)
  
  ##################################
  # STEP 2:
  # Alternative estimated ONLY on D1
  ##################################
  
  fit_left <- mle_segment(
    D1[1:khat]
  )
  
  fit_right <- mle_segment(
    D1[(khat+1):m]
  )
  
  ###########################
  # STEP 3:
  # Null estimated ONLY on D2
  ###########################
  
  fit_null <- mle_segment(D2)
  
  ###############################
  # STEP 4:
  # Transfer relative CP location
  ###############################
  
  rel_cp <- khat / m
  
  n2 <- length(D2)
  
  kD2 <- round(rel_cp * n2)
  
  kD2 <- max(
    floor(0.10*n2),
    min(
      n2-floor(0.10*n2),
      kD2
    )
  )
  
  ################################
  # STEP 5:
  # Evaluate D2 using D1 estimates
  ################################
  
  ll_alt <-
    
    loglik_segment(
      D2[1:kD2],
      fit_left$par
    ) +
    
    loglik_segment(
      D2[(kD2+1):n2],
      fit_right$par
    )
  
  #######################
  # Null likelihood on D2
  #######################
  
  ll_null <- fit_null$ll
  
  ##################
  # Valid UI E-value
  ##################
  
  if(!is.finite(ll_alt) ||
     !is.finite(ll_null))
    return(NA)
  
  EUI <- exp(
    ll_alt - ll_null
  )
  
  return(EUI)
}

####################
# 7. SIGNAL SETTINGS
####################

params1 <- c(
  5.6,
  3.8,
  -0.36,
  1.12
)

signal_list <- list(
  
  Weak = c(
    6.2,
    4.0,
    -0.40,
    1.20
  ),
  
  Moderate = c(
    8.4,
    5.7,
    -0.54,
    1.68
  ),
  
  Strong = c(
    12.0,
    8.0,
    -0.70,
    1.90
  )
)
#######################
# 8. SIMULATION
#######################

simulate_UI <- function(
    n,
    signal_par,
    case_type = "Balanced",
    R = 2000){
  
  E <- numeric(R)
  
  ########################
  # Change-point locations
  ########################
  
  if(case_type == "Balanced"){
    
    cp <- floor(n/2)
    
  } else if(case_type == "Unbalanced"){
    
    cp <- floor(0.30*n)
    
  } else if(case_type == "H0"){
    
    cp <- floor(n/2)
  }
  
  #################
  # Simulation loop
  #################
  
  for(r in 1:R){
    
    #########
    # H0 CASE
    #########
    
    if(case_type == "H0"){
      
      y <- IIT_fast(
        n,
        params1[1],
        params1[2],
        params1[3],
        params1[4]
      )
      
    } else {
      
      ##################################
      # Balanced / Unbalanced alternatives
      #################################
      
      y <- c(
        
        IIT_fast(
          cp,
          params1[1],
          params1[2],
          params1[3],
          params1[4]
        ),
        
        IIT_fast(
          n-cp,
          signal_par[1],
          signal_par[2],
          signal_par[3],
          signal_par[4]
        )
      )
    }
    
    E[r] <- UI_evalue(y)
  }
  
  ######################
  # Remove invalid values
  ######################
  
  E <- E[is.finite(E)]
  
  #################
  # Summary measures
  ##################
  
  c(
    
    Mean_EUI =
      mean(E),
    
    Median_EUI =
      median(E),
    
    P_EUI_20 =
      mean(E >= 20),
    
    P_EUI_100 =
      mean(E >= 100)
  )
}

#########
# Results
#########

ALL_RESULTS <- data.frame()

for(n_current in c(50,100,200)){
  
  cat("\n")
  cat("=========================================\n")
  cat("RUNNING SAMPLE SIZE =", n_current, "\n")
  cat("=========================================\n")
  
  RESULTS <- data.frame()
  
  for(sig in names(signal_list)){
    
    for(case_now in c(
      "H0",
      "Balanced",
      "Unbalanced"
    )){
      
      cat(
        "n =", n_current,
        "| Signal =", sig,
        "| Case =", case_now,
        "\n"
      )
      
      out <- simulate_UI(
        n = n_current,
        signal_par = signal_list[[sig]],
        case_type = case_now,
        R = 1000
      )
      
      RESULTS <- rbind(
        RESULTS,
        data.frame(
          Sample_Size = n_current,
          Signal = sig,
          Case = case_now,
          Mean_EUI   = round(out[1],3),
          Median_EUI = round(out[2],3),
          P_EUI_20   = round(out[3],3),
          P_EUI_100  = round(out[4],3)
        )
      )
    }
  }
  
  #############################
  # PRINT RESULTS FOR CURRENT n
  #############################
  
  cat("\n")
  cat("=========================================\n")
  cat("RESULTS FOR n =", n_current, "\n")
  cat("=========================================\n")
  
  print(RESULTS, row.names = FALSE)
  
  ######################
  # SAVE INDIVIDUAL FILE
  ######################
  
  write.csv(
    RESULTS,
    paste0(
      "UI_Evalues_n_",
      n_current,
      ".csv"
    ),
    row.names = FALSE
  )
  
  ###########################
  # APPEND TO OVERALL RESULTS
  ###########################
  
  ALL_RESULTS <- rbind(
    ALL_RESULTS,
    RESULTS
  )
}

########################
# FINAL COMBINED RESULTS
########################

cat("\n")
cat("=========================================\n")
cat("FINAL COMBINED RESULTS\n")
cat("=========================================\n")

print(ALL_RESULTS,row.names=FALSE)

write.csv(
  ALL_RESULTS,
  "Revised_Table7_UI_Evalues.csv",
  row.names = FALSE
)
#############
# 10. TABLE 
############

cat("\n\n")
cat("=========================================\n")
cat("REVISED TABLE 7\n")
cat("UNIVERSAL INFERENCE E-VALUES\n")
cat("=========================================\n\n")

print(RESULTS,row.names=FALSE)

write.csv(
  RESULTS,
  "Revised_Table7_UI_Evalues.csv",
  row.names=FALSE
)

cat(
  "\n\nCSV file saved as:\n",
  "Revised_Table7_UI_Evalues.csv\n"
)