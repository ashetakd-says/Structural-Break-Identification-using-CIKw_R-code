# ============================================================
# FINAL UNIFIED SIMULATION
# CIKw PARAMETER-SHIFT + VARIANCE-SHIFT BENCHMARK
#
# Methods:
#   GLRT + CUSUM + BIC + MIC
#
# Scenarios:
#   1. CIKw parameter shifts:
#      Weak / Moderate / Strong
#      Balanced and Unbalanced change-points
#
#   2. Variance-shift benchmark:
#      Variance ratios = 1, 1.25, 1.5, 2, 3, 4
#
# Sample sizes:
#   n = 50, 100, 200
# ============================================================

rm(list = ls())
set.seed(2025)


# ============================================================
# 1) CDF + SAMPLER
# ============================================================

CDF_y <- function(y, a, b, l1, l2){
  
  A <- 1 - l1 * (l2 - 1)
  H <- l1 * (2 * l2 - 1)
  K <- l1 * l2
  
  t <- pmax(1 - y^(-a), 0)
  
  A * t^b +
    H * t^(2 * b) -
    K * t^(3 * b)
}


IIT_fast <- function(n, a, b, l1, l2){
  
  grid <- c(
    seq(1.0001, 10, length = 300),
    seq(10, 100, length = 300),
    exp(seq(log(100), log(1e4), length = 300))
  )
  
  F <- CDF_y(
    grid,
    a,
    b,
    l1,
    l2
  )
  
  F <- pmin(cummax(F), 1)
  
  u <- runif(n)
  
  approx(
    F,
    grid,
    u
  )$y
}


# ============================================================
# 2) LOG-LIKELIHOOD
# ============================================================

loglik <- function(y, p){
  
  a  <- p[1]
  b  <- p[2]
  l1 <- p[3]
  l2 <- p[4]
  
  if(a <= 0 || b <= 0)
    return(-Inf)
  
  t <- 1 - y^(-a)
  
  if(any(t <= 0))
    return(-Inf)
  
  A <- 1 - l1 * (l2 - 1)
  H <- l1 * (2 * l2 - 1)
  K <- l1 * l2
  
  term <-
    A * t^(b - 1) +
    2 * H * t^(2 * b - 1) -
    3 * K * t^(3 * b - 1)
  
  if(any(term <= 0))
    return(-Inf)
  
  sum(
    log(a * b) -
      (a + 1) * log(y) +
      log(term)
  )
}


# ============================================================
# 3) MLE
# ============================================================

mle <- function(y){
  
  starts <- list(
    c(3, 2, -0.3, 1.1),
    c(5, 2, -0.5, 1.5)
  )
  
  best <- -Inf
  
  for(s in starts){
    
    res <- try(
      optim(
        s,
        function(p) -loglik(y, p),
        method = "L-BFGS-B",
        lower = c(
          1e-6,
          1e-6,
          -0.99,
          0.001
        ),
        upper = c(
          Inf,
          Inf,
          0.99,
          2
        )
      ),
      silent = TRUE
    )
    
    if(class(res) != "try-error"){
      
      ll <- loglik(
        y,
        res$par
      )
      
      if(
        is.finite(ll) &&
        ll > best
      ){
        best <- ll
      }
    }
  }
  
  best
}


# ============================================================
# 4) GLRT STATISTIC
# ============================================================

GLRT_stat <- function(
    y,
    lmin = 10
){
  
  n <- length(y)
  
  ll0 <- mle(y)
  
  if(!is.finite(ll0))
    return(0)
  
  Tmax <- -Inf
  
  for(k in lmin:(n - lmin)){
    
    ll1 <-
      mle(y[1:k]) +
      mle(y[(k + 1):n])
    
    if(is.finite(ll1)){
      
      Tmax <- max(
        Tmax,
        2 * (ll1 - ll0)
      )
    }
  }
  
  if(!is.finite(Tmax))
    Tmax <- 0
  
  Tmax
}


# ============================================================
# OPTIONAL GLRT LAMBDA
# ============================================================

GLRT_lambda <- function(
    y,
    lmin = 10
){
  
  n <- length(y)
  
  ll0 <- mle(y)
  
  if(!is.finite(ll0))
    return(0)
  
  best <- -Inf
  
  for(k in lmin:(n - lmin)){
    
    ll1 <-
      mle(y[1:k]) +
      mle(y[(k + 1):n])
    
    if(is.finite(ll1)){
      
      best <- max(
        best,
        ll1 - ll0
      )
    }
  }
  
  if(!is.finite(best))
    best <- 0
  
  best
}


# ============================================================
# 5) CUSUM STATISTIC
# ============================================================

CUSUM_stat <- function(y){
  
  n <- length(y)
  
  mu <- mean(y)
  
  s <- cumsum(
    y - mu
  )
  
  max(abs(s)) /
    (sd(y) * sqrt(n))
}


# ============================================================
# 6) BIC / MIC
# ============================================================

BIC_MIC <- function(
    y,
    c = 0.1
){
  
  n <- length(y)
  
  p0 <- 4
  p1 <- 8
  
  ll0 <- mle(y)
  
  B0 <-
    -2 * ll0 +
    p0 * log(n)
  
  M0 <- B0
  
  bestB <- Inf
  bestM <- Inf
  
  for(k in 10:(n - 10)){
    
    ll1 <-
      mle(y[1:k]) +
      mle(y[(k + 1):n])
    
    if(!is.finite(ll1))
      next
    
    B1 <-
      -2 * ll1 +
      p1 * log(n)
    
    w <-
      c *
      p1 *
      (
        1/k +
          1/(n-k)
      )
    
    M1 <-
      -2 * ll1 +
      (p1 + w) * log(n)
    
    bestB <- min(
      bestB,
      B1
    )
    
    bestM <- min(
      bestM,
      M1
    )
  }
  
  c(
    bestB < B0,
    bestM < M0
  )
}


# ============================================================
# 7) BOOTSTRAP CALIBRATION
# ============================================================

calibrate_glrt <- function(
    n,
    B = 200
){
  
  vals <- numeric(B)
  
  for(i in 1:B){
    
    y <- IIT_fast(
      n,
      3,
      2,
      0,
      1
    )
    
    vals[i] <-
      GLRT_stat(y)
  }
  
  quantile(
    vals,
    0.95
  )
}


calibrate_cusum <- function(
    n,
    B = 200
){
  
  vals <- numeric(B)
  
  for(i in 1:B){
    
    y <- IIT_fast(
      n,
      3,
      2,
      0,
      1
    )
    
    vals[i] <-
      CUSUM_stat(y)
  }
  
  quantile(
    vals,
    0.95
  )
}


# ============================================================
# 8) MCP DETECTION
# ============================================================

detect_cp <- function(
    y,
    C_glrt
){
  
  T <- GLRT_stat(y)
  
  T > C_glrt
}


# ============================================================
# 8A) CIKw SIGNAL PARAMETERS
# ============================================================

params1 <- c(
  5.6,
  3.8,
  -0.36,
  1.12
)


params_weak <- c(
  6.2,
  4.0,
  -0.40,
  1.20
)


params_mod <- c(
  8.4,
  5.7,
  -0.54,
  1.68
)


params_str <- c(
  12.0,
  8.0,
  -0.70,
  1.90
)


signal_list <- list(
  
  Weak = params_weak,
  
  Moderate = params_mod,
  
  Strong = params_str
)


# ============================================================
# 8B) VARIANCE-SHIFT BENCHMARK GENERATOR
#
# Mean is unchanged.
# Only the variance changes after tau.
# ============================================================

generate_variance_shift <- function(
    n,
    tau = floor(n / 2),
    mu = 5,
    sigma1 = 0.5,
    variance_ratio = 2
){
  
  sigma2 <-
    sigma1 *
    sqrt(variance_ratio)
  
  y1 <- rnorm(
    tau,
    mean = mu,
    sd = sigma1
  )
  
  y2 <- rnorm(
    n - tau,
    mean = mu,
    sd = sigma2
  )
  
  c(
    y1,
    y2
  )
}


# ============================================================
# 8C) VARIANCE-SHIFT NULL CALIBRATION
# ============================================================

calibrate_variance_benchmark <- function(
    n,
    R = 200
){
  
  glrt_vals <- numeric(R)
  cusum_vals <- numeric(R)
  
  for(i in 1:R){
    
    # H0: NO variance change
    y <- generate_variance_shift(
      n = n,
      tau = floor(n / 2),
      variance_ratio = 1
    )
    
    glrt_vals[i] <-
      GLRT_stat(y)
    
    cusum_vals[i] <-
      CUSUM_stat(y)
  }
  
  list(
    
    C_glrt =
      quantile(
        glrt_vals,
        0.95
      ),
    
    C_cusum =
      quantile(
        cusum_vals,
        0.95
      )
  )
}


# ============================================================
# 9) UNIFIED SIMULATION
#
# CIKw PARAMETER-SHIFT
# +
# VARIANCE-SHIFT BENCHMARK
# ============================================================

simulate_all <- function(
    n,
    R = 200,
    B = 200,
    variance_ratios =
      c(
        1,
        1.25,
        1.5,
        2,
        3,
        4
      )
){
  
  cat(
    "\n====================================================\n"
  )
  
  cat(
    "SAMPLE SIZE:",
    n,
    "\n"
  )
  
  cat(
    "====================================================\n"
  )
  
  
  # ==========================================================
  # A. CALIBRATE CIKw GLRT AND CUSUM
  # ==========================================================
  
  cat(
    "\nCalibrating GLRT threshold...\n"
  )
  
  C_glrt <-
    calibrate_glrt(
      n,
      B
    )
  
  
  cat(
    "Calibrating CUSUM threshold...\n"
  )
  
  C_cusum <-
    calibrate_cusum(
      n,
      B
    )
  
  
  cat(
    "\nCritical values:\n"
  )
  
  cat(
    "GLRT  =",
    C_glrt,
    "\n"
  )
  
  cat(
    "CUSUM =",
    C_cusum,
    "\n"
  )
  
  
  # ==========================================================
  # B. CIKw PARAMETER-SHIFT SIMULATION
  # ==========================================================
  
  signal_results <- data.frame()
  
  evalue_results <- data.frame()
  
  
  for(signal_name in names(signal_list)){
    
    cat(
      "\n----------------------------------------------------\n"
    )
    
    cat(
      "Signal:",
      signal_name,
      "\n"
    )
    
    cat(
      "----------------------------------------------------\n"
    )
    
    
    params2 <-
      signal_list[[signal_name]]
    
    
    # --------------------------------------------------------
    # Storage
    # --------------------------------------------------------
    
    GLRT_H0  <- numeric(R)
    GLRT_BAL <- numeric(R)
    GLRT_UNB <- numeric(R)
    
    CUSUM_H0  <- numeric(R)
    CUSUM_BAL <- numeric(R)
    CUSUM_UNB <- numeric(R)
    
    BIC_H0  <- numeric(R)
    BIC_BAL <- numeric(R)
    BIC_UNB <- numeric(R)
    
    MIC_H0  <- numeric(R)
    MIC_BAL <- numeric(R)
    MIC_UNB <- numeric(R)
    
    E_H0  <- numeric(R)
    E_BAL <- numeric(R)
    E_UNB <- numeric(R)
    
    
    # ========================================================
    # MONTE CARLO LOOP
    # ========================================================
    
    for(i in 1:R){
      
      
      # ------------------------------------------------------
      # H0: NO CHANGE
      # ------------------------------------------------------
      
      y0 <- IIT_fast(
        n,
        3,
        2,
        0,
        1
      )
      
      
      # ------------------------------------------------------
      # BALANCED CHANGE
      # tau = 0.5n
      # ------------------------------------------------------
      
      k_bal <-
        floor(n / 2)
      
      y_bal <- c(
        
        IIT_fast(
          k_bal,
          params1[1],
          params1[2],
          params1[3],
          params1[4]
        ),
        
        IIT_fast(
          n - k_bal,
          params2[1],
          params2[2],
          params2[3],
          params2[4]
        )
      )
      
      
      # ------------------------------------------------------
      # UNBALANCED CHANGE
      # tau = 0.1n
      # ------------------------------------------------------
      
      k_unb <-
        max(
          floor(0.1 * n),
          10
        )
      
      y_unb <- c(
        
        IIT_fast(
          k_unb,
          params1[1],
          params1[2],
          params1[3],
          params1[4]
        ),
        
        IIT_fast(
          n - k_unb,
          params2[1],
          params2[2],
          params2[3],
          params2[4]
        )
      )
      
      
      # ======================================================
      # H0
      # ======================================================
      
      T0 <-
        GLRT_stat(y0)
      
      C0 <-
        CUSUM_stat(y0)
      
      BM0 <-
        BIC_MIC(y0)
      
      
      GLRT_H0[i] <-
        as.numeric(
          T0 > C_glrt
        )
      
      CUSUM_H0[i] <-
        as.numeric(
          C0 > C_cusum
        )
      
      BIC_H0[i] <-
        BM0[1]
      
      MIC_H0[i] <-
        BM0[2]
      
      E_H0[i] <-
        exp(
          min(
            T0 / 2,
            700
          )
        )
      
      
      # ======================================================
      # BALANCED
      # ======================================================
      
      TB <-
        GLRT_stat(y_bal)
      
      CB <-
        CUSUM_stat(y_bal)
      
      BMB <-
        BIC_MIC(y_bal)
      
      
      GLRT_BAL[i] <-
        as.numeric(
          TB > C_glrt
        )
      
      CUSUM_BAL[i] <-
        as.numeric(
          CB > C_cusum
        )
      
      BIC_BAL[i] <-
        BMB[1]
      
      MIC_BAL[i] <-
        BMB[2]
      
      E_BAL[i] <-
        exp(
          min(
            TB / 2,
            700
          )
        )
      
      
      # ======================================================
      # UNBALANCED
      # ======================================================
      
      TU <-
        GLRT_stat(y_unb)
      
      CU <-
        CUSUM_stat(y_unb)
      
      BMU <-
        BIC_MIC(y_unb)
      
      
      GLRT_UNB[i] <-
        as.numeric(
          TU > C_glrt
        )
      
      CUSUM_UNB[i] <-
        as.numeric(
          CU > C_cusum
        )
      
      BIC_UNB[i] <-
        BMU[1]
      
      MIC_UNB[i] <-
        BMU[2]
      
      E_UNB[i] <-
        exp(
          min(
            TU / 2,
            700
          )
        )
    }
    
    
    # ========================================================
    # STORE CIKw RESULTS
    # ========================================================
    
    signal_results <-
      rbind(
        
        signal_results,
        
        data.frame(
          n = n,
          Signal = signal_name,
          Scenario = "H0",
          GLRT = mean(GLRT_H0),
          CUSUM = mean(CUSUM_H0),
          BIC = mean(BIC_H0),
          MIC = mean(MIC_H0)
        ),
        
        data.frame(
          n = n,
          Signal = signal_name,
          Scenario = "Balanced",
          GLRT = mean(GLRT_BAL),
          CUSUM = mean(CUSUM_BAL),
          BIC = mean(BIC_BAL),
          MIC = mean(MIC_BAL)
        ),
        
        data.frame(
          n = n,
          Signal = signal_name,
          Scenario = "Unbalanced",
          GLRT = mean(GLRT_UNB),
          CUSUM = mean(CUSUM_UNB),
          BIC = mean(BIC_UNB),
          MIC = mean(MIC_UNB)
        )
      )
    
    
    # ========================================================
    # E-VALUE SUMMARY
    # ========================================================
    
    evalue_results <-
      rbind(
        
        evalue_results,
        
        data.frame(
          n = n,
          Signal = signal_name,
          Scenario = "H0",
          Mean_E = mean(E_H0),
          Median_E = median(E_H0),
          Q95_E =
            quantile(
              E_H0,
              0.95
            )
        ),
        
        data.frame(
          n = n,
          Signal = signal_name,
          Scenario = "Balanced",
          Mean_E = mean(E_BAL),
          Median_E = median(E_BAL),
          Q95_E =
            quantile(
              E_BAL,
              0.95
            )
        ),
        
        data.frame(
          n = n,
          Signal = signal_name,
          Scenario = "Unbalanced",
          Mean_E = mean(E_UNB),
          Median_E = median(E_UNB),
          Q95_E =
            quantile(
              E_UNB,
              0.95
            )
        )
      )
    
    
    # ========================================================
    # PRINT CIKw RESULTS
    # ========================================================
    
    cat(
      "\nResults:",
      signal_name,
      "\n"
    )
    
    cat(
      "Type I error:",
      round(
        mean(GLRT_H0),
        3
      ),
      "\n"
    )
    
    cat(
      "GLRT Power (balanced):",
      round(
        mean(GLRT_BAL),
        3
      ),
      "\n"
    )
    
    cat(
      "GLRT Power (unbalanced):",
      round(
        mean(GLRT_UNB),
        3
      ),
      "\n"
    )
    
    cat(
      "CUSUM Power (balanced):",
      round(
        mean(CUSUM_BAL),
        3
      ),
      "\n"
    )
    
    cat(
      "CUSUM Power (unbalanced):",
      round(
        mean(CUSUM_UNB),
        3
      ),
      "\n"
    )
  }
  
  
  # ==========================================================
  # C. VARIANCE-SHIFT BENCHMARK
  # ==========================================================
  
  cat(
    "\n====================================================\n"
  )
  
  cat(
    "VARIANCE-SHIFT BENCHMARK | n =",
    n,
    "\n"
  )
  
  cat(
    "====================================================\n"
  )
  
  
  # ----------------------------------------------------------
  # Calibrate variance benchmark
  # ----------------------------------------------------------
  
  var_thresholds <-
    calibrate_variance_benchmark(
      n = n,
      R = R
    )
  
  C_glrt_var <-
    var_thresholds$C_glrt
  
  C_cusum_var <-
    var_thresholds$C_cusum
  
  
  variance_results <-
    data.frame()
  
  
  # ----------------------------------------------------------
  # Variance ratios
  # ----------------------------------------------------------
  
  for(vr in variance_ratios){
    
    GLRT_rej <- 0
    CUSUM_rej <- 0
    BIC_rej <- 0
    MIC_rej <- 0
    
    
    for(i in 1:R){
      
      y <-
        generate_variance_shift(
          n = n,
          tau = floor(n / 2),
          variance_ratio = vr
        )
      
      
      # GLRT
      T <-
        GLRT_stat(y)
      
      if(T > C_glrt_var)
        GLRT_rej <-
        GLRT_rej + 1
      
      
      # CUSUM
      C <-
        CUSUM_stat(y)
      
      if(C > C_cusum_var)
        CUSUM_rej <-
        CUSUM_rej + 1
      
      
      # BIC / MIC
      BM <-
        BIC_MIC(y)
      
      BIC_rej <-
        BIC_rej + BM[1]
      
      MIC_rej <-
        MIC_rej + BM[2]
    }
    
    
    variance_results <-
      rbind(
        
        variance_results,
        
        data.frame(
          n = n,
          Variance_Ratio = vr,
          GLRT = GLRT_rej / R,
          CUSUM = CUSUM_rej / R,
          BIC = BIC_rej / R,
          MIC = MIC_rej / R
        )
      )
    
    
    cat(
      "Variance ratio =",
      vr,
      "| GLRT =",
      round(
        GLRT_rej / R,
        3
      ),
      "| CUSUM =",
      round(
        CUSUM_rej / R,
        3
      ),
      "| BIC =",
      round(
        BIC_rej / R,
        3
      ),
      "| MIC =",
      round(
        MIC_rej / R,
        3
      ),
      "\n"
    )
  }
  
  
  # ==========================================================
  # RETURN ALL RESULTS
  # ==========================================================
  
  list(
    
    signal_results =
      signal_results,
    
    evalue_results =
      evalue_results,
    
    variance_results =
      variance_results,
    
    C_glrt =
      C_glrt,
    
    C_cusum =
      C_cusum,
    
    C_glrt_variance =
      C_glrt_var,
    
    C_cusum_variance =
      C_cusum_var
  )
}


# ============================================================
# 10) RUN FOR ALL SAMPLE SIZES
# ============================================================

sample_sizes <-
  c(
    50,
    100,
    200
  )


all_results <-
  lapply(
    sample_sizes,
    function(n){
      
      simulate_all(
        n = n,
        R = 200,
        B = 200,
        variance_ratios =
          c(
            1,
            1.25,
            1.5,
            2,
            3,
            4
          )
      )
    }
  )


names(all_results) <-
  paste0(
    "n_",
    sample_sizes
  )


# ============================================================
# 11) COMBINE CIKw SIGNAL RESULTS
# ============================================================

signal_results_all <-
  do.call(
    rbind,
    lapply(
      all_results,
      function(x)
        x$signal_results
    )
  )


# ============================================================
# 12) COMBINE E-VALUE RESULTS
# ============================================================

evalue_results_all <-
  do.call(
    rbind,
    lapply(
      all_results,
      function(x)
        x$evalue_results
    )
  )


# ============================================================
# 13) COMBINE VARIANCE RESULTS
# ============================================================

variance_results_all <-
  do.call(
    rbind,
    lapply(
      all_results,
      function(x)
        x$variance_results
    )
  )


# ============================================================
# 14) DISPLAY RESULTS
# ============================================================

cat("\n\n====================================================\n")
cat("CIKw PARAMETER-SHIFT RESULTS\n")
cat("====================================================\n")

print(
  signal_results_all
)


cat("\n\n====================================================\n")
cat("E-VALUE RESULTS\n")
cat("====================================================\n")

print(
  evalue_results_all
)


cat("\n\n====================================================\n")
cat("VARIANCE-SHIFT BENCHMARK RESULTS\n")
cat("====================================================\n")

print(
  variance_results_all
)


# ============================================================
# 15) SAVE RESULTS
# ============================================================

write.csv(
  signal_results_all,
  "CIKw_signal_results.csv",
  row.names = FALSE
)


write.csv(
  evalue_results_all,
  "CIKw_evalue_results.csv",
  row.names = FALSE
)


write.csv(
  variance_results_all,
  "CIKw_variance_benchmark_results.csv",
  row.names = FALSE
)


# ============================================================
# 16) FINAL CRITICAL VALUES
# ============================================================

critical_values <- do.call(
  rbind,
  lapply(
    seq_along(all_results),
    function(i){
      
      x <- all_results[[i]]
      
      data.frame(
        n = sample_sizes[i],
        C_GLRT_CIKw = x$C_glrt,
        C_CUSUM_CIKw = x$C_cusum,
        C_GLRT_Variance = x$C_glrt_variance,
        C_CUSUM_Variance = x$C_cusum_variance
      )
    }
  )
)


cat("\n\n====================================================\n")
cat("CRITICAL VALUES\n")
cat("====================================================\n")

print(
  critical_values
)


write.csv(
  critical_values,
  "CIKw_critical_values.csv",
  row.names = FALSE
)