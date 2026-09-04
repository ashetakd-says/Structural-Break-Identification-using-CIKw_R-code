# ============================================
# FINAL SCRIPT: MCP + TYPE I ERROR + E-VALUES
# GLRT + CUSUM + BIC + MIC
# =============================================

rm(list = ls())
set.seed(2025)

# ================
# 1) CDF + SAMPLER
# ================
CDF_y <- function(y,a,b,l1,l2){
  A <- 1 - l1*(l2 - 1)
  H <- l1*(2*l2 - 1)
  K <- l1*l2
  t <- pmax(1 - y^(-a),0)
  A*t^b + H*t^(2*b) - K*t^(3*b)
}

IIT_fast <- function(n,a,b,l1,l2){
  grid <- c(seq(1.0001,10,length=300),
            seq(10,100,length=300),
            exp(seq(log(100),log(1e4),length=300)))
  
  F <- CDF_y(grid,a,b,l1,l2)
  F <- pmin(cummax(F),1)
  
  approx(F,grid,runif(n))$y
}

# =================
# 2) LOG-LIKELIHOOD
# =================
loglik <- function(y,p){
  a=p[1]; b=p[2]; l1=p[3]; l2=p[4]
  
  if(a<=0 || b<=0) return(-Inf)
  
  t <- 1 - y^(-a)
  if(any(t<=0)) return(-Inf)
  
  A <- 1 - l1*(l2 - 1)
  H <- l1*(2*l2 - 1)
  K <- l1*l2
  
  term <- A*t^(b-1) + 2*H*t^(2*b-1) - 3*K*t^(3*b-1)
  if(any(term<=0)) return(-Inf)
  
  sum(log(a*b) - (a+1)*log(y) + log(term))
}

# ======
# 3) MLE
# ======
mle <- function(y){
  starts <- list(c(3,2,-0.3,1.1), c(5,2,-0.5,1.5))
  best <- -Inf
  
  for(s in starts){
    res <- try(optim(s, function(p) -loglik(y,p),
                     method="L-BFGS-B",
                     lower=c(1e-6,1e-6,-0.99,0.001),
                     upper=c(Inf,Inf,0.99,2)), silent=TRUE)
    
    if(class(res)!="try-error"){
      ll <- loglik(y,res$par)
      if(is.finite(ll) && ll > best){
        best <- ll
      }
    }
  }
  best
}

# =====================
# 4) GLRT STAT + LAMBDA
# =====================
GLRT_stat <- function(y, lmin=10){
  n <- length(y)
  ll0 <- mle(y)
  Tmax <- -Inf
  
  for(k in lmin:(n-lmin)){
    ll1 <- mle(y[1:k]) + mle(y[(k+1):n])
    if(is.finite(ll1)){
      Tmax <- max(Tmax, 2*(ll1 - ll0))
    }
  }
  if(!is.finite(Tmax)) Tmax <- 0
  Tmax
}

GLRT_lambda <- function(y, lmin=10){
  n <- length(y)
  ll0 <- mle(y)
  best <- -Inf
  
  for(k in lmin:(n-lmin)){
    ll1 <- mle(y[1:k]) + mle(y[(k+1):n])
    if(is.finite(ll1)){
      best <- max(best, ll1 - ll0)
    }
  }
  if(!is.finite(best)) best <- 0
  best
}

# ========
# 5) CUSUM
# ========
CUSUM_stat <- function(y){
  n <- length(y)
  mu <- mean(y)
  s <- cumsum(y - mu)
  max(abs(s)) / (sd(y) * sqrt(n))
}

# ============
# 6) BIC / MIC
# ============
BIC_MIC <- function(y,c=0.1){
  n <- length(y)
  p0 <- 4; p1 <- 8
  
  ll0 <- mle(y)
  B0 <- -2*ll0 + p0*log(n)
  M0 <- B0
  
  bestB <- Inf
  bestM <- Inf
  
  for(k in 10:(n-10)){
    ll1 <- mle(y[1:k]) + mle(y[(k+1):n])
    if(!is.finite(ll1)) next
    
    B1 <- -2*ll1 + p1*log(n)
    w <- c*p1*(1/k + 1/(n-k))
    M1 <- -2*ll1 + (p1+w)*log(n)
    
    bestB <- min(bestB,B1)
    bestM <- min(bestM,M1)
  }
  
  c(bestB < B0, bestM < M0)
}

# ========================
# 7) BOOTSTRAP CALIBRATION
# ========================
calibrate_glrt <- function(n,B=200){
  vals <- replicate(B,{
    y <- IIT_fast(n,3,2,0,1)
    GLRT_stat(y)
  })
  quantile(vals,0.95)
}

calibrate_cusum <- function(n,B=200){
  vals <- replicate(B,{
    y <- IIT_fast(n,3,2,0,1)
    CUSUM_stat(y)
  })
  quantile(vals,0.95)
}

# ======================================
# 8) MCP DETECTION (Binary Segmentation)
# ======================================
detect_cp <- function(y, C_glrt){
  T <- GLRT_stat(y)
  return(T > C_glrt)
}

# ===============================
# 9) TYPE I ERROR SIMULATION (MCP)
# ================================
simulate_type1 <- function(n,R=200,B=200){
  
  cat("\nCalibrating thresholds for n =", n, "...\n")
  C_glrt <- calibrate_glrt(n,B)
  C_cusum <- calibrate_cusum(n,B)
  
  cat("GLRT Threshold =", C_glrt, "\n")
  cat("CUSUM Threshold =", C_cusum, "\n")
  
  GLRT_rej <- CUSUM_rej <- BIC_rej <- MIC_rej <- 0
  
  for(i in 1:R){
    
    y <- IIT_fast(n,3,2,0,1)  # H0
    
    # GLRT
    if(GLRT_stat(y) > C_glrt) GLRT_rej <- GLRT_rej + 1
    
    # CUSUM
    if(CUSUM_stat(y) > C_cusum) CUSUM_rej <- CUSUM_rej + 1
    
    # BIC / MIC
    bm <- BIC_MIC(y)
    BIC_rej <- BIC_rej + bm[1]
    MIC_rej <- MIC_rej + bm[2]
  }
  
  cat("\n--- TYPE I ERROR (MCP) ---\n")
  cat("GLRT :", GLRT_rej/R, "\n")
  cat("CUSUM:", CUSUM_rej/R, "\n")
  cat("BIC  :", BIC_rej/R, "\n")
  cat("MIC  :", MIC_rej/R, "\n")
}

# =======
# 10) RUN
# =======
for(n in c(50,100,200)){
  cat("\n====================================\n")
  cat("RUNNING FOR n =", n, "\n")
  simulate_type1(n)
}



# ============================================================
# 11) VARIANCE-SHIFT BENCHMARK
#     Mean remains unchanged; only variance changes
# ============================================================

generate_variance_shift <- function(
    n,
    tau = floor(n/2),
    mu = 5,
    sigma1 = 0.5,
    variance_ratio = 2
){
  
  sigma2 <- sigma1 * sqrt(variance_ratio)
  
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
  
  c(y1, y2)
}


# ============================================================
# 12) Calibrate benchmark thresholds under NO CHANGE
# ============================================================

calibrate_variance_benchmark <- function(
    n,
    R = 200
){
  
  tau <- floor(n / 2)
  
  glrt_vals <- numeric(R)
  cusum_vals <- numeric(R)
  bic_vals <- numeric(R)
  mic_vals <- numeric(R)
  
  for(i in 1:R){
    
    # No variance change
    y <- generate_variance_shift(
      n = n,
      tau = tau,
      variance_ratio = 1
    )
    
    glrt_vals[i] <- GLRT_stat(y)
    cusum_vals[i] <- CUSUM_stat(y)
    
    bm <- BIC_MIC(y)
    
    bic_vals[i] <- bm[1]
    mic_vals[i] <- bm[2]
  }
  
  list(
    C_glrt = quantile(glrt_vals, 0.95),
    C_cusum = quantile(cusum_vals, 0.95)
  )
}

# ============================================================
# 13) UNIVARIATE VARIANCE-SHIFT POWER SIMULATION
#     Mean remains unchanged; only variance changes
# ============================================================

simulate_variance_power <- function(
    n,
    variance_ratios = c(1, 1.25, 1.5, 2, 3, 4),
    R = 200,
    B = 200
){
  
  tau <- floor(n / 2)
  
  # ----------------------------------------------------------
  # Calibrate critical values under H0:
  # no variance change (variance ratio = 1)
  # ----------------------------------------------------------
  
  cat(
    "\nCalibrating univariate variance-shift benchmark for n =",
    n, "...\n"
  )
  
  thresholds <- calibrate_variance_benchmark(
    n = n,
    B = B
  )
  
  C_glrt <- thresholds$C_glrt
  C_cusum <- thresholds$C_cusum
  
  cat("GLRT critical value  =", C_glrt, "\n")
  cat("CUSUM critical value =", C_cusum, "\n")
  
  results <- data.frame()
  
  # ----------------------------------------------------------
  # Evaluate empirical power for each variance ratio
  # ----------------------------------------------------------
  
  for(vr in variance_ratios){
    
    GLRT_rej <- 0
    CUSUM_rej <- 0
    BIC_rej <- 0
    MIC_rej <- 0
    
    for(i in seq_len(R)){
      
      # ------------------------------------------------------
      # Generate UNIVARIATE variance-shift series
      # ------------------------------------------------------
      
      y <- generate_variance_shift(
        n = n,
        tau = tau,
        variance_ratio = vr
      )
      
      # ------------------------------------------------------
      # GLRT
      # ------------------------------------------------------
      
      T_glrt <- GLRT_stat(y)
      
      if(T_glrt > C_glrt){
        GLRT_rej <- GLRT_rej + 1
      }
      
      # ------------------------------------------------------
      # CUSUM
      # ------------------------------------------------------
      
      T_cusum <- CUSUM_stat(y)
      
      if(T_cusum > C_cusum){
        CUSUM_rej <- CUSUM_rej + 1
      }
      
      # ------------------------------------------------------
      # BIC / MIC
      # ------------------------------------------------------
      
      bm <- BIC_MIC(y)
      
      BIC_rej <- BIC_rej + bm[1]
      MIC_rej <- MIC_rej + bm[2]
    }
    
    # --------------------------------------------------------
    # Store empirical power
    # --------------------------------------------------------
    
    results <- rbind(
      results,
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
      "n =", n,
      "| Variance ratio =", vr,
      "| GLRT =", round(GLRT_rej / R, 3),
      "| CUSUM =", round(CUSUM_rej / R, 3),
      "| BIC =", round(BIC_rej / R, 3),
      "| MIC =", round(MIC_rej / R, 3),
      "\n"
    )
  }
  
  return(results)
}


# ============================================================
# 14) RUN UNIVARIATE VARIANCE-SHIFT BENCHMARK
# ============================================================

variance_results <- do.call(
  rbind,
  lapply(
    c(50, 100, 200),
    function(n){
      
      simulate_variance_power(
        n = n,
        variance_ratios =
          c(1, 1.25, 1.5, 2, 3, 4),
        R = 200,
        B = 200
      )
    }
  )
)


# ============================================================
# 15) PRINT AND SAVE RESULTS
# ============================================================

print(variance_results)

write.csv(
  variance_results,
  "MVCIKw_Univariate_Variance_Shift_Benchmark.csv",
  row.names = FALSE
)


# ============================================================
# 16) POWER CURVES FOR UNIVARIATE VARIANCE SHIFTS
# ============================================================

library(ggplot2)
library(tidyr)

plot_data <- variance_results |>
  pivot_longer(
    cols = c(
      GLRT,
      CUSUM,
      BIC,
      MIC
    ),
    names_to = "Method",
    values_to = "Power"
  )


p_variance <- ggplot(
  plot_data,
  aes(
    x = Variance_Ratio,
    y = Power,
    colour = Method,
    group = Method
  )
) +
  
  geom_line(
    linewidth = 0.9
  ) +
  
  geom_point(
    size = 2
  ) +
  
  facet_wrap(
    ~ n,
    nrow = 1
  ) +
  
  labs(
    title = "Univariate Benchmark Under Variance Shifts",
    
    x = expression(
      "Variance ratio " ~
        sigma[2]^2 / sigma[1]^2
    ),
    
    y = "Empirical power",
    
    colour = "Method"
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "bottom"
  )


# Display plot
print(p_variance)


# ------------------------------------------------------------
# Save manuscript-quality figures
# ------------------------------------------------------------

ggsave(
  "MVCIKw_Univariate_Variance_Shift_Benchmark.pdf",
  plot = p_variance,
  width = 9,
  height = 3.8,
  units = "in",
  device = "pdf"
)

ggsave(
  "MVCIKw_Univariate_Variance_Shift_Benchmark.png",
  plot = p_variance,
  width = 9,
  height = 3.8,
  units = "in",
  dpi = 600
)