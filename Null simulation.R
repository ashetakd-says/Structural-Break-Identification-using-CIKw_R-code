rm(list = ls())
library(stats)
library(randtests)

# ---------------------------
# 0) DATA
# ---------------------------

covid <- c(
  8.826,6.105,10.383,7.26,13.220,6.015,10.855,6.122,10.685,10.035,
  5.242,7.630,14.604,7.903,6.327,9.391,14.962,4.730,3.215,16.498,
  11.665,9.284,12.878,6.656,3.440,5.854,8.813,10.043,7.260,5.985,
  4.424,4.344,5.143,9.935,7.840,9.550,6.968,6.370,3.537,3.286,
  10.158,8.108,6.697,7.151,6.560,2.988,3.336,6.814,8.325,7.854,
  8.551,3.228,3.499,3.751,7.486,6.625,6.140,4.909,4.661,1.867,
  2.838,5.392,12.042,8.696,6.412,3.395,1.815,3.327,5.406,6.182,
  4.949,4.089,3.359,2.070,3.298,5.317,5.442,4.557,4.292,2.500,
  6.535,4.648,4.697,5.459,4.120,3.922,3.219,1.402,2.438,3.257,
  3.632,3.233,3.027,2.352,1.205,2.077,3.778,3.218,2.926,2.601,
  2.065,1.041,1.800,3.029,2.058,2.326,2.506,1.923
)

india_revenue <- c(
  8.1811,9.0245,9.1818,8.5971,9.3023,9.6890,9.0261,9.2137,9.2691,9.2685,
  9.3065,10.0754,10.4273,10.4083,10.3567,10.4683,9.9935,10.1709,9.8052,
  8.5352,8.9817,9.2268,9.2314,9.0085,8.1131,8.6383,8.8101,8.0794,
  8.6762,9.1080,9.5708,10.0809,11.1293,12.1083,10.9771,9.8097,
  10.3880,10.1773,10.8367,11.0016,9.9846,10.5697,11.1477,11.3874,12.0174
)

eps <- 1e-12

# ---------------------------
# 1) RUNS TEST
# ---------------------------
runs_test <- function(y) {
  med <- median(y)
  signs <- as.numeric(y >= med)
  
  runs <- 1 + sum(signs[-1] != signs[-length(signs)])
  n1 <- sum(signs == 1)
  n2 <- sum(signs == 0)
  
  expected <- 1 + 2*n1*n2/(n1+n2)
  var <- (2*n1*n2*(2*n1*n2 - n1 - n2)) / ((n1+n2)^2*(n1+n2-1))
  
  z <- (runs - expected)/sqrt(var)
  p <- 2*(1 - pnorm(abs(z)))
  
  return(p)
}

# ---------------------------
# 2) PREPROCESSING
# ---------------------------
preprocess_series <- function(y) {
  
  p_runs <- runs_test(log(y + eps))
  
  if (p_runs < 0.05) {
    yd <- diff(log(y + eps))
    y_new <- exp(yd)
    
    shift <- max(0, 1.001 - min(y_new))
    y_new <- y_new + shift
    
    cat(sprintf("Transformation applied (p=%.4f), shift=%.4f\n", p_runs, shift))
    return(y_new)
  } else {
    cat(sprintf("No transformation (p=%.4f)\n", p_runs))
    return(y)
  }
}

# ---------------------------
# 3) CIKw PDF & CDF
# ---------------------------
pdf_CIKw <- function(x, a, b, l1, l2) {
  A <- 1 - l1*(l2 - 1)
  H <- l1*(2*l2 - 1)
  K <- l1*l2
  om <- pmax(pmin(1 - x^(-a), 1 - 1e-12), 1e-12)
  
  val <- a*b*x^(-a-1)*om^(b-1)*(A + 2*H*om^b - 3*K*om^(2*b))
  return(pmax(val, 1e-300))
}

cdf_CIKw <- function(x, a, b, l1, l2) {
  A <- 1 - l1*(l2 - 1)
  H <- l1*(2*l2 - 1)
  K <- l1*l2
  om <- pmax(pmin(1 - x^(-a), 1 - 1e-12), 1e-12)
  
  val <- A*om^b + H*om^(2*b) - K*om^(3*b)
  return(pmax(pmin(val,1),0))
}

# ---------------------------
# 4) SIMPLE FIT (placeholder optimization)
# ---------------------------
neglog <- function(p, y) {
  a=p[1]; b=p[2]; l1=p[3]; l2=p[4]
  if (a<=0 || b<=0 || l2<=0) return(Inf)
  f <- pdf_CIKw(y,a,b,l1,l2)
  return(-sum(log(f)))
}

fit_CIKw <- function(y) {
  res <- optim(c(2,2,0.1,1), neglog, y=y,
               lower=c(1e-6,1e-6,-1,1e-6),
               upper=c(Inf,Inf,1,2),
               method="L-BFGS-B")
  return(list(par=res$par, ll=-res$value))
}

# ---------------------------
# 5) GLR SCAN
# ---------------------------
glr_scan <- function(y, min_size=8) {
  fit_all <- fit_CIKw(y)
  ll_all <- fit_all$ll
  n <- length(y)
  
  Tvals <- c()
  taus <- c()
  
  for (t in min_size:(n-min_size)) {
    ll1 <- fit_CIKw(y[1:t])$ll
    ll2 <- fit_CIKw(y[(t+1):n])$ll
    
    Tvals <- c(Tvals, 2*((ll1+ll2)-ll_all))
    taus <- c(taus, t)
  }
  
  return(list(taus=taus, T=Tvals))
}

# ---------------------------
# 6) SAMPLE FROM CIKw
# ---------------------------
sample_CIKw <- function(n, params) {
  a=params[1]; b=params[2]; l1=params[3]; l2=params[4]
  
  x_grid <- seq(1.001, 50, length.out=5000)
  cdf_vals <- cdf_CIKw(x_grid,a,b,l1,l2)
  cdf_vals <- cummax(cdf_vals)
  
  u <- runif(n)
  return(approx(cdf_vals, x_grid, u)$y)
}

# ---------------------------
# 7) BOOTSTRAP C_alpha
# ---------------------------
bootstrap_Calpha <- function(y, B=200) {
  n <- length(y)
  Tmax <- c()
  
  for (b in 1:B) {
    yb <- sample(y, n, replace=TRUE)
    res <- glr_scan(yb)
    
    if (length(res$T) > 0)
      Tmax <- c(Tmax, max(res$T))
  }
  
  return(quantile(Tmax, 0.95))
}

# ---------------------------
# 8) ONE NULL RUN
# ---------------------------
one_null_run <- function(y_ref) {
  n <- length(y_ref)
  
  fit <- fit_CIKw(y_ref)
  y_sim <- sample_CIKw(n, fit$par)
  
  res <- glr_scan(y_sim)
  if (length(res$T) == 0) return(FALSE)
  
  T_best <- max(res$T)
  
  chi_cut <- qchisq(0.95,4)
  C_alpha <- bootstrap_Calpha(y_sim)
  
  return( (T_best > chi_cut) || (T_best > C_alpha) )
}

# ---------------------------
# 9) FULL NULL SIMULATION
# ---------------------------
null_simulation <- function(y_raw, R=100) {
  
  cat("\nPreprocessing...\n")
  y_ref <- preprocess_series(y_raw)
  
  rejections <- 0
  
  for (r in 1:R) {
    if (one_null_run(y_ref)) {
      rejections <- rejections + 1
    }
    
    if (r %% 20 == 0) {
      cat(sprintf("Progress: %d/%d\n", r, R))
    }
  }
  
  type1 <- rejections / R
  
  cat("\n====================================\n")
  cat(sprintf("Estimated Type I error: %.4f\n", type1))
  cat("====================================\n")
  
  return(type1)
}

# ---------------------------
# 10) RUN
# ---------------------------
cat("\n===== COVID DATA =====\n")
type1_covid <- null_simulation(covid, R=100)

cat("\n===== INDIA REVENUE =====\n")
type1_india <- null_simulation(india_revenue, R=100)



# ==========================================================
# DATA
# ==========================================================
datasets <- list(Floyd=floyd, COVID=covid, India=india)

# ==========================================================
# FIXED SEGMENTS (USER-DEFINED)
# ==========================================================
fixed_segments <- list(
  Floyd = list(c(1,20), c(21,39)),
  COVID = list(c(1,65), c(66,87), c(88,108)),
  India = list(c(1,30), c(31,45))
)

# ==========================================================
# CIKw PDF + LOG-LIK
# ==========================================================
pdf_CIKw <- function(x,a,b,l1,l2){
  om <- pmax(pmin(1 - x^(-a),1),1e-10)
  A <- 1 - l1*(l2-1)
  H <- l1*(2*l2-1)
  K <- l1*l2
  d <- a*b*x^(-a-1)*(om^(b-1))*(A + 2*H*om^b - 3*K*om^(2*b))
  pmax(d,1e-300)
}

loglik <- function(par,y){
  if(par[1]<=0 || par[2]<=0 || par[4]<=0) return(-1e12)
  sum(log(pdf_CIKw(y,par[1],par[2],par[3],par[4])))
}

# ==========================================================
# MULTI-START MLE (PSO-LIKE DIAGNOSTIC)
# ==========================================================
fit_MLE <- function(y){
  
  starts <- list(
    c(2,2,0.1,1),
    c(3,2,-0.3,1.2),
    c(5,3,-0.5,1.5),
    c(7,4,-0.7,1.8)
  )
  
  best_ll <- -Inf
  best_par <- NULL
  converged <- 0
  
  for(s in starts){
    res <- try(optim(s, fn=function(p) -loglik(p,y),
                     method="L-BFGS-B",
                     lower=c(1e-3,1e-3,-1,1e-3),
                     upper=c(15,15,1,2)),
               silent=TRUE)
    
    if(class(res)!="try-error"){
      ll <- loglik(res$par,y)
      if(is.finite(ll)){
        converged <- converged + 1
        if(ll > best_ll){
          best_ll <- ll
          best_par <- res$par
        }
      }
    }
  }
  
  list(par=best_par,ll=best_ll,conv_rate=converged/length(starts))
}

# ==========================================================
# DIAGNOSTICS (GLRT, BIC, MIC, CUSUM)
# ==========================================================
diagnostics_cp <- function(y, c_mic=0.26, min_size=8){
  
  n <- length(y)
  if(n < 2*min_size) return(NULL)
  
  ll0 <- fit_MLE(y)$ll
  
  res <- data.frame()
  
  for(tau in min_size:(n-min_size)){
    
    ll1 <- fit_MLE(y[1:tau])$ll + fit_MLE(y[(tau+1):n])$ll
    
    GLRT <- 2*(ll1 - ll0)
    BIC  <- (-2*ll1 + 8*log(n)) - (-2*ll0 + 4*log(n))
    MIC  <- (-2*ll1 + (8 + c_mic*(1/tau + 1/(n-tau)))*log(n)) - (-2*ll0 + 4*log(n))
    
    res <- rbind(res, data.frame(tau=tau, GLRT=GLRT, BIC=BIC, MIC=MIC))
  }
  
  # best candidates
  best_glrt <- res$tau[which.max(res$GLRT)]
  best_bic  <- res$tau[which.min(res$BIC)]
  best_mic  <- res$tau[which.min(res$MIC)]
  
  list(table=res,
       best=c(GLRT=best_glrt, BIC=best_bic, MIC=best_mic))
}

cusum_cp <- function(y){
  s <- cumsum(y - mean(y))
  which.max(abs(s))
}

# ==========================================================
# MAIN ANALYSIS
# ==========================================================
results <- data.frame()
agreement <- data.frame()

for(nm in names(datasets)){
  
  cat("\n====================\n",nm,"\n")
  
  y <- datasets[[nm]]
  segs <- fixed_segments[[nm]]
  
  for(i in seq_along(segs)){
    
    s <- segs[[i]][1]
    e <- segs[[i]][2]
    seg <- y[s:e]
    
    cat("\nSegment",i,"[",s,":",e,"]\n")
    
    # ---- MLE ----
    fit <- fit_MLE(seg)
    
    # ---- diagnostics ----
    diag <- diagnostics_cp(seg)
    cus <- cusum_cp(seg)
    
    if(!is.null(diag)){
      best <- diag$best
      
      agreement <- rbind(agreement, data.frame(
        Dataset=nm,
        Segment=i,
        GLRT_CP=best["GLRT"] + s - 1,
        BIC_CP=best["BIC"] + s - 1,
        MIC_CP=best["MIC"] + s - 1,
        CUSUM_CP=cus + s - 1
      ))
    }
    
    # ---- store results ----
    results <- rbind(results, data.frame(
      Dataset=nm,
      Segment=i,
      Start=s,
      End=e,
      alpha=fit$par[1],
      beta=fit$par[2],
      lambda1=fit$par[3],
      lambda2=fit$par[4],
      Mean=mean(seg),
      SD=sd(seg),
      Length=length(seg),
      ConvRate=fit$conv_rate
    ))
  }
}

# ==========================================================
# OUTPUT
# ==========================================================
cat("\n====================\nSEGMENT VALIDITY\n")
print(results)

cat("\n====================\nCP ROBUSTNESS (AGREEMENT)\n")
print(agreement)