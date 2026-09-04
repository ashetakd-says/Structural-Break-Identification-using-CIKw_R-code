# =========================================
# FULL STABLE R SCRIPT (ERROR-FREE VERSION)
# =========================================

library(stats)

# -----------------------------------------
# DATA
# -----------------------------------------
covid <- c(8.826,6.105,10.383,7.26,13.220,6.015,10.855,6.122,10.685,10.035,5.242,7.630,14.604,7.903,6.327,9.391,14.962,4.730,3.215,16.498,
           11.665,9.284,12.878,6.656,3.440,5.854,8.813,10.043,7.260,5.985,4.424,4.344,5.143,9.935,7.840,9.550,6.968,6.370,3.537,3.286,10.158,
           8.108,6.697,7.151,6.560,2.988,3.336,6.814,8.325,7.854,8.551,3.228,3.499,3.751,7.486,6.625,6.140,4.909,4.661,1.867,2.838,5.392,
           12.042,8.696,6.412,3.395,1.815,3.327,5.406,6.182,4.949,4.089,3.359,2.070,3.298,5.317,5.442,4.557,4.292,2.500,6.535,4.648,4.697,
           5.459,4.120,3.922,3.219,1.402,2.438,3.257,3.632,3.233,3.027,2.352,1.205,2.077,3.778,3.218,2.926,2.601,2.065,1.041,1.800,3.029,
           2.058,2.326,2.506,1.923)

india_revenue <- c(8.1811,9.0245,9.1818,8.5971,9.3023,9.6890,9.0261,9.2137,9.2691,9.2685,9.3065,10.0754,10.4273,10.4083,
                   10.3567,10.4683,9.9935,10.1709,9.8052,8.5352,8.9817,9.2268,9.2314,9.0085,8.1131,8.6383,8.8101,8.0794,
                   8.6762,9.1080,9.5708,10.0809,11.1293,12.1083,10.9771,9.8097,10.3880,10.1773,10.8367,11.0016,9.9846,
                   10.5697,11.1477,11.3874,12.0174)

datasets <- list(COVID_Mexico=covid, India_Revenue=india_revenue)


pdf_CIKw <- function(x,a,b,l1,l2){
  
  x <- pmax(x, 1e-8)
  
  A <- 1 - l1*(l2-1)
  H <- l1*(2*l2-1)
  K <- l1*l2
  
  logx <- log(x)
  
  om <- 1 - exp(-a * logx)
  om <- pmax(om, 1e-12)
  
  val <- a*b * exp(-(a+1)*logx) *
    om^(b-1) *
    (A + 2*H*om^b - 3*K*om^(2*b))
  
  val <- pmax(val, 1e-300)
  val[!is.finite(val)] <- 1e-300
  
  return(val)
}

# ---
# CDF
# ---
cdf_CIKw <- function(x,a,b,l1,l2){
  A <- 1 - l1*(l2-1)
  H <- l1*(2*l2-1)
  K <- l1*l2
  om <- pmax(1 - x^(-a),0)
  pmin(pmax(A*om^b + H*om^(2*b) - K*om^(3*b),0),1)
}

# -----------------------
# SAFE NEG LOG-LIKELIHOOD
# -----------------------
neglog <- function(par,y){
  
  a<-par[1]; b<-par[2]; l1<-par[3]; l2<-par[4]
  
  if(any(!is.finite(par)) || a<=0 || b<=0 || l2<=0)
    return(1e10)
  
  dens <- pdf_CIKw(y,a,b,l1,l2)
  
  if(any(!is.finite(dens)) || any(dens<=0,na.rm=TRUE))
    return(1e10)
  
  -sum(log(dens))
}

# ---
# FIT
# --
fit_CIKw <- function(y){
  
  starts <- list(c(2,2,0.1,1), c(3,1.5,0.2,1), c(6,5,0.8,0.5))
  
  results <- list()
  
  for(s in starts){
    res <- tryCatch(
      optim(s,neglog,y=y,method="L-BFGS-B",
            lower=c(1e-6,1e-6,-1,1e-6),
            upper=c(50,50,1,2)),
      error=function(e) NULL)
    
    if(!is.null(res) && res$convergence==0){
      results[[length(results)+1]] <- list(par=res$par,ll=-res$value)
    }
  }
  
  if(length(results)==0) return(NULL)
  
  best <- results[[which.max(sapply(results,function(x)x$ll))]]
  
  list(par=best$par,ll=best$ll)
}

# ----
# GLRT
# ----
glrt_scan <- function(y,min_size=8){
  
  base <- fit_CIKw(y)
  if(is.null(base)) return(NULL)
  
  bestT <- -Inf; best_tau <- NA
  n <- length(y)
  
  for(tau in min_size:(n-min_size)){
    
    l <- fit_CIKw(y[1:tau])
    r <- fit_CIKw(y[(tau+1):n])
    
    if(is.null(l)||is.null(r)) next
    
    T <- 2*((l$ll+r$ll)-base$ll)
    
    if(is.finite(T) && T>bestT){
      bestT <- T; best_tau <- tau
    }
  }
  
  list(tau=best_tau,T=bestT)
}

# ------------
# SEGMENTATION
# ------------
segment <- function(y,min_size=8){
  
  segs <- list(list(s=1,e=length(y)))
  cps <- c()
  
  iter <- 0
  max_iter <- 50
  
  while(length(segs)>0 && iter<max_iter){
    
    iter <- iter + 1
    
    seg <- segs[[1]]
    segs <- segs[-1]
    
    s<-seg$s; e<-seg$e
    segy <- y[s:e]
    
    cat("\nSegment",s,":",e,"\n")
    
    if(length(segy)<2*min_size) next
    
    res <- glrt_scan(segy,min_size)
    if(is.null(res)) next
    
    if(!is.finite(res$T)) next
    
    crit <- qchisq(0.95,4)
    
    if(res$T > crit){
      cp <- s + res$tau - 1
      cps <- c(cps,cp)
      
      segs <- append(list(list(s=s,e=cp)),segs)
      segs <- append(list(list(s=cp+1,e=e)),segs)
    }
  }
  
  sort(cps)
}

# ---
# RUN
# ----
for(name in names(datasets)){
  
  cat("\n============================\n",name,"\n")
  
  y <- datasets[[name]]
  
  cps <- segment(y)
  
  cat("Detected CPs:",cps,"\n")
  
  all_cuts <- c(0,cps,length(y))
  
  for(i in 1:(length(all_cuts)-1)){
    
    segy <- y[(all_cuts[i]+1):all_cuts[i+1]]
    
    fit <- fit_CIKw(segy)
    if(is.null(fit)) next
    
    ks <- tryCatch(
      ks.test(segy,function(x)cdf_CIKw(x,fit$par[1],fit$par[2],fit$par[3],fit$par[4])),
      error=function(e) list(p.value=NA)
    )
    
    cat("\nSegment",i,"KS p:",ks$p.value,"\n")
  }
  
  plot(y,type="b",main=name)
  abline(v=cps,col="red",lty=2)
}

# ==========
# SIMULATION
# ==========

simulate_india_like <- function(n=45,case="H0"){
  
  if(case=="H0") return(runif(n,2,6))
  
  if(case=="CP_balanced"){
    k<-31
    return(c(runif(k,3,6),runif(n-k,1,3)))
  }
  
  if(case=="CP_unbalanced"){
    k<-35
    return(c(runif(k,3,6),runif(n-k,1,3)))
  }
}

cat("\n=== SIMULATION ===\n")

for(case in c("H0","CP_balanced","CP_unbalanced")){
  
  detect<-0
  
  for(b in 1:50){
    y<-simulate_india_like(45,case)
    cps<-segment(y)
    if(length(cps)>0) detect<-detect+1
  }
  
  cat(case,"| Detection =",detect/50,"\n")
}