rm(list = ls())
set.seed(123)

library(stats)
library(goftest)

# ======================================================
# Utilities
# ======================================================
scale_to_unit <- function(y){
  eps <- 1e-6
  (y - min(y))/(max(y)-min(y))*(1 - 2*eps) + eps
}

jitter_y <- function(y){
  jitter(y, factor = 1e-8)
}

safe_log <- function(x){
  log(pmax(x, 1e-12))
}

# ======================================================
# Quadratic Rank Transmutation
# ======================================================
QRT_CDF <- function(G, lambda){
  (1 + lambda)*G - lambda*G^2
}

QRT_PDF <- function(g, G, lambda){
  g*((1 + lambda) - 2*lambda*G)
}

# ======================================================
# 1. Transmuted Pareto (TP)
# ======================================================
pTP <- function(y, xm, a, l){
  G <- 1 - (xm/y)^a
  QRT_CDF(G, l)
}

dTP <- function(y, xm, a, l){
  g <- a*xm^a/y^(a+1)
  G <- 1 - (xm/y)^a
  QRT_PDF(g, G, l)
}

ll_TP <- function(p, y){
  xm <- exp(p[1])
  a  <- exp(p[2])
  l  <- 2*plogis(p[3]) - 1
  if(any(y <= xm)) return(-Inf)
  sum(safe_log(dTP(y, xm, a, l)))
}

# ======================================================
# 2. Transmuted Exponentiated Lomax (TEL)
# ======================================================
pTEL <- function(y, a, b, t, l){
  G0 <- 1 - (1 + y/b)^(-a)
  G  <- G0^t
  QRT_CDF(G, l)
}

dTEL <- function(y, a, b, t, l){
  g0 <- (a/b)*(1 + y/b)^(-(a+1))
  G0 <- 1 - (1 + y/b)^(-a)
  g  <- t*g0*G0^(t-1)
  G  <- G0^t
  QRT_PDF(g, G, l)
}

ll_TEL <- function(p, y){
  a <- exp(p[1]); b <- exp(p[2]); t <- exp(p[3])
  l <- 2*plogis(p[4]) - 1
  sum(safe_log(dTEL(y, a, b, t, l)))
}

# ======================================================
# 3. Transmuted Lomax (TL)
# ======================================================
pTL <- function(y, a, s, l){
  G <- 1 - (1 + y/s)^(-a)
  QRT_CDF(G, l)
}

dTL <- function(y, a, s, l){
  g <- (a/s)*(1 + y/s)^(-(a+1))
  G <- 1 - (1 + y/s)^(-a)
  QRT_PDF(g, G, l)
}

ll_TL <- function(p, y){
  a <- exp(p[1]); s <- exp(p[2])
  l <- 2*plogis(p[3]) - 1
  sum(safe_log(dTL(y, a, s, l)))
}

# ======================================================
# 4. Transmuted Weibull (TW)
# ======================================================
pTW <- function(y, k, s, l){
  G <- 1 - exp(-(y/s)^k)
  QRT_CDF(G, l)
}

dTW <- function(y, k, s, l){
  g <- (k/s)*(y/s)^(k-1)*exp(-(y/s)^k)
  G <- 1 - exp(-(y/s)^k)
  QRT_PDF(g, G, l)
}

ll_TW <- function(p, y){
  k <- exp(p[1]); s <- exp(p[2])
  l <- 2*plogis(p[3]) - 1
  sum(safe_log(dTW(y, k, s, l)))
}

# ======================================================
# 5. Transmuted Inverted Kumaraswamy (TIKw)
# ======================================================
pTIKw <- function(y, a, b, l){
  G <- (1 - y^(-a))^b
  QRT_CDF(G, l)
}

dTIKw <- function(y, a, b, l){
  g <- a*b*y^(-(a+1))*(1 - y^(-a))^(b-1)
  G <- (1 - y^(-a))^b
  QRT_PDF(g, G, l)
}

ll_TIKw <- function(p, y){
  a <- exp(p[1]); b <- exp(p[2])
  l <- 2*plogis(p[3]) - 1
  if(any(y <= 1)) return(-Inf)
  sum(safe_log(dTIKw(y, a, b, l)))
}

# ======================================================
# 6. Cubic Inverted Kumaraswamy (CIKw) – Stabilized
# ======================================================
pCIKw <- function(y, th1, th2, th3, th4){
  a  <- exp(th1)
  b  <- exp(th2)
  l1 <- 2*plogis(th3) - 1
  l2 <- 2*plogis(th4)
  
  G <- (1 - y^(-a))^b
  A <- 1 - l1*(l2 - 1)
  B <- l1*(2*l2 - 1)
  C <- l1*l2
  
  pmin(pmax(A*G + B*G^2 - C*G^3, 0), 1)
}

dCIKw <- function(y, th1, th2, th3, th4){
  a  <- exp(th1)
  b  <- exp(th2)
  l1 <- 2*plogis(th3) - 1
  l2 <- 2*plogis(th4)
  
  G <- (1 - y^(-a))^b
  g <- a*b*y^(-(a+1))*(1 - y^(-a))^(b-1)
  
  A <- 1 - l1*(l2 - 1)
  B <- l1*(2*l2 - 1)
  C <- l1*l2
  
  pmax(g*(A + 2*B*G - 3*C*G^2), 1e-12)
}

ll_CIKw <- function(p, y){
  if(any(y <= 1)) return(-Inf)
  sum(log(dCIKw(y, p[1], p[2], p[3], p[4])))
}

# ======================================================
# 7–8. IKw1 and IKw2
# ======================================================
# IKw1
pIKw1 <- function(y, a, b){
  if(any(y <= 1)) return(rep(NA, length(y)))
  1 - (1 - y^(-a))^b
}

dIKw1 <- function(y, a, b){
  if(any(y <= 1)) return(rep(0, length(y)))
  a*b*y^(-(a+1))*(1 - y^(-a))^(b-1)
}

ll_IKw1 <- function(p, y){
  a <- exp(p[1]); b <- exp(p[2])
  if(any(y <= 1)) return(-Inf)
  sum(log(dIKw1(y, a, b)))
}


# IKw2
pIKw2 <- function(y, a, b){
  if(any(y >= 0)) return(rep(NA, length(y)))
  1 - (1 - (1 - y)^(-a))^b
}

dIKw2 <- function(y, a, b){
  if(any(y >= 0)) return(rep(0, length(y)))
  a*b*(1-y)^(-(a+1))*(1 - (1 - y)^(-a))^(b-1)
}

ll_IKw2 <- function(p, y){
  a <- exp(p[1]); b <- exp(p[2])
  if(any(y >= 0)) return(-Inf)
  sum(log(dIKw2(y, a, b)))
}

# ======================================================
# Lognormal Distribution (CORRECTED)
# ======================================================
pLN <- function(y, mu, log_sigma){
  sigma <- exp(log_sigma)
  plnorm(y, meanlog = mu, sdlog = sigma)
}

dLN <- function(y, mu, log_sigma){
  sigma <- exp(log_sigma)
  dlnorm(y, meanlog = mu, sdlog = sigma)
}

ll_LN <- function(p, y){
  mu        <- p[1]
  log_sigma <- p[2]
  sigma     <- exp(log_sigma)
  
  if(any(y <= 0) || sigma <= 0) return(-Inf)
  
  sum(dlnorm(y, meanlog = mu, sdlog = sigma, log = TRUE))
}
# ======================================================
# Unified fitting function (ALL MODELS – CORRECTED)
# ======================================================
fit_all_models <- function(y){
  
  yj <- jitter_y(y)
  n  <- length(y)
  out <- data.frame()
  
  models <- list(
    TP   = list(ll=ll_TP,   p=pTP,   k=3, start=c(0,0,0)),
    TL   = list(ll=ll_TL,   p=pTL,   k=3, start=c(0,0,0)),
    TEL  = list(ll=ll_TEL,  p=pTEL,  k=4, start=c(0,0,0,0)),
    TW   = list(ll=ll_TW,   p=pTW,   k=3, start=c(0,0,0)),
    TIKw = list(ll=ll_TIKw, p=pTIKw, k=3, start=c(0,0,0)),
    IKw1 = list(ll=ll_IKw1, p=pIKw1, k=2, start=c(0,0)),
    IKw2 = list(ll=ll_IKw2, p=pIKw2, k=2, start=c(0,0)),
    CIKw = list(ll=ll_CIKw, p=pCIKw, k=4, start=c(0,0,0,0)),
    LN   = list(
      ll = ll_LN,
      p  = pLN,
      k  = 2,
      start = c(mean(log(y)), log(sd(log(y))))
    )
  )
  
  fitted_models <- list()
  
  for(m in names(models)){
    
    fit <- try(
      optim(
        par     = models[[m]]$start,
        fn      = models[[m]]$ll,
        y       = y,
        method  = "BFGS",
        control = list(fnscale = -1)
      ),
      silent = TRUE
    )
    
    if(inherits(fit, "try-error")) next
    
    fitted_models[[m]] <- fit$par
    
    # Safe KS test
    ks <- try(
      ks.test(
        yj,
        function(x){
          p <- do.call(models[[m]]$p, c(list(x), as.list(fit$par)))
          pmin(pmax(p, 0), 1)
        }
      ),
      silent = TRUE
    )
    
    if(inherits(ks, "try-error")) next
    
    out <- rbind(
      out,
      data.frame(
        Model  = m,
        LogLik = fit$value,
        AIC    = -2*fit$value + 2*models[[m]]$k,
        BIC    = -2*fit$value + models[[m]]$k * log(n),
        KS     = as.numeric(ks$statistic),
        pvalue = ks$p.value
      )
    )
  }
  
  attr(out, "fits") <- fitted_models
  out[order(out$AIC), ]
}





# ======================================================
# Data
# ======================================================
IndGini <- read.csv("D:/Training Material/India Gini.csv")
y_gini  <- IndGini$Gini

CF <- c(2.247,2.64,2.908,3.099,3.126,3.245,3.328,3.355,3.383,
        3.572,3.581,3.681,3.726,3.727,3.728,3.783,3.785,3.786,
        3.896,3.912,3.964,4.05,4.063,4.082,4.111,4.118,4.141,
        4.246,4.251,4.262,4.326,4.402,4.457,4.466,4.519,4.542,
        4.555,4.614,4.632,4.634,4.636,4.678,4.698,4.738,4.832,
        4.924,5.043,5.099,5.134,5.359,5.473,5.571,5.684,5.721,
        5.998,6.06)

RCC <- c(3.96,4.41,4.14,4.11,4.45,4.10,4.31,4.42,4.30,4.51,
         4.71,4.62,4.35,4.26,4.63,4.36,3.91,4.51,4.37,4.90,
         4.46,3.95,4.46,5.02,4.26,4.46,4.16,4.49,4.21,4.57,
         4.87,4.44,4.45,4.41,4.87,4.56,4.15,4.16,4.32,4.06,
         4.12,4.17,3.80,3.96,4.44,4.27,3.90,4.02,4.39,4.52,
         4.25,4.46,4.40,4.83,4.23,4.24,3.95,4.03,4.36,4.07,
         4.17,4.23,4.46,4.38,4.31,4.51,4.13,4.48,5.31,4.58,
         4.81,4.51,4.77,5.33,4.75,4.11,4.76,4.27,4.44,4.20,
         4.71,4.09,4.24,3.90)

ccf<-c(1.901, 2.132, 2.203, 2.228, 2.257, 2.350, 2.361, 2.396, 2.397, 2.445, 2.454, 2.474, 2.518, 2.522,
2.525, 2.532, 2.575, 2.614, 2.616, 2.618, 2.624, 2.659, 2.675, 2.738, 2.740, 2.856, 2.917, 2.928,
2.937, 2.937, 2.977, 2.996, 3.030, 3.125, 3.139, 3.145, 3.220, 3.223, 3.235, 3.243, 3.264, 3.272,
3.294, 3.332, 3.346, 3.377, 3.408, 3.435, 3.493, 3.501, 3.537, 3.554, 3.562, 3.628, 3.852, 3.871,
3.886, 3.971, 4.024, 4.027, 4.225, 4.395, 5.020)
fit_gini <- fit_all_models(y_gini)
fit_CF   <- fit_all_models(CF)
fit_RCC  <- fit_all_models(RCC)
############################################################
# MVCIKw Distribution Example: Wine Data
############################################################

## ---------------------------------------------------------
## STEP 1: Read Wine data
## ---------------------------------------------------------

wine <- read.csv(
  "C:/Users/hp/Desktop/Ph.D Activities/MVCIKw/Wine.csv",
  header = FALSE
)

## Assign column names
colnames(wine) <- c(
  "Class",
  "Alcohol",
  "MalicAcid",
  "Ash",
  "AlcalinityAsh",
  "Magnesium",
  "TotalPhenols",
  "Flavanoids",
  "NonflavanoidPhenols",
  "Proanthocyanins",
  "ColorIntensity",
  "Hue",
  "OD280_OD315",
  "Proline"
)

## Remove cultivar/class label
wine_num <- wine[, -1]


############################################################
## STEP 2: Select the three variables
############################################################

Phenols <- wine_num$TotalPhenols
Flavanoids <- wine_num$Flavanoids
Proanthocyanins <- wine_num$Proanthocyanins


############################################################
## STEP 3: Descriptive statistics
############################################################

summary(Phenols)
summary(Flavanoids)
summary(Proanthocyanins)


############################################################
## STEP 4: Shift variables to satisfy CIKw support y > 1
############################################################

c_shift <- 1

Phenols_shifted <- Phenols + c_shift
Flavanoids_shifted <- Flavanoids + c_shift
Proanthocyanins_shifted <- Proanthocyanins + c_shift


## Check the support
min(Phenols_shifted)
min(Flavanoids_shifted)
min(Proanthocyanins_shifted)


############################################################
## STEP 5: Fit competing marginal distributions
############################################################

fit_phenol <- fit_all_models(Phenols_shifted)

fit_flav <- fit_all_models(Flavanoids_shifted)

fit_pro <- fit_all_models(Proanthocyanins_shifted)


############################################################
## STEP 6: Display fitting results
############################################################

fit_phenol
fit_flav
fit_pro
############################################################
## Plotting functions
############################################################

plot_hist_with_pdfs <- function(y, fit_object, main_title) {
  
  par(mar = c(4, 4, 2, 1))
  
  pars <- attr(fit_object, "fits")
  models <- names(pars)
  
  ## Fixed colors and line types
  model_cols <- c(
    TP   = "black",
    TL   = "red3",
    TEL  = "blue3",
    TW   = "darkgreen",
    TIKw = "orange3",
    IKw1 = "purple4",
    IKw2 = "brown3",
    CIKw = "deepskyblue4",
    LN   = "gray25"
  )
  
  model_lty <- c(
    TP   = 1,
    TL   = 1,
    TEL  = 2,
    TW   = 1,
    TIKw = 2,
    IKw1 = 2,   # dashed
    IKw2 = 3,
    CIKw = 1,   # solid
    LN   = 4
  )
  
  model_lwd <- c(
    TP   = 2,
    TL   = 2,
    TEL  = 2,
    TW   = 2,
    TIKw = 2,
    IKw1 = 3,   # thicker dashed
    IKw2 = 2,
    CIKw = 3,   # thicker solid
    LN   = 2
  )
  
  ## Histogram
  hist(
    y,
    probability = TRUE,
    col = "gray85",
    border = "white",
    main = main_title,
    xlab = "y"
  )
  
  ## Grid
  xg <- seq(
    min(y),
    max(y),
    length.out = 500
  )
  
  ## Add fitted densities
  for (i in seq_along(models)) {
    
    m <- models[i]
    p <- pars[[m]]
    
    dens <- switch(
      m,
      
      TP = dTP(
        xg,
        exp(p[1]),
        exp(p[2]),
        2 * plogis(p[3]) - 1
      ),
      
      TL = dTL(
        xg,
        exp(p[1]),
        exp(p[2]),
        2 * plogis(p[3]) - 1
      ),
      
      TEL = dTEL(
        xg,
        exp(p[1]),
        exp(p[2]),
        exp(p[3]),
        2 * plogis(p[4]) - 1
      ),
      
      TW = dTW(
        xg,
        exp(p[1]),
        exp(p[2]),
        2 * plogis(p[3]) - 1
      ),
      
      TIKw = dTIKw(
        xg,
        exp(p[1]),
        exp(p[2]),
        2 * plogis(p[3]) - 1
      ),
      
      IKw1 = dIKw1(
        xg,
        exp(p[1]),
        exp(p[2])
      ),
      
      IKw2 = dIKw2(
        xg,
        exp(p[1]),
        exp(p[2])
      ),
      
      CIKw = dCIKw(
        xg,
        p[1],
        p[2],
        p[3],
        p[4]
      ),
      
      LN = dLN(
        xg,
        p[1],
        p[2]
      )
    )
    
    lines(
      xg,
      dens,
      col = model_cols[m],
      lty = model_lty[m],
      lwd = model_lwd[m]
    )
  }
  
  legend(
    "topright",
    legend = models,
    col = model_cols[models],
    lty = model_lty[models],
    lwd = model_lwd[models],
    cex = 0.75,
    bty = "n"
  )
}


############################################################
## ECDF + fitted CDFs
############################################################

plot_ecdf_with_cdfs <- function(y, fit_object, main_title) {
  
  par(mar = c(4, 4, 2, 1))
  
  pars <- attr(fit_object, "fits")
  models <- names(pars)
  
  ## Same fixed colors and line types
  model_cols <- c(
    TP   = "black",
    TL   = "red3",
    TEL  = "blue3",
    TW   = "darkgreen",
    TIKw = "orange3",
    IKw1 = "purple4",
    IKw2 = "brown3",
    CIKw = "deepskyblue4",
    LN   = "gray25"
  )
  
  model_lty <- c(
    TP   = 1,
    TL   = 1,
    TEL  = 2,
    TW   = 1,
    TIKw = 2,
    IKw1 = 2,
    IKw2 = 3,
    CIKw = 1,
    LN   = 4
  )
  
  model_lwd <- c(
    TP   = 2,
    TL   = 2,
    TEL  = 2,
    TW   = 2,
    TIKw = 2,
    IKw1 = 3,
    IKw2 = 2,
    CIKw = 3,
    LN   = 2
  )
  
  ## Empirical CDF
  plot(
    ecdf(y),
    main = main_title,
    xlab = "y",
    ylab = "CDF",
    lwd = 2
  )
  
  ## Grid
  xg <- seq(
    min(y),
    max(y),
    length.out = 500
  )
  
  ## Add fitted CDFs
  for (i in seq_along(models)) {
    
    m <- models[i]
    p <- pars[[m]]
    
    Fhat <- switch(
      m,
      
      TP = pTP(
        xg,
        exp(p[1]),
        exp(p[2]),
        2 * plogis(p[3]) - 1
      ),
      
      TL = pTL(
        xg,
        exp(p[1]),
        exp(p[2]),
        2 * plogis(p[3]) - 1
      ),
      
      TEL = pTEL(
        xg,
        exp(p[1]),
        exp(p[2]),
        exp(p[3]),
        2 * plogis(p[4]) - 1
      ),
      
      TW = pTW(
        xg,
        exp(p[1]),
        exp(p[2]),
        2 * plogis(p[3]) - 1
      ),
      
      TIKw = pTIKw(
        xg,
        exp(p[1]),
        exp(p[2]),
        2 * plogis(p[3]) - 1
      ),
      
      IKw1 = pIKw1(
        xg,
        exp(p[1]),
        exp(p[2])
      ),
      
      IKw2 = pIKw2(
        xg,
        exp(p[1]),
        exp(p[2])
      ),
      
      CIKw = pCIKw(
        xg,
        p[1],
        p[2],
        p[3],
        p[4]
      ),
      
      LN = pLN(
        xg,
        p[1],
        p[2]
      )
    )
    
    lines(
      xg,
      Fhat,
      col = model_cols[m],
      lty = model_lty[m],
      lwd = model_lwd[m]
    )
  }
  
  legend(
    "bottomright",
    legend = models,
    col = model_cols[models],
    lty = model_lty[models],
    lwd = model_lwd[models],
    cex = 0.75,
    bty = "n"
  )
}


############################################################
## Dataset list
############################################################

plot_sets <- list(
  
  "Total Phenols" = list(
    y = Phenols_shifted,
    fit = fit_phenol
  ),
  
  "Flavanoids" = list(
    y = Flavanoids_shifted,
    fit = fit_flav
  ),
  
  "Proanthocyanins" = list(
    y = Proanthocyanins_shifted,
    fit = fit_pro
  )
)


############################################################
## Create PDF
############################################################

graphics.off()

pdf(
  "Wine_MVCIKw_Marginal_Fits.pdf",
  width = 7,
  height = 5
)

for (name in names(plot_sets)) {
  
  y <- plot_sets[[name]]$y
  fit <- plot_sets[[name]]$fit
  
  ## Histogram + fitted PDFs
  plot_hist_with_pdfs(
    y,
    fit,
    paste(name, "- Histogram and Fitted PDFs")
  )
  
  ## ECDF + fitted CDFs
  plot_ecdf_with_cdfs(
    y,
    fit,
    paste(name, "- Empirical and Fitted CDFs")
  )
}

dev.off()

############################################################
## MVCIKw REAL-DATA APPLICATION: WINE DATA
## Three-variable analysis
## TotalPhenols, Flavanoids, Proanthocyanins
############################################################

## Required packages
library(copula)
library(stats)
library(xtable)

############################################################
## STEP 1: LOAD AND PREPARE WINE DATA
############################################################

wine <- read.csv(
  "C:/Users/hp/Desktop/Ph.D Activities/MVCIKw/Wine.csv",
  header = FALSE
)

colnames(wine) <- c(
  "Class",
  "Alcohol",
  "MalicAcid",
  "Ash",
  "AlcalinityAsh",
  "Magnesium",
  "TotalPhenols",
  "Flavanoids",
  "NonflavanoidPhenols",
  "Proanthocyanins",
  "ColorIntensity",
  "Hue",
  "OD280_OD315",
  "Proline"
)

## Select the three variables
Y_raw <- wine[, c(
  "TotalPhenols",
  "Flavanoids",
  "Proanthocyanins"
)]

############################################################
## CIKw support transformation
############################################################

## CIKw requires y > 1
c_shift <- 1

Y <- data.frame(
  TotalPhenols =
    Y_raw$TotalPhenols + c_shift,
  
  Flavanoids =
    Y_raw$Flavanoids + c_shift,
  
  Proanthocyanins =
    Y_raw$Proanthocyanins + c_shift
)

Y <- as.matrix(Y)

colnames(Y) <- c(
  "TotalPhenols",
  "Flavanoids",
  "Proanthocyanins"
)

summary(Y)


############################################################
## STEP 2: MARGINAL CIKw FITTING
############################################################

fit_phenol <- fit_all_models(Y[, "TotalPhenols"])
fit_flav   <- fit_all_models(Y[, "Flavanoids"])
fit_pro    <- fit_all_models(Y[, "Proanthocyanins"])

print(fit_phenol)
print(fit_flav)
print(fit_pro)


############################################################
## Extract CIKw parameter estimates
############################################################

get_cikw_parameters <- function(fit_object) {
  
  pars <- attr(fit_object, "fits")
  
  if (is.null(pars$CIKw)) {
    stop("CIKw fit was not found.")
  }
  
  p <- pars$CIKw
  
  names(p) <- c(
    "alpha",
    "beta",
    "lambda1",
    "lambda2"
  )
  
  return(p)
}

theta1 <- get_cikw_parameters(fit_phenol)
theta2 <- get_cikw_parameters(fit_flav)
theta3 <- get_cikw_parameters(fit_pro)

theta1
theta2
theta3


############################################################
## STEP 3: KENDALL DEPENDENCE
############################################################

kendall_tau <- cor(
  Y,
  method = "kendall"
)

cat("\nKendall's tau matrix:\n")
print(round(kendall_tau, 4))


############################################################
## Convert Kendall's tau to Gaussian copula correlation
############################################################

R_tau <- sin(
  (pi / 2) * kendall_tau
)

## Force exact diagonal = 1
diag(R_tau) <- 1

cat("\nInitial Gaussian-copula correlation matrix:\n")
print(round(R_tau, 4))


############################################################
## Check positive definiteness
############################################################

eigen(R_tau)$values

if (any(eigen(R_tau)$values <= 0)) {
  
  warning(
    "R_tau is not positive definite. ",
    "A near-positive-definite correction may be required."
  )
  
}


############################################################
## STEP 4: GAUSSIAN COPULA FIT
############################################################

## Convert observations to pseudo-observations
## using fitted CIKw marginal CDFs

U <- matrix(
  NA_real_,
  nrow = nrow(Y),
  ncol = 3
)

colnames(U) <- colnames(Y)

U[, 1] <- pCIKw(
  Y[, 1],
  theta1[1],
  theta1[2],
  theta1[3],
  theta1[4]
)

U[, 2] <- pCIKw(
  Y[, 2],
  theta2[1],
  theta2[2],
  theta2[3],
  theta2[4]
)

U[, 3] <- pCIKw(
  Y[, 3],
  theta3[1],
  theta3[2],
  theta3[3],
  theta3[4]
)

## Numerical protection
U <- pmin(
  pmax(U, 1e-10),
  1 - 1e-10
)

head(U)


############################################################
## Fit Gaussian copula
############################################################

gaussian_cop <- normalCopula(
  param = c(
    R_tau[1, 2],
    R_tau[1, 3],
    R_tau[2, 3]
  ),
  dim = 3,
  dispstr = "un"
)

fit_gaussian <- fitCopula(
  gaussian_cop,
  data = U,
  method = "ml"
)

summary(fit_gaussian)


############################################################
## Gaussian copula fitted correlation matrix
############################################################

R_hat <- getSigma(
  fit_gaussian@copula
)

cat("\nEstimated Gaussian copula correlation matrix:\n")
print(round(R_hat, 4))


############################################################
## STEP 5: OPTIONAL COPULA COMPARISON
############################################################
## Compare Gaussian, Clayton and Frank copulas
## where the selected copula is appropriate for the data.

## Gaussian
loglik_gaussian <- as.numeric(logLik(fit_gaussian))

k_gaussian <- length(
  coef(fit_gaussian)
)

n <- nrow(U)

AIC_gaussian <- -2 * loglik_gaussian +
  2 * k_gaussian

BIC_gaussian <- -2 * loglik_gaussian +
  k_gaussian * log(n)


############################################################
## Clayton copula
############################################################

fit_clayton <- tryCatch({
  
  clayton_cop <- claytonCopula(
    dim = 3
  )
  
  fitCopula(
    clayton_cop,
    data = U,
    method = "ml"
  )
  
}, error = function(e) NULL)


if (!is.null(fit_clayton)) {
  
  loglik_clayton <- as.numeric(
    logLik(fit_clayton)
  )
  
  k_clayton <- length(
    coef(fit_clayton)
  )
  
  AIC_clayton <- -2 * loglik_clayton +
    2 * k_clayton
  
  BIC_clayton <- -2 * loglik_clayton +
    k_clayton * log(n)
  
} else {
  
  loglik_clayton <- NA
  AIC_clayton <- NA
  BIC_clayton <- NA
}


############################################################
## Frank copula
############################################################

fit_frank <- tryCatch({
  
  frank_cop <- frankCopula(
    dim = 3
  )
  
  fitCopula(
    frank_cop,
    data = U,
    method = "ml"
  )
  
}, error = function(e) NULL)


if (!is.null(fit_frank)) {
  
  loglik_frank <- as.numeric(
    logLik(fit_frank)
  )
  
  k_frank <- length(
    coef(fit_frank)
  )
  
  AIC_frank <- -2 * loglik_frank +
    2 * k_frank
  
  BIC_frank <- -2 * loglik_frank +
    k_frank * log(n)
  
} else {
  
  loglik_frank <- NA
  AIC_frank <- NA
  BIC_frank <- NA
}


############################################################
## Copula comparison table
############################################################

copula_comparison <- data.frame(
  Copula = c(
    "Gaussian",
    "Clayton",
    "Frank"
  ),
  
  LogLik = c(
    loglik_gaussian,
    loglik_clayton,
    loglik_frank
  ),
  
  AIC = c(
    AIC_gaussian,
    AIC_clayton,
    AIC_frank
  ),
  
  BIC = c(
    BIC_gaussian,
    BIC_clayton,
    BIC_frank
  )
)

print(
    copula_comparison)


############################################################
## STEP 6: FULL MVCIKw JOINT LOG-LIKELIHOOD
############################################################

## This is the important joint likelihood:
##
## log c_R(U_i) +
## log f_1(Y_i1) +
## log f_2(Y_i2) +
## log f_3(Y_i3)

mvcikw_loglik <- function(
    Y,
    theta_list,
    R
) {
  
  n <- nrow(Y)
  
  ## Marginal parameters
  th1 <- theta_list[[1]]
  th2 <- theta_list[[2]]
  th3 <- theta_list[[3]]
  
  ## Marginal CDFs
  U <- cbind(
    
    pCIKw(
      Y[, 1],
      th1[1],
      th1[2],
      th1[3],
      th1[4]
    ),
    
    pCIKw(
      Y[, 2],
      th2[1],
      th2[2],
      th2[3],
      th2[4]
    ),
    
    pCIKw(
      Y[, 3],
      th3[1],
      th3[2],
      th3[3],
      th3[4]
    )
  )
  
  ## Numerical protection
  U <- pmin(
    pmax(U, 1e-10),
    1 - 1e-10
  )
  
  ## Gaussian copula
  cop <- normalCopula(
    param = c(
      R[1, 2],
      R[1, 3],
      R[2, 3]
    ),
    dim = 3,
    dispstr = "un"
  )
  
  ## Copula log-density
  log_cop <- dCopula(
    U,
    copula = cop,
    log = TRUE
  )
  
  ## Marginal log-densities
  log_f1 <- log(
    pmax(
      dCIKw(
        Y[, 1],
        th1[1],
        th1[2],
        th1[3],
        th1[4]
      ),
      1e-300
    )
  )
  
  log_f2 <- log(
    pmax(
      dCIKw(
        Y[, 2],
        th2[1],
        th2[2],
        th2[3],
        th2[4]
      ),
      1e-300
    )
  )
  
  log_f3 <- log(
    pmax(
      dCIKw(
        Y[, 3],
        th3[1],
        th3[2],
        th3[3],
        th3[4]
      ),
      1e-300
    )
  )
  
  sum(
    log_cop +
      log_f1 +
      log_f2 +
      log_f3
  )
}


############################################################
## Full-data MVCIKw likelihood
############################################################

theta_list <- list(
  theta1,
  theta2,
  theta3
)

full_loglik <- mvcikw_loglik(
  Y,
  theta_list,
  R_hat
)

cat(
  "\nFull MVCIKw log-likelihood =",
  full_loglik,
  "\n"
)


############################################################
## STEP 7: GLRT CHANGE-POINT SCAN
############################################################

## Minimum segment size
minseg <- 20

## Candidate change-points
tau_grid <- seq(
  minseg,
  n - minseg,
  by = 5
)


############################################################
## Function for fitting one segment
############################################################

fit_mvcikw_segment <- function(Yseg) {
  
  fit1 <- fit_all_models(
    Yseg[, 1]
  )
  
  fit2 <- fit_all_models(
    Yseg[, 2]
  )
  
  fit3 <- fit_all_models(
    Yseg[, 3]
  )
  
  th1 <- get_cikw_parameters(fit1)
  th2 <- get_cikw_parameters(fit2)
  th3 <- get_cikw_parameters(fit3)
  
  theta <- list(
    th1,
    th2,
    th3
  )
  
  ## Pseudo-observations
  U <- cbind(
    
    pCIKw(
      Yseg[, 1],
      th1[1],
      th1[2],
      th1[3],
      th1[4]
    ),
    
    pCIKw(
      Yseg[, 2],
      th2[1],
      th2[2],
      th2[3],
      th2[4]
    ),
    
    pCIKw(
      Yseg[, 3],
      th3[1],
      th3[2],
      th3[3],
      th3[4]
    )
  )
  
  U <- pmin(
    pmax(U, 1e-10),
    1 - 1e-10
  )
  
  ## Estimate Gaussian copula
  cop <- normalCopula(
    param = c(
      0.2,
      0.2,
      0.2
    ),
    dim = 3,
    dispstr = "un"
  )
  
  copfit <- tryCatch(
    
    fitCopula(
      cop,
      U,
      method = "ml"
    ),
    
    error = function(e) NULL
  )
  
  if (is.null(copfit)) {
    return(NULL)
  }
  
  Rseg <- getSigma(
    copfit@copula
  )
  
  ## Joint likelihood
  ll <- mvcikw_loglik(
    Yseg,
    theta,
    Rseg
  )
  
  list(
    theta = theta,
    R = Rseg,
    loglik = ll
  )
}


############################################################
## GLRT scan
############################################################

scan_results <- data.frame(
  tau = tau_grid,
  GLRT = NA_real_
)

for (i in seq_along(tau_grid)) {
  
  tau <- tau_grid[i]
  
  cat(
    "Processing tau =",
    tau,
    "\n"
  )
  
  Y1 <- Y[
    1:tau,
    ,
    drop = FALSE
  ]
  
  Y2 <- Y[
    (tau + 1):n,
    ,
    drop = FALSE
  ]
  
  fit_seg1 <- tryCatch(
    fit_mvcikw_segment(Y1),
    error = function(e) NULL
  )
  
  fit_seg2 <- tryCatch(
    fit_mvcikw_segment(Y2),
    error = function(e) NULL
  )
  
  if (
    !is.null(fit_seg1) &&
    !is.null(fit_seg2)
  ) {
    
    scan_results$GLRT[i] <-
      2 * (
        fit_seg1$loglik +
          fit_seg2$loglik -
          full_loglik
      )
  }
}


############################################################
## Best change-point
############################################################

best_index <- which.max(
  scan_results$GLRT
)

best_tau <- scan_results$tau[
  best_index
]

best_GLRT <- scan_results$GLRT[
  best_index
]

cat(
  "\nEstimated change-point =",
  best_tau,
  "\n"
)

cat(
  "Maximum GLRT =",
  best_GLRT,
  "\n"
)


############################################################
## STEP 8: GLRT PROFILE PLOT
############################################################

pdf(
  "Wine_GLRT_Profile.pdf",
  width = 7,
  height = 5
)

plot(
  scan_results$tau,
  scan_results$GLRT,
  type = "b",
  pch = 19,
  lwd = 2,
  xlab = "Candidate change-point",
  ylab = expression(G[n](tau)),
  main = "MVCIKw GLRT Change-point Profile"
)

abline(
  v = best_tau,
  lty = 2,
  lwd = 2
)

text(
  best_tau,
  best_GLRT,
  labels = paste0(
    "CP = ",
    best_tau
  ),
  pos = 4
)

dev.off()


############################################################
## STEP 9: SEGMENT PROFILE PLOT
############################################################

## Standardize variables only for visualization
## (not for model fitting)

Y_std <- scale(Y)

pdf(
  "Wine_MVCIKw_Detected_Segments.pdf",
  width = 8,
  height = 6
)

matplot(
  Y_std,
  type = "l",
  lty = 1,
  lwd = 1.5,
  xlab = "Observation order",
  ylab = "Standardized value",
  main = "MVCIKw GLRT-BS Detected Change-point"
)

abline(
  v = best_tau,
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = colnames(Y),
  lty = 1,
  lwd = 2,
  bty = "n"
)

dev.off()


############################################################
## Segment-specific plots
############################################################

pdf(
  "Wine_MVCIKw_Segment_Profiles.pdf",
  width = 8,
  height = 6
)

plot(
  Y_std[, 1],
  type = "l",
  lwd = 2,
  xlab = "Observation order",
  ylab = "Standardized value",
  main = "Wine MVCIKw Segment Profile"
)

lines(
  Y_std[, 2],
  lwd = 2,
  lty = 2
)

lines(
  Y_std[, 3],
  lwd = 2,
  lty = 3
)

abline(
  v = best_tau,
  lty = 2,
  lwd = 2
)

legend(
  "topright",
  legend = colnames(Y),
  lty = 1:3,
  lwd = 2,
  bty = "n"
)

dev.off()


############################################################
## STEP 10: BIC AND MIC
############################################################

## Function to calculate BIC and MIC
## q = number of parameters per marginal
## p = number of variables
## m = number of change-points

calculate_model_criteria <- function(
    loglik,
    n,
    p = 3,
    m = 0,
    q = 4,
    beta_MIC = 2
) {
  
  ## Marginal parameters
  marginal_parameters <-
    (m + 1) * p * q
  
  ## Gaussian copula parameters
  copula_parameters <-
    (m + 1) * p * (p - 1) / 2
  
  k <-
    marginal_parameters +
    copula_parameters
  
  BIC <- -2 * loglik +
    k * log(n)
  
  MIC <- -2 * loglik +
    k * log(n) +
    beta_MIC * m * log(n)
  
  data.frame(
    CPs = m,
    LogLik = loglik,
    Parameters = k,
    BIC = BIC,
    MIC = MIC
  )
}


############################################################
## No-change model
############################################################

criteria_0 <- calculate_model_criteria(
  loglik = full_loglik,
  n = n,
  p = 3,
  m = 0
)


############################################################
## One-change-point model
############################################################

Y1 <- Y[
  1:best_tau,
  ,
  drop = FALSE
]

Y2 <- Y[
  (best_tau + 1):n,
  ,
  drop = FALSE
]

fit1 <- fit_mvcikw_segment(Y1)
fit2 <- fit_mvcikw_segment(Y2)

two_segment_loglik <-
  fit1$loglik +
  fit2$loglik


criteria_1 <- calculate_model_criteria(
  loglik = two_segment_loglik,
  n = n,
  p = 3,
  m = 1
)


############################################################
## Compare model complexity
############################################################

criteria_table <- rbind(
  criteria_0,
  criteria_1
)

print(
  round(
    criteria_table,
    4
  )
)


############################################################
## Plot BIC and MIC
############################################################

pdf(
  "Wine_BIC_MIC.pdf",
  width = 7,
  height = 5
)

plot(
  criteria_table$CPs,
  criteria_table$BIC,
  type = "b",
  pch = 19,
  lwd = 2,
  xlab = "Number of change-points",
  ylab = "Information criterion",
  ylim = range(
    c(
      criteria_table$BIC,
      criteria_table$MIC
    )
  ),
  main = "BIC and MIC for Wine Data"
)

lines(
  criteria_table$CPs,
  criteria_table$MIC,
  type = "b",
  pch = 17,
  lwd = 2,
  lty = 2
)

legend(
  "topright",
  legend = c(
    "BIC",
    "MIC"
  ),
  lty = c(1, 2),
  pch = c(19, 17),
  lwd = 2,
  bty = "n"
)

dev.off()


############################################################
## UI E-VALUE
############################################################

## IMPORTANT:
## This is a simple GLRT-based diagnostic.
## For the manuscript, use your formally defined
## split-sample/cross-fit UI e-value implementation.

log_e_value <- best_GLRT / 2

UI_E_value <- exp(
  min(
    log_e_value,
    700
  )
)

cat(
  "\nGLRT =",
  best_GLRT,
  "\n"
)

cat(
  "Log UI-e diagnostic =",
  log_e_value,
  "\n"
)

cat(
  "UI-e diagnostic =",
  UI_E_value,
  "\n"
)


############################################################
## FINAL SUMMARY
############################################################

cat("\n====================================\n")
cat("      MVCIKw WINE DATA RESULTS\n")
cat("====================================\n")

cat(
  "Sample size:",
  n,
  "\n"
)

cat(
  "Number of variables:",
  ncol(Y),
  "\n"
)

cat(
  "Detected change-point:",
  best_tau,
  "\n"
)

cat(
  "Maximum GLRT:",
  best_GLRT,
  "\n"
)

cat(
  "Full-data log-likelihood:",
  full_loglik,
  "\n"
)

cat(
  "Two-segment log-likelihood:",
  two_segment_loglik,
  "\n"
)

cat(
  "BIC (no CP):",
  criteria_0$BIC,
  "\n"
)

cat(
  "BIC (1 CP):",
  criteria_1$BIC,
  "\n"
)

cat(
  "MIC (no CP):",
  criteria_0$MIC,
  "\n"
)

cat(
  "MIC (1 CP):",
  criteria_1$MIC,
  "\n"
)

cat(
  "====================================\n"
)





############################################################
## Wine data: detected change-point
############################################################

tau_hat <- 130

Y <- cbind(
  TotalPhenols      = Phenols_shifted,
  Flavanoids        = Flavanoids_shifted,
  Proanthocyanins   = Proanthocyanins_shifted
)

n <- nrow(Y)

Y1 <- Y[1:tau_hat, , drop = FALSE]
Y2 <- Y[(tau_hat + 1):n, , drop = FALSE]

cat("Segment 1:", nrow(Y1), "observations\n")
cat("Segment 2:", nrow(Y2), "observations\n")

############################################################
## Segment-wise CIKw MLE
############################################################

fit_segment_cikw <- function(Yseg) {
  
  results <- list()
  
  for (j in 1:ncol(Yseg)) {
    
    fit <- fit_all_models(Yseg[, j])
    
    cikw_row <- fit[fit$Model == "CIKw", ]
    
    results[[colnames(Yseg)[j]]] <- cikw_row
  }
  
  results
}

fit_seg1 <- fit_segment_cikw(Y1)
fit_seg2 <- fit_segment_cikw(Y2)

fit_seg1
fit_seg2
############################################################
## Extract CIKw parameter estimates
############################################################

extract_cikw_par <- function(y) {
  
  fit <- fit_all_models(y)
  
  pars <- attr(fit, "fits")$CIKw
  
  data.frame(
    alpha   = pars[1],
    beta    = pars[2],
    lambda1 = pars[3],
    lambda2 = pars[4]
  )
}

MLE_seg1 <- do.call(
  rbind,
  lapply(1:ncol(Y1), function(j)
    extract_cikw_par(Y1[, j])
  )
)

MLE_seg2 <- do.call(
  rbind,
  lapply(1:ncol(Y2), function(j)
    extract_cikw_par(Y2[, j])
  )
)

rownames(MLE_seg1) <- colnames(Y1)
rownames(MLE_seg2) <- colnames(Y2)

MLE_seg1
MLE_seg2

############################################################
## GLRT at tau = 130
############################################################
############################################################
## Segment log-likelihood for the MVCIKw marginal model
############################################################

segment_loglik <- function(Y) {
  
  Y <- as.matrix(Y)
  
  total_logLik <- 0
  
  for (j in 1:ncol(Y)) {
    
    y <- Y[, j]
    
    ## Fit all candidate marginal models
    fit <- fit_all_models(y)
    
    ## Extract CIKw log-likelihood
    cikw_row <- fit[fit$Model == "CIKw", ]
    
    if (nrow(cikw_row) == 0) {
      stop("CIKw model was not found in fit_all_models().")
    }
    
    total_logLik <- total_logLik + cikw_row$LogLik
  }
  
  return(total_logLik)
}

ll_full <- segment_loglik(Y)

ll_seg1 <- segment_loglik(Y1)
ll_seg2 <- segment_loglik(Y2)
############################################################
## Two segments: tau = 130
############################################################

tau_hat <- 130

Y1 <- Y[1:130, , drop = FALSE]
Y2 <- Y[131:nrow(Y), , drop = FALSE]

ll_seg1 <- segment_loglik(Y1)
ll_seg2 <- segment_loglik(Y2)

cat("Segment 1 log-likelihood =", ll_seg1, "\n")
cat("Segment 2 log-likelihood =", ll_seg2, "\n")
GLRT_130 <- 2 * (ll_seg1 + ll_seg2 - ll_full)

cat("Full log-likelihood =", ll_full, "\n")
cat("Segment 1 log-likelihood =", ll_seg1, "\n")
cat("Segment 2 log-likelihood =", ll_seg2, "\n")
cat("GLRT at tau = 130 =", GLRT_130, "\n")

############################################################
## Full-sample CIKw parameters
############################################################

full_cikw_par <- lapply(
  1:ncol(Y),
  function(j) extract_cikw_par(Y[, j])
)

names(full_cikw_par) <- colnames(Y)

full_cikw_par
############################################################
## Parametric bootstrap for C_alpha
############################################################

set.seed(123)

B_boot <- 1000
alpha_level <- 0.05

n <- nrow(Y)
p <- ncol(Y)

minseg <- 20
step <- 5

tau_grid <- seq(
  minseg,
  n - minseg,
  by = step
)

GLRT_boot <- numeric(B_boot)

for (b in 1:B_boot) {
  
  ## Generate data under H0: no change
  Y_boot <- rmvcikw(
    n = n,
    p = p,
    theta_list = full_cikw_par,
    psi = R_hat
  )
  
  ## Calculate full likelihood
  ll0 <- segment_loglik(Y_boot)
  
  ## GLRT scan
  glrt_values <- numeric(length(tau_grid))
  
  for (i in seq_along(tau_grid)) {
    
    tau <- tau_grid[i]
    
    Y_left <- Y_boot[1:tau, , drop = FALSE]
    Y_right <- Y_boot[(tau + 1):n, , drop = FALSE]
    
    ll_left <- segment_loglik(Y_left)
    ll_right <- segment_loglik(Y_right)
    
    glrt_values[i] <-
      2 * (ll_left + ll_right - ll0)
  }
  
  GLRT_boot[b] <- max(glrt_values)
  
  cat(
    "\rBootstrap replication:",
    b,
    "of",
    B_boot
  )
}

cat("\n")

############################################################
## Critical value C_alpha
############################################################

C_alpha <- quantile(
  GLRT_boot,
  probs = 1 - alpha_level,
  na.rm = TRUE
)

C_alpha

############################################################
## GLRT decision
############################################################

GLRT_obs <- GLRT_130

decision <- ifelse(
  GLRT_obs > C_alpha,
  "Reject H0: change-point detected",
  "Do not reject H0"
)

cat("Observed GLRT =", GLRT_obs, "\n")
cat("C_alpha =", C_alpha, "\n")
cat("Decision:", decision, "\n")

#UI e-value
############################################################
## UI E-value
############################################################

log_E <- GLRT_obs / 2

E_value <- exp(
  min(log_E, 700)
)

cat("Log E-value =", log_E, "\n")
cat("E-value =", E_value, "\n")
############################################################
## Bootstrap UI E-values
############################################################

log_E_boot <- GLRT_boot / 2

P_E_gt_1 <- mean(log_E_boot > log(1), na.rm = TRUE)

P_E_gt_20 <- mean(log_E_boot > log(20), na.rm = TRUE)

P_E_gt_100 <- mean(log_E_boot > log(100), na.rm = TRUE)

E_summary <- data.frame(
  Measure = c(
    "Expected E-value",
    "P(E > 1)",
    "P(E > 20)",
    "P(E > 100)"
  ),
  Value = c(
    mean(exp(pmin(log_E_boot, 700)), na.rm = TRUE),
    P_E_gt_1,
    P_E_gt_20,
    P_E_gt_100
  )
)

print(E_summary)
############################################################
## Final Wine change-point summary
############################################################

wine_cp_summary <- data.frame(
  ChangePoint = tau_hat,
  Segment1 = paste0("1-", tau_hat),
  Segment2 = paste0(tau_hat + 1, "-", n),
  Full_LogLik = ll_full,
  Segmented_LogLik = ll_seg1 + ll_seg2,
  GLRT = GLRT_obs,
  C_alpha = as.numeric(C_alpha),
  E_value = E_value,
  P_E_gt_1 = P_E_gt_1,
  P_E_gt_20 = P_E_gt_20,
  P_E_gt_100 = P_E_gt_100
)

print(wine_cp_summary)

############################################################
## MVCIKw random generator
## Gaussian copula + CIKw margins
############################################################

rmvcikw <- function(n, theta_list, R) {
  
  ## Number of dimensions
  p <- length(theta_list)
  
  ## Gaussian copula
  Z <- MASS::mvrnorm(
    n = n,
    mu = rep(0, p),
    Sigma = R
  )
  
  ## Convert to uniform variables
  U <- pnorm(Z)
  
  ## Generate CIKw margins
  Y <- matrix(NA, nrow = n, ncol = p)
  
  for (j in 1:p) {
    
    theta <- theta_list[[j]]
    
    Y[, j] <- qcikw(
      U[, j],
      alpha   = theta$alpha,
      beta    = theta$beta,
      lambda1 = theta$lambda1,
      lambda2 = theta$lambda2
    )
  }
  
  colnames(Y) <- names(theta_list)
  
  return(Y)
}
