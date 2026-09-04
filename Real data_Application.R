# ============================================================
# SUPPLEMENTARY ANALYSIS
# CIKw-GLRT Binary Segmentation with PSO
# and Parametric Bootstrap Critical Values
# Datasets:
#   1. COVID-19 mortality
#   2. India Revenue
# Procedure:
#   CIKw likelihood -> PSO estimation -> GLRT
#   -> parametric bootstrap -> Binary Segmentation
#   -> segment summaries and distributional plots
# ============================================================

library(ggplot2)

set.seed(2025)

# ============================================================
# 1. SETTINGS
# ============================================================

n_min <- 20
alpha_glrt <- 0.05
B_boot <- 120

# ============================================================
# 2. CIKw LOG-LIKELIHOOD
# ============================================================

ll_cikw <- function(par, y) {
  
  alpha <- par[1]
  beta  <- par[2]
  lam1  <- par[3]
  lam2  <- par[4]
  
  if (alpha <= 2 || beta <= 0 ||
      lam1 < -1 || lam1 > 1 ||
      lam2 < 0 || lam2 > 2)
    return(-Inf)
  
  z <- y + lam1 * sin(y)
  
  if (any(z <= 0)) return(-Inf)
  
  logC <- log(alpha) +
    (alpha + 1) * log(beta) -
    log(2) - lgamma(alpha)
  
  sum(logC +
        (alpha - 1) * log(z) -
        (beta + lam2) * z)
}

# ============================================================
# 3. PSO PARAMETER ESTIMATION
# ============================================================

pso_fit <- function(y, particles = 20, iterations = 30) {
  
  lb <- c(2.01, 1e-3, -1, 0)
  ub <- c(10, 20, 1, 2)
  
  swarm <- matrix(runif(particles * 4),
                  particles, 4)
  
  for (j in 1:4)
    swarm[, j] <- lb[j] +
    (ub[j] - lb[j]) * swarm[, j]
  
  velocity <- matrix(0, particles, 4)
  
  pbest <- swarm
  pval <- apply(swarm, 1, ll_cikw, y = y)
  
  if (all(!is.finite(pval)))
    return((lb + ub) / 2)
  
  g <- which.max(pval)
  gbest <- pbest[g, ]
  gval <- pval[g]
  
  w <- 0.7
  c1 <- 1.5
  c2 <- 1.5
  
  for (it in 1:iterations) {
    
    for (i in 1:particles) {
      
      r1 <- runif(1)
      r2 <- runif(1)
      
      velocity[i, ] <-
        w * velocity[i, ] +
        c1 * r1 * (pbest[i, ] - swarm[i, ]) +
        c2 * r2 * (gbest - swarm[i, ])
      
      swarm[i, ] <- swarm[i, ] + velocity[i, ]
      
      swarm[i, ] <- pmax(pmin(swarm[i, ], ub), lb)
      
      val <- ll_cikw(swarm[i, ], y)
      
      if (is.finite(val) && val > pval[i]) {
        pbest[i, ] <- swarm[i, ]
        pval[i] <- val
      }
    }
    
    gnew <- which.max(pval)
    
    if (pval[gnew] > gval) {
      gval <- pval[gnew]
      gbest <- pbest[gnew, ]
    }
  }
  
  gbest
}

# ============================================================
# 4. PARAMETRIC BOOTSTRAP CRITICAL VALUE
# ============================================================

bootstrap_critical <- function(y, B = B_boot) {
  
  n <- length(y)
  
  if (n < 2 * n_min)
    return(qchisq(0.95, df = 4))
  
  theta <- pso_fit(y)
  
  alpha <- theta[1]
  beta  <- theta[2]
  lam1  <- theta[3]
  
  Tboot <- numeric(B)
  
  for (b in 1:B) {
    
    # Generate bootstrap data
    z <- rgamma(n, shape = alpha, rate = beta)
    
    x <- numeric(n)
    
    for (i in 1:n) {
      xi <- z[i]
      for (iter in 1:30)
        xi <- z[i] - lam1 * sin(xi)
      x[i] <- xi
    }
    
    theta0 <- pso_fit(x)
    ll0 <- ll_cikw(theta0, x)
    
    Tvals <- numeric(n - 2 * n_min + 1)
    
    for (j in seq_along(Tvals)) {
      
      k <- n_min + j - 1
      
      thL <- pso_fit(x[1:k])
      thR <- pso_fit(x[(k + 1):n])
      
      llL <- ll_cikw(thL, x[1:k])
      llR <- ll_cikw(thR, x[(k + 1):n])
      
      Tvals[j] <- if (
        is.finite(llL) && is.finite(llR)
      ) {
        2 * (llL + llR - ll0)
      } else 0
    }
    
    Tboot[b] <- max(Tvals)
    
    if (b %% 20 == 0)
      cat("Bootstrap:", b, "/", B, "\n")
  }
  
  as.numeric(quantile(
    Tboot,
    1 - alpha_glrt,
    na.rm = TRUE
  ))
}

# ============================================================
# 5. GLRT FOR ONE SEGMENT
# ============================================================

glrt_segment <- function(y, start = 1) {
  
  n <- length(y)
  
  if (n < 2 * n_min)
    return(NULL)
  
  theta0 <- pso_fit(y)
  ll0 <- ll_cikw(theta0, y)
  
  candidates <- n_min:(n - n_min)
  Tvals <- numeric(length(candidates))
  
  for (j in seq_along(candidates)) {
    
    k <- candidates[j]
    
    left <- y[1:k]
    right <- y[(k + 1):n]
    
    thL <- pso_fit(left)
    thR <- pso_fit(right)
    
    llL <- ll_cikw(thL, left)
    llR <- ll_cikw(thR, right)
    
    Tvals[j] <- if (
      is.finite(llL) && is.finite(llR)
    ) {
      2 * (llL + llR - ll0)
    } else -Inf
  }
  
  if (!any(is.finite(Tvals)))
    return(NULL)
  
  i <- which.max(Tvals)
  
  list(
    k = candidates[i],
    k_abs = start + candidates[i] - 1,
    T = Tvals[i]
  )
}

# ============================================================
# 6. BINARY SEGMENTATION
# ============================================================

binary_segmentation <- function(y) {
  
  results <- list()
  cps <- numeric(0)
  
  queue <- list(
    list(y = y, start = 1)
  )
  
  while (length(queue) > 0) {
    
    current <- queue[[1]]
    queue <- queue[-1]
    
    seg <- current$y
    start <- current$start
    
    cat("\nSegment:",
        start, "-", start + length(seg) - 1,
        " n =", length(seg), "\n")
    
    if (length(seg) < 2 * n_min)
      next
    
    glr <- glrt_segment(seg, start)
    
    if (is.null(glr))
      next
    
    cat("Best CP:", glr$k_abs,
        " GLRT =", round(glr$T, 4), "\n")
    
    cat("Computing bootstrap critical value...\n")
    
    ca <- bootstrap_critical(seg)
    
    cat("c_alpha =", round(ca, 4), "\n")
    
    if (glr$T > ca) {
      
      cp <- glr$k_abs
      
      cat("SIGNIFICANT SPLIT at", cp, "\n")
      
      cps <- c(cps, cp)
      
      queue <- append(
        queue,
        list(
          list(
            y = seg[1:glr$k],
            start = start
          ),
          list(
            y = seg[(glr$k + 1):length(seg)],
            start = start + glr$k
          )
        )
      )
      
    } else {
      
      cat("No significant change-point.\n")
    }
  }
  
  sort(unique(cps))
}

# ============================================================
# 7. DATA
# ============================================================

covid <- c(
  8.826,6.105,10.383,7.26,13.220,6.015,10.855,6.122,
  10.685,10.035,5.242,7.630,14.604,7.903,6.327,9.391,
  14.962,4.730,3.215,16.498,11.665,9.284,12.878,6.656,
  3.440,5.854,8.813,10.043,7.260,5.985,4.424,4.344,
  5.143,9.935,7.840,9.550,6.968,6.370,3.537,3.286,
  10.158,8.108,6.697,7.151,6.560,2.988,3.336,6.814,
  8.325,7.854,8.551,3.228,3.499,3.751,7.486,6.625,
  6.140,4.909,4.661,1.867,2.838,5.392,12.042,8.696,
  6.412,3.395,1.815,3.327,5.406,6.182,4.949,4.089,
  3.359,2.070,3.298,5.317,5.442,4.557,4.292,2.500,
  6.535,4.648,4.697,5.459,4.120,3.922,3.219,1.402,
  2.438,3.257,3.632,3.233,3.027,2.352,1.205,2.077,
  3.778,3.218,2.926,2.601,2.065,1.041,1.800,3.029,
  2.058,2.326,2.506,1.923
)

india_revenue <- c(
  8.1811,9.0245,9.1818,8.5971,9.3023,9.6890,9.0261,
  9.2137,9.2691,9.2685,9.3065,10.0754,10.4273,10.4083,
  10.3567,10.4683,9.9935,10.1709,9.8052,8.5352,9.9817,
  9.2268,9.2314,9.0085,8.1131,8.6383,8.8101,8.0794,
  8.6762,9.1080,9.5708,10.0809,11.1293,12.1083,10.9771,
  9.8097,10.3880,10.1773,10.8367,11.0016,9.9846,10.5697,
  11.1477,11.3874,12.0174
)

datasets <- list(
  COVID_19 = covid,
  India_Revenue = india_revenue
)

# ============================================================
# 8. RUN THE SAME PSO-GLRT PROCEDURE FOR BOTH DATASETS
# ============================================================

all_cps <- list()

for (nm in names(datasets)) {
  
  cat("\n========================================\n")
  cat("DATASET:", nm, "\n")
  cat("========================================\n")
  
  cps <- binary_segmentation(datasets[[nm]])
  
  all_cps[[nm]] <- cps
  
  cat("\nDetected change-points:", cps, "\n")
}

# ============================================================
# 9. PLOT DETECTED CHANGE-POINTS
# ============================================================

for (nm in names(datasets)) {
  
  y <- datasets[[nm]]
  cps <- all_cps[[nm]]
  
  df <- data.frame(
    time = seq_along(y),
    value = y
  )
  
  p <- ggplot(df, aes(time, value)) +
    geom_line() +
    geom_point(size = 1.5) +
    geom_vline(
      xintercept = cps,
      linetype = "dashed"
    ) +
    labs(
      title = paste(
        "CIKw-GLRT Binary Segmentation:",
        gsub("_", " ", nm)
      ),
      x = "Observation",
      y = "Value"
    ) +
    theme_minimal()
  
  print(p)
}

# ============================================================
# 10. CREATE SEGMENTS
# ============================================================

make_segments <- function(y, cps) {
  
  br <- c(0, cps, length(y))
  
  lapply(
    seq_len(length(br) - 1),
    function(i)
      y[(br[i] + 1):br[i + 1]]
  )
}

# ============================================================
# 11. SEGMENT SUMMARY
# ============================================================

segment_summary <- function(y, cps) {
  
  br <- c(0, cps, length(y))
  
  do.call(
    rbind,
    lapply(
      seq_len(length(br) - 1),
      function(i) {
        
        x <- y[(br[i] + 1):br[i + 1]]
        
        data.frame(
          Segment = paste0(
            br[i] + 1, "-", br[i + 1]
          ),
          Start = br[i] + 1,
          End = br[i + 1],
          n = length(x),
          Mean = mean(x),
          SD = sd(x),
          Skewness =
            mean((x - mean(x))^3) /
            sd(x)^3
        )
      }
    )
  )
}

# ============================================================
# 12. DISTRIBUTIONAL DENSITY PLOTS
# ============================================================

for (nm in names(datasets)) {
  
  y <- datasets[[nm]]
  cps <- all_cps[[nm]]
  
  segs <- make_segments(y, cps)
  
  d <- do.call(
    rbind,
    lapply(
      seq_along(segs),
      function(i)
        data.frame(
          value = segs[[i]],
          Segment = paste0(
            "Segment ",
            i,
            " (",
            ifelse(i == 1, 1, cps[i - 1] + 1),
            "-",
            ifelse(
              i <= length(cps),
              cps[i],
              length(y)
            ),
            ")"
          )
        )
    )
  )
  
  print(
    ggplot(
      d,
      aes(value, colour = Segment,
          fill = Segment)
    ) +
      geom_density(
        alpha = 0.20,
        linewidth = 0.9
      ) +
      labs(
        title = paste(
          "Distributional Density Across",
          gsub("_", " ", nm),
          "Segments"
        ),
        x = "Value",
        y = "Density"
      ) +
      theme_minimal() +
      theme(
        legend.position = "bottom"
      )
  )
  
  cat("\n\n", nm, "SEGMENT SUMMARY\n")
  print(segment_summary(y, cps))
}

# =============================
# END OF SUPPLEMENTARY ANALYSIS
# =============================