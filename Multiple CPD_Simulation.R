# ============================================================
# CIKw-GLRT Binary Segmentation
# R implementation of the supplied Python code
# ============================================================

# Required packages
library(ggplot2)

# ============================================================
# Simulation: true parameters
# ============================================================

true_params <- list(
  c(5.6, 4.8, -0.36, 1.2),
  c(9.0, 5.0,  0.50, 1.5),
  c(12.0, 7.0, -0.8, 1.0)
)

# ============================================================
# CIKw data generation
# ============================================================

simulate_cikw <- function(n, params, seed = NULL) {
  
  if (!is.null(seed))
    set.seed(seed)
  
  a  <- params[1]
  b  <- params[2]
  l1 <- params[3]
  l2 <- params[4]
  
  u <- runif(n)
  
  inner <- (1 - (1 - u)^(1 / b))^(-1 / a) - 1
  
  inner <- pmax(inner, 1e-12)
  
  # Same transformation as supplied Python code
  x <- (inner^(1 / l2)) / (l1 + 1e-9)
  
  return(x)
}


# ============================================================
# Generate data
# ============================================================

n1 <- 120
n2 <- 150
n3 <- 130

set.seed(123)

segment1 <- simulate_cikw(n1, true_params[[1]])
segment2 <- simulate_cikw(n2, true_params[[2]])
segment3 <- simulate_cikw(n3, true_params[[3]])

data <- c(segment1, segment2, segment3)

n <- length(data)


# ============================================================
# Stable log-likelihood
# ============================================================

EPS <- 1e-12

ll_cikw <- function(params, x) {
  
  a  <- params[1]
  b  <- params[2]
  l1 <- params[3]
  l2 <- params[4]
  
  # Parameter constraints
  if (a <= 2 ||
      b <= 0 ||
      l2 <= 0 ||
      l2 > 2 ||
      l1 < -1 ||
      l1 > 1) {
    return(-1e20)
  }
  
  # Safe x
  x_safe <- pmax(x, EPS)
  
  # z = 1 + lambda1 * x^lambda2
  x_pow_l2 <- x_safe^l2
  
  z <- 1 + l1 * x_pow_l2
  
  # Domain checks
  if (any(!is.finite(z)) || any(z <= EPS))
    return(-1e20)
  
  # 1 - z^(-a)
  z_pow_minus_a <- z^(-a)
  
  one_minus_z <- 1 - z_pow_minus_a
  
  if (any(!is.finite(one_minus_z)))
    return(-1e20)
  
  # For beta - 1 < 0
  if ((b - 1) < 0 && any(one_minus_z <= EPS))
    return(-1e20)
  
  if (any(one_minus_z <= 0))
    return(-1e20)
  
  # Reject lambda1 = 0
  if (isTRUE(all.equal(l1, 0)))
    return(-1e20)
  
  # Same restriction used in supplied Python code
  if (l1 < 0)
    return(-1e20)
  
  # Log likelihood components
  log_const <-
    log(a) +
    log(b) +
    log(l2) +
    log(abs(l1) + EPS)
  
  log_x_part <- (l2 - 1) * log(x_safe)
  
  log_z_part <- (-a - 1) * log(z)
  
  log_f1 <- log_const + log_x_part + log_z_part
  
  log_f2 <- (b - 1) * log(one_minus_z)
  
  logf <- log_f1 + log_f2
  
  if (any(!is.finite(logf)))
    return(-1e20)
  
  ll <- sum(logf)
  
  return(ll)
}


# ============================================================
# MLE: constrained and safe
# ============================================================

fit_mle <- function(x, p0 = NULL, verbose = FALSE) {
  
  neg_ll <- function(p) {
    
    v <- ll_cikw(p, x)
    
    if (!is.finite(v) || v < -1e19) {
      
      penalty <- 1e12 +
        sum(pmax(
          0,
          c(2, 0, -1, 0) - p
        ) * 1e3)
      
      return(penalty)
    }
    
    return(-v)
  }
  
  
  # Initial values
  if (is.null(p0))
    p0 <- c(4, 3, 0.2, 1.0)
  
  
  # Bounds
  lower <- c(
    2.0001,  # alpha
    1e-3,    # beta
    1e-3,    # lambda1
    1e-2     # lambda2
  )
  
  upper <- c(
    25.0,
    25.0,
    1.0,
    2.0
  )
  
  
  # L-BFGS-B optimization
  res <- tryCatch(
    optim(
      par = p0,
      fn = neg_ll,
      method = "L-BFGS-B",
      lower = lower,
      upper = upper,
      control = list(maxit = 2000)
    ),
    error = function(e) NULL
  )
  
  
  # If optimization fails, try another starting point
  if (is.null(res) ||
      res$convergence != 0) {
    
    p0 <- c(6.0, 2.0, 0.5, 0.8)
    
    res <- tryCatch(
      optim(
        par = p0,
        fn = neg_ll,
        method = "L-BFGS-B",
        lower = lower,
        upper = upper,
        control = list(maxit = 2000)
      ),
      error = function(e) NULL
    )
  }
  
  
  # Final fallback
  if (is.null(res))
    return(list(
      params = p0,
      loglik = -1e20
    ))
  
  
  params_hat <- res$par
  ll_hat <- -res$value
  
  
  if (!is.finite(ll_hat) || ll_hat < -1e19)
    ll_hat <- -1e20
  
  
  # Check whether parameters hit bounds
  hits <- list()
  
  for (j in 1:length(params_hat)) {
    
    if (isTRUE(all.equal(params_hat[j], lower[j])) ||
        isTRUE(all.equal(params_hat[j], upper[j]))) {
      
      hits[[length(hits) + 1]] <-
        c(params_hat[j], lower[j], upper[j])
    }
  }
  
  
  if (verbose && length(hits) > 0) {
    
    print("Warning: MLE hit bounds for parameters:")
    
    print(hits)
  }
  
  
  return(list(
    params = params_hat,
    loglik = ll_hat
  ))
}


# ============================================================
# GLRT statistic across candidate change-points
# ============================================================

glrt_stat <- function(x, min_seg = 30) {
  
  n <- length(x)
  
  full_fit <- fit_mle(x)
  
  params_full <- full_fit$params
  ll_full <- full_fit$loglik
  
  stats <- list()
  
  ks <- min_seg:(n - min_seg - 1)
  
  for (k in ks) {
    
    x1 <- x[1:k]
    x2 <- x[(k + 1):n]
    
    fit1 <- fit_mle(x1)
    fit2 <- fit_mle(x2)
    
    ll1 <- fit1$loglik
    ll2 <- fit2$loglik
    
    
    if (!is.finite(ll1) ||
        !is.finite(ll2) ||
        !is.finite(ll_full)) {
      
      T <- 0
      
    } else {
      
      T <- 2 * (ll1 + ll2 - ll_full)
      
      if (T < 0)
        T <- 0
    }
    
    stats[[length(stats) + 1]] <-
      c(k = k, T = T)
  }
  
  
  stats <- do.call(rbind, stats)
  
  return(list(
    stats = stats,
    params_full = params_full,
    ll_full = ll_full
  ))
}


# ============================================================
# Binary Segmentation
# ============================================================

binary_segmentation <- function(
    x,
    start = 0,
    min_seg = 30,
    threshold = NULL) {
  
  n <- length(x)
  
  if (n < 2 * min_seg)
    return(integer(0))
  
  
  result <- glrt_stat(
    x,
    min_seg = min_seg
  )
  
  stats <- result$stats
  
  
  if (nrow(stats) == 0)
    return(integer(0))
  
  
  ks <- stats[, "k"]
  Ts <- stats[, "T"]
  
  
  max_index <- which.max(Ts)
  
  rel_k <- ks[max_index]
  
  T_max <- max(Ts)
  
  abs_k <- start + rel_k
  
  
  # Chi-square threshold
  if (is.null(threshold))
    threshold <- qchisq(
      0.95,
      df = 4
    )
  
  
  if (T_max <= threshold)
    return(integer(0))
  
  
  # Left segment
  left_cps <- binary_segmentation(
    x = x[1:rel_k],
    start = start,
    min_seg = min_seg,
    threshold = threshold
  )
  
  
  # Right segment
  right_cps <- binary_segmentation(
    x = x[(rel_k + 1):n],
    start = abs_k,
    min_seg = min_seg,
    threshold = threshold
  )
  
  
  return(
    c(
      left_cps,
      abs_k,
      right_cps
    )
  )
}


# ============================================================
# Estimate parameters for final segments
# ============================================================

estimate_segments <- function(data, cps) {
  
  edges <- c(
    0,
    sort(cps),
    length(data)
  )
  
  seg_params <- list()
  
  
  for (i in 1:(length(edges) - 1)) {
    
    s <- edges[i]
    e <- edges[i + 1]
    
    # Convert Python-style [s:e] to R indexing
    segment_data <-
      data[(s + 1):e]
    
    fit <- fit_mle(
      segment_data,
      verbose = TRUE
    )
    
    seg_params[[i]] <- list(
      range = c(s, e),
      params = fit$params,
      loglik = fit$loglik
    )
  }
  
  
  return(seg_params)
}


# ============================================================
# Empirical threshold using bootstrap
# ============================================================

empirical_threshold <- function(
    data,
    B = 120,
    min_seg = 30,
    seed = 999) {
  
  set.seed(seed)
  
  boot_max <- numeric(B)
  
  n <- length(data)
  
  
  for (i in 1:B) {
    
    # Same as Python permutation bootstrap
    perm <- sample(data, size = n, replace = FALSE)
    
    result <- glrt_stat(
      perm,
      min_seg = min_seg
    )
    
    stats <- result$stats
    
    
    if (nrow(stats) == 0) {
      
      boot_max[i] <- 0
      
    } else {
      
      boot_max[i] <- max(
        stats[, "T"]
      )
    }
  }
  
  
  C_alpha <- as.numeric(
    quantile(
      boot_max,
      probs = 0.95
    )
  )
  
  
  return(list(
    threshold = C_alpha,
    boot_max = boot_max
  ))
}


# ============================================================
# Run detection
# ============================================================

min_seg <- 30

cat("Computing GLRT across k for plotting...\n")

start_time <- Sys.time()

result <- glrt_stat(
  data,
  min_seg = min_seg
)

elapsed <- Sys.time() - start_time

stats <- result$stats

params_full <- result$params_full
ll_full <- result$ll_full

ks <- stats[, "k"]
Ts <- stats[, "T"]

Tmax <- max(Ts)

k_hat_rel <- ks[which.max(Ts)]


cat(
  "Elapsed (GLRT pass):",
  round(
    as.numeric(elapsed),
    2
  ),
  "seconds\n"
)

cat(
  "Tmax:",
  Tmax,
  " k_hat (relative):",
  k_hat_rel,
  "\n"
)


# ============================================================
# Chi-square threshold
# ============================================================

chi_threshold <- qchisq(
  0.95,
  df = 4
)

cat(
  "Chi-square 95% cutoff (df=4):",
  chi_threshold,
  "\n"
)


# ============================================================
# Bootstrap threshold
# ============================================================

cat(
  "Computing empirical threshold with bootstrap ",
  "(this may take a while)...\n"
)

boot_result <- empirical_threshold(
  data,
  B = 120,
  min_seg = min_seg,
  seed = 999
)

C_alpha <- boot_result$threshold
boot_vals <- boot_result$boot_max


cat(
  "Empirical C_alpha (95%):",
  C_alpha,
  "\n"
)


# ============================================================
# Binary Segmentation
# ============================================================

cps_emp <- binary_segmentation(
  data,
  start = 0,
  min_seg = min_seg,
  threshold = C_alpha
)

cps_chi <- binary_segmentation(
  data,
  start = 0,
  min_seg = min_seg,
  threshold = chi_threshold
)


cat(
  "\nDetected change-points ",
  "(empirical threshold): ",
  paste(cps_emp, collapse = ", "),
  "\n"
)

cat(
  "Detected change-points ",
  "(chi-square threshold): ",
  paste(cps_chi, collapse = ", "),
  "\n"
)


# ============================================================
# Estimate segment parameters
# ============================================================

seg_params_emp <- estimate_segments(
  data,
  cps_emp
)


cat(
  "\nEstimated parameters for segments ",
  "(empirical threshold):\n"
)


for (i in seq_along(seg_params_emp)) {
  
  seg <- seg_params_emp[[i]]
  
  s <- seg$range[1]
  e <- seg$range[2]
  
  cat(
    "Segment ",
    s,
    ":",
    e,
    " length=",
    e - s,
    " -> params (a,b,l1,l2) = ",
    paste(
      round(seg$params, 6),
      collapse = ", "
    ),
    ", log-lik=",
    round(seg$loglik, 2),
    "\n",
    sep = ""
  )
}


# ============================================================
# Plot 1: GLRT statistic across k
# ============================================================

plot_data <- data.frame(
  k = ks,
  T = Ts
)

p1 <- ggplot(
  plot_data,
  aes(x = k, y = T)
) +
  geom_line(linewidth = 0.8) +
  
  geom_hline(
    yintercept = chi_threshold,
    linetype = "dashed"
  ) +
  
  geom_hline(
    yintercept = C_alpha,
    linetype = "dashed"
  ) +
  
  geom_point(
    data = data.frame(
      k = k_hat_rel,
      T = Tmax
    ),
    size = 3
  ) +
  
  geom_vline(
    xintercept = n1,
    linetype = "dotted"
  ) +
  
  geom_vline(
    xintercept = n1 + n2,
    linetype = "dotted"
  ) +
  
  geom_vline(
    xintercept = cps_emp,
    linetype = "solid",
    alpha = 0.7
  ) +
  
  geom_vline(
    xintercept = cps_chi,
    linetype = "dashed",
    alpha = 0.7
  ) +
  
  labs(
    title = "GLRT statistic across k with detected change-points",
    x = "k (split index)",
    y = "T*(k)"
  ) +
  
  theme_minimal()

print(p1)


# ============================================================
# Plot 2: Data with detected change-points
# ============================================================

data_df <- data.frame(
  Index = seq_along(data),
  Value = data
)

p2 <- ggplot(
  data_df,
  aes(
    x = Index,
    y = Value
  )
) +
  
  geom_line(
    linewidth = 0.5
  ) +
  
  geom_point(
    size = 0.8
  ) +
  
  geom_vline(
    xintercept = cps_emp,
    linetype = "solid",
    linewidth = 0.8,
    alpha = 0.8
  ) +
  
  geom_vline(
    xintercept = cps_chi,
    linetype = "dashed",
    linewidth = 0.8,
    alpha = 0.8
  ) +
  
  labs(
    title = "Data with detected change-points",
    x = "Index",
    y = "Value"
  ) +
  
  theme_minimal()

print(p2)