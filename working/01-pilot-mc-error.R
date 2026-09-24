set.seed(1)
lam_ml  <- function(x) 1 / mean(x)
lam_iqr <- function(x, type = 8) {
  q <- quantile(x, c(.25, .75), type = type, names = FALSE)
  log(3) / (q[2] - q[1])
}
# vectorised bootstrap: B resamples at once
lam_bc <- function(x, B = 200, type = 8) {
  n <- length(x); hat <- lam_iqr(x, type)
  m <- matrix(sample(x, n * B, replace = TRUE), nrow = B)
  q <- t(apply(m, 1, quantile, probs = c(.25, .75), type = type, names = FALSE))
  star <- log(3) / (q[, 2] - q[, 1])
  ok <- is.finite(star)
  c(est = 2 * hat - mean(star[ok]), nbad = sum(!ok))
}

cat("=== timing one replication ===\n")
for (n in c(5, 30, 100)) {
  x <- rexp(n, 1)
  t0 <- system.time(for (i in 1:20) lam_bc(x))[["elapsed"]] / 20
  cat(sprintf("n=%3d  %.1f ms/rep  -> %.1f s per 1000 reps\n", n, t0 * 1000, t0 * 1000))
}

cat("\n=== pilot MSE, R=2000, lambda=1 ===\n")
R <- 2000
for (n in c(5, 30, 100)) {
  ml <- iqr <- bc <- numeric(R); bad <- 0
  for (r in 1:R) {
    x <- rexp(n, 1)
    ml[r] <- lam_ml(x); iqr[r] <- lam_iqr(x)
    o <- lam_bc(x); bc[r] <- o[1]; bad <- bad + o[2]
  }
  f <- function(e) { se <- (e - 1)^2; c(mse = mean(se), mcse = sd(se)/sqrt(R), rel = sd(se)/sqrt(R)/mean(se)) }
  cat(sprintf("n=%3d  ML  mse=%8.4f mcse=%7.4f rel=%.3f\n", n, f(ml)[1], f(ml)[2], f(ml)[3]))
  cat(sprintf("       IQR mse=%8.4f mcse=%7.4f rel=%.3f   max=%.1f\n", f(iqr)[1], f(iqr)[2], f(iqr)[3], max(iqr)))
  cat(sprintf("       BC  mse=%8.4f mcse=%7.4f rel=%.3f   finite=%d/%d  bad boot=%.2f%%\n",
      f(bc[is.finite(bc)])[1], f(bc[is.finite(bc)])[2], f(bc[is.finite(bc)])[3],
      sum(is.finite(bc)), R, 100*bad/(R*200)))
}
