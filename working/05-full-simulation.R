q8i <- function(n, p, m) { h <- switch(m, "7" = (n-1)*p + 1, "8" = (n + 1/3)*p + 1/3)
                           j <- floor(h); list(j = j, g = h - j) }
lam_ml  <- function(xs) 1/mean(xs)
lam_iqr <- function(xs, m = "8") { n <- length(xs); a <- q8i(n,.25,m); b <- q8i(n,.75,m)
  log(3)/((xs[b$j] + b$g*(xs[b$j+1]-xs[b$j])) - (xs[a$j] + a$g*(xs[a$j+1]-xs[a$j]))) }
boot_star <- function(xs, B, m = "8") {
  n <- length(xs); idx <- sample.int(n, n*B, TRUE)
  cum <- t(apply(matrix(tabulate(rep((0:(B-1))*n, each=n) + idx, nbins = n*B),
                        nrow = B, byrow = TRUE), 1, cumsum))
  a <- q8i(n,.25,m); b <- q8i(n,.75,m)
  pk <- function(k) xs[max.col(cum >= k, ties.method = "first")]
  log(3)/((pk(b$j) + b$g*(pk(b$j+1)-pk(b$j))) - (pk(a$j) + a$g*(pk(a$j+1)-pk(a$j))))
}
run_cell <- function(lambda, n, R, B = 200) {
  ml <- i8 <- i7 <- bc <- numeric(R); nbad <- 0L
  for (r in seq_len(R)) {
    xs <- sort(rexp(n, lambda))
    ml[r] <- lam_ml(xs); i8[r] <- lam_iqr(xs, "8"); i7[r] <- lam_iqr(xs, "7")
    s <- boot_star(xs, B); ok <- is.finite(s); nbad <- nbad + sum(!ok)
    bc[r] <- 2*i8[r] - mean(s[ok])
  }
  stat <- function(e) { se <- (e - lambda)^2
    c(mse = mean(se), mcse = sd(se)/sqrt(R), top1 = max(se)/sum(se)) }
  list(lambda = lambda, n = n, R = R,
       ML = stat(ml), IQR8 = stat(i8), IQR7 = stat(i7), BC = stat(bc),
       p_neg_bc = mean(bc < 0), bad_rate = nbad/(R*B), med_bc = median(bc))
}
args <- commandArgs(TRUE); R <- as.integer(args[1])
t0 <- Sys.time(); set.seed(20260922)
res <- list()
for (lambda in c(.1, 1, 10)) for (n in c(5, 30, 100))
  res[[sprintf("%g_%d", lambda, n)]] <- run_cell(lambda, n, R)
cat(sprintf("R=%d  elapsed %.1f s\n\n", R, as.numeric(difftime(Sys.time(), t0, units="secs"))))
for (k in names(res)) { z <- res[[k]]
  cat(sprintf("lam=%-4g n=%3d | ML %8.4g (%.1f%%) | IQR8 %8.4g (%.1f%%) | IQR7 %8.4g | BC %10.4g (%.1f%%) top1=%.0f%% Pneg=%.2f\n",
    z$lambda, z$n, z$ML["mse"], 100*z$ML["mcse"]/z$ML["mse"], z$IQR8["mse"], 100*z$IQR8["mcse"]/z$IQR8["mse"],
    z$IQR7["mse"], z$BC["mse"], 100*z$BC["mcse"]/z$BC["mse"], 100*z$BC["top1"], z$p_neg_bc)) }
saveRDS(res, sprintf("res_%d.rds", R))
