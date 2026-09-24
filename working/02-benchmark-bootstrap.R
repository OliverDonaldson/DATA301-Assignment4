set.seed(2)
# direct type-8 quantile from a sorted vector / sorted matrix rows
q8_idx <- function(n, p) { h <- (n + 1/3) * p + 1/3; j <- floor(h); list(j = j, g = h - j) }

lam_iqr_fast <- function(xs) {          # xs already sorted
  n <- length(xs)
  a <- q8_idx(n, .25); b <- q8_idx(n, .75)
  lo <- xs[a$j] + a$g * (xs[a$j + 1] - xs[a$j])
  hi <- xs[b$j] + b$g * (xs[b$j + 1] - xs[b$j])
  log(3) / (hi - lo)
}

boot_apply_quantile <- function(x, B) {
  n <- length(x); m <- matrix(sample(x, n*B, TRUE), nrow = B)
  q <- t(apply(m, 1, quantile, probs = c(.25,.75), type = 8, names = FALSE))
  log(3)/(q[,2]-q[,1])
}
boot_sort_direct <- function(x, B) {
  n <- length(x); m <- matrix(sample(x, n*B, TRUE), nrow = B)
  ms <- t(apply(m, 1, sort))
  a <- q8_idx(n,.25); b <- q8_idx(n,.75)
  lo <- ms[,a$j] + a$g*(ms[,a$j+1]-ms[,a$j])
  hi <- ms[,b$j] + b$g*(ms[,b$j+1]-ms[,b$j])
  log(3)/(hi-lo)
}
# column-sort trick: sort whole matrix by row using order on a keyed vector is slow;
# instead resample INDICES into sorted x, then use counting via tabulate
boot_counting <- function(x, B) {
  n <- length(x); xs <- sort(x)
  idx <- sample.int(n, n*B, TRUE)
  off <- rep((0:(B-1))*n, each = n)
  cnt <- matrix(tabulate(off + idx, nbins = n*B), nrow = B, byrow = TRUE)
  cum <- t(apply(cnt, 1, cumsum))
  a <- q8_idx(n,.25); b <- q8_idx(n,.75)
  pick <- function(k) xs[max.col(cum >= k, ties.method = "first")]
  lo <- pick(a$j) + a$g*(pick(a$j+1) - pick(a$j))
  hi <- pick(b$j) + b$g*(pick(b$j+1) - pick(b$j))
  log(3)/(hi-lo)
}
for (n in c(5,30,100)) {
  x <- rexp(n,1); B <- 200
  set.seed(9); t1 <- system.time(for(i in 1:30) boot_apply_quantile(x,B))[["elapsed"]]/30
  set.seed(9); t2 <- system.time(for(i in 1:30) boot_sort_direct(x,B))[["elapsed"]]/30
  set.seed(9); t3 <- system.time(for(i in 1:30) boot_counting(x,B))[["elapsed"]]/30
  set.seed(9); v2 <- boot_sort_direct(x,B); set.seed(9); v3 <- boot_counting(x,B)
  cat(sprintf("n=%3d  quantile()=%.2fms  sort+direct=%.2fms  counting=%.2fms   agree(sort,count)=%s\n",
      n, t1*1000, t2*1000, t3*1000, isTRUE(all.equal(sort(v2), sort(v3)))))
}
# correctness of direct formula vs stats::quantile
x <- rexp(37,1); xs <- sort(x)
cat("\ndirect vs quantile(type=8): ", lam_iqr_fast(xs),
    log(3)/diff(quantile(x, c(.25,.75), type=8, names=FALSE)), "\n")
