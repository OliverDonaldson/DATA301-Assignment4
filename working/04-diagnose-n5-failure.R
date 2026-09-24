set.seed(20260922)
q8 <- function(n,p){h<-(n+1/3)*p+1/3;j<-floor(h);list(j=j,g=h-j)}
boot_star <- function(xs, B) {
  n <- length(xs); idx <- sample.int(n, n*B, TRUE)
  cnt <- matrix(tabulate(rep((0:(B-1))*n, each=n)+idx, nbins=n*B), nrow=B, byrow=TRUE)
  cum <- t(apply(cnt,1,cumsum)); a<-q8(n,.25); b<-q8(n,.75)
  pick <- function(k) xs[max.col(cum>=k, ties.method="first")]
  log(3)/((pick(b$j)+b$g*(pick(b$j+1)-pick(b$j))) - (pick(a$j)+a$g*(pick(a$j+1)-pick(a$j))))
}
lam_iqr <- function(xs){n<-length(xs);a<-q8(n,.25);b<-q8(n,.75)
  log(3)/((xs[b$j]+b$g*(xs[b$j+1]-xs[b$j]))-(xs[a$j]+a$g*(xs[a$j+1]-xs[a$j])))}

cat("=== n=5, lambda=1: anatomy of the bootstrap correction (R=5000) ===\n")
R <- 5000; n <- 5
hat <- bcor <- mstar <- numeric(R)
for (r in 1:R) {
  xs <- sort(rexp(n,1)); hat[r] <- lam_iqr(xs)
  s <- boot_star(xs, 200); mstar[r] <- mean(s[is.finite(s)])
  bcor[r] <- 2*hat[r] - mstar[r]
}
cat(sprintf("lam_iqr : mean=%.3f  median=%.3f  max=%.1f  P(>5)=%.3f\n", mean(hat), median(hat), max(hat), mean(hat>5)))
cat(sprintf("mean(*) : mean=%.3f  median=%.3f  max=%.1f\n", mean(mstar), median(mstar), max(mstar)))
cat(sprintf("lam_tilde: mean=%.2f median=%.3f  min=%.1f  P(<0)=%.3f  P(<-1)=%.3f\n",
    mean(bcor), median(bcor), min(bcor), mean(bcor<0), mean(bcor< -1)))
cat("\nMSE contributions (lambda=1):\n")
se <- (bcor-1)^2; o <- order(se, decreasing=TRUE)
cat(sprintf("  total MSE=%.1f ; top 1 rep contributes %.1f%% ; top 10 contribute %.1f%%\n",
    mean(se), 100*se[o[1]]/sum(se), 100*sum(se[o[1:10]])/sum(se)))
cat(sprintf("  trimmed (drop top 1%%) MSE=%.3f   median sq err=%.3f\n",
    mean(se[o[-(1:ceiling(R*.01))]]), median(se)))

cat("\n=== is MSE finite? MSE vs R, lambda=1, n=5 ===\n")
for (Rk in c(1000,5000,20000,80000)) {
  set.seed(7)
  v <- replicate(Rk, lam_iqr(sort(rexp(5,1))))
  cat(sprintf("  R=%6d  IQR: mse=%8.3f  max=%8.1f | ML: mse=%.4f (analytic 0.5833)\n",
      Rk, mean((v-1)^2), max(v), { set.seed(7); mean((replicate(Rk, 1/mean(rexp(5,1)))-1)^2) }))
}
