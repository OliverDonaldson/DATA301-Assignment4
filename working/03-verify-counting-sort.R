set.seed(11)
q8_idx <- function(n,p){h<-(n+1/3)*p+1/3;j<-floor(h);list(j=j,g=h-j)}
# same index set -> both routes must agree exactly
for (n in c(5,6,7,30,100)) {
  x <- rexp(n,1); xs <- sort(x); B <- 500
  idx <- sample.int(n, n*B, TRUE)
  # route A: materialise, sort each row, direct formula
  m <- matrix(xs[idx], nrow=B, byrow=TRUE)
  ms <- t(apply(m,1,sort))
  a <- q8_idx(n,.25); b <- q8_idx(n,.75)
  loA <- ms[,a$j]+a$g*(ms[,a$j+1]-ms[,a$j]); hiA <- ms[,b$j]+b$g*(ms[,b$j+1]-ms[,b$j])
  A <- log(3)/(hiA-loA)
  # route B: counting sort
  off <- rep((0:(B-1))*n, each=n)
  cnt <- matrix(tabulate(off+idx, nbins=n*B), nrow=B, byrow=TRUE)
  cum <- t(apply(cnt,1,cumsum))
  pick <- function(k) xs[max.col(cum>=k, ties.method="first")]
  loB <- pick(a$j)+a$g*(pick(a$j+1)-pick(a$j)); hiB <- pick(b$j)+b$g*(pick(b$j+1)-pick(b$j))
  B_ <- log(3)/(hiB-loB)
  # route C: stats::quantile on the materialised rows
  C <- apply(m,1,function(r) log(3)/diff(quantile(r,c(.25,.75),type=8,names=FALSE)))
  cat(sprintf("n=%3d  A==B %s   A==C %s   #Inf=%d\n", n,
      isTRUE(all.equal(A,B_)), isTRUE(all.equal(A,C)), sum(!is.finite(A))))
}
