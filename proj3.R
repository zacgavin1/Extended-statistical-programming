library(splines)

modelling <- function(t, K) {
  range_t <- range(t)
  internal_knots <- seq(range_t[1], range_t[2], length.out = K - 2)
  all_knots <- c(seq(range_t[1] - 3, range_t[1], length.out = 2),internal_knots, seq(range_t[2], range_t[2] + 3, length.out = 2))
  
  X_tilde <- splineDesign(all_knots, t, ord = 4, outer.ok = TRUE)
  
  d <- 1:80; edur <- 3.151; sdur <- .469
  pd <- dlnorm(d, edur, sdur); pd <- pd / sum(pd)
  
  n <- length(t)
  X <- matrix(0, n, K)
  for (i in 1:n) {
    for (j in 1:min(29 + i, 80)) {
      if (30 + i - j > 0 && 30 + i - j <= n) {
        X[i, ] <- X[i, ] + X_tilde[30 + i - j, ] * pd[j]
      }
    }
  }
  
  S <- crossprod(diff(diag(K), diff = 2))
  list(X_tilde = X_tilde, X = X, S = S, pd = pd)
}

t <- 1:200
K <- 80
out <- modelling(t, K)
dim(out$X_tilde)
print(out)

##Q2


