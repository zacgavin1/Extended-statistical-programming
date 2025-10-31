library(splines)

### Question1 

modelling <- function(t, K) {

  knots <- seq(min(t), max(t), length.out = K - 2)
  all_knots <- seq(min(t), max(t), length.out = K + 4)
  X_tilde <- splineDesign(all_knots, t, outer.ok=TRUE) ##need outer.ok = TRUE to allow x to be outside the inner knots
  
  
  d <- 1:80; edur <- 3.151; sdur <- .469
  pd <- dlnorm(d, edur, sdur);
  pd <- pd/sum(pd)
  
  n <- length(t)
  X <- matrix(0, n, K)
  for (i in 1:n) {
    for (j in 1:min(29+i, 80)){
      if (30 + i - j > 0 && 30 + i - j <= n) {
        X[i,] <- X[i,] + X_tilde[30 + i - j, ] * pd[j]
      }
    }
  }
  
  S <- crossprod(diff(diag(K), diff = 2))
  list(X_tilde = X_tilde, X = X, S = S)
}

### Question 2

