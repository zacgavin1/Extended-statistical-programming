library(splines)

data <- read.table("engcov.txt", header=T, stringsAsFactor=T)


# t should
modelling <- function(t, K) {

  range_t <- range(t)
  internal_knots <- seq(range_t[1]-30, range_t[2], length.out = K - 2)
  diff <- internal_knots[2]-internal_knots[1]
  all_knots <- c(range_t[1]-30-3*diff, range_t[1]-30-2*diff, range_t[1]-30-diff, internal_knots,
                 range_t[2]+diff, range_t[2]+2*diff, range_t[2]+3*diff)

  X_tilde <- splineDesign(all_knots, (min(t)-30):max(t) , ord = 4)
  
  d <- 1:80; edur <- 3.151; sdur <- .469
  pd <- dlnorm(d, edur, sdur); pd <- pd / sum(pd)
  
  # make sure to find a way to put no. days in here rather than hard coding 150
  n <- length(t)
  X <- matrix(0, 150, K)
  for (i in 1:150) {
    for (j in 1:min(29 + i, 80)) {
      if (30 + i - j > 0 && 30 + i - j <= 150) {
        X[i, ] <- X[i, ] + X_tilde[30 + i - j, ] * pd[j]
      }
    }
  }
  
  S <- crossprod(diff(diag(K), diff = 2))
  list(X_tilde = X_tilde, X = X, S = S, pd = pd)
}

# t contains the start and end days of the data we have
t <- c(min(data$julian),max(data$julian))
K <- 80
out <- modelling(t, K)
dim(out$X_tilde)
print(out)

#### ----- Question 2 ----- ####

pnll <- function(gamma, y, X, lambda, S){
  beta <- exp(gamma)
  nlogl <- ( - t(y) %*% log(X %*% beta) + sum(X %*% beta )
           + .5*lambda * t(beta) %*% S %*% beta )

  as.numeric(nlogl)
}

# test
pnll(rep(0, 80), y, X, lambda, S)

d_nll <- function(gamma, y, X, lambda, S){
  beta <- exp(gamma)
  mu <- X %*% beta
  F <- diag(drop(y/mu-1)) %*% X %*% diag(beta)  # each term here is a matrix   
  deriv <- (apply(F, MARGIN=2, FUN=sum)
            + diag(beta) %*% S %*% beta)
  deriv
}

y <- data$nhs
pnll(gamma, y, X, lambda, S)
d_nll(gamma, y, X, lambda, S)

## plotting to get an idea of what pen log likelihood looks like for constant gamma
y_plt <- rep(0,100)
for (i in 1:100){
  y_plt[i] <- pnll(rep(i/8-4, 80), y, X, lambda, S)
}
plot(1:100/8-4, y_plt, type="l")


# Finite differencing check - to do!


#### ----- Question 3 ----- ####

# Fit the model - ie use optim to find the mle for gamma

g_mle <- optim(par=rep(0,80), fn=pnll, gr = d_nll, y=y, X=X, 
               lambda=lambda, S=S, method='BFGS')

g_mle

# estimate for beta
b_hat = exp(g_mle$par)


# use this to get fitted values (mu) for deaths on each of the days
mu <- X %*% b_hat

# we find f from beta using X_tilde
f <- X_tilde %*% b_hat

# plotting overlaid graphs
plot(data$julian, data$nhs, cex=.5, pch=19, col='blue', xlab='time', ylab='', 
     xlim=c(30,220))
points(data$julian, mu, cex=.5, pch=19, col='red') # this is not currently correct
# as the above is not correct, this must also be wrong
points((min(data$julian)-30):max(data$julian), f, cex=.5, pch=19, col='green')

# Remark: this did not pass the sanity check, either in time or in th expected 
# kinds of fitted values



#### --- Question 4 --- ####
 
# Choosing lambda - once previous parts sorted, this should
# be straightforwards to get working properly

H <- function(lambda, X, b_hat, S, y){
  mu_hat <- X %*% b_hat
  W <- diag(y/(mu_hat)^2)
  t(X) %*%W %*%X + lambda*S
}
EDF <- function(H0, H_l){
  # this could well be quite expensive
  trace((H_l^-1)%*%H0)
}

for (lambda in some_range){
  # need to get g_mle working for this to work
  # compute the estimator of beta for that lambda
  g_mle <- optim(par=rep(0,80), fn=pnll, gr = d_nll, y=y, X=X, 
                 lambda=lambda, S=S, method='BFGS')
  b_hat <- exp(g_mle)
  
  # Compute its BIC score
  EDF <- EDF(H(0,X,b_hat,S,y), H(lambda, X,b_hat, S,y))
  BIC <- 2*pnll(log(b_hat), y,X, lambda, S) + 2*log(n)*EDF
  
  # we'll want to add this to some vector
  # we'll want to minimise BIC - ie higher EDF will be a penalty
}



#### ------ Questioon 5 ------ ####













