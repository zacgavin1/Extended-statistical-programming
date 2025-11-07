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
  
  # I think this is wrong and that mine works. I haven't deleted incase you disagree
  #n <- length(t)
  #X <- matrix(0, 150, K)
  #for (i in 1:150) {
    #for (j in 1:min(29 + i, 80)) {
      #if (30 + i - j > 0 && 30 + i - j <= 150) {
      #  X[i, ] <- X[i, ] + X_tilde[30 + i - j, ] * pd[j]
     # }
    #}
  #}
  
  # note: the last column of X_tilde does not affect X. That is because the 
  # new infections on the last day do not affect the deaths on the last day
  # (the model assumes you cannot die of the disease on the same day as infection)
  
  # rows of X give the effect of each spline on that fitted value
  # columns give the effect that spline has on each fitted value
  
  X <- matrix(0, max(t)-min(t)+1, K)
  for (i in 1:(max(t)-min(t)+1)){
    if (i<=51){ # if min(29+i, 80)= 29+i
      X[i,] = t(pd[(29+i):1]) %*% X_tilde[1:(29+i),] 
      
    } else if(i>51){ # if min(29+i, 80)=80
      X[i,] = t(pd[80:1]) %*% X_tilde[(-50+i):(29+i), ]
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

X_tilde <- out$X_tilde
X <- out$X
S <- out$S
pd <- out$pd
lambda <- 10^-5




#### ----- Question 2 ----- ####

pnll <- function(gamma, y, X, lambda, S, w=rep(1,150)){
  beta <- exp(gamma)
  nlogl <- ( - t(w*y) %*% log(X %*% beta) + t(w) %*% X %*% beta
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
plot(1:100/8-4, y_plt, type="l", xlab='const gamma')


# Finite differencing check - using package
library(numDeriv)

x <- rep(0,80)
num_deriv <- grad(function(g) pnll(g, y, X, lambda, S), g_mle$par)
anal_deriv <- d_nll(g_mle$par, y, X, lambda, S)






#### ----- Question 3 ----- ####

# Fit the model - ie use optim to find the mle for gamma

# try without grad first (optim will finite diff the derivs)
g_mle <- optim(par=rep(0,80), fn=pnll,  y=y, X=X, 
               lambda=lambda, S=S, method='BFGS',  control = list(maxit = 5000))
# this gives a pretty good min (checking using numDeriv.grad)

#### BIGGEST REMAINING TASK - get this working with analytic derivative. 
# Need to redo d_nll in Q2


# estimate beta, mu and f from the mle for gamma, and matrices X and X_tilde
b_hat <- exp(g_mle$par)
mu <- X %*% b_hat       
f <- X_tilde %*% b_hat


# plotting overlaid graphs
plot(data$julian, data$nhs, cex=.5, pch=19, col='blue', xlab='time', ylab='', 
     xlim=c(30,220))
points(data$julian, mu, cex=.5, pch=19, col='red') 
points((min(data$julian)-30):max(data$julian), f, cex=.5, pch=19, col='green')

# Remarks: this has now passed the sanity check (after code setting X in Q1 redone)
# f is very 'wiggly'. A stronger penalty is likely required



#### --- Question 4 --- ####
 
# Choosing lambda - once previous parts sorted, this should
# be straightforwards to get working properly

H <- function(lambda, X, mu_hat, S, y){
  W <- diag(drop(y/(mu_hat)^2))
  t(X) %*%W %*%X + lambda*S
}

EDF <- function(H0, H_l){
  # note that H_l is symmetric, as t(X)WX is symm, and so is S. CHOLESKY!
  A <- chol(H_l)
  ATI <- forwardsolve(t(A), diag(rep(1, 80)))
  Hl_inv <- backsolve(A, ATI)
  sum(diag(Hl_inv%*%H0))
}

test_range <- exp(seq(-13,-7,length=50))
BIC <- rep(0, 50)
i<-1
for (lambda in test_range){
  
  # compute the estimator of gamma, thus beta, thus mu for that lambda
  g_mle <- optim(par=rep(0,80), fn=pnll, y=y, X=X, 
                 lambda=lambda, S=S, method='BFGS')
  b_hat <- exp(g_mle$par)
  mu_hat <- X %*% b_hat
  
  ## Computing the H matrices
  # The fitted values mu_hat are under penalty par. lambda, even for H0
  Hl <- H(lambda, X, mu_hat, S, y)
  H0 <- H(0, X, mu_hat, S, y)
  
  
  # Compute its BIC score
  # pnll has penalty zero here
  BIC[i] <- 2*pnll(g_mle$par, y,X, 0, S) + log(length(y))*EDF(H0, Hl)
  
  # we'll want to add this to some vector
  # we'll want to minimise BIC - ie higher EDF will be a penalty
  
  print(i)
  i <- i+1
}

plot(log(test_range), BIC, xlab="log(lambda)")
bi <- which(BIC==min(BIC))
lambda_opt <- test_range[bi]

lambda_opt
# note: this is larger than the initial guess of 10^-5, suggesting
# I might have been right to say I thought the wiggly f meant we 
# needed a bigger penalty

# we do still get a wiggly f here tho so not sure whats going on there. It is
# mildly less wiggly than initial guess 

#### ------ Question 5 ------ ####

n <- length(y)

n_rep <- 200
g_mle <- matrix(0, n_rep, 80)
for (i in 1:n_rep){
  wb <- tabulate(sample(n,replace=TRUE), n)
  
  # calculate the sample mle
  min <- optim(par=rep(0,80), fn=pnll,  y=y, X=X, 
                 lambda=lambda_opt, S=S, w=wb, method='BFGS',  control = list(maxit = 5000))
  g_mle[i,] <- min$par
  print(i)
}

beta_boot <- exp(g_mle)
f_boot <- X_tilde %*% t(beta_boot)

CI <- apply(f_boot, MARGIN=1, FUN=quantile, probs=c(0.025, 0.975))


##### ---- Question 6 ---- #######
# Now plot all this information on a graph

# finding the actual prediction using the actual data, lambda=lambda_opt
g_mle <- optim(par=rep(0,80), fn=pnll,  y=y, X=X, 
               lambda=lambda_opt, S=S, method='BFGS',  control = list(maxit = 5000))
b_hat <- exp(g_mle$par)
mu <- X %*% b_hat
f <- X_tilde %*% b_hat

# plot actual deaths, and fitted deaths 
plot(data$julian, data$nhs, cex=.5, pch=19, col='blue', xlab='time', ylab='', 
     xlim=c(30,220), ylim=c(0, 1700))
points(data$julian, mu, cex=.5, pch=19, col='red')

# now plot infection predictions with empirical CI
points((min(data$julian)-30):max(data$julian), f, cex=.5, pch=19, col='green')
lines((min(data$julian)-30):max(data$julian), CI[1,])
lines((min(data$julian)-30):max(data$julian), CI[2,])


# note this currently runs much too slowly. I expect this to be fixed once 
# optim is implemented with an analytic derivative, but will be worth 
# rechecking in on the speed when this is fixed.

# Want to double check lambda_opt as this has changed for a reason I can't 
# determine

# It's a bit weird that the infection predictions are wavy - this might
# be something to investigate
