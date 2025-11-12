## GENERAL DESCRIPTION

# The aim of this project is to use data from the year 2020 on daily deaths from 
# COVID-19 in English hospitals in order to estimate the (unobserved) daily counts 
# of new infections that ultimately resulted in these deaths through fitting a 
# deconvolution model. This model is fit through B-splines (smoothly connected 
# curve segments) and a smoothing penalty (to avoid overfitting the variation).
# After choosing an appropriate smoothing parameter for the fit (based on the BIC 
# criterion), we assess the uncertainty of the fit through non-parametric bootstrap 
# confidence intervals. Finally, we visualize our findings in a plot.


library(splines) # splineDesign()
library(numDeriv) # grad() (finite differencing)
library(ggplot2) # visuals


data <- read.table("engcov.txt", header=T, stringsAsFactor=T)


#################################################################
####### ---------- BUILDING THE MODEL MATRICES ---------- #######
#################################################################

# We first want to build the matrices that will be used for fitting the model.
# In particular, these are:
# * X_tilde - the B-spline basis for the infections,
# * X - the model matrix for the deaths, and
# * S - the penalty matrix


# * t = start and end dates (calendar day) 
# * K = the number of basis functions
# modelling() returns X_tilde, X, S, and pd (infection-to-death time probability 
# distribution)
modelling <- function(t, K) {
  
  ## --- X_tilde --- ##
  
  # first and last day
  range_t <- range(t)
  # evenly spaced time intervals (knots) starting 30 days before first death
  # -- accounts for delay between infection and death
  # -- defines B-splines
  internal_knots <- seq(range_t[1] - 30, range_t[2], length.out = K - 2)
  # space between intervals
  diff <- internal_knots[2] - internal_knots[1]
  # extra knots to make proper 4th-order (cubic) spline
  all_knots <- c(range_t[1]-30-3*diff, range_t[1]-30-2*diff, range_t[1]-30-diff, 
                 internal_knots,
                 range_t[2]+diff, range_t[2]+2*diff, range_t[2]+3*diff)
  # spline basis matrix (rows = days, columns = spline basis functions)
  X_tilde <- splineDesign(all_knots, (min(t)-30):max(t) , ord = 4)
  
  ##### should this part be outside of the function and be an input?
  d <- 1:80; edur <- 3.151; sdur <- .469
  pd <- dlnorm(d, edur, sdur); pd <- pd / sum(pd)
  
  
  ## --- X --- ##
  
  # note: the last column of X_tilde does not affect X. That is because the 
  # new infections on the last day do not affect the deaths on the last day
  # (the model assumes you cannot die of the disease on the same day as infection)
  
  # rows of X give the effect of each spline on that fitted value
  # columns give the effect that spline has on each fitted value
  
  X <- matrix(0, max(t) - min(t) + 1, K)
  # X_i = sum (from j = 1 to min(29 + i, 80)) X_tilde_{30 + i - j} * pd(j)
  for (i in 1:(max(t) - min(t) + 1)) {
    if (i <= 51) { # if min(29+i, 80) = 29+i
      X[i,] = t(pd[(29 + i):1]) %*% X_tilde[1:(29 + i),] 
      
    } else if(i > 51) { # if min(29+i, 80)=80
      X[i,] = t(pd[80:1]) %*% X_tilde[(-50 + i):(29 + i),]
    }
  }
  
  ## --- S --- ##
  
  # for calculating smoothing penalty
  S <- crossprod(diff(diag(K), diff = 2))
  
  # output
  list(X_tilde = X_tilde, X = X, S = S, pd = pd)
}



#################################################################
######### ---------- SMOOTHING PENALIZATION ---------- ##########
#################################################################

# We now want to be able to evaluate the penalized negative log-likelihood and
# its derivative vector for the future application of finding an optimal 
# smoothing parameter.


# pnll() returns the penalized negative log-likelihood
# * gamma = some set of K parameters
# * y = deaths
# * X = model matrix
# * lambda = smoothing parameter
# * S = penalty matrix
# * w = weights (used in bootstrapping later)
pnll <- function(gamma, y, X, lambda, S, w = rep(1, 150)) {
  
  # beta = exp(gamma)
  # -- ensures model parameters, beta, are positive
  # -- therefore ensures infection curve f is positive
  beta <- exp(gamma) 
  
  # penalized negative log-likelihood = NLL + P
  # -- NLL = - sum(w*y*log(mu)) + sum(w*mu); mu = X*beta, the death model fit
  # -- P = 0.5*lambda*(beta^T*S*beta)
  nlogl <- - t(w*y) %*% log(X %*% beta) + t(w) %*% X %*% beta + 
    .5*lambda * t(beta) %*% S %*% beta
  
  # convert resulting 1x1 matrix to scalar
  # return PNLL
  as.numeric(nlogl)
}

# d_nll() returns gradient of penalized negative log-likelihood w.r.t. gamma
# * gamma, y, X, lambda, S, w defined as before in pnll()
d_nll <- function(gamma, y, X, lambda, S, w = rep(1, 150)) {
  
  # beta = exp(gamma); mu = X*beta (as in pnll())
  beta <- exp(gamma)
  mu <- X %*% beta
  
  # dNLL/dgamma = diag(y*w/mu - w)*X*diag(beta)
  d_likelihood <- beta * t(y*w/mu - w) %*% X
  #d_likelihood <- apply(as.vector(y/mu - 1)*t(beta*t(X))*w, 2, sum)
  # dP/dgamma = lambda*diag(beta)*S*beta
  d_penalty <- lambda * diag(beta) %*% S %*% beta
  # dPNLL/dgamma = dNLL/dgamma + dP/dgamma
  deriv <- -t(d_likelihood) + d_penalty
  #deriv <- -d_likelihood + d_penalty
  
  # return derivative vector
  deriv
}


##### -------- Testing the Derivative Function -------- #####

# We make sure the derivative function is correct by comparing the output to 
# finite differencing, an approximation of the derivative.

# start and end days of the data
t <- c(min(data$julian), max(data$julian))
# number of basis functions
K <- 80
# X_tilde, X, S, pd
out <- modelling(t, K)
X_tilde <- out$X_tilde; X <- out$X; S <- out$S; pd <- out$pd
# deaths
y <- data$nhs

# lambda for testing purposes
lambda_test <- 1e-05

# find gammas that minimize PNLL using the BFGS method
g_mle_test <- optim(par = rep(0, 80), fn = pnll, gr = d_nll, 
                    y = y, X = X, lambda = lambda_test, S = S, 
                    method = 'BFGS', control = list(maxit = 5000))

# approximate gradient of PNLL
num_deriv <- grad(function(g) pnll(g, y, X, lambda_test, S), g_mle_test$par)
# exact gradient of PNLL
anal_deriv <- d_nll(g_mle_test$par, y, X, lambda_test, S)

# gamma0 <- rep(0, 80)
# num_deriv <- numeric(80)
# eps <- 5e-07
# for (i in 1:length(gamma0)) {
#   gamma1 <- gamma0
#   gamma1[i] <- gamma0[i] + eps
#   pen_nll0 <- pnll(gamma0, y, X, lambda_test, S)
#   pen_nll1 <- pnll(gamma1, y, X, lambda_test, S)
#   num_deriv[i] <- (pen_nll1 - pen_nll0)/eps
# }


# find the most that the resulting gradients differ from each other
max(abs(num_deriv - anal_deriv)) # very small ~ 10^-5, so should be correct 


######## -------- Preliminary Sanity Check -------- #########

# Before proceeding, and as part of finding sane starting values for gamma, we
# fit the model using lambda = 5 * 10^-5 and plot the actual and fitted deaths,
# as well as the fitted infection curve, f.

# for sanity check
lambda_sanity <- 5e-05
# find minimizing gammas for PNLL based on sanity check lambda
g_mle_sanity <- optim(par = rep(0, 80), fn = pnll, gr = d_nll,
                      y = y, X = X, lambda = lambda_sanity, S = S,
                      method = 'BFGS',  control = list(maxit = 5000))


# estimate beta, mu, and f (= X_tilde*beta) from gamma MLE
b_hat_sanity <- exp(g_mle_sanity$par)
mu_sanity <- X %*% b_hat_sanity       
f_sanity <- X_tilde %*% b_hat_sanity


# turn range of potential infection days into data frame for plotting purposes
range_dt <- data.frame(range = (min(data$julian) - 30):max(data$julian))

# plot the sanity check fitted deaths and fitted infections
data |> ggplot(aes(x = julian, y = nhs)) + 
  geom_point() + # actual deaths
  geom_line(aes(x = julian, y = mu_sanity, color = 'blue')) + # fitted deaths
  geom_line(data = range_dt, aes(x = range, y = f_sanity, color = 'red')) + # fitted infs
  theme_bw() + 
  scale_color_discrete(labels = c("Fitted Deaths", "Estimated New Infections")) +
  labs(title = "Sanity Check", subtitle = "(lambda = 5x10^-5)",
       x = "Day of the Year", y = "Counts") +
  theme(plot.subtitle = element_text(size = 10),
        legend.title = element_blank(),
        legend.position = c(0.7, 0.8),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = 'black'))

# Remarks: this has now passed the sanity check
# -- f is very 'wiggly'; a stronger penalty is likely required



#################################################################
###### ---------- CHOICE OF SMOOTHING PARAMETER ---------- ######
#################################################################

# We now want to choose an appropriate smoothing parameter, lambda. This will
# be chosen to minimize the BIC.
 
# calculate H_lambda, the Hessian w.r.t. beta of NLL at beta_hat + lambda*S
# --- H_lambda = X^T*W*X + lambda*S, where W = diag(y/mu_hat^2)
H <- function(lambda, X, mu_hat, S, y){
  W <- diag(drop(y/(mu_hat)^2))
  t(X) %*% W %*% X + lambda*S
}

# calculate effective degrees of freedom (EDF)
# --- EDF = tr(H_lambda^-1)
EDF <- function(H0, H_l){
  # note that H_l is symmetric, as t(X)WX and S are symm -> CHOLESKY!
  A <- chol(H_l)
  ATI <- forwardsolve(t(A), diag(rep(1, 80)))
  Hl_inv <- backsolve(A, ATI)
  sum(diag(Hl_inv %*% H0))
}


# find_lambda_opt() finds lambda that minimizes the BIC
# --- BIC = 2*PNLL + log(n)*EDF
# * test_range = range of lambdas to look over for a minimizer
# * par = set of gamma parameters to obtain gamma MLEs
# * pnll = PNLL evaluating function
# * d_nll = derivative of PNLL evaluating function
# * y, X, S defined as before
find_lambda_opt <- function(test_range, par, pnll, d_nll, y, X, S) {
  
  # initialize list of BIC values
  BIC <- rep(0, length(test_range))
  converges <- c()
  i <- 1
  # calculate BIC for each potential lambda in range
  for (lambda in test_range) {
  
    # compute the gamma MLE -> beta -> mu for that lambda
    g_mle <- optim(par = par, fn = pnll, gr = d_nll, 
                   y = y, X = X, lambda = lambda, S = S, method = 'BFGS')#,
                   #control = list(abstol = 1e-30, maxit = 100000))
    converges <- append(converges, g_mle$convergence)
    b_hat <- exp(g_mle$par)
    mu_hat <- X %*% b_hat
  
    # compute H matrices
    # -- fitted values mu_hat are under penalty par. lambda, even for H0
    Hl <- H(lambda, X, mu_hat, S, y)
    H0 <- H(0, X, mu_hat, S, y)
  
    # compute BIC
    # -- pnll has penalty zero here
    # -- higher EDF will be a penalty
    BIC[i] <- 2*pnll(g_mle$par, y,X, 0, S) + log(length(y))*EDF(H0, Hl)
  
    i <- i+1
  }
  
  # find the lambda that minimizes the BIC
  lambda_opt_index <- BIC == min(BIC)
  lambda_opt <- test_range[lambda_opt_index]
  
  plot(test_range, BIC, type = 'l')
  
  print(converges)
  
  # output
  lambda_opt
}

# search over log(lambda) values
test_range <- exp(seq(-13,-7,length=50))
# get the optimal lambda over this range
lambda_opt <- find_lambda_opt(test_range, par = g_mle_sanity$par, pnll, d_nll, y, X, S)

lambda_opt

# note: this is larger than the initial guess of 5*10^-5, validating claim in 
# sanity check that we needed a bigger penalty



#################################################################
####### ---------- ASSESSING MODEL UNCERTAINTY ---------- #######
#################################################################


#### ------ Question 5 ------ ####

n <- length(y)

n_rep <- 200
g_mle <- matrix(0, n_rep, 80)
for (i in 1:n_rep){
  wb <- tabulate(sample(n,replace=TRUE), n)
  
  # calculate the sample mle
  min <- optim(par=rep(0,80), fn=pnll, gr=d_nll,  y=y, X=X, 
                 lambda=lambda_opt, S=S, w=wb, method='BFGS')
  g_mle[i,] <- min$par
  #print(i)
}

beta_boot <- exp(g_mle)
f_boot <- X_tilde %*% t(beta_boot)

CI <- apply(f_boot, MARGIN=1, FUN=quantile, probs=c(0.025, 0.975))


#potential optimization suggestion?
wb <- matrix(0, n, n_rep)
for (i in 1:n_rep){
   wb[, i] <- tabulate(sample(n, replace = TRUE), n)
}
 
g_mle <- apply(wb, MARGIN = 2, FUN = function (wb) {
  min <- optim(par = rep(0, 80), fn = pnll, gr = d_nll,
                y = y, X = X, lambda = lambda_opt, S = S, w = wb,
                method = 'BFGS')
  min$par
})
 
beta_boot <- exp(g_mle)
f_boot <- X_tilde %*% beta_boot
CI <- apply(f_boot, MARGIN = 1, FUN = quantile, probs = c(0.025, 0.975))



#################################################################
######## ---------- VISUALIZING THE MODEL FIT ---------- ########
#################################################################


##### ---- Question 6 ---- #######
# Now plot all this information on a graph

# finding the actual prediction using the actual data, lambda=lambda_opt
g_mle <- optim(par=rep(0,80), fn=pnll, gr = d_nll,  y=y, X=X, 
               lambda=lambda_opt, S=S, method='BFGS')
b_hat <- exp(g_mle$par)
mu <- X %*% b_hat
f <- X_tilde %*% b_hat


CI_dt <- data.frame(lower = CI[1,], upper = CI[2,])

data |> ggplot(aes(x = julian, y = nhs)) +
  geom_point() + 
  #geom_line() + 
  geom_line(data = data, aes(x = julian, y = mu, color = 'blue')) +
  geom_line(data = range_dt, aes(x = range, y = f, color = 'red')) +
  geom_ribbon(data = range_dt, aes(x = range, y = f, 
                                   ymin = CI_dt$lower, ymax = CI_dt$upper), 
              alpha = 0.2) + 
  theme_bw() + 
  labs(x = "Day of the Year", y = "Counts", 
       title = "Daily Infections and Deaths from COVID-19", color = NULL) +
  scale_color_discrete(labels = c("Fitted Deaths", "Estimated New Infections")) + 
  theme(legend.position = c(0.7, 0.8),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = 'black'),
        plot.title = element_text(face = "bold", size = 15))
  



