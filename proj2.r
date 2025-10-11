## GENERAL DESCRIPTION


################################################################################
######### ------- START OF MODEL SETUP AND POPULATION GENERATION ------- #########
################################################################################

# This section defines the global parameters for the model and generates the
# base population, assigning each of the 'n' individuals to a household.


n <- 1000
people <- 1:n
h_max <- 5

set.seed(13)

#### Putting people in households
# this does use more computation than necessary, but its such a tiny prop. of 
# the total computation here it doesn't really make any difference overall.

#sizes = sample(1:h_max, 1000, replace=TRUE) # generates 1000 household sizes uniform from {1,2,3,4,5}
#h <- rep(1:length(sizes), sizes) # repeats i sizes[i] times. 
#h <- h[1:n] # takes only the first n of these

h <- rep(1:n, sample(1:h_max, n, replace = TRUE))[1:n] #Above code written in one line. replace 1000 with n as n can change.


links <- cbind(person=people, household=h)

# system.time seems to suggest this is roughly O(n^2). About 8s for n=10000
### get.net function
beta <- runif(n,0,1)
nc=15 # get rid of this
get.net <- function(beta, h, nc=15){
  n <- length(beta)
  
  # this is a bit of a hack to make a "household matrix" from h
  # the H[i,j] is perfect sqaure only if h[i]==h[j], ie if i,j from same hshld
  H <- outer(h,h, function(x1, x2) x1 == x2) # household connections (1/0)
  diag(H) <- FALSE # remove self
  #H <- matrix(rep(0,n^2),n,n)
  #H[sqrt(outer(h,h))%%1==0]<-1
  #diag(H) <- 0
  H <- outer(h,h, "==")
  diag(H) <- FALSE

  
  b_bar <- sum(beta)/n
  M_prob <-  nc*outer(beta, beta)/(b_bar^2 *(n-1))
  v_prob <- M_prob[upper.tri(M_prob)] # getting vector of upper triangular els of M_prob
  v_cons <- rbinom(length(v_prob), 1, v_prob) # generate bernoullis 
  
  
  M_cons <- matrix(0,nrow=n, ncol=n) # creating new matrix to be filled with bern realisations
  M_cons[upper.tri(M_cons)] <- v_cons
  
  M_cons <- M_cons + t(M_cons)
  M <- (!H)*M_cons  #moving family connections by multiplying by negated household matrix
  conns_list <- as.list(as.data.frame(matrix(M_cons==1, n, n)))
  i_list <- lapply(conns_list, which)
  
}

get.net2 <- function(beta, h, nc=15) {
  
  n <- length(h)
  b_bar <- sum(beta)/length(beta) 
  conns_init <- list()
  conns <- vector(mode="list", length=n)
  for (i in 1:n) {
    b <- rbinom(n-i, 1, nc*beta[i]*beta[i+1:n]/(b_bar^2*(n-1)) )
    conns_init[[i]] <- which(b==1) + i 
  }
  
  pairs <- cbind(rep(1:n, lapply(conns_init, length)), unlist(conns_init)) # might want this as a list
  
  pairs <- pairs[h[pairs[,1]] != h[pairs[,2]] , ] # removing family connections.
    
  for (i in 1:nrow(pairs)){ # putting into a list
    conns[[pairs[i,1]]] <- append(conns[[pairs[i,1]]],  pairs[i,2] )
    conns[[pairs[i,2]]] <- append(conns[[pairs[i,2]]],  pairs[i,1] )
  }
  
  conns
    
}


alink <- get.net(beta, h) # alink is the variable we need to put into nseir
alink2 <- get.net2(beta,h)


nseir <- function(beta, h, alink, 
                  alpha=c(0.1, 0.01, 0.01), 
                  delta=.2, 
                  gamma=.4, 
                  nc=15, nt=100, pinf=.005) {
  
  n <- length(beta)
  x <- rep(0,n)  # everyone starts out in susceptible state
  num_to_infect <- max(1, round(pinf * n))
  x[sample(1:n, num_to_infect)] <- 2 # random subset of specified size is infected
  
  H <- outer(h,h, function(x1, x2) x1 == x2) # household connections (1/0)
  diag(H) <- FALSE # remove self
  #H_prob <- alpha[1]*H # household connections (alpha/0's)
  H_conns_list <- as.list(as.data.frame(H))
  H_i_list <- lapply(H_conns_list, which) # list of indices for household conns
  
  b_bar <- sum(beta)/n
  M_prob <-  nc*outer(beta, beta)/(b_bar^2 *(n-1)) # this matrix is for general mixing daily prob of infection
  v_prob <- M_prob[upper.tri(M_prob)] # getting vector of upper triangular els of M_prob
  
  t <- S <- E <- I <- R <- rep(0, nt)
  t[1] <- 1; S[1] <- n*(1-pinf); I[1] <- n*pinf 
  #daily_counts <- list(cbind(t = 0, S = n*(1-pinf), E = 0, I = pinf*n, R = 0))
  for(i in 2:nt){
    u <- runif(n) # random variables to compare against

    ## random mixing set up
    v_mix <- rbinom(length(v_prob), 1, v_prob) # generate bernoullis 
    
    M_mix <- matrix(0, nrow = n, ncol = n) # creating new matrix to be filled with bern realisations
    M_mix[upper.tri(M_mix)] <- v_mix
    M_mix <- M_mix + t(M_mix)
    
    mix_list <- as.list(as.data.frame(matrix(M_mix == 1, n, n)))
    i_mix_list <- lapply(mix_list, which) # who interacted with each other on each day
    
    ## who was in contact with someone in state I
    h_exposed <- unique(unlist(H_i_list[which(x == 2)])) # everyone in household where there's an infection
    net_exposed <- unique(unlist(alink[which(x == 2)])) # everyone in network where there's an infection
    mix_infected <- unique(unlist(i_mix_list[which(x == 2)])) # everyone who got in contact with someone infected
    
    # if in state I (2), prob = delta -> R
    # if in state E (1), prob = gamma -> I
    x[x==2 & u<delta] <- 3
    x[x==1 & u<gamma] <- 2
    
    # if in state S (0),
    # -- prob = a_h -> I (2), if in household
    # -- prob = a_c -> I (2), if in network
    # -- prob = a_r*P(mixing) -> I (2), regardless
    
    #x[h_exposed][x[h_exposed] == 0 & u[h_exposed] < alpha[1]] <- 1 # if S in hh w/ I 
    #x[net_exposed][x[net_exposed] == 0 & u[net_exposed] < alpha[2]] <- 1 # if S in net w/ I
    #x[mix_infected][x[mix_infected] == 0 & u[mix_infected] < alpha[3]] <- 1
    
    x_h_exposed <- h_exposed[x[h_exposed] == 0 & u[h_exposed] < alpha[1]] # if S in hh w/ I 
    x_net_exposed <- net_exposed[x[net_exposed] == 0 & u[net_exposed] < alpha[2]] # if S in net w/ I
    x_mix_exposed <- mix_infected[x[mix_infected] == 0 & u[mix_infected] < alpha[3]]
    
    exposed_list <- list(x_h_exposed, x_net_exposed, x_mix_exposed)
    
    # randomize order of whether infected from household, network, or random mixing
    expose_order <- sample(1:3, 3)
    
    x[exposed_list[[expose_order[1]]]] <- 1
    x[exposed_list[[expose_order[2]]]] <- 1
    x[exposed_list[[expose_order[3]]]] <- 1
    
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
    t[i] <- i
    
  }
  return(list(t = t, S = S, E = E, I = I, R = R))
}



###### optimized version? ########
nseir <- function(beta, h, alink, 
                  alpha=c(0.1, 0.01, 0.01), 
                  delta=.2, 
                  gamma=.4, 
                  nc=15, nt=100, pinf=.005) {
  
  n <- length(beta)
  x <- rep(0,n)  # everyone starts out in susceptible state
  num_to_infect <- max(1, round(pinf*n))
  x[sample(1:n, num_to_infect)] <- 2 # random subset of specified size is infected
  
  # list of who shares a household with each individual
  hh_index_list <- split(1:n, h) 
  # remove self from each item in list
  hh_infect <- vector(mode = "list", length = n)
  for (i in 1:n) {
    hh_infect[[i]] <- setdiff(hh_index_list[[h[i]]], i)
  }
  
  b_bar <- sum(beta)/n
  # general mixing daily prob of infection
  mix_prob <-  nc*(beta %o% beta)/(b_bar^2 *(n-1)) # %o% upper triangular mult
  
  # vector of upper triangular els of mix_prob
  v_prob <- mix_prob[upper.tri(mix_prob)] 
  #rm(mix_prob) # save memory
  
  # initialize i,j pairs of upper triangular matrix
  ij <- which(upper.tri(matrix(1, n, n)), arr.ind = TRUE)
  
  t <- S <- E <- I <- R <- rep(0, nt)
  t[1] <- 1
  S[1] <- n - num_to_infect
  I[1] <- num_to_infect
  
  for(i in 2:nt) {
    
    u <- runif(n) # random variables to compare against
    
    ## random mixing set up
    v_mix <- rbinom(length(v_prob), 1, v_prob) # generate bernoullis 
    
    infected_indices <- which(x == 2)
    
    # who interacted with someone infected through random mixing
    if (length(infected_indices) > 0) {
      mix_pairs <- ij[v_mix == 1, , drop = FALSE]
      mix_infected <- unique(c(
        mix_pairs[mix_pairs[ , 1] %in% infected_indices, 2],
        mix_pairs[mix_pairs[ , 2] %in% infected_indices, 1]
      ))
    } else {
      mix_infected <- integer(0)
    }
    
    # if in state I (2), prob = delta -> R
    # if in state E (1), prob = gamma -> I
    x[x == 2 & u < delta] <- 3
    x[x == 1 & u < gamma] <- 2
    
    # if in state S (0),
    # -- prob = a_h -> I (2), if in household
    # -- prob = a_c -> I (2), if in network
    # -- prob = a_r*P(mixing) -> I (2), regardless
    newly_exposed <- integer(0)
    
    if (length(infected_indices) > 0) {
      # indices of those in household of someone infected
      hh_exposed <- unique(unlist(hh_infect[infected_indices]))
      hh_exposed <- hh_exposed[x[hh_exposed] == 0 & u[hh_exposed] < alpha[1]]
      
      # indices of those in regular network of someone infected
      net_exposed <- unique(unlist(alink[infected_indices]))
      net_exposed <- net_exposed[x[net_exposed] == 0 & u[net_exposed] < alpha[2]]
      
      # indices of those in contact with someone infected via random mixing
      mix_exposed <- mix_infected[x[mix_infected] == 0 & u[mix_infected] < alpha[3]]
      
      # all (unique) indices of those that went from S to E
      newly_exposed <- unique(c(hh_exposed, net_exposed, mix_exposed))  
    }
    
    x[newly_exposed] <- 1
    
    
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
    t[i] <- i
    
  }
  return(list(t = t, S = S, E = E, I = I, R = R))
}






# default parameters
epi <- nseir(beta, h, alink)

# adjusted parameters
epi1 <- nseir(beta, h, alink, 
              alpha = c(0.5, 0.3, 0.1), gamma = 0.5, delta = 0.1)

# plotting
epi_plot <- function(beta, h, alink, 
                     alpha = c(0.1, 0.01, 0.01),
                     delta = 0.2,
                     gamma = 0.4,
                     nc = 15, nt = 100, pinf = 0.005,
                     title = "SEIR Dynamics") {
  
  epi <- nseir(beta, h, alink, alpha, delta, gamma, nc, nt, pinf)
  
  plot(x = epi$t, y = epi$S, ylim = c(0, max(epi$S)), 
       xlab = "day", ylab = "N", main = title)
  points(epi$E, col = 4); points(epi$I, col = 2); points(epi$R, col = 3)
  
  #op <- par(cex = 0.7)
  
  legend(x = "topright", 
         legend = c("Susceptible", "Exposed",
                    "Infected", "Recovered"),
         fill = c("black", "blue",
                  "red", "green"),
         bty = "n")
  
}

epi_plot(beta, h, alink, title = "Epidemic with Default Parameters")
epi_plot(beta, h, alink,
         alpha = c(0.5, 0.3, 0.1), gamma = 0.5, delta = 0.1,
         title = "Epidemic with Adjusted Parameters")
