n <- 1000
people <- 1:n
h_max <- 5

set.seed(13)

#### Putting people in households
# this does use more computation than necessary, but its such a tiny prop. of 
# the total computation here it doesn't really make any difference overall.
sizes = sample(1:h_max, 1000, replace=TRUE) # generates 1000 household sizes uniform from {1,2,3,4,5}
h <- rep(1:length(sizes), sizes) # repeats i sizes[i] times. 
h <- h[1:n] # takes only the first n of these

links <- cbind(person=people, household=h)

# system.time seems to suggest this is roughly O(n^2). About 8s for n=10000
### get.net function
beta <- runif(n,0,1)
nc=15 # get rid of this
get.net <- function(beta, h, nc=15){
  n <- length(beta)
  
  # this is a bit of a hack to make a "household matrix" from h
  # the H[i,j] is perfect sqaure only if h[i]==h[j], ie if i,j from same hshld
  H <- matrix(rep(0,n^2),n,n)
  H[sqrt(outer(h,h))%%1==0]<-1
  diag(H) <- 0
  
  b_bar <- sum(beta)/n
  M_prob <-  nc*outer(beta, beta)/(b_bar^2 *(n-1))
  v_prob <- M_prob[upper.tri(M_prob)] # getting vector of upper triangular els of M_prob
  v_cons <- rbinom(length(v_prob), 1, v_prob) # generate bernoullis 
  
  M_cons <- matrix(0,nrow=n, ncol=n) # creating new matrix to be filled with bern realisations
  M_cons[upper.tri(M_cons)] <- v_cons
  
  M_cons <- M_cons + t(M_cons)
  M <- (!H)*M_cons  #moving family connections by multiplying by negated household matrix
  conns_list <- as.list(as.data.frame(matrix(M_cons==1, n,n)))
  i_list <- lapply(conns_list, which )
  
}

alink <- get.net(beta, h) # alink is the variable we need to put into nseir



nseir <- function(beta, h, alink, alpha=c(.1,.01,.01), delta=.2, gamma=.4, nc=15, nt=100, pinf=.005){
  n <- length(beta)
  x <- rep(0,n)  # everyone starts out in susceptible state
  x[sample(1:n, pinf*n)] <- 2 # random subset of specified size is infected
  
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

# default parameters
epi <- nseir(beta, h, alink)

# adjusted parameters
epi1 <- nseir(beta, h, alink, 
              alpha = c(0.5, 0.3, 0.1), gamma = 0.5, delta = 0.1)

# plotting
plot(x = epi$t, y = epi$S, ylim = c(0, max(epi$S)), 
     xlab = "day", ylab = "N")
points(epi$E, col = 4); points(epi$I, col = 2); points(epi$R, col = 3)
legend(x = "right", y = "center", 
       legend = c("Susceptible", "Exposed",
                  "Infected", "Recovered"),
       fill = c("black", "blue",
                "red", "green"))
