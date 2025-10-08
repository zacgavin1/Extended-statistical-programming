n <- 10000
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


nseir <- function(beta, h, alink, alpha=c(.1,.01,.01), delta=.2, gamma=.4, nc=15, nt=100, pinf=.005){
  n <- length(beta)
  x <- rep(0,n)  # everyone starts out in susceptible state
  x[sample(1:n, pinf*n)] <- 2 # random subset of specified size is infected
  H <- outer(h,h)
  
  b_bar <- sum(beta)/n
  M_prob <-  nc*outer(beta, beta)/(b_bar^2 *(n-1)) # this matrix is for general mixing daily prob of infection
  
  for(i in 1:nt){
    u <- runif(n) # random variables to compare against
    x[x==2&u<delta] <- 3
    x[x==1&u<gamma] <- 2
    
    
    
  }
}


alink <- get.net(beta, h) # alink is the variable we need to put into nseir


