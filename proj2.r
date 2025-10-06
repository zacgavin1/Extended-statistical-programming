n <- 100
people <- 1:n
h_max <- 5

set.seed(13)

#### Putting people in households
sizes = sample(1:h_max, 1000, replace=TRUE) # generates 1000 household sizes uniform from{1,2,3,4,5}
h <- rep(1:length(sizes), sizes) # repeats i sizes[i] times. 
h <- h[1:n] # takes only the first n of these

links <- cbind(person=people, household=h)

### get.net function
beta <- runif(n,0,1)
nc=15 # get rid of his
get.net <- function(beta, nc=15){
  n <- length(beta)
  
  b_bar <- sum(beta)/n
  M_prob <-  nc*outer(beta, beta)/(b_bar^2 *(n-1))
  v_prob <- M_prob[upper.tri(M_prob)] # getting vector of upper triangular els of M_prob
  v_cons <- rbinom(length(v_prob), 1, v_prob) # generate bernoullis 
  
  M_cons <- matrix(0,nrow=n, ncol=n) # creating new matrix to be filled with bern realisations
  M_cons[upper.tri(M_cons)] <- v_cons
  
  M_cons <- M_cons + t(M_cons)
  conns_list <- as.list(as.data.frame(M_cons))
}




connections <- get.net(beta)


