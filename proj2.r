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
h <- rep(1:n, sample(1:h_max, n, replace = TRUE))[1:n] 

links <- cbind(person=people, household=h)

beta <- runif(n,0,1)

get.net <- function(beta, h, nc=15){
  n <- length(beta)
  
  H <- outer(h,h, "==")
  diag(H) <- FALSE
  
  b_bar <- sum(beta)/n
  M_prob <-  nc*outer(beta, beta)/(b_bar^2 *(n-1))
  v_prob <- M_prob[upper.tri(M_prob)] 
  v_cons <- rbinom(length(v_prob), 1, v_prob) 
  
  M_cons <- matrix(0,nrow=n, ncol=n) 
  M_cons[upper.tri(M_cons)] <- v_cons
  
  M_cons <- M_cons + t(M_cons)
  M <- (!H)*M_cons  
  conns_list <- as.list(as.data.frame(matrix(M_cons==1, n, n)))
  i_list <- lapply(conns_list, which)
  
  return(i_list) 
}

get.net2 <- function(beta, h, nc=15) {
  n <- length(h)
  b_bar <- sum(beta)/length(beta) 
  conns_init <- vector("list", n)
  conns <- vector(mode="list", length=n)
  for (i in 1:n) {
    if (i < n) {
      b <- rbinom(n-i, 1, nc*beta[i]*beta[(i+1):n]/(b_bar^2*(n-1)) )
      conns_init[[i]] <- which(b==1) + i 
    }
  }
  
  # <-- FIX IS HERE: Changed lapply to sapply
  pairs <- cbind(rep(1:n, sapply(conns_init, length)), unlist(conns_init)) 
  
  if (nrow(pairs) > 0) {
    pairs <- pairs[h[pairs[,1]] != h[pairs[,2]], , drop = FALSE] 
    
    if(nrow(pairs) > 0) {
      for (i in 1:nrow(pairs)){
        conns[[pairs[i,1]]] <- append(conns[[pairs[i,1]]],  pairs[i,2] )
        conns[[pairs[i,2]]] <- append(conns[[pairs[i,2]]],  pairs[i,1] )
      }
    }
  }
  
  return(conns)
}


alink <- get.net(beta, h) 
alink2 <- get.net2(beta,h)


nseir <- function(beta, h, alink, 
                  alpha=c(0.1, 0.01, 0.01), 
                  delta=.2, 
                  gamma=.4, 
                  nc=15, nt=100, pinf=.005) {
  
  n <- length(beta)
  x <- rep(0,n)  
  num_to_infect <- max(1, round(pinf * n))
  x[sample(1:n, num_to_infect)] <- 2 
  
  H <- outer(h,h, "==") 
  diag(H) <- FALSE 
  H_conns_list <- as.list(as.data.frame(H))
  H_i_list <- lapply(H_conns_list, which)
  
  b_bar <- sum(beta)/n
  M_prob <-  nc*outer(beta, beta)/(b_bar^2 *(n-1)) 
  v_prob <- M_prob[upper.tri(M_prob)] 
  
  t <- S <- E <- I <- R <- rep(0, nt)
  t[1] <- 1; S[1] <- n*(1-pinf); I[1] <- n*pinf 
  
  for(i in 2:nt){
    u <- runif(n) 
    v_mix <- rbinom(length(v_prob), 1, v_prob) 
    
    M_mix <- matrix(0, nrow = n, ncol = n) 
    M_mix[upper.tri(M_mix)] <- v_mix
    M_mix <- M_mix + t(M_mix)
    
    mix_list <- as.list(as.data.frame(matrix(M_mix == 1, n, n)))
    i_mix_list <- lapply(mix_list, which) 
    
    infected_indices <- which(x == 2)
    h_exposed <- unique(unlist(H_i_list[infected_indices])) 
    net_exposed <- unique(unlist(alink[infected_indices])) 
    mix_infected <- unique(unlist(i_mix_list[infected_indices])) 
    
    x[x==2 & u<delta] <- 3
    x[x==1 & u<gamma] <- 2
    
    x_h_exposed <- h_exposed[x[h_exposed] == 0 & u[h_exposed] < alpha[1]] 
    x_net_exposed <- net_exposed[x[net_exposed] == 0 & u[net_exposed] < alpha[2]] 
    x_mix_exposed <- mix_infected[x[mix_infected] == 0 & u[mix_infected] < alpha[3]]
    
    exposed_list <- list(x_h_exposed, x_net_exposed, x_mix_exposed)
    
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
  x <- rep(0,n)  
  num_to_infect <- max(1, round(pinf*n))
  x[sample(1:n, num_to_infect)] <- 2 
  
  hh_index_list <- split(1:n, h) 
  hh_infect <- vector(mode = "list", length = n)
  for (i in 1:n) {
    hh_infect[[i]] <- setdiff(hh_index_list[[h[i]]], i)
  }
  
  b_bar <- sum(beta)/n
  mix_prob <-  nc*(beta %o% beta)/(b_bar^2 *(n-1)) 
  
  v_prob <- mix_prob[upper.tri(mix_prob)] 
  
  ij <- which(upper.tri(matrix(1, n, n)), arr.ind = TRUE)
  
  t <- S <- E <- I <- R <- rep(0, nt)
  t[1] <- 1
  S[1] <- n - num_to_infect
  I[1] <- num_to_infect
  
  for(i in 2:nt) {
    
    u <- runif(n) 
    v_mix <- rbinom(length(v_prob), 1, v_prob) 
    
    infected_indices <- which(x == 2)
    
    if (length(infected_indices) > 0) {
      mix_pairs <- ij[v_mix == 1, , drop = FALSE]
      mix_infected <- unique(c(
        mix_pairs[mix_pairs[ , 1] %in% infected_indices, 2],
        mix_pairs[mix_pairs[ , 2] %in% infected_indices, 1]
      ))
    } else {
      mix_infected <- integer(0)
    }
    
    x[x == 2 & u < delta] <- 3
    x[x == 1 & u < gamma] <- 2
    
    newly_exposed <- integer(0)
    
    if (length(infected_indices) > 0) {
      hh_exposed <- unique(unlist(hh_infect[infected_indices]))
      hh_exposed <- hh_exposed[x[hh_exposed] == 0 & u[hh_exposed] < alpha[1]]
      
      net_exposed <- unique(unlist(alink[infected_indices]))
      net_exposed <- net_exposed[x[net_exposed] == 0 & u[net_exposed] < alpha[2]]
      
      mix_exposed <- mix_infected[x[mix_infected] == 0 & u[mix_infected] < alpha[3]]
      
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


plot_original_dynamics <- function(epi_data, title = "SEIR Dynamics"){
  # plotting
  plot(x = epi_data$t, y = epi_data$S, ylim = c(0, max(epi_data$S)), 
       xlab = "day", ylab = "N", main = title, las = 1)
  points(epi_data$E, col = 4); points(epi_data$I, col = 2); points(epi_data$R, col = 3)
}

# plotting
epi_plot <- function(beta, h, alink, 
                     alpha = c(0.1, 0.01, 0.01),
                     delta = 0.2,
                     gamma = 0.4,
                     nc = 15, nt = 100, pinf = 0.005,
                     title = "SEIR Dynamics") {
  
  epi <- nseir(beta, h, alink, alpha, delta, gamma, nc, nt, pinf)
  
  plot(x = epi$t, y = epi$S, ylim = c(0, max(epi$S)), 
       xlab = "day", ylab = "N", main = title, las = 1)
  points(epi$E, col = 4); points(epi$I, col = 2); points(epi$R, col = 3)
  
  legend(x = "topright", 
         legend = c("Susceptible", "Exposed",
                    "Infected", "Recovered"),
         fill = c("black", "blue",
                  "red", "green"),
         bty = "n")
  
}

par(mfrow = c(1,1))
plot_original_dynamics(epi, title = "Epidemic with Default Parameters")
plot_original_dynamics(epi1, title = "Epidemic with Adjusted Parameters")

epi_plot(beta, h, alink, title = "Epidemic with Default Parameters")
epi_plot(beta, h, alink,
         alpha = c(0.5, 0.3, 0.1), gamma = 0.5, delta = 0.1,
         title = "Epidemic with Adjusted Parameters")

##q5

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1)) #2x2 grid, mar allows good spacing

adjacencyList_variable_beta <- get.net2(beta, h)

#scenario 1: full model with default parameters
epi_full <- nseir(beta, h, adjacencyList_variable_beta)
plot_original_dynamics(epi_full, title = "1. Full Model")


#scenario 2: random mixing only
epi_random <- nseir(beta, h, adjacencyList_variable_beta, alpha = c(0, 0, 0.04))
plot_original_dynamics(epi_random, title = "2. Random Mixing Only")

#scenario 3: full model with constant beta value
beta_const <- rep(mean(beta), n)
adjacencyList_constant_beta <- get.net2(beta_const, h)
epi_const_beta <- nseir(beta_const, h, adjacencyList_constant_beta)
plot_original_dynamics(epi_const_beta, title = "3. Full Model + Constant Beta Value")

#scenario 4: random mixing with constant beta
epi_random_const_beta <- nseir(beta_const, h, adjacencyList_constant_beta, alpha = c(0, 0, 0.04))
plot_original_dynamics(epi_random_const_beta, title = "4. Random Mixing + Constant Beta Value")

par(mfrow = c(1, 1))