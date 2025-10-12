## GENERAL DESCRIPTION


################################################################################
######### ------- START OF MODEL SETUP AND POPULATION GENERATION ------- #########
################################################################################

# This section defines the global parameters for the model and generates the
# base population, assigning each of the 'n' individuals to a household.


n <- 10000
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



#### higher weights for multiple infections in hh, network, random mix ######
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
    
      # how many people are infected within each household
      hh_exposed <- c()
      hh_infected_counts <- tabulate(h[infected_indices], nbins = max(h))
      # who got infected (household)
      for (hh in which(hh_infected_counts > 0)) {
        hh_members <- hh_infect[[hh]]
        sus <- hh_members[x[hh_members] == 0]
        if (length(sus) == 0) next
        prob_hh_infected <- 1 - (1-alpha[1])^hh_infected_counts[hh]
        hh_exposed_indices <- sus[u[sus] < prob_hh_infected]
        hh_exposed <- c(hh_exposed, hh_exposed_indices) 
      }
      
      # how many people are infected within each network
      net_exposed <- c()
      net_infected_counts <- mapply(function(lst, index) {
                                      ifelse(index %in% infected_indices, 
                                             sum(!(lst %in% which(x == 0))) + 1, 0)
                                    }, 
                                    alink2, seq_along(alink2))
      # who got infected (network)
      for (net in which(net_infected_counts > 0)) {
        net_members <- alink[[net]]
        sus <- net_members[x[net_members] == 0]
        if (length(sus) == 0) next
        prob_net_infected <- 1 - (1-alpha[2])^net_infected_counts[net]
        net_exposed_indices <- sus[u[sus] < prob_net_infected]
        net_exposed <- c(net_exposed, net_exposed_indices)
      }
      
      # set-up of random mixing interactions with infected people
      mix_pairs <- ij[v_mix == 1, , drop = FALSE]
      mix_infected <- c(
        mix_pairs[mix_pairs[ , 1] %in% infected_indices, 2],
        mix_pairs[mix_pairs[ , 2] %in% infected_indices, 1]
      )
      
      # how many infected people did each person randomly interact with
      mix_exposed <- c()
      mix_infected_counts <- tabulate(mix_infected, nbins = n)
      # who got infected (mixing)
      for (mix in which(mix_infected_counts > 0)) {
        sus <- mix[x[mix] == 0]
        if (length(sus) == 0) next
        # need to think about this; might not be exactly what the probability is
        prob_mix_infected <- 1 - (1-alpha[3])^mix_infected_counts[mix]
        mix_exposed_indices <- sus[u[sus] < prob_mix_infected]
        mix_exposed <- c(mix_exposed, mix_exposed_indices)
      }
    }
    
    x[x == 2 & u < delta] <- 3
    x[x == 1 & u < gamma] <- 2
    
    newly_exposed <- integer(0)
    
    if (length(infected_indices) > 0) {
      
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



###### optimized v2.0 ##########
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
  
  # remove to save memory
  rm(mix_prob)
  
  ij <- which(upper.tri(matrix(1, n, n)), arr.ind = TRUE)
  
  t <- S <- E <- I <- R <- rep(0, nt)
  t[1] <- 1
  S[1] <- n - num_to_infect
  I[1] <- num_to_infect
  
  for(i in 2:nt) {
    
    u <- runif(n) 
    v_mix <- runif(length(v_prob)) < v_prob # maybe faster?
    
    infected_indices <- which(x == 2)
    
    if (length(infected_indices) > 0) {
      
      #maybe faster?
      # how many people are infected within each household
      hh_infected_counts <- tabulate(h[infected_indices], nbins = max(h))
      # who got infected (household)
      infected_hh <- which(hh_infected_counts > 0)
      hh_probs <- 1 - (1-alpha[1])^hh_infected_counts[infected_hh]
      hh_exposed <- unlist(
        Map(function(hh, prob) {
          hh_members <- hh_infect[[hh]]
          sus <- hh_members[x[hh_members] == 0]
          sus[u[sus] < prob]
        }, hh = infected_hh, prob = hh_probs)
      )
      
      
      # maybe faster?
      # how many people are infected within each network
      net_infected_counts <- mapply(function(lst, index) {
        ifelse(index %in% infected_indices, 
               sum(!(lst %in% which(x == 0))) + 1, 0)
      }, alink, seq_along(alink))
      # who got infected (network)
      infected_net <- which(net_infected_counts > 0)
      net_probs <- 1 - (1-alpha[2])^net_infected_counts[infected_net]
      net_exposed <- unlist(
        Map(function(net, prob) {
          net_members <- alink[[net]]
          sus <- net_members[x[net_members] == 0]
          sus[u[sus] < prob]
        }, net = infected_net, prob = net_probs)
      )
      
      
      # maybe faster?
      # set-up of random mixing interactions with infected people
      mix_pairs <- ij[v_mix, , drop = FALSE]
      
      is_infected <- logical(n)
      is_infected[infected_indices] <- TRUE
      
      mix_infected <- c(
        mix_pairs[is_infected[mix_pairs[, 1]], 2],
        mix_pairs[is_infected[mix_pairs[, 2]], 1]
      )
      
      
      # maybe faster?
      # how many infected people did each person randomly interact with
      mix_infected_counts <- tabulate(mix_infected, nbins = n)
      # who got infected (mixing)
      infected_mix <- which(mix_infected_counts > 0)
      # need to think about this; might not be exactly what the probability is
      mix_probs <- 1 - (1-alpha[3])^mix_infected_counts[infected_mix]
      mix_exposed <- unlist(
        Map(function(mix, prob) {
          sus <- mix[x[mix] == 0]
          sus[u[sus] < prob]
        }, mix = infected_mix, prob = mix_probs)
      )
    }
    
    
    x[x == 2 & u < delta] <- 3
    x[x == 1 & u < gamma] <- 2
    
    newly_exposed <- integer(0)
    
    if (length(infected_indices) > 0) {
      
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




######## EVEN FASTER ##############
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
  
  t <- S <- E <- I <- R <- rep(0, nt)
  t[1] <- 1
  S[1] <- n - num_to_infect
  I[1] <- num_to_infect
  
  for(i in 2:nt) {
    
    u <- runif(n) 
    
    infected_indices <- which(x == 2)
    n_infected <- length(infected_indices)
    
    if (n_infected > 0) {
      
      #maybe faster?
      # how many people are infected within each household
      hh_infected_counts <- tabulate(h[infected_indices], nbins = max(h))
      # who got infected (household)
      infected_hh <- which(hh_infected_counts > 0)
      hh_probs <- 1 - (1-alpha[1])^hh_infected_counts[infected_hh]
      hh_exposed <- unlist(
        Map(function(hh, prob) {
          hh_members <- hh_infect[[hh]]
          sus <- hh_members[x[hh_members] == 0]
          sus[u[sus] < prob]
        }, hh = infected_hh, prob = hh_probs)
      )
      
      # maybe faster?
      # how many people are infected within each network
      net_infected_counts <- mapply(function(lst, index) {
        ifelse(index %in% infected_indices, 
               sum(!(lst %in% which(x == 0))) + 1, 0)
      }, alink, seq_along(alink))
      # who got infected (network)
      infected_net <- which(net_infected_counts > 0)
      net_probs <- 1 - (1-alpha[2])^net_infected_counts[infected_net]
      net_exposed <- unlist(
        Map(function(net, prob) {
          net_members <- alink[[net]]
          sus <- net_members[x[net_members] == 0]
          sus[u[sus] < prob]
        }, net = infected_net, prob = net_probs)
      )
      
      # maybe faster? also more correct probabilities
      # who got infected (mixing)
      v_probs <- c()
      v_mix_id <- c()
      mix_exposed <- c()
      # randomly connect infected people to nc people
      for (infected in infected_indices) {
        infect_interact_rand <- unique(sample((1:n)[-infected], nc, 
                                                replace = TRUE))
        mix_probs <- beta[infected]*beta[infect_interact_rand]/(b_bar^2)
        # store those nc people
        v_mix_id <- c(v_mix_id, infect_interact_rand)
        # store their corresponding interaction probs with infecteds
        v_probs <- c(v_probs, mix_probs)
      }
      # for each of those nc people, get their corresponding probs
      for (id in unique(v_mix_id)) {
        id_mix_prob <- v_probs[v_mix_id == id]
        # calculate prob of infection from at least one of their infected contacts
        prob_not_infected <- 1
        for (prob in id_mix_prob) {
          prob_not_infected <- prob_not_infected*(1-alpha[3]*prob)
        }
        prob_infected <- 1 - prob_not_infected
        
        exposed <- x[id] == 0 & u[id] < prob_infected
        mix_exposed <- c(mix_exposed, id[exposed])
      }
    
    }
    
    
    x[x == 2 & u < delta] <- 3
    x[x == 1 & u < gamma] <- 2
    
    newly_exposed <- integer(0)
    
    if (length(infected_indices) > 0) {
      
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
epi <- nseir(beta, h, alink2)

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
       xlab = "day", ylab = "N", main = title)
  points(epi$E, col = 4); points(epi$I, col = 2); points(epi$R, col = 3)
  
  legend(x = "topright", 
         legend = c("Susceptible", "Exposed",
                    "Infected", "Recovered"),
         fill = c("black", "blue",
                  "red", "green"),
         bty = "n",
         cex = 0.5)
  
}

par(mfrow = c(1,1))
plot_original_dynamics(epi, title = "Epidemic with Default Parameters")
plot_original_dynamics(epi1, title = "Epidemic with Adjusted Parameters")

epi_plot(beta, h, alink2, title = "Epidemic with Default Parameters")
epi_plot(beta, h, alink2,
         alpha = c(0.5, 0.3, 0.1), gamma = 0.5, delta = 0.1,
         title = "Epidemic with Adjusted Parameters")

##q5

par(mfrow = c(2, 2), mar = c(4, 4, 2, 1)) #2x2 grid, mar allows good spacing

adjacencyList_variable_beta <- get.net2(beta, h)

#scenario 1: full model with default parameters
epi_plot(beta, h, adjacencyList_variable_beta, title = "1. Full Model")


#scenario 2: random mixing only
epi_plot(beta, h, adjacencyList_variable_beta, alpha = c(0, 0, 0.04),
         title = "2. Random Mixing Only")

#scenario 3: full model with constant beta value
beta_const <- rep(mean(beta), n)
adjacencyList_constant_beta <- get.net2(beta_const, h)
epi_plot(beta_const, h, adjacencyList_constant_beta, 
         title = "3. Full Model + Constant Beta Value")

#scenario 4: random mixing with constant beta
epi_plot(beta_const, h, adjacencyList_constant_beta, alpha = c(0, 0, 0.04),
         title = "4. Random Mixing + Constant Beta Value")

par(mfrow = c(1, 1))
