## GENERAL DESCRIPTION
# 


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

get.net <- function(beta, h, nc=15) {
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


# alink <- get.net(beta,h)


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
        prob_infected <- 1 - prod(1-alpha[3]*id_mix_prob)
        
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



# plotting
epi_plot <- function(beta, h, alink, 
                     alpha = c(0.1, 0.01, 0.01),
                     delta = 0.2,
                     gamma = 0.4,
                     nc = 15, nt = 100, pinf = 0.005,
                     title = "SEIR Dynamics") {
  
  epi <- nseir(beta, h, alink, alpha, delta, gamma, nc, nt, pinf)
  
  plot(x = epi$t, y = epi$S, ylim = c(0, max(epi$S)), 
       xlab = "", ylab = "", main = title, las = 1)
  
  title(xlab = "Day", mgp = c(2.2, 0.7, 0)) #Spacing for x axis and title
  title(ylab = "N", mgp = c(2.5, 0.7, 0)) #Spacing for y axis and title
  
  points(epi$E, col = 4); points(epi$I, col = 2); points(epi$R, col = 3)

  legend(x = "right", 
         legend = c("Susceptible", "Exposed",
                    "Infected", "Recovered"),
         col = c("black", "blue",
                  "red", "green"),
         pch = 16,
         inset = 0.005,
         bty = "n",
         cex = 0.9)
}


##q5

par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.2, 0.7, 0)) #2x2 grid, mar allows good spacing, mgp moves axis

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
