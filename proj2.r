# Zachary Gavin (s2222962): Wrote get.net, wrote household allocation, created initial nseir function, debugging and commenting
# Shaehroz Khalid (s2869421): All plotting and visualisation, turned household allocation vector to one line, debugging and commenting
# Brandon Causing (s2901457): Completion and performance optimisation of the main nseir simulation function, debugging and commenting
# division of work was close to 1/3 each

################################################################################


## GENERAL DESCRIPTION

# This script implements an agent-based SEIR (Susceptible-Exposed-Infected-Recovered)
# model to simulate the spread of an epidemic in a population of n individuals.
# The model incorporates three modes of transmission: within-household, through a
# pre-defined social network, and via random mixing.
# The script defines functions to generate the population and social network,
# run the SEIR simulation, and plot the results under different scenarios.
# The final section compares four scenarios to analyse the impact of social
# structure and individual heterogeneity (beta) on the epidemic dynamics.



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

######################################################################
######### ------- NETWORK GENERATION FUNCTIONS ------- ###############
######################################################################

# These functions create the social network connections between individuals
# who are not in the same household.


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

###########################################################
######### ------- SEIR SIMULATION FUNCTION ------- ##########
###########################################################

# This function is the core of the model. It simulates the epidemic over 
# nt time steps, tracking the number of individuals in each SEIR compartment.

nseir <- function(beta, h, alink, 
                  alpha=c(0.1, 0.01, 0.01), 
                  delta=.2, 
                  gamma=.4, 
                  nc=15, nt=100, pinf=.005) {
  
  n <- length(beta)
  x <- rep(0,n) # everyone starts off susceptible 
  num_to_infect <- max(1, round(pinf*n))
  x[sample(1:n, num_to_infect)] <- 2 # pick a few people to start off infected
  
  # set up households (who's in each household)
  hh_index_list <- split(1:n, h) 
  
  # for future calculation of mixing probs
  b_bar <- sum(beta)/n
  
  # initialize counts of statuses
  t <- S <- E <- I <- R <- rep(0, nt)
  t[1] <- 1
  S[1] <- n - num_to_infect
  I[1] <- num_to_infect
  
  # daily simulation
  for(i in 2:nt) {
    
    # random numbers to compare with probabilities
    u <- runif(n) 
    
    # store who is infected
    infected_indices <- which(x == 2)
    n_infected <- length(infected_indices)
    
    # initialize lists of those exposed through infected people 
    hh_exposed <- net_exposed <- mix_exposed <- integer(0)
    
    # find who transitions from susceptible to exposed
    if (n_infected > 0) {
      
      ## -- find who gets infected from household contacts -- ##
      # count how many people are infected in each household
      hh_infected_counts <- tabulate(h[infected_indices], nbins = max(h))
      # which households have someone infected in them
      infected_hh <- which(hh_infected_counts > 0)
      # prob of infection is higher for households with multiple people infected
      if (length(infected_hh) > 0) {
        hh_probs <- 1 - (1-alpha[1])^hh_infected_counts[infected_hh]
        # who transitions to infected in the household
        hh_exposed <- unlist(
          Map(function(hh, prob) { # this defines a function and then applys it to each pair
            hh_members <- hh_index_list[[as.character(hh)]]
            sus <- hh_members[x[hh_members] == 0] 
            sus[u[sus] < prob]
          }, hh = infected_hh, prob = hh_probs), use.names = FALSE
        )
      }
      
      ## -- find who gets infected from regular network contacts -- ##
      # who is in a regular network with someone infected
      neighbors_of_infected <- unlist(alink[infected_indices], use.names = FALSE)
      if(length(neighbors_of_infected) > 0) {
        # count how many people are infected in each regular network 
        net_infected_counts <- tabulate(neighbors_of_infected, nbins = n)
        # which networks have someone infected in them
        sus_exposed_to_net <- which(net_infected_counts > 0 & x == 0)
        # prob of infection is higher for networks with multiple people infected
        if (length(sus_exposed_to_net) > 0) {
          net_probs <- 1 - (1-alpha[2])^net_infected_counts[sus_exposed_to_net]
          # who transitions to infected in network
          net_exposed <- sus_exposed_to_net[u[sus_exposed_to_net] < net_probs]
        }
      }
      
<<<<<<< HEAD
      if (n > 1) {
        total_contacts <- n_infected * nc
        contacts <- sample(n, total_contacts, replace = TRUE)
        infected_repeats <- rep(infected_indices, each = nc)
        
        log_prob_no_infection <- log1p(-alpha[3] * beta[infected_repeats] * beta[contacts] / (b_bar^2))
        
        sum_log_probs <- tapply(log_prob_no_infection, contacts, sum)
        ids_contacted <- as.integer(names(sum_log_probs))
        prob_infection <- 1 - exp(sum_log_probs)
        
        
        
        sus_contacted <- ids_contacted[x[ids_contacted] == 0]
        if(length(sus_contacted) > 0){
          u_mix <- runif(length(sus_contacted),0,1) # new unifs to avoid dep. of mixing infection and net/hh infection
          prob_idx <- match(sus_contacted, ids_contacted)
          mix_exposed <- sus_contacted[u_mix < prob_infection[prob_idx]]
        }
=======
      ## -- find who gets infected from random mixing -- ##
      # total number of contacts between all infected people
      total_contacts <- n_infected * nc
      # randomly select who those contacts are
      contacts <- sample(n, total_contacts, replace = TRUE)
      # set up for calculating mixing probs for infected people with corresponding contacts
      infected_repeats <- rep(infected_indices, each = nc)
      
      # exp^(sum(log(x_i))) faster than product(x_i)
      # prob of infection is higher for those that mixed with multiple infected people
      log_prob_no_infection <- log1p(-alpha[3] * beta[infected_repeats] * beta[contacts] / (b_bar^2))
      sum_log_probs <- tapply(log_prob_no_infection, contacts, sum)
      prob_infection <- 1 - exp(sum_log_probs)
      
      # store who came into contact with those infected
      ids_contacted <- as.integer(names(sum_log_probs))
      # subset to those that are susceptible
      sus_contacted <- ids_contacted[x[ids_contacted] == 0]
      # who transitions to infected from mixing
      if(length(sus_contacted) > 0){
        prob_idx <- match(sus_contacted, ids_contacted)
        mix_exposed <- sus_contacted[u[sus_contacted] < prob_infection[prob_idx]]
>>>>>>> a98c33b6ed19042094fc23103d9ab16f88e8ca06
      }
    }
    
    # transition from infected to recover with probability = delta
    x[x == 2 & u < delta] <- 3
    # transition from exposed to infected with probability = gamma
    x[x == 1 & u < gamma] <- 2
    
    # gather everyone who transitions from susceptible to exposed
    if (length(hh_exposed) > 0 || length(net_exposed) > 0 || length(mix_exposed) > 0) {
      newly_exposed <- unique(c(hh_exposed, net_exposed, mix_exposed))
      x[newly_exposed] <- 1
    }
    
    # get daily counts of statuses and store it in a list
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
    t[i] <- i
  }
  return(list(t = t, S = S, E = E, I = I, R = R))
}


###########################################################
######### ------- PLOTTING HELPER FUNCTION ------- ##########
###########################################################

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


###########################################################
######### ------- MODEL SCENARIO COMPARISON (Q5) ------- ####
###########################################################

par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.2, 0.7, 0)) #2x2 grid, mar allows good spacing, mgp moves axis

adjacencyList_variable_beta <- get.net(beta, h)

#scenario 1: full model with default parameters
epi_plot(beta, h, adjacencyList_variable_beta, title = "1. Full Model")


#scenario 2: random mixing only
epi_plot(beta, h, adjacencyList_variable_beta, alpha = c(0, 0, 0.04),
         title = "2. Random Mixing Only")

#scenario 3: full model with constant beta value
beta_const <- rep(mean(beta), n)
adjacencyList_constant_beta <- get.net(beta_const, h)
epi_plot(beta_const, h, adjacencyList_constant_beta, 
         title = "3. Full Model + Constant Beta Value")

#scenario 4: random mixing with constant beta
epi_plot(beta_const, h, adjacencyList_constant_beta, alpha = c(0, 0, 0.04),
         title = "4. Random Mixing + Constant Beta Value")

par(mfrow = c(1, 1))
