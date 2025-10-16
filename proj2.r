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




######################################################################
### ------- START OF MODEL SETUP AND POPULATION GENERATION ------- ###
######################################################################

# We first define the global parameters for the model and generate the base 
# population, setting the maximum size of a household and assigning each of the 
# 'n' individuals a "sociability" parameter.

n <- 10000
people <- 1:n
h_max <- 5
beta <- runif(n, 0, 1)


######################################################################
######## ---------- HOUSEHOLD GENERATION FUNCTION ---------- #########
######################################################################

# We then generate n hh sizes uniform on {1,2,3,4,5}, make a vector repeating
# the "household number" the size of the hh times, then cut off at length n.
# This assigns each of the n individuals to a household.

h <- rep(1:n, sample(1:h_max, n, replace = TRUE))[1:n] 


######################################################################
############ ------- NETWORK GENERATION FUNCTION ------- #############
######################################################################

# We proceed to make a function that creates regular social network connections 
# between individuals who are not in the same household.

# * beta = length n vector of sociability parameters
# * h = household connections
# * nc = average number of network connections someone has
# get.net initially generates a vector of Bernoullis for each person i, containing
# their connections with people (i+1):n. For each person, this is immediately 
# turned into a list of indices to save computation time, then into a list
# of connection pairs for inversion and removing hh connections.

get.net <- function(beta, h, nc=15) {
  
  n <- length(h)
  b_bar <- sum(beta)/length(beta)
  # initialising two lists to store connections
  conns_init <- vector(mode = "list", length = n) 
  conns <- vector(mode = "list", length = n)
  for (i in 1:n) {
    if (i < n) {
      b <- rbinom(n-i, 1, nc*beta[i]*beta[(i+1):n]/(b_bar^2*(n-1)) )
      # which(b==1) returns values in 1:n-i, so +i to get indices in (i+1):n
      conns_init[[i]] <- which(b == 1) + i  
    }
  }
  
  # vectorised rep, sapply returns vector. pairs is all unique net connections
  pairs <- cbind(rep(1:n, sapply(conns_init, length)), unlist(conns_init)) 
  
  if (nrow(pairs) > 0) { 
    # keep only rows of pairs where members are not in the same household
    pairs <- pairs[h[pairs[,1]] != h[pairs[,2]], , drop = FALSE] 
    
    if(nrow(pairs) > 0) {
      for (i in 1:nrow(pairs)){
        # add person 1 in pair to connections list of person 2, and vice-versa
        conns[[pairs[i,1]]] <- append(conns[[pairs[i,1]]],  pairs[i,2] )
        conns[[pairs[i,2]]] <- append(conns[[pairs[i,2]]],  pairs[i,1] )
      }
    }
  }
  
  return(conns)
}


######################################################################
############## ------- SEIR SIMULATION FUNCTION ------- ##############
######################################################################

# We now create a function that simulates the epidemic over nt time steps, 
# incorporating into the model the household and regular network contacts that 
# can be established using the previous functions, and tracking the number of 
# individuals in each SEIR compartment.

# * beta = sociability vector 
# * h = list of what household each person is in
# * alink = nested list of the net connections of each person 
# * alpha = parameters governing the prob. of exposure from each pathway 
#  -- alpha[1] concerns hh infection, alpha[2] network infections, alpha[3] gen. mixing 
# * delta =  daily prob of I->R 
# * gamma = daily prob of E->R 
# * nc = measure of how much "social interaction" happens in the population
# * nt = number of days the simulation runs over
# * pinf is prop of population initially infected

nseir <- function(beta, h, alink, 
                  alpha=c(0.1, 0.01, 0.01), 
                  delta=.2, 
                  gamma=.4, 
                  nc=15, nt=100, pinf=.005) {
  
  ### --- 1. Set-Up --- ###
  n <- length(beta)
  x <- rep(0, n) # everyone starts off susceptible 
  num_to_infect <- max(1, round(pinf*n))
  x[sample(1:n, num_to_infect)] <- 2 # pick a few people to start off infected
  
  # set up households (who's in each household)
  hh_index_list <- split(1:n, h) 
  
  # for future calculation of mixing probs
  b_bar <- sum(beta)/n
  
  # initialize counts of statuses for day 1
  t <- S <- E <- I <- R <- rep(0, nt)
  t[1] <- 1
  S[1] <- n - num_to_infect
  I[1] <- num_to_infect
  
  ### --- 2. Daily Simulation --- ###
  # each day:
  # -- find who transitions to a different status 
  # -- count how many people are in each status
  for(i in 2:nt) {
    
    # random numbers assigned to each person to compare with transition probs
    u <- runif(n) 
    
    # keep track of infected people
    infected_indices <- which(x == 2)
    susceptible_indices <- which(x==0)
    n_infected <- length(infected_indices)
    
    ## -- 2a. Infectious to Recovered -- ##
    # transition from infected to recover with probability = delta
    x[x == 2 & u < delta] <- 3
    
    ## -- 2b. Exposed to Infectious -- ##
    # transition from exposed to infected with probability = gamma
    x[x == 1 & u < gamma] <- 2
    
    ## -- 2c. Set-Up of Susceptible to Exposed -- ##
    # find who transitions from susceptible to exposed
    # initialize lists of those exposed through infected people 
    hh_exposed <- net_exposed <- mix_exposed <- integer(0)
    
    # if nobody's infected, nobody can become exposed
    if (n_infected > 0) {
      
      ## -- 2c(i). Household Exposure -- ##
      # find who gets infected from household contacts
      # ** prob of getting infected by a household member = alpha[1]
      # -----
      # count how many people are infected in each household
      hh_infected_counts <- tabulate(h[infected_indices], nbins = max(h))
      # which households have someone infected in them
      infected_hh <- which(hh_infected_counts > 0)
      # prob of infection by at least one infected household member
      if (length(infected_hh) > 0) {
        hh_probs <- 1 - (1-alpha[1])^hh_infected_counts[infected_hh]
        # who transitions to infected in the household
        hh_exposed <- unlist(
          Map(function(hh, prob) { # this applies function to each pair
            hh_members <- hh_index_list[[as.character(hh)]]
            sus <- hh_members[x[hh_members] == 0] 
            sus[u[sus] < prob]
          }, hh = infected_hh, prob = hh_probs), use.names = FALSE
        )
      }
      
      ## -- 2c(ii). Regular Network Exposure -- ##
      # find who gets infected from regular network contacts
      # ** prob of getting infected by regular network contact = alpha[2]
      # -----
      # who is in a regular network with someone infected
      neighbors_of_infected <- unlist(alink[infected_indices], use.names = FALSE)
      if(length(neighbors_of_infected) > 0) {
        # count how many people are infected in each regular network 
        net_infected_counts <- tabulate(neighbors_of_infected, nbins = n)
        # which networks have someone infected in them
        sus_exposed_to_net <- which(net_infected_counts > 0 & x == 0)
        # prob of infection by at least one infected regular network contact
        if (length(sus_exposed_to_net) > 0) {
          net_probs <- 1 - (1-alpha[2])^net_infected_counts[sus_exposed_to_net]
          # who transitions to infected in network
          net_exposed <- sus_exposed_to_net[u[sus_exposed_to_net] < net_probs]
        }
      }
      
      ## -- 2c(iii). Random Mixing Exposure -- ##
      # find who gets infected from random mixing
      # ** prob of infection by random mixing = a[3]*b[i]*b[j]/(b_bar^2*(n-1))
      # ** b[i], b[j] are "sociability" parameters for people i, j
      # -----
      u <- runif(length(susceptible_indices)) # generate new unifs so no undue corr. between rand. mix. infs & other inf. pathways
      # this is shape #infected x #susceptible
      
      trans_probs <- nc*alpha[3]*outer(beta[infected_indices],beta[susceptible_indices])/(b_bar^2 * (n-1))
      
      #log_trans_probs <- 
      
      # gives infection prob for each susceptible person
      inf_probs <- 1 - apply(1-trans_probs, MARGIN=2, FUN=prod) # sum over columns (possible infectors)
      
      mix_exposed <- susceptible_indices[u<inf_probs]
      
    }
    
    ## -- 2d. Susceptible to Exposed -- ##
    # gather everyone who transitions from susceptible to exposed
    if (length(hh_exposed) > 0 || length(net_exposed) > 0 || length(mix_exposed) > 0) {
      newly_exposed <- unique(c(hh_exposed, net_exposed, mix_exposed))
      x[newly_exposed] <- 1
    }
    
    ## -- 2e. Daily Counts -- ##
    # keep track of daily counts of statuses
    S[i] <- sum(x == 0)
    E[i] <- sum(x == 1)
    I[i] <- sum(x == 2)
    R[i] <- sum(x == 3)
    t[i] <- i
  }
  
  ### --- 3. Output --- ###
  # each component of list is an SEIR status with corresponding daily counts
  return(list(t = t, S = S, E = E, I = I, R = R))
}


######################################################################
############## ------- PLOTTING HELPER FUNCTION ------- ##############
######################################################################

# We next create a function to cleanly plot SEIR simulations. The function allows
# for users to easily change parameters of the SEIR model to visualize changes
# in behavior by taking the same inputs as the nseir() function.

epi_plot <- function(beta, h, alink, 
                     alpha = c(0.1, 0.01, 0.01),
                     delta = 0.2,
                     gamma = 0.4,
                     nc = 15, nt = 100, pinf = 0.005,
                     title = "SEIR Dynamics") {
  
  # simulate the model with the chosen parameters from the epi_plot arguments
  epi <- nseir(beta, h, alink, alpha, delta, gamma, nc, nt, pinf)
  
  # set up the plot limits, title, and axis label orientations
  plot(x = NA, y = NA, xlim = c(0, max(epi$t)), ylim = c(0, max(epi$S)), 
       xlab = "", ylab = "", main = title, las = 1)
  
  # set up dashed grid lines
  grid(nx = NULL, ny = NULL, lty = 2, col = "gray", lwd = 1)
  
  # set axis titles, spacing, and title orientation
  title(xlab = "Day", mgp = c(2.2, 0.7, 0)) # x axis
  mtext("N", side = 2, line = 2.8, las = 1) # y axis
  
  # plot the simulation with each day connected via a line
  lines(x = epi$t, y = epi$S, lty = 1, lwd = 2)
  lines(epi$E, col = 4, lty = 1, lwd = 2) 
  lines(epi$I, col = 2, lty = 1, lwd = 2)
  lines(epi$R, col = 3, lty = 1, lwd = 2)
  
  # create the legend and set aesthetics
    legend(x = "right", 
         legend = c("Susceptible", "Exposed",
                    "Infected", "Recovered"),
         col = c("black", "blue", "red", "green"),
         lty = 1, lwd = 2,
         # --- Key Adjustments for Spacing ---
         y.intersp = 0.9,    # Increase vertical spacing between lines
         seg.len = .5,      # Make the line segments longer and more visible
         x.intersp = 0.8,    # Adjust horizontal space between line and text
         cex = 0.85,       # Reduce the overall font and symbol size
         bty = "n")         # No box
}


######################################################################
############# ------- MODEL SCENARIO COMPARISON ------- ##############
######################################################################

# Finally, we plot four different scenarios for the model: (1) the full model with
# default parameters, (2) random mixing as the only method of infection
# (i.e. removing spread via household and network contacts), (3) the full model 
# but setting a constant sociability parameter across the population (set as the 
# mean of the betas), and (4) a mix of the previous two scenarios that considers
# both random mixing only and constant sociability.

# 2x2 grid, mar allows good spacing, mgp moves axis
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), mgp = c(2.2, 0.7, 0)) 

# regular network contacts with random betas (scenarios 1, 2)
alink <- get.net(beta, h)

# regular network contacts with constant beta (scenarios 3, 4)
beta_const <- rep(mean(beta), n)
alink_const_beta <- get.net(beta_const, h)

# ------------------

### --- SCENARIO 1: full model with default parameters --- ###
epi_plot(beta, h, alink, title = "1. Full Model")


### --- SCENARIO 2: random mixing only --- ###
epi_plot(beta, h, alink, alpha = c(0, 0, 0.04),
         title = "2. Random Mixing Only")

### --- SCENARIO 3: full model with constant beta value --- ###
epi_plot(beta_const, h, alink_const_beta, 
         title = "3. Full Model + Constant Beta Value")

### --- SCENARIO 4: random mixing with constant beta --- ###
epi_plot(beta_const, h, alink_const_beta, alpha = c(0, 0, 0.04),
         title = "4. Random Mixing + Constant Beta Value")

# ------------------
# reset plot settings
par(mfrow = c(1, 1))


######################################################################
################ ------- DISCUSSION OF PLOTS -------- ################
######################################################################

# When we remove household and network effects, while holding constant the
# total no. of daily contacts, we see the epidemic proceeding noticeably
# faster. This occurs both in the case of non-constant and constant beta.
# In the full model, although the total number of contacts in a day
# is the same as the mixing only model, these contacts will be
# concentrated locally. This introduces the possibility, for instance,
# for household groups to all get infected, at which point there can 
# be no disease spread in this pathway. These types of effects will
# combine to slow the spread when compared to the mixing only model. 
