
rm(list = ls())  # Delete everything that is in R's memory

#install.packages(c("reshape2"), dependencies = TRUE)

source("Functions/PSA_functions.R") # load custom made functions for PSA

######################### INPUT PARAMETERS ####################################
v.n   <- c("Healthy", "Sick", "Dead")     # state names
n.s   <- length(v.n)                      # number of states
n.t   <- 60                               # number of cycles
d.r   <- 0.03                             # discount rate 
v.dw  <- 1 / (1 + d.r) ^ (0:n.t)          # calculate discount weight for each cycle based on discount rate d.r

####### probabilistic analysis ########################################### 
####### Create function that generate random samples ###################################

## Function that generates random sample for PSA
gen_psa <- function(n.sim = 1000, seed = 1){
  set.seed(seed)              # set a seed to be able to reproduce the same results
  
  df.psa <- data.frame(
    p.HD = rbeta(n.sim, 20,  980),                # probability from healthy to dead
    p.HS = rbeta(n.sim, 50,  950),                # probability from healthy to sick
    p.SD = rbeta(n.sim, 100, 900),                # probability from sick to dead
    c.H  = rnorm(n.sim, 400, 50) ,                # cost of being in healthy state
    c.S  = rnorm(n.sim, 100, 80) ,                # cost of being in sick state
    c.D  = 0,
    
    u.H  = rnorm(n.sim, 0.8, 0.02),               # utility of being in healthy state
    u.S  = rnorm(n.sim, 0.5, 0.02),               # utility of being in the sick state
    u.D  = 0                                     # utility of being dead 
  )
  return(df.psa)
}
# Try it
gen_psa(n.sim = 10) # Works!

# wrap the code 3-state Markov model in a function
MM.3state <- function(params) {
  with(as.list(params),
       {
         v.dw <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weight for each cycle based on discount rate d.r
         
         ####### INITIALIZATION ##########################################
         # create the cohort trace
         m.TR <- matrix(NA, nrow = n.t + 1 , ncol = n.s, 
                        dimnames = list(0:n.t, v.n))     # create Markov trace (n.t + 1 because R doesn't understand  Cycle 0)
         
         m.TR[1, ] <- c(1, 0, 0)                     # initialize Markov trace
         
         # create the transition probability matrix
         m.P  <- matrix(0,
                        nrow = n.s, ncol = n.s,
                        dimnames = list(v.n, v.n), byrow = TRUE) # name the columns and rows of the transition probability matrix
         
         # fill in the transition probability matrix
         ### From Healthy
         m.P["Healthy", "Healthy"] <- 1 - p.HD - p.HS
         m.P["Healthy", "Sick"]    <- p.HS
         m.P["Healthy", "Dead"]    <- p.HD
         
         ### From Sick
         m.P["Sick", "Sick"] <- 1 - p.SD
         m.P["Sick", "Dead"] <- p.SD
         
         ### From Dead
         m.P["Dead", "Dead"] <- 1
         
         ####### PROCESS ###########################################
         
         for (t in 1:n.t){                              # throughout the number of cycles
           m.TR[t + 1, ] <- m.TR[t, ] %*% m.P           # estimate the Markov trace for cycle the next cycle (t + 1)
         }
         
         ####### OUTPUT  ###########################################
         
         # mean cost and QALYs per cycle
         v.c <- m.TR %*% c(c.H, c.S, c.D)  # calculate expected costs by multiplying m.TR with the cost vector for the different health states   
         v.u <- m.TR %*% c(u.H, u.S, u.D)  # calculate expected QALYs by multiplying m.TR with the utilities for the different health states   
         
         # discounted QALYs and costs
         
         TC <- t(v.c) %*% v.dw    # Discount costs by multiplying the cost vector with discount weights (v.dw) 
         TE <- t(v.u) %*% v.dw    # Discount QALYS by multiplying the QALYs vector with discount weights (v.dw)
         
         results <- c(Cost = TC, QALY = TE)
         return(results)
       }
  )
}
MM.3state(gen_psa(1))

############################### Run PSA  ###########################
## Number of simulations
n.sim <- 10000                             
## Draw random sample for PSA
df.psa <- gen_psa(n.sim = n.sim)
## Initialize matrix of outcomes
df.out <- matrix(NaN, nrow = n.sim, ncol = 2)
colnames(df.out) <- c("TE", "TC")

## Run PSA
p <- Sys.time()   # save system time 
for(i in 1:n.sim){
  df.out[i,] <- MM.3state(df.psa[i, ])
  cat('\r', paste(round(i/n.sim * 100), "% done", sep = " "))       # display the progress of the simulation
}
t.psa <- Sys.time() - p  # calculate time to run the analusis by extracting p from the current system time 

t.psa
source("Functions/PSA_functions.R")
ScatterCE(strategies = "treatment", m.e = data.frame(df.out[, "TE"]), m.c = data.frame(df.out[, "TC"])) # cost-effectiveness plane


