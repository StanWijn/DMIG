##############################################################################
####  Simple individual level 3-state transition model with fixed         ####
####    probabilities over time                                           ####    

# Developed by the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (1)	
# M.G. Myriam Hunink, MD, PhD (2,3)
# Hawre J. Jalal, MD, PhD (4) 
# Eline M. Krijkamp, MSc (2)	6t
# Petros Pechlivanoglou, PhD (5) 

# In collaboration of: 		
# 1 University of Minnesota School of Public Health, Minneapolis, MN, USA
# 2 Erasmus MC, Rotterdam, The Netherlands
# 3 Harvard T.H. Chan School of Public Health, Boston, USA
# 4 University of Pittsburgh Graduate School of Public Health, Pittsburgh, PA, USA
# 5 The Hospital for Sick Children, Toronto and University of Toronto, Toronto ON, Canada

#####################################################################################
# Please cite our publications when using this code
# Jalal H, et al. An Overview of R in Health Decision Sciences. Med. Decis. Making. 2017; 37(3): 735-746. 
# Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: a tutorial. Med. Decis. Making. 2018; (in press). 

#####################################################################################
# Copyright 2017, THE HOSPITAL FOR Sick CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property are exclusively owned by THE HOSPITAL FOR Sick CHILDREN and the 
# collaborating institutions and may not be used, reproduced, modified, distributed or adapted in any way without written permission.

#####################################################################################
rm(list =ls())

# setwd("your working directory")  #set working directory as the folder where the course material is stored

# INPUT PARAMETERS
set.seed(1987)                            # set the seed 

# Model structure
v.n       <- c("Healthy", "Sick", "Dead")    # state names
n.s       <- length(v.n)                     # number of states
n.i       <- 10000                           # number of individuals
n.sim     <- 100                             # number of PSA simulations
v.M_Init  <- rep("Healthy",n.i)              # initial state for all individuals
n.t       <- 60                              # number of cycles

# Costs and utilities  
p.HD <- rbeta(n.sim, 20,  980)                # probability from healthy to dead
p.HS <- rbeta(n.sim, 50,  950)                # probability from healthy to sick
p.SD <- rbeta(n.sim, 100, 900)                # probability from sick to dead
c.H  <- rnorm(n.sim, 400, 50)                 # cost of being in healthy state
c.S  <- rnorm(n.sim, 100, 80)                 # cost of being in sick state
c.D  <- 0

u.H  <- rnorm(n.sim, 0.8, 0.02)               # utility of being in healthy state
u.S  <- rnorm(n.sim, 0.5, 0.02)               # utility of being in the sick state
u.D  <- 0
d.r  <- 0.03                    # discount rate per cycle
v.dw <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for each cycle based on discount rate d.r

########################################## samplev solution #######################3
# samplev : efficient implementation of the rMultinom() function of the Hmisc package 
# probs: the matrix of probabilities for each individual for each state at each time point
# m: the number of values that need to be sampled at a time per individual

samplev <- function (probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}

# CONSTRUCT PROBABILITY FUNCTION
# Probs function outputs transition probabilities for next cycle

# Application of the samplev() function requires slight modification of the function updating the transition probabilities. This modification is necessary for the matrix of probabilities to be updated for a cohort of individuals rather than a specific individual. Hence, the output of the Probs() function is now a matrix of n.sim x n.s. Since the samplev() function samples the state individuals conditional on their probability, state names are required as columnames 
set.seed(1987)                            # set the seed 


Probs <- function(M_it, k) { 
  # M_it:    current health state
  # m.n :    the matrix of possible states
  p.it <- matrix(0, nrow = n.s, ncol = n.i)             # create matrix of state transition probabilities
  rownames(p.it) <-  v.n                                # give the state names to the rows
  
  # update p.it with the appropriate probabilities   
  p.it [, M_it == "Healthy"] <- c(1 - p.HD[k] - p.HS[k], p.HS[k], p.HD[k])     # transition probabilities when Healthy 
  p.it [, M_it == "Sick"]    <- c(0,     1 - p.SD[k],       p.SD[k])     # transition probabilities when Sick 
  p.it [, M_it == "Dead"]    <- c(0,            0,          1)     # transition probabilities when Dead      
  return(t(p.it))    # return the transition probabilities
}  

# CONSTRUCT COST FUNCTION
# Costs function calculates the cost accrued by an individual this cycle
Costs <- function (M_it, k) {
  # M_it: current health state
  c.it <- c()
  c.it[M_it == "Dead"]    <- c.D        # costs at Dead state
  c.it[M_it == "Healthy"] <- c.H[k]     # costs accrued by being Healthy this cycle
  c.it[M_it == "Sick"]    <- c.S[k]     # costs accrued by being Sick this cycle
  
  return(c.it)  # return costs accrued this cycle
}

samplev <- function (probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}

# CONSTRUCT OUTCOME (EFFECTIVENESS) FUNCTION
# Effs function outputs QALYs accrued by an individual for this cycle
Effs <- function (M_it, k) {
  # M_it: current health state
  q.it <- c() 
  q.it[M_it == "Dead"]    <- u.D           # QALYs at Dead state
  q.it[M_it == "Healthy"] <- u.H[k]        # QALYs accrued by being Healthy this cycle
  q.it[M_it == "Sick"]    <- u.S[k]        # QALYs accrued by being Sick this cycle
  
  return(q.it)  # return the QALYs accrued this cycle
}

# start PSA
p <- Sys.time()
tc_PSA <- te_PSA <- matrix(NA, nrow = n.i, ncol = n.sim)


for (k in  1: n.sim)
    {# CREATE RESULT MATRICES
  # m.M: health state for each patient at each cycle
  # m.E: outcomes (e.g. QALYs) accrued by each patient at each cycle
  # m.C: costs accrued by each individual at each cycle
  m.M  <- m.E <- m.C <- 
    matrix(nrow = n.i, ncol = n.t + 1, 
           dimnames = list(paste("ind",   1:n.i, sep = " "),   # name the rows ind1, ind2, ind3, etc.
                           paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle1, cycle2, cycle3, etc.
  
  # Start Microsimulation

  # we will simulate all individuals simultaneously (no looping over individuals)
  m.M[, 1] <- v.M_Init         # initial health state for all individuals
  m.C[, 1] <- Costs(m.M[, 1], k)  # costs accrued by each individual during cycle 0
  m.E[, 1] <- Effs( m.M[, 1] , k)  # QALYs accrued by each individual during cycle 0
  
  # cycles 1 through n.t
  for (t in 1:n.t) { # open time loop
    
    # get transition probabilities based on health state at t
    m.p <- Probs(m.M[, t], k)
    
    # sample the current health state based on transition probabilities v.p 
    m.M[, t+1] <- samplev(m.p, 1)      # health states for all individuals during cycle t + 1
    m.C[, t+1] <- Costs(m.M[, t + 1], k)  # costs accrued by each individual during cycle  t + 1
    m.E[, t+1] <- Effs( m.M[, t + 1], k)  # QALYs accrued by each individual during cycle  t + 1
    
  } 
  
  
  # CACLULATE TOTAL COSTS AND QALYS
  v.dw <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weight for each cycle based on discount rate d.r
  
  tc_PSA[, k] <- m.C %*% v.dw    # total discounted cost per individual
  te_PSA[, k] <- m.E %*% v.dw    # total discounted QALYs per individual 
}

tc_avg2 <- colMeans(tc_PSA)    # average discounted cost 
te_avg2 <- colMeans(te_PSA)    # average discounted QALYs

hist(tc_avg2)

print(Sys.time() - p)

