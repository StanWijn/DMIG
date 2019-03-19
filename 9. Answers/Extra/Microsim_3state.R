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
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))  #set working directory as the folder where the course material is stored

# INPUT PARAMETERS
set.seed(1987)                            # set the seed 

# Model structure
v.n       <- c("Healthy", "Sick", "Dead")    # state names
n.s       <- length(v.n)                     # number of states
n.i       <- 10000                            # number of individuals
v.M_Init  <- rep("Healthy",n.i)              # initial state for all individuals
n.t       <- 60                              # number of cycles
p.HD      <- 0.02                            # probability to die when Healthy
p.HS      <- 0.05                            # probability to become Sick when Healthy
p.SD      <- 0.1                             # probability to die when Sick

# Costs and utilities  
c.H  <- 400                     # cost of remaining one cycle Healthy
c.S  <- 100                     # cost of remaining one cycle Sick
c.D  <- 0                       # cost of remaining one cycle Dead
u.H  <- 0.8                     # utility when Healthy 
u.S  <- 0.5                     # utility when Sick
u.D  <- 0                       # utility when Dead
d.r  <- 0.03                    # discount rate per cycle
v.dw <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for each cycle based on discount rate d.r

# CONSTRUCT PROBABILITY FUNCTION
# Probs function outputs transition probabilities for next cycle
Probs <- function(M_it) { 
  # M_it:    current health state
  
  p.it <- vector("numeric", length(v.n))     # intiatilize vector of state transition probabilities
  
  # update p.it with the appropriate probabilities   
  p.it[M_it == "Healthy"] <- c(1 - p.HD - p.HS,     p.HS, p.HD)     # transition probabilities when Healthy 
  p.it[M_it == "Sick"]    <- c(              0, 1 - p.SD, p.SD)     # transition probabilities when Sick 
  p.it[M_it == "Dead"]    <- c(              0,        0,    1)     # transition probabilities when Dead      
  
  return(p.it)  # return transition probability
}  


# CONSTRUCT COST FUNCTION
# Costs function calculates the cost accrued by an individual this cycle
Costs <- function (M_it) {
  # M_it: current health state
  c.it <- c()
  c.it[M_it == "Dead"]    <- c.D     # costs at Dead state
  c.it[M_it == "Healthy"] <- c.H     # costs accrued by being Healthy this cycle
  c.it[M_it == "Sick"]    <- c.S     # costs accrued by being Sick this cycle
  
  return(c.it)  # return costs accrued this cycle
}


# CONSTRUCT OUTCOME (EFFECTIVENESS) FUNCTION
# Effs function outputs QALYs accrued by an individual for this cycle
Effs <- function (M_it) {
  # M_it: current health state
  q.it <- c() 
  q.it[M_it == "Dead"]    <- u.D        # QALYs at Dead state
  q.it[M_it == "Healthy"] <- u.H        # QALYs accrued by being Healthy this cycle
  q.it[M_it == "Sick"]    <- u.S        # QALYs accrued by being Sick this cycle
  
  return(q.it)  # return the QALYs accrued this cycle
}


################################################################################


# CREATE RESULT MATRICES
# m.M: health state for each patient at each cycle
# m.E: outcomes (e.g. QALYs) accrued by each patient at each cycle
# m.C: costs accrued by each individual at each cycle
m.M  <- m.E <- m.C <- 
  matrix(nrow = n.i, ncol = n.t + 1, 
         dimnames = list(paste(  "ind", 1:n.i, sep = " "),     # name the rows ind1, ind2, ind3, etc.
                         paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle0, cycle1, cycle2, cycle3, etc.


# START SIMULATION
p <- Sys.time()

# individuals 1 thru n.i
for (i in 1:n.i) { # open individuals loop
  
  m.M[i, 1] <- v.M_Init[i]       # initial health state for individual i
  m.C[i, 1] <- Costs(m.M[i, 1])  # costs accrued individual i during cycle 0
  m.E[i, 1] <- Effs( m.M[i, 1])  # QALYs accrued individual i during cycle 0
  
  # cycles 1 thru n.t
  for (t in 1:n.t) { # open time loop
    
    # get transition probabilities based on health state at t
    v.p <- Probs(m.M[i, t])
    
    # sample the next health state based on transition probabilities v.p 
    m.M[i, t + 1] <- sample(x = v.n, prob = v.p, size = 1)
    
    # calculate costs and effects accrued by individual i for cycle t + 1
    m.C[i, t + 1] <- Costs(m.M[i, t + 1])    # costs
    m.E[i, t + 1] <- Effs( m.M[i, t + 1])    # QALYs
    
  } # close time loop
} # close individuals loop
comp.time = Sys.time() - p


# CACLULATE TOTAL COSTS AND QALYS

tc <- m.C %*% v.dw    # total discounted cost per individual
te <- m.E %*% v.dw    # total discounted QALYs per individual 


tc_avg <- mean(tc)    # average discounted cost 
te_avg <- mean(te)    # average discounted QALYs

# VISUALIZE RESULTS
source("Microsim_visualize.R")



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


Probs <- function(M_it) { 
  # M_it:    current health state
  # m.n :    the matrix of possible states
  p.it <- matrix(0, nrow = n.s, ncol = n.i)             # create matrix of state transition probabilities
  rownames(p.it) <-  v.n                                # give the state names to the rows
  
  # update p.it with the appropriate probabilities   
  p.it [, M_it == "Healthy"] <- c(1 - p.HD - p.HS, p.HS, p.HD)     # transition probabilities when Healthy 
  p.it [, M_it == "Sick"]    <- c(0,     1 - p.SD,       p.SD)     # transition probabilities when Sick 
  p.it [, M_it == "Dead"]    <- c(0,            0,          1)     # transition probabilities when Dead      
  return(t(p.it))    # return the transition probabilities
}  

# CONSTRUCT COST FUNCTION
# Costs function calculates the cost accrued by an individual this cycle
Costs <- function (M_it) {
  # M_it: current health state
  c.it <- c()
  c.it[M_it == "Dead"]    <- c.D     # costs at Dead state
  c.it[M_it == "Healthy"] <- c.H     # costs accrued by being Healthy this cycle
  c.it[M_it == "Sick"]    <- c.S     # costs accrued by being Sick this cycle
  
  return(c.it)  # return costs accrued this cycle
}


# CONSTRUCT OUTCOME (EFFECTIVENESS) FUNCTION
# Effs function outputs QALYs accrued by an individual for this cycle
Effs <- function (M_it) {
  # M_it: current health state
  q.it <- c() 
  q.it[M_it == "Dead"]    <- u.D        # QALYs at Dead state
  q.it[M_it == "Healthy"] <- u.H        # QALYs accrued by being Healthy this cycle
  q.it[M_it == "Sick"]    <- u.S        # QALYs accrued by being Sick this cycle
  
  return(q.it)  # return the QALYs accrued this cycle
}



# CREATE RESULT MATRICES
# m.M: health state for each patient at each cycle
# m.E: outcomes (e.g. QALYs) accrued by each patient at each cycle
# m.C: costs accrued by each individual at each cycle
m.M  <- m.E <- m.C <- 
  matrix(nrow = n.i, ncol = n.t + 1, 
         dimnames = list(paste("ind",   1:n.i, sep = " "),   # name the rows ind1, ind2, ind3, etc.
                         paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle1, cycle2, cycle3, etc.


# START SIMULATION
p = Sys.time()

# we will simulate all individuals simultaneously (no looping over individuals)
m.M[, 1] <- v.M_Init         # initial health state for all individuals
m.C[, 1] <- Costs(m.M[, 1])  # costs accrued by each individual during cycle 0
m.E[, 1] <- Effs( m.M[, 1])  # QALYs accrued by each individual during cycle 0

# cycles 1 through n.t
for (t in 1:n.t) { # open time loop
  
  # get transition probabilities based on health state at t
  m.p <- Probs(m.M[, t])
  
  # sample the current health state based on transition probabilities v.p 
  m.M[, t+1] <- samplev(m.p, 1)      # health states for all individuals during cycle t + 1
  m.C[, t+1] <- Costs(m.M[, t + 1])  # costs accrued by each individual during cycle  t + 1
  m.E[, t+1] <- Effs( m.M[, t + 1])  # QALYs accrued by each individual during cycle  t + 1
  
} # close time loop
comp.time = Sys.time() - p


# CACLULATE TOTAL COSTS AND QALYS
v.dw <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weight for each cycle based on discount rate d.r

tc <- m.C %*% v.dw    # total discounted cost per individual
te <- m.E %*% v.dw    # total discounted QALYs per individual 

tc_avg2 <- mean(tc)    # average discounted cost 
te_avg2 <- mean(te)    # average discounted QALYs


# VISUALIZE RESULTS
source("Microsim_visualize.R")
