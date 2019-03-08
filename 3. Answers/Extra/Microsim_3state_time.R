  ##################################################################################
####    Simple individual level 3-state transition model with time varying    ####
####    probabilities                                                         ####
##################################################################################

# Developed by the Decision Analysis in R for Technologies in Health (DARTH) group
# Fernando Alarid-Escudero, PhD (1) 
# Eva A. Enns, MS, PhD (1)	
# M.G. Myriam Hunink, MD, PhD (2,3)
# Hawre J. Jalal, MD, PhD (4) 
# Eline M. Krijkamp, MSc (2)	
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

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#setwd("")

#####################################################################################
# Model structure
v.n       <- c("Healthy", "Sick", "Dead")    # state names
n.s       <- length(v.n)                     # number of states
n.i       <- 10                              # number of individuals
n.t       <- 60                              # number of cycles
p.HD      <- seq(0.001, 0.05, 
                 length.out = (n.t + 1))     # probability Healthy -> Dead (increases over time, e.g. age-specific mortality) 
p.HS      <- 0.05                            # probability to become Sick when Healthy

r.SD5 <- seq(0, 10, length.out = (n.t + 1 - 5)) # increasing rate death when sick for more than 5 cycles
p.SD5 <- 1 - exp(- r.SD5)                     # Convert rate to a probability  
p.SD  <- c(rep(0.1, 5), 0.1 + 0.2 * p.SD5 )   # probability Sick -> Dead (depends on time since onset)  , stable at 10% the first 5 years and increasing thereafter                       
p.SD  <- matrix(p.SD, ncol = n.t + 1, nrow = n.i, byrow = T) # create a matrix of all possible probabilities per time per individual.

v.M_Init  <- rep("Healthy", n.i)             # initial state for all individuals
v.Ts_Init <- rep(0, times = n.i)             # no illness onset at start of model

# Costs and utilities  
c.H  <- 400                     # cost of remaining one cycle Healthy
c.S  <- 100                     # cost of remaining one cycle Sick
c.D  <- 0                       # cost of remaining one cycle Dead
u.H  <- 0.8                     # utility when Healthy 
u.S  <- 0.5                     # utility when Sick
u.D  <- 0                       # utility when Dead
d.r  <- 0.03                    # discount rate per cycle
v.dw <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for each cycle based on discount rate d.r
set.seed(1987)                            # set the seed 


# CONSTRUCT PROBABILITY FUNCTION
# Probs function outputs transition probabilities for next cycle
Probs <- function(M_it, tSick, t) { 
  # M_it:    current health state
  # tSick:   number of cycles since illness onset
  # t:       current cycle
  
  p.it <- matrix(0, nrow = n.s, ncol = n.i)             # create matrix of state transition probabilities
  rownames(p.it) <-  v.n                                # give the state names to the rows
  pos <- cbind(M_it == "Sick", tSick)              # the time points for the sick individuals
  # update p.it with the appropriate probabilities   
  p.it[, M_it == "Healthy"] <-     c(1 - p.HD[t] - p.HS, p.HS,                p.HD[t])  # trans. probabilities when Healthy 
  p.it[, M_it == "Sick"]    <- rbind(0,                  1 - p.SD[pos],     p.SD[pos])  # trans. probabilities when Sick 
  p.it[, M_it == "Dead"]    <-     c(0,                  0,                         1)  # trans. probabilities when Dead      
  
  return(t(p.it))  # return transition probability
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
         dimnames = list(paste("ind", 1:n.i, sep = " "), # name the rows ind1, ind2, ind3, etc.
                         paste("cycle", 0:n.t, sep = " ")))  # name the columns cycle1, cycle2, cycle3, etc.

# CREATE ATTRIBUTE MATRICES
# m.TSick: track time since illness onset
m.Ts = m.M


# START SIMULATION
p = Sys.time()

  m.M[, 1]  <- v.M_Init    # initial health state for individual i
  m.Ts[, 1] <- v.Ts_Init   # initialize time since illnses onset for individual i
  
  m.C[, 1] <- Costs(m.M[, 1])  # costs accrued individual i during cycle 0
  m.E[, 1] <- Effs( m.M[, 1])   # QALYs accrued individual i during cycle 0

  # cycles 1 throuugh n.t
  for (t in 1:n.t) { # open time loop
    
    # get transition probabilities based on current health state and time since illness
    m.p <- Probs(m.M[, t], m.Ts[, t], t)

    # sample the  health state at t + 1 based on transition probabilities v.p 
    m.M[, t + 1] <- samplev(probs = m.p, m = 1) 
    # update time since illness onset for t + 1 
    m.Ts[, t + 1] = ifelse(m.M[, t + 1] == "Sick", m.Ts[, t] + 1, 0)
    
    # calculate costs and effects accrued by individual i for cycle t + 1
    m.C[, t + 1] <- Costs(m.M[, t + 1])   # costs
    m.E[, t + 1] <-  Effs(m.M[, t + 1])    # QALYs
     
  } # close the loop for the time

comp.time = Sys.time() - p



# CACLULATE TOTAL COSTS AND QALYS
v.dw <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weight for each cycle based on discount rate d.r

tc <- m.C %*% v.dw    # total discounted cost per individual
te <- m.E %*% v.dw    # total discounted QALYs per individual 

tc_avg2 <- mean(tc)    # average discounted cost 
te_avg2 <- mean(te)    # average discounted QALYs


# VISUALIZE RESULTS
source("Microsim_visualize.R")

