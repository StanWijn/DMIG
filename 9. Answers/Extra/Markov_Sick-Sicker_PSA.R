
#########              Sick-Sicker Markov model                     #################
#########             With Sensitivity Analysis                     #################

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
# Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: a tutorial. Med. Decis. Making. 2018;38(3):400-422.

#####################################################################################
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the 
# collaborating institutions and may not be used, reproduced, modified, distributed or adapted in any way without written permission.

rm(list = ls())  # Delete everything that is in R's memory

###########################################################################################
# set working directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# First run the Excersice Markov Sick Sicker model. The next step is to define the model input
source("Functions/PSA_functions.R") # load custom made functions to create: 
# Cost-effectiveness plane, 
# Cost-effectiveness acceptability curves/frontiers and
# Expected value of Perfect Information

source("Functions/functions.R") 
# install.packages("scales")
# install.packages("reshape2")
# install.packages("ggplot2")
# install.packages("ellipse")
# install.packages("plyr")
# install.packages("rpgm", dependencies = TRUE)
# install.packages("truncnorm", dependencies = TRUE)  # install package for truncated normal
library(truncnorm)  # load the package truncnorm

##################################### Model input #########################################

# Model input
v.Strategies <- c("No Treatment", "Treatment") # strategy names 
n.sim   <- 10000                               # number of simulations 
age     <- 25                                  # age at baseline
max.age <- 55                                  # maximum age of follow up
n.t  <- max.age - age                          # time horizon
v.n  <- c("H", "S1", "S2", "D")                # the 4 states of the model: Healthy, Sick, Sicker, Dead              
n.s <- length(v.n)                             # number of health states
d.r  <- d.r <- 0.03                            # equal discount of costs and QALYs by 3%

####### probabilistic analysis ########################################### 
####### Create function that generate random samples ###################################

## Function that generates random sample for PSA
gen_psa <- function(n.sim = 1000, seed = 1){
  set.seed(seed)              # set a seed to be able to reproduce the same results
  
  df.psa <- data.frame(
    # Transition probabilities (per cycle)
    p.HS1   = rbeta(n.sim, 30, 170),        # probability to become sick when healthy
    p.S1H   = rbeta(n.sim, 60, 60) ,        # probability to become healthy when sick
    p.S1S2  = rbeta(n.sim, 84, 716),        # probability to become sicker when sick
    
    p.HD    = rbeta(n.sim, 10, 1990)      ,  # probability to die when healthy
    rr.S1   = rlnorm(n.sim, log(3),  0.01),  # rate ratio of death in S1 vs healthy
    rr.S2   = rlnorm(n.sim, log(10), 0.02),  # rate ratio of death in S2 vs healthy 
    # Cost vectors with length n.sim
    c.H   = rgamma(n.sim, shape = 100, scale = 20)    , # cost of remaining one cycle in state H
    c.S1  = rgamma(n.sim, shape = 177.8, scale = 22.5), # cost of remaining one cycle in state S1
    c.S2  = rgamma(n.sim, shape = 225, scale = 66.7)  , # cost of remaining one cycle in state S2
    c.Trt = rgamma(n.sim, shape = 73.5, scale = 163.3), # cost of treatment (per cycle)
    c.D   = 0                                         , # cost of being in the death state
    # Utility vectors with length n.sim 
    u.H   = rtruncnorm(n.sim, mean =    1, sd = 0.01, b = 1), # utility when healthy
    u.S1  = rtruncnorm(n.sim, mean = 0.75, sd = 0.02, b = 1), # utility when sick
    u.S2  = rtruncnorm(n.sim, mean = 0.50, sd = 0.03, b = 1), # utility when sicker
    u.D   = 0                                               , # utility when dead
    u.Trt = rtruncnorm(n.sim, mean = 0.95, sd = 0.02, b = 1) # utility when being treated
  )
  return(df.psa)
}
# Try it
gen_psa(10) # Works!


############################### Markov Model  ###########################
MM.SickSicker <- function(params) {
  with(as.list(params), {
    # compute internal paramters as a function of external parameter
    r.HD    = - log(1 - p.HD) # rate of death in healthy
    r.S1D   = rr.S1 * r.HD 	  # rate of death in sick
    r.S2D   = rr.S2 * r.HD  	# rate of death in sicker
    p.S1D   = 1 - exp(-r.S1D) # probability to die in sick
    p.S2D   = 1 - exp(-r.S2D) # probability to die in sicker
    
    v.dwe <- v.dwc <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weight for each cycle based on discount rate d.r
    
    # create transition probability matrix for NO treatment
    m.P <- matrix(0,
                  nrow = n.s, ncol = n.s,
                  dimnames = list(v.n, v.n))
    # fill in the transition probability array
    ### From Healthy
    m.P["H", "H"]  <- 1 - (p.HS1 + p.HD)
    m.P["H", "S1"] <- p.HS1
    m.P["H", "D"]  <- p.HD
    ### From Sick
    m.P["S1", "H"]  <- p.S1H
    m.P["S1", "S1"] <- 1 - (p.S1H + p.S1S2 + p.S1D)
    m.P["S1", "S2"] <- p.S1S2
    m.P["S1", "D"]  <- p.S1D
    ### From Sicker
    m.P["S2", "S2"] <- 1 - p.S2D
    m.P["S2", "D"]  <- p.S2D
    ### From Dead
    m.P["D", "D"] <- 1
    
    m.TR <- matrix(NA, nrow = n.t + 1 , ncol = n.s, 
                   dimnames = list(0:n.t, v.n))     # create Markov trace (n.t + 1 because R doesn't understand  Cycle 0)
    
    m.TR[1, ] <- c(1, 0, 0, 0)                      # initialize Markov trace
    
    ############# PROCESS ###########################################
    
    for (t in 1:n.t){                              # throughout the number of cycles
      m.TR[t + 1, ] <- m.TR[t, ] %*% m.P           # estimate the Markov trace for cycle the next cycle (t + 1)
    }
    
    ############ OUTPUT  ###########################################
    # create vectors of utility and costs for each state
    v.u_trt    <- c(u.H, u.Trt, u.S2, u.D)
    v.u_no_trt <- c(u.H, u.S1, u.S2, u.D)
    
    v.c_trt    <- c(c.H, c.S1 + c.Trt, c.S2 + c.Trt, c.D)
    v.c_no_trt <- c(c.H, c.S1, c.S2, c.D)
    
    # estimate mean QALys and costs
    v.E_no_trt <- m.TR %*% v.u_no_trt
    v.E_trt    <- m.TR %*% v.u_trt
    
    v.C_no_trt <- m.TR %*% v.c_no_trt
    v.C_trt    <- m.TR %*% v.c_trt
    
    ### discount costs and QALYs
    te_no_trt <- t(v.E_no_trt) %*% v.dwe  # 1x31 %*% 31x1 -> 1x1
    te_trt    <- t(v.E_trt) %*% v.dwe
    
    tc_no_trt <- t(v.C_no_trt) %*% v.dwc
    tc_trt    <- t(v.C_trt)    %*% v.dwc
    
    results <- c("Cost_NoTrt" = tc_no_trt, 
                 "Cost_Trt"   = tc_trt, 
                 "QALY_NoTrt" = te_no_trt, 
                 "QALY_Trt"   = te_trt,
                 "ICER"       = (tc_trt - tc_no_trt)/(te_trt - te_no_trt))
    
    return(results)
  }
  )
}
MM.SickSicker(params = gen_psa(1))

############################### Run PSA  ###########################
## Number of simulations
n.sim <- 1000                             
## Draw random sample for PSA
df.psa <- gen_psa(n.sim = n.sim)
## Initialize matrix of outcomes
df.out <- matrix(NaN, nrow = n.sim, ncol = 5)
colnames(df.out) <- c("Cost_NoTrt", "Cost_Trt",
                      "QALY_NoTrt", "QALY_Trt",
                      "ICER")
# start the clock and run the model
## Run PSA
p <- Sys.time()   # save system time 
for(i in 1:n.sim){
  df.out[i,] <- MM.SickSicker(df.psa[i, ])
  cat('\r', paste(round(i/n.sim * 100), "% done", sep = " "))       # display the progress of the simulation
}
t.psa <- Sys.time() - p  # calculate time to run the analusis by extracting p from the current system time 

t.psa

tc_no_trt <- df.out[, 1]
te_no_trt <- df.out[, 3]
tc_trt    <- df.out[, 2]
te_trt    <- df.out[, 4]
icer      <- df.out[, 5]

hist(tc_no_trt)   # histogram of total discounted cost for control 
hist(te_no_trt)   # histogram of total discounted QALYs for control 

hist(tc_trt)      # histogram of total discounted cost for treatment
hist(te_trt)      # histogram of total discounted QALYs for treatment

hist(icer)      # histogram of treatment's ICER

v.wtp <- c(1, seq(5000, 200000, length.out = 31)) # create vector with willingness to pay values

m.c <- data.frame(df.out[, 1:2])          # create data frame of costs for both strategies
m.e <- data.frame(df.out[, 3:4])          # create data frame of effectiveness for both strategies

#install.packages("ellipse")
library("ellipse")

ScatterCE(strategies = v.Strategies, m.c = m.c, m.e = m.e)           # create cost-effectiveness plane
ceaf(v.wtp = v.wtp, strategies = v.Strategies, 
     m.e = m.e, m.c = m.c) # create cost-effectiveness acceptability curve
evpi(v.wtp = v.wtp, m.e = m.e, m.c = m.c)                            # create expected value of information plot

