
#########              Sick-Sicker Markov model                 #####################

#####################################################################################
## This code forms the basis for the cohort model of the article:                  ## 
## 'Microsimulation modeling for health decision sciences using R: a tutorial'     ##
## Authors: Eline Krijkamp, Fernando Alarid-Escudero,                              ##
##          Eva Enns, Hawre Jalal, Myriam Hunink and  Petros Pechlivanoglou        ##
## citation below                                                                  ##
#####################################################################################

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
# Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: a tutorial. Med. Decis. Making. 2018;. 

#####################################################################################
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the 
# collaborating institutions and may not be used, reproduced, modified, distributed or adapted in any way without written permission.

#####################################################################################

rm(list = ls())  # delete everything that is in R's memory

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
##################################### Model input #########################################

Strategies <- c("No Treatment", "Treatment")  # strategy names 
age     <- 25                                 # age at baseline
max.age <- 55                                 # maximum age of follow up
n.t  <- max.age - age                         # time horizon, number of cycles
v.n  <- c("H", "S1", "S2", "D")               # the 4 states of the model: Healthy (H), Sick (S1), Sicker (S2), Dead (D)
n.s <- length(v.n)                            # number of health states 
d.r <- 0.03                                   # equal discount of costs and QALYs by 3%

# Transition probabilities (per cycle)
p.HD    <- 0.005           # probability to die when healthy
p.HS1   <- 0.15          	 # probability to become sick when healthy
p.S1H   <- 0.5           	 # probability to become healthy when sick
p.S1S2  <- 0.105         	 # probability to become sicker when sick
rr.S1   <- 3             	 # rate ratio of death in sick vs healthy
rr.S2   <- 10            	 # rate ratio of death in sicker vs healthy 
r.HD    <- - log(1 - p.HD) # rate of death in healthy
r.S1D   <- rr.S1 * r.HD  	 # rate of death in sick
r.S2D   <- rr.S2 * r.HD  	 # rate of death in sicker
p.S1D   <- 1 - exp(-r.S1D) # probability to die in sick
p.S2D   <- 1 - exp(-r.S2D) # probability to die in sicker

# Cost and utility inputs 
c.H     <- 2000            # cost of remaining one cycle in the healthy state
c.S1    <- 4000            # cost of remaining one cycle in the sick state
c.S2    <- 15000           # cost of remaining one cycle in the sicker state
c.Trt   <- 12000           # cost of treatment(per cycle)
c.D     <- 0               # cost of being in the death state
u.H     <- 1               # utility when healthy
u.S1    <- 0.75            # utility when sick
u.S2    <- 0.5             # utility when sicker
u.D     <- 0               # utility when dead
u.Trt   <- 0.95            # utility when being treated


############################### Markov Model  ###########################

v.dwe <- v.dwc <- 1 / ((1 + d.r) ^ (0:n.t))  # discount weight (equal discounting is assumed for costs and effects)

# create transition probability matrix for NO treatment
m.P_notrt <- matrix(0,
                  nrow = n.s, ncol = n.s,
                  dimnames = list(v.n, v.n))
# fill in the transition probability array
### From Healthy
m.P_notrt["H", "H"]  <- 1 - (p.HS1 + p.HD)
m.P_notrt["H", "S1"] <- p.HS1
m.P_notrt["H", "D"]  <- p.HD
### From Sick
m.P_notrt["S1", "H"]  <- p.S1H
m.P_notrt["S1", "S1"] <- 1 - (p.S1H + p.S1S2 + p.S1D)
m.P_notrt["S1", "S2"] <- p.S1S2
m.P_notrt["S1", "D"]  <- p.S1D
### From Sicker
m.P_notrt["S2", "S2"] <- 1 - p.S2D
m.P_notrt["S2", "D"]  <- p.S2D
### From Dead
m.P_notrt["D", "D"] <- 1

# check rows add up to 1
rowSums(m.P_notrt)

# create transition probability matrix for treatment same as NO treatment
m.P_trt <- m.P_notrt

# create the markov trace matrix M capturing the proportion of the cohort in each state at each cycle
m.M_no_trt <- m.M_trt <- matrix(NA, 
                                nrow = n.t + 1, ncol = n.s,
                                dimnames = list(paste("cycle", 0:n.t, sep = " "), v.n))

head(m.M_no_trt) # show first 6 rows of the matrix 

# The cohort starts as healthy
m.M_no_trt[1, ] <- m.M_trt[1, ] <- c(1, 0, 0, 0) # initiate the Markov trace 

for (t in 1:n.t){
  ######### using transition matrices ###########
  # calculate the proportion of the cohort in each state at time t
  m.M_no_trt[t + 1, ] <- t(m.M_no_trt[t, ]) %*% m.P_notrt
     m.M_trt[t + 1, ] <- t(m.M_trt[t, ])    %*% m.P_trt
} # close the loop

head(m.M_no_trt)  # show the first 6 lines of the matrix


####### Survival ############################################

matplot(0:n.t, m.M_no_trt, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Markov Trace")                 # create a plot of the Markov trace
legend("topright", v.n, col = 1:n.s,lty = 1:n.s, bty = "n")  # add a legend to the graph

v.os_no_trt <- 1 - m.M_no_trt[, 4]                    # calculate the overall survival (OS) probability for no treatment
v.os_trt    <- 1 -    m.M_trt[, 4]                    # calculate the overall survival (OS) probability for no treatment

plot(0:n.t, v.os_trt, type='l',
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")         # create a simple plot showing the OS probability
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

LE <- sum(v.os_trt)                     # summing probability of OS across time  (i.e. life expectancy)

####### Prevalence ############################################
# Calculate prevalence of out of Markov trace.
prev.curve <- rowSums(m.M_no_trt[, 2:3])/v.os_no_trt
plot(0:n.t, prev.curve, xlab = "Cycle", ylab = "Prevalence", col = "black", main = "Prevalence of sick and sicker", type = "l")

## ratio of sick(S1) vs sicker(S2)
ratio.SS <- m.M_no_trt[, 2] / m.M_no_trt[, 3]
plot(0:n.t, ratio.SS, xlab = "Cycle", ylab = "Ratio S2 vs S2", col = "black", main = "Ratio of sick and sicker", type = "l")

plot(0:n.t, m.M_no_trt[, "H"], type = "l", ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycles", main = "Sick-Sicker") # create a plot including the healthy line
lines(0:n.t, m.M_no_trt[, "S1"], col = "green")  # add the line for sick
lines(0:n.t, m.M_no_trt[, "S2"], col = "blue")   # add the line for sicker
lines(0:n.t, m.M_no_trt[, "D"],  col = "red")    # add the line for dead
legend("topright", v.n, col = c("black", "blue", "green", "red"), lty = rep(1, 4), bty = "n")  # add a legend to the graph


# create vectors of utility and costs for each state
v.u_trt    <- c(u.H, u.Trt, u.S2, u.D)
v.u_no_trt <- c(u.H, u.S1, u.S2, u.D)

v.c_trt    <- c(c.H, c.S1 + c.Trt, c.S2 + c.Trt, c.D)
v.c_no_trt <- c(c.H, c.S1, c.S2, c.D)

# estimate mean QALys and costs
v.E_no_trt <- m.M_no_trt %*% v.u_no_trt
v.E_trt    <- m.M_trt %*% v.u_trt

v.C_no_trt <- m.M_no_trt %*% v.c_no_trt
v.C_trt    <- m.M_trt %*% v.c_trt

### discount costs and QALYs
te_no_trt <- t(v.E_no_trt) %*% v.dwe  # 1x31 %*% 31x1 -> 1x1
te_trt    <- t(v.E_trt) %*% v.dwe

tc_no_trt <- t(v.C_no_trt) %*% v.dwc
tc_trt    <- t(v.C_trt)    %*% v.dwc

# calculate lifelong per patient discounted cost and QALYs.
DC          <- tc_trt - tc_no_trt      # calculate the difference in discounted costs between the two strategies 
names(DC)   <- "Incremental costs"
DE          <- te_trt - te_no_trt      # calculate the difference in discounted effects between the two strategies 
names(DE)   <- "QALYs gained"
ICER        <- DC / DE                 # calculate the ICER
names(ICER) <- "ICER" 
results     <- c(DC, DE, ICER)         # combine the results 

# create full incremental cost-effectiveness analysis table 
C <- round(c(tc_no_trt, tc_trt), 2)  # bind and round the total costs of the two strategies
E <- round(c(te_no_trt, te_trt), 2)  # bind and round the total effects of the two strategies

DC   <- c("", as.character(round(DC, 2)))   # round the delta of the costs (No Treatment is reference)
DE   <- c("", as.character(round(DE, 2)))   # round the delta of the effects (No Treatment is reference)
ICER <- c("", as.character(round(ICER, 2))) # round the ICER 

table_Markov_SickSicker <- cbind(Strategies, C, E, DC, DE, ICER)     # combine all data in a table
table_Markov_SickSicker <- as.data.frame(table_Markov_SickSicker)    # create a data frame 
table_Markov_SickSicker                                              # print the table    
write.table(table_Markov_SickSicker, "table_Markov_SickSicker.txt")  # export the results to a .txt file 


############ One Way Sensitivity Analysis #############################

# wrap the code 3-state Markov model in a function
input <- data.frame(
  p.HD    = 0.005,            # probability to die when healthy
  p.HS1   = 0.15,        	    # probability to become sick when healthy
  p.S1H   = 0.5,           	  # probability to become healthy when sick
  p.S1S2  = 0.105,         	  # probability to become sicker when sick
  rr.S1   = 3,             	  # rate ratio of death in sick vs healthy
  rr.S2   = 10,            	  # rate ratio of death in sicker vs healthy
  c.H     = 2000,             # cost of remaining one cycle in the healthy state
  c.S1    = 4000,             # cost of remaining one cycle in the sick state
  c.S2    = 15000,            # cost of remaining one cycle in the sicker state
  c.Trt   = 12000,            # cost of treatment(per cycle)
  c.D     = 0,                # cost of being in the death state
  u.H     = 1,                # utility when healthy
  u.S1    = 0.75,             # utility when sick
  u.S2    = 0.5,              # utility when sicker
  u.D     = 0,                # utility when dead
  u.Trt   = 0.95              # utility when treated
)

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
                 "Cost_Trt" = tc_trt, 
                 "QALY_NoTrt" = te_no_trt, 
                 "QALY_Trt" = te_trt,
                 "ICER" = (tc_trt - tc_no_trt)/(te_trt - te_no_trt))
    
    return(results)
    }
  )
}

MM.SickSicker(params = input)

p.HD_range  <- c(BaseCase = 0.005, Low = 0.002, High = 0.01)
c.Trt_range <- c(BaseCase = 12000, Low = 6000, High = 18000)
u.S2_range  <- c(BaseCase = 0.5,  Low = 0.40, High = 0.70)

paramNames <- c("p.HD", "c.Trt", "u.S2")

## List of inputs
l.tor.in <- vector("list", 3)
names(l.tor.in) <- paramNames
l.tor.in$p.HD  <- cbind(p.HD = p.HD_range, input[-input$p.HD])
l.tor.in$c.Trt <- cbind(c.Trt = c.Trt_range, input[-input$c.Trt])
l.tor.in$u.S2  <- cbind(u.S2 = u.S2_range, input[-input$u.S2])

## List of outputs
l.tor.out <- vector("list", 3)
names(l.tor.out) <- paramNames

## Run model on different parameters
# Run for costs
for (i in 1:3) { # i <- 2
  l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, MM.SickSicker))[, 2] # select the cost column
}

## Data structure: ymean	ymin	ymax
m.tor <- matrix(unlist(l.tor.out), nrow = 3, ncol = 3, byrow = TRUE,
                dimnames = list(paramNames, c("basecase", "low", "high")))

m.tor      # print the data structure 

# Plot tornado costs
source("Functions/tornado_diagram_code.R")
TornadoPlot(Parms = paramNames, Outcomes = m.tor, 
            titleName = "Tornado Plot", outcomeName = "Total costs")

# Different values for u.S2 do not have an effect on the costs
 