

rm(list = ls())  # Delete everything that is in R's memory

####### INPUT PARAMETERS ####################################

v.n  <- c("Healthy", "Sick", "Dead")    # state names
n.s  <- length(v.n)                     # number of states
n.t  <- 60                              # number of cycles

p.HD <- 0.02                    # probability to die when healthy
p.HS <- 0.05                    # probability to become sick when healthy
p.SD <- 0.1                     # probability to die when sick

# Costs and utilities  
c.H  <- 400                     # cost of remaining one cycle healthy
c.S  <- 100                     # cost of remaining one cycle sick
c.D  <- 0                       # cost of remaining one cycle dead
u.H  <- 0.8                     # utility when healthy 
u.S  <- 0.5                     # utility when sick
u.D  <- 0                       # utility when dead
d.r  <- 0.03                    # discount rate per cycle
v.dw <- 1 / (1 + d.r) ^ (0:n.t) # calculate discount weights for each cycle based on discount rate d.r


####### INITIALIZATION ##########################################
# create the cohort trace
m.TR <- matrix(NA, nrow = n.t + 1 , ncol = n.s, 
               dimnames = list(0:n.t, v.n))     # create Markov trace (n.t + 1 because R doesn't understand  Cycle 0)

m.TR[1, ] <- c(1, 0, 0)                     # initialize Markov trace

# create the transition probability matrix
m.P  <- matrix(0,
               nrow = n.s, ncol = n.s,
               dimnames = list(v.n, v.n)) # name the columns and rows of the transition probability matrix

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

####### PROCESS  ###########################################
# Run Markov model
for (t in 1:n.t){                           # throughout the number of cycles
  m.TR[t + 1, ] <- m.TR[t, ] %*% m.P        # estimate the Markov trace for cycle the next cycle (t + 1)
}

####### Survival ############################################
matplot(m.TR, type = 'l', 
        ylab = "Probability of state occupancy",
        xlab = "Cycle",
        main = "Markov Trace")              # create a plot of the data
legend("right", v.n, col = c("black", "red", "green"), lty = 1:3, bty = "n")  # add a legend to the graph

v.os <- 1 - m.TR[, "Dead"]                  # calculate the overall survival (OS) probability
v.os <- rowSums(m.TR[, 1:2])                # alternative way of calculating the OS probability   

plot(v.os, type = 'l', 
     ylim = c(0, 1),
     ylab = "Survival probability",
     xlab = "Cycle",
     main = "Overall Survival")             # create a simple plot showing the OS
grid(nx = n.t, ny = 10, col = "lightgray", lty = "dotted", lwd = par("lwd"), equilogs = TRUE) # add grid 

LE <- sum(v.os)                             # summing probablity of OS across time  (i.e. life expectancy)

v.prev <- m.TR[, "Sick"]/v.os
plot(v.prev)
####### OUTPUT  ###########################################

# mean cost and QALYs per cycle
v.c <- m.TR %*% c(c.H, c.S, c.D)  # calculate expected costs by multiplying m.TR with the cost vector for the different health states   
v.u <- m.TR %*% c(u.H, u.S, u.D)  # calculate expected QALYs by multiplying m.TR with the utilities for the different health states   

# discounted QALYs and costs

TC <-  t(v.c) %*% v.dw   # Discount costs  by multiplying the cost vector with discount weights (v.dw) 
TE <-  t(v.u) %*% v.dw   # Discount QALYS  by multiplying the QALYs vector with discount weights (v.dw)

results <- data.frame( "Total Cost" = TC, "Life Expectancy" = LE, "Total QALYs" = TE )
results



####### One Way Sensitivity Analysis #############################

# Create data frame "input" with all input probabilities, cost and utilities for the 3-state model.
input <- data.frame(
  p.HD = 0.02,                    # probability to die when sick
  p.HS = 0.05,                    # probability to become healthy when sick
  p.SD = 0.1,                     # probability to die when sick
  c.H  = 400,                     # cost of remaining one cycle healthy
  c.S  = 100,                     # cost of remaining one cycle sick
  c.D  = 0,                       # cost of remaining one cycle dead
  u.H  = 0.8,                     # utility when healthy 
  u.S  = 0.5,                     # utility when sick
  u.D  = 0                        # utility when dead
)

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


#######  one-way sensitivity analysis #################

# create a range of low,basecase and high for the one-way sensitivity parameters
p.HD_range <- c(BaseCase = 0.02, Low = 0.01, High = 0.1)  
p.HS_range <- c(BaseCase = 0.05, Low = 0.01, High = 0.1)
p.SD_range <- c(BaseCase = 0.1,  Low = 0.05, High = 0.2)

paramNames <- c("p.HD", "p.HS", "p.SD")  # names of parameters in sensitivity

####### create lists with input parameters for the sensitivity analaysis 

l.tor.in <- vector("list", 3)
names(l.tor.in) <- paramNames
l.tor.in$p.HD <- cbind(p.HD = p.HD_range, input[-1])
l.tor.in$p.HS <- cbind(p.HS = p.HS_range, input[-2])
l.tor.in$p.SD <- cbind(p.SD = p.SD_range, input[-3])

## List of outputs
l.tor.out <- vector("list", 3)
names(l.tor.out) <- paramNames

## Run model on different parameters
for (i in 1:3) {
  l.tor.out[[i]] <- t(apply(l.tor.in[[i]], 1, MM.3state))[, 1]  # select the cost column 
}

## Data structure: ymean	ymin	ymax
m.tor <- matrix(unlist(l.tor.out), nrow = 3, ncol = 3, byrow = TRUE,
                dimnames = list(paramNames, c("basecase", "low", "high"	)))
m.tor # plot 

####### Plot tornado  ########################################
# file needed for tornado diagram
source("/Functions/tornado_diagram_code.R")

TornadoPlot(Parms = paramNames, Outcomes = m.tor, 
            titleName = "Tornado Plot", outcomeName = "Total costs")

