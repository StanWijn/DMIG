

rm(list = ls())      # clear memory (removes all the variables from the workspace)

####################### Input Model Parameters   ############################################

Strategies <- c("Routine Practice", "Primary Care", "Hospital Care")  

p.PCed <- 0.40   # Probability of early detection PC
p.HCed <- 0.45   # Probability of early detection HC
p.RPed <- 0.35   # Probability of early detection RP 

le.ed  <- 7      # Life expectancy after early detection
le.ld  <- 1      # Life expectancy after  late detection 

c.PCed <-  3900  # Total costs after early detection with PC
c.HCed <-  6200  # Total costs after early detection with HC
c.RPed <-  3030  # Total costs after early detection with RP
c.PCld <- 12800  # Total costs after  late detection with PC
c.HCld <- 14400  # Total costs after  late detection with HC
c.RPld <- 12020  # Total costs after  late detection with RP
wtp    <- 10000  # Define WTP


#################### Estimation of the Decision  Tree  ########################################

# the solution of the tree is the sum of the weights (probabilities of each leaf) times the reward of every leaf.

# for costs
c.RP <- p.RPed * c.RPed + (1 - p.RPed) * c.RPld # RP cost
c.PC <- p.PCed * c.PCed + (1 - p.PCed) * c.PCld # PC cost
c.HC <- p.HCed * c.HCed + (1 - p.HCed) * c.HCld # HC cost

# ...and for effects
e.RP <- p.RPed * le.ed + (1 - p.RPed) * le.ld   # RP life expectancy
e.PC <- p.PCed * le.ed + (1 - p.PCed) * le.ld   # PC life expectancy
e.HC <- p.HCed * le.ed + (1 - p.HCed) * le.ld   # HC life expectancy

LE <- c(e.RP, e.PC, e.HC)                       # vector of life expectancies
C  <- c(c.RP, c.PC, c.HC)                       # vector of total costs
names(LE) <- names(C) <- c("RP", "PC", "HC")    # names for the elements of the two vectors


# estimating pairwise incremental costs and effects
DC      <-  C  - C[c(1,1,1)]      # incremental costs
DE      <- LE - LE[c(1,1,1)]      # incremental effectiveness
ICER    <- DC / DE                # Incremental cost-effectiveness ratios
ICER[1] <- NA
NMB     <- LE * wtp - C           # calculate net monetary benefit
# create full incremental cost-effectiveness analysis table 
table           <- cbind(C, LE, DC, DE, round(ICER, 2))  # bind together the results into a table
table           <- as.data.frame(table)                  # as the table has both text and numbers, define as a data frame
colnames(table) <- c("Costs", "LE", "Inc.Cost", "Inc.Effects", "ICER")  # give column names
rownames(table) <- Strategies                                           # give row names
write.table(table, "table.txt")                                         # store the table 
table                                                                   # present the results

# load functions needed for plotting and CEA Frontier
source("Functions/CEA_functions.R")

### Get CEA frontier 
ce.mat <- cbind(Strategy = 1:3, 
                Cost = C, 
                Effectiveness = LE)
ce.front <- getFrontier(ce.mat, plot = FALSE)

### Plot frontier
## Using basic `plot` function
plot(ce.mat[ce.front, 3], ce.mat[ce.front, 2], 
     col = 1:3, pch = 1:3,
     xlab = "Effectiveness", ylab = "Cost")
lines(ce.mat[ce.front, 3], ce.mat[ce.front, 2])
legend("bottomright", Strategies, pch = 1:3, col = 1:3, bty = "n", cex = 0.8)
## Using `ggplot2` function
library(ggplot2)
ce.df <- data.frame(Strategy = Strategies,
                    Cost = C,
                    Effectiveness = LE)
ce.df
ce.front <- getFrontier(ce.mat, plot = F)
plotFrontier(CEmat = ce.df, frontier = ce.front)


##############################################################################
############# PSA ############################################################


#################### Estimation of the Decision  Tree  ########################################

sim <- 1000 # number of simulations

psa_dt <- function(){

Strategies <- c("Routine Practice", "Primary Care", "Hospital Care")  
p.PCed <- rbeta(1, 320, 480)   # Probability of early detection PC
p.HCed <- rbeta(1, 360, 440)   # Probability of early detection HC
p.RPed <- rbeta(1, 280, 520)   # Probability of early detection RP 

le.ed  <- rnorm(1, 7, 1.2)      # Life expectancy after early detection
le.ld  <- rnorm(1, 1, 0.03)      # Life expectancy after  late detection 

c.PCed <-  rgamma(1, shape = (3900^2)/(1000^2), rate = 3900/1000^2)        # Total costs after early detection with PC
c.HCed <-  rgamma(1, shape = (6200^2)/(1500^2), rate = 6200/1500^2)        # Total costs after early detection with HC
c.RPed <-  rgamma(1, shape = (6200^2)/(1500^2), rate = 6200/1500^2)        # Total costs after early detection with RP
c.PCld <-  rgamma(1, shape = (6200^2)/(1500^2), rate = 6200/1500^2)        # Total costs after  late detection with PC
c.HCld <-  rgamma(1, shape = (6200^2)/(1500^2), rate = 6200/1500^2)        # Total costs after  late detection with HC
c.RPld <-  rgamma(1, shape = (6200^2)/(1500^2), rate = 6200/1500^2)        # Total costs after  late detection with RP
wtp <- 80000
  
# for costs
c.RP <- p.RPed * c.RPed + (1 - p.RPed) * c.RPld # RP cost
c.PC <- p.PCed * c.PCed + (1 - p.PCed) * c.PCld # PC cost
c.HC <- p.HCed * c.HCed + (1 - p.HCed) * c.HCld # HC cost

# ...and for effects
e.RP <- p.RPed * le.ed + (1 - p.RPed) * le.ld   # RP life expectancy
e.PC <- p.PCed * le.ed + (1 - p.PCed) * le.ld   # PC life expectancy
e.HC <- p.HCed * le.ed + (1 - p.HCed) * le.ld   # HC life expectancy

LE <- c(e.RP, e.PC, e.HC)                       # vector of life expectancies
C  <- c(c.RP, c.PC, c.HC)                       # vector of total costs
names(LE) <- names(C) <- c("RP", "PC", "HC")    # names for the elements of the two vectors


# estimating pairwise incremental costs and effects
DC      <-  C  - C[c(1,1,2)]      # incremental costs
DE      <- LE - LE[c(1,1,2)]      # incremental effectiveness
ICER    <- DC / DE                # Incremental cost-effectiveness ratios
ICER[1] <- 0
ICER <- round(ICER)
NMB     <- round(LE * wtp - C)           # calculate net monetary benefit
return(NMB)

}


psa_output <- data.frame(RP = 1:sim, PC= NA, HC = NA)

for(i in 1:sim){
psa_output[i,] <- psa_dt()

}

# count treatments most cost-effective

psa_output$best <- as.factor(ifelse(psa_output$RP > psa_output$PC & psa_output$RP > psa_output$HC, "RP",
                          ifelse(psa_output$PC > psa_output$RP & psa_output$PC > psa_output$HC,"PC",
                                ifelse(psa_output$HC > psa_output$PC & psa_output$HC > psa_output$RP, "HC", NA))))
library(plyr)
count(psa_output$best)

# conclusions: 
# given a WTP of 80.000 euro's: HC is costeffective in 97% of the simulations.


















