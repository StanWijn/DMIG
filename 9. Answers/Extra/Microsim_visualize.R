#####################################################################################
## Function to make a simple plot from the data generated in a Microsimulation model ##
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
# Krijkamp EM, et al. Microsimulation modeling for health decision sciences using R: a tutorial. Med. Decis. Making. 2018; (in press). 

#####################################################################################
# ? Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the 
# collaborating institutions and may not be used, reproduced, modified, distributed or adapted in any way without written permission.

#####################################################################################
# Histogram showing variability in individual total costs
plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")

# Histogram showing variability in individual total QALYs
plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")


# PLOT THE DISTRIBUTION OF THE POPULATON ACROSS HEALTH STATES OVER TIME (TRACE)
# count the number of individuals in each health state at each cycle
m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
m.TR <- m.TR / n.i                                       # calculate the proportion of individuals 
colnames(m.TR) <- v.n                                    # name the rows of the matrix
rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")       # name the columns of the matrix

# Plot trace of first health state
plot(0:n.t, m.TR[, 1], type = "l", main = "Health state trace", 
     ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
# add a line for each additional state
for (n.s in 2:length(v.n)) {
  lines(m.TR[, n.s], col = n.s)       # adds a line to current plot
}
legend("topright", v.n, col = 1:3,    # add a legend to current plot
       lty = rep(1, 3), bty = "n", cex = 0.65)


