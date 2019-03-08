
######  Functions to perform a one way sensitifity analysis in R   ###########

##############################################################################
## Credits for the R Code: Fernando Alarid- Escudero and Hawre Jalal        ##
## Applied Methods of Cost-effectiveness Analysis in Healthcare             ##
##############################################################################


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
# Copyright 2017, THE HOSPITAL FOR SICK CHILDREN AND THE COLLABORATING INSTITUTIONS. 
# All rights reserved in Canada, the United States and worldwide.  
# Copyright, trademarks, trade names and any and all associated intellectual property are exclusively owned by THE HOSPITAL FOR SICK CHILDREN and the 
# collaborating institutions and may not be used, reproduced, modified, distributed or adapted in any way without written permission.

#####################################################################################


### ====================================================
###     Function for plotting One-way SA Diagrams
### ====================================================
owsa.plot.det <- function(param, outcomes,
                          paramName,
                          strategyNames, 
                          outcomeName = "Outcome"){
  ## Load required packages
  library(ggplot2)
  library(reshape2)
  library(scales)
  
  owsa.df <- data.frame(outcomes)
  colnames(owsa.df) <- strategyNames
  owsa.df$param <- param
  
  print(
    ggplot(data = melt(owsa.df, id.vars = "param", 
                       variable.name = "Strategy"), 
           aes(x = param, y = value, color = Strategy)) +
      geom_line() +
      xlab(paramName) +
      ggtitle("One-way sensitivity analysis", subtitle = outcomeName)+
      theme_bw(base_size = 14) +
      theme(legend.position = "bottom")
  )
}
