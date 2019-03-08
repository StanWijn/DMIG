
######   Function for plotting Tornado Diagrams                    ###########

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


#### =========================
### Tornado Plots
## Test code:
# paramNames <-  c( "Param 1 [Low/High CI]",
#                   "Param 2 [Low/High CI]",
#                   "Param 3 [-/+ 15%]",
#                   "Param 4 [-/+ 15%]"
# )
# 
# # data structure: ymean	ymin	ymax
# data <- matrix(c(100,	80,	120,
#                  100,	25,	150,
#                  100, 95,	120,
#                  100, 75, 160), nrow = 4, ncol = 3, byrow = TRUE)
# 
# data
# TornadoPlot(Parms = paramNames, Outcomes = data, titleName = "Tornado Plot")
# Parms = paramNames
# Outcomes = data
# titleName = "Tornado Plot"

### ====================================================
###     Function for plotting Tornado Diagrams
### ====================================================
TornadoPlot <- function(Parms, Outcomes, titleName, outcomeName){
  # Parm:        vector with parameter names  
  # Outcomes:    matrix including parameter specific outcomes (number of Parm x 3)
  # titleName:   title of the plot (e.g Tornado Plot)
  # outcomeName: name of the outcome shown in the Tornado plot
  
  library(ggplot2)
  library(reshape2)
  library(scales)
  
  # Grouped Bar Plot
  # Determine the overall optimal strategy
  paramNames2 <- Parms
  
  # Combine the parameter list with the data
  ymean <- Outcomes[1, 1]
  
  yMin <- Outcomes[, 2] - ymean
  yMax <- Outcomes[, 3] - ymean
  ySize <- abs(yMax - yMin)  # High value - Low value
  
  rankY<- order(ySize)
  nParams <- length(paramNames2)
  
  Tor <- data.frame(
    Parameter = c(paramNames2[rankY], paramNames2[rankY]),  
    Level = c(rep("Low", nParams), rep("High", nParams)),
    value = ymean + c(yMin[rankY], yMax[rankY]),
    sort = seq(1, nParams)
  )
  
  #re-order the levels in the order of appearance in the data.frame
  Tor$Parameter2 <- ordered(Tor$Parameter, Tor$Parameter[1:(length(Tor$Parameter) / 2)])
  # Tor$Parameter2 <- factor(Tor$Parameter, as.character(Tor$Parameter))
  #Define offset as a new axis transformation. Source: http://blog.ggplot2.org/post/25938265813/defining-a-new-transformation-for-ggplot2-scales  
  offset_trans <- function(offset = 0) {
    trans_new(paste0("offset-", format(offset)), function(x) x-offset, function(x) x + offset)
  }
  #Plot the Tornado diagram.
  txtsize <- 12
  print(
    ggplot(Tor[Tor$Level == "Low", ], aes(x = Parameter2, y = value, fill = level)) +
      geom_bar(stat = "identity", fill = "blue") +
      ggtitle("Tornado Plot", subtitle = outcomeName) +
      scale_fill_discrete("Parameter Level: ", l = 50) +
      scale_y_continuous(name = "$", trans = offset_trans(offset = ymean)) +
      scale_x_discrete(name = "Parameter") +
      geom_bar(data = Tor[Tor$Level == "High", ], aes(x = Parameter2, y = value, fill = level), stat = "identity", fill = "red", alpha = 0.5) +
      geom_hline(yintercept = ymean, linetype = "dotted", size = 0.5) +
      theme_bw(base_size = 14) +
      coord_flip() +
      theme(legend.position = "bottom",
            legend.title = element_text(size = txtsize, angle = 0, hjust = 1),
            legend.key = element_rect(colour = "black"),
            legend.text = element_text(size = txtsize),
            title = element_text(face="bold", size=15),
            axis.title.x = element_text(face = "bold", size = txtsize),
            axis.title.y = element_text(face = "bold", size = txtsize),
            axis.text.y = element_text(size = txtsize),
            axis.text.x = element_text(size = txtsize),
            axis.ticks.y = element_blank())
  )
  # ggsave(paste("results/", titleName,".png"))
}
