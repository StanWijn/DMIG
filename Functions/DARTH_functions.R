######################
#### CEA Frontier ####
######################
getFrontier <- function(CEmat, maxWTP = Inf, plot = TRUE){
  # Name: getFrontier.R
  # Goal: Find the CEA frontier, up to a given WTP level, by 
  #       identifying strategies with the highest NMB
  # Originally written by: Sze Suen on Feb 25, 2015
  # Modified by: Fernando Alarid-Escudero on July 20, 2015 
  # Notes: 
  #    ~Frontier strategies are displaced on the R output screen and 
  #      plotted in red on the scatter plot.  
  #
  #    ~User needs to provide a csv file of costs and QALYs 
  #	    (CostQalyInputFile_online_supp.csv) inside the folder specified    
  #	    below (inputFolder). The CSV should have three columns (labeled 
  #     in first row) in this order: 
  #      Strategy number, costs, and QALYs.
  #
  #    ~User can specify the maximum willingness-to-pay level to 
  #      consider (maxWTP).  Can be Inf for infinity.
  #
  #    ~QALY-reducing strategies will be on the frontier if they save
  #      enough money (script assumes maximum willingness to save money 
  #      per QALY lost is equivalent to maximum willingness to pay per QALY
  #      gained). If the user does not wish to consider such policies as
  #      being on the frontier, do not include strategies with negative 
  #      QALYs in the input csv file.
  #
  #    ~Script does not use the one-line code cited in the text
  #      as the max function is slow. This implementation is
  #      faster and methodologically does the same thing.
  #
  #    ~May take a few minutes if thousands of strategies and 
  #       processing resources are low.  Please be patient.
  #
  #    Please cite article if this code is used.
  #
  # USER INPUTS:
  #inputFolder <- "CostEffectivenessFrontier_MDM/"
  #maxWTP <- Inf        # any positive value or Inf
  
  ## Clean everythng from workspace
  #rm(list=ls())
  ####################################################################
  ####################################################################
  
  # check for duplicated strategies
  dups <- CEmat[c(duplicated(CEmat[,2:3]) | duplicated(CEmat[,2:3], fromLast = TRUE)),1]
  
  # initialize some variables
  costsCol <- 2; qalyCol <- 3
  numStrat <- nrow(CEmat)
  
  # find WTP levels to test so that all strategies on frontier will be captured
  # this means testing on either side of all NMB intersections, which are just all the pairwise ICERs
  ICERmat <- matrix(1, numStrat, numStrat)
  for (i in 1:numStrat ) {
    indexStrat <- matrix(1, numStrat, 3)
    indexStrat[,costsCol] <- indexStrat[,costsCol]*CEmat[i,costsCol]
    indexStrat[,qalyCol] <- indexStrat[,qalyCol]*CEmat[i,qalyCol]
    delCostQalys <- CEmat - indexStrat
    ICERmat[,i] <- delCostQalys[,costsCol] / delCostQalys[,qalyCol]
  }  
  intersections <- sort(unique(c(ICERmat)))
  intersections <- intersections[is.finite(intersections)]
  WTPtestPoints <- c(0, intersections [intersections >= 0 & intersections <= maxWTP ], maxWTP)
  
  # Find the strategy with the max NMB at each of the WTP test points
  indiciesOfMax <- vector()
  NMBmat <- matrix(0, numStrat, length(WTPtestPoints))
  for (i in 1:length(WTPtestPoints) ) {
    NMBmat[,i] <- (WTPtestPoints[i]*CEmat[,qalyCol]) - CEmat[,costsCol]
  }
  if (is.infinite(maxWTP)) {
    #WTP of infinity means costs are not considered
    NMBmat[,length(WTPtestPoints)] = CEmat[,qalyCol] - (0*CEmat[,costsCol]); 
  }
  maxVals <- apply(NMBmat, 2, max)  #find strategy that maximizes NMB at each WTP
  for (i in 1:length(WTPtestPoints) ) {  #find all strategies that match max at each WTP
    indiciesOfMax <- c(indiciesOfMax,which( NMBmat[,i] == maxVals[i]))
  }
  frontier <- unique(indiciesOfMax)  #find strategy that maximizes NMB at each WTP
  
  if (plot == TRUE){
    # display out: make plot and print to output screen
    plot(CEmat[frontier,qalyCol], CEmat[frontier,costsCol], col = 'red', pch = 16, 
         xlab = "Effectiveness", ylab = "Cost")
    points(CEmat[,qalyCol], CEmat[,costsCol])
    lines(CEmat[frontier,qalyCol], CEmat[frontier,costsCol])
    if (length(dups)>0){
      warning("Strategies have the same costs and benefits (displayed above)")
      print(dups)
    }
  }
  sprintf("Frontier is formed by strategies: %s", paste( sort(CEmat[frontier,1]), collapse=" "))
  
  return(frontier)
}

plotFrontier <- function(CEmat, frontier, 
                         ncol = 1,
                         coord.flip = F,
                         txtsize = 12)
{
  # A function to plot CE frontier
  # USER INPUTS:
  #   CEmat: A CE matrix arranged as: Col1: Strategy; Col2: Cost; Col3: Effectiveness
  # Create a dataframe from matrix
  CEmat.df <- data.frame(CEmat)
  colnames(CEmat.df)[3] <- "Effectiveness"
  n.strategies <- nrow(CEmat.df)
  # Make Strategies as factor
  CEmat.df$Strategy <- as.factor(CEmat.df$Strategy)
  #
  if (coord.flip == T){
    ggplot(CEmat.df, aes(Effectiveness, Cost)) +
      geom_point(aes(color = Strategy, shape = Strategy), size = 4) + 
      coord_flip() +
      ggtitle("Cost-Effectiveness Frontier") +
      geom_point(data = CEmat.df[frontier,], 
                 aes(Effectiveness, Cost, shape = Strategy, color = Strategy), size = 4) +
      geom_line(data = CEmat.df[frontier,], aes(Effectiveness, Cost)) +
      scale_shape_manual(values = 0:(n.strategies-1)) +
      guides(shape = guide_legend(ncol = ncol)) +
      theme_bw() +
      theme(title = element_text(face = "bold", size = txtsize+2),
            axis.title.x = element_text(face = "bold", size = txtsize),
            axis.title.y = element_text(face = "bold", size = txtsize),
            axis.text.y = element_text(size = txtsize),
            axis.text.x = element_text(size = txtsize))
  } else {
    ggplot(CEmat.df, aes(Effectiveness, Cost)) +
      geom_point(aes(color = Strategy, shape = Strategy), size = 4) + 
      ggtitle("Cost-Effectiveness Frontier") +
      geom_point(data = CEmat.df[frontier,], 
                 aes(Effectiveness, Cost, shape = Strategy, color = Strategy), size = 4) +
      geom_line(data = CEmat.df[frontier,], aes(Effectiveness, Cost)) +
      scale_shape_manual(values = 0:(n.strategies-1)) +
      guides(shape = guide_legend(ncol = ncol)) +
      theme_bw() +
      theme(title = element_text(face = "bold", size = txtsize+2),
            axis.title.x = element_text(face = "bold", size = txtsize),
            axis.title.y = element_text(face = "bold", size = txtsize),
            axis.text.y = element_text(size = txtsize),
            axis.text.x = element_text(size = txtsize))
  }
}

################################
#### Metamodeling Functions ####
################################

#######################
#### PSA Functions ####
#######################
TornadoOpt <-function(Parms,Outcomes){
  # Grouped Bar Plot
  # Determine the overall optimal strategy
  opt<-which.max(colMeans(Outcomes))
  # calculate min and max vectors of the parameters (e.g., lower 2.5% and 97.5%)
  X <- as.matrix(Parms)
  y <- as.matrix(Outcomes[,opt])
  ymean <- mean(y)
  n <- nrow(Parms)
  nParams <- ncol(Parms)
  #paramNames <- Names[seq(8,7+nParams)]
  paramNames <- colnames(Parms)
  Parms.sorted <- apply(Parms,2,sort,decreasing=F)#Sort in increasing order each column of Parms
  lb <- 2.5
  ub <- 97.5 
  Xmean <- rep(1,nParams) %*% t(colMeans(X))
  XMin <- Xmean
  XMax <- Xmean
  paramMin <- as.vector(Parms.sorted[round(lb*n/100),])
  paramMax <- as.vector(Parms.sorted[round(ub*n/100),])
  paramNames2 <- paste(paramNames, "[", round(paramMin,2), ",", round(paramMax,2), "]")
  
  diag(XMin) <- paramMin
  diag(XMax) <- paramMax
  
  XMin <- cbind(1, XMin)
  XMax <- cbind(1, XMax)
  
  X <- cbind(1,X)
  B <- solve(t(X) %*% X) %*% t(X) %*% y
  yMin <- XMin %*% B - ymean
  yMax <- XMax %*% B - ymean
  ySize <- abs(yMax - yMin) 
  
  rankY<- order(ySize)
  xmin <- min(c(yMin, yMax)) + ymean
  xmax <- max(c(yMin, yMax)) + ymean
  
  Tor <- data.frame(
    Parameter=c(paramNames2[rankY],paramNames2[rankY]),  
    Level=c(rep("Low",nParams),rep("High",nParams)),
    value=ymean+c(yMin[rankY],yMax[rankY]),
    sort=seq(1,nParams)
  )
  #re-order the levels in the order of appearance in the data.frame
  Tor$Parameter2 <- factor(Tor$Parameter, as.character(Tor$Parameter))
  #Define offset as a new axis transformation. Source: http://blog.ggplot2.org/post/25938265813/defining-a-new-transformation-for-ggplot2-scales  
  offset_trans <- function(offset=0) {
    trans_new(paste0("offset-", format(offset)), function(x) x-offset, function(x) x+offset)
  }
  #Plot the Tornado diagram.
  txtsize<-12
  ggplot(Tor[Tor$Level=="Low",], aes(x=Parameter2,y=value, fill=Level)) +
    geom_bar(stat="identity") +
    ggtitle("Tornado Diagram")+
    scale_fill_discrete("Parameter Level: ", l=50)+
    scale_y_continuous(name="Net Benefit",trans=offset_trans(offset=ymean)) +
    scale_x_discrete(name="Parameter") +
    geom_bar(data=Tor[Tor$Level=="High",], aes(x=Parameter2,y=value, fill=Level), stat="identity") +
    geom_hline(yintercept = ymean, linetype = "dotted", size=0.5) +
    theme_bw()+
    theme(legend.position="bottom",legend.title=element_text(size = txtsize,angle = 0, hjust = 1),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize),
          axis.ticks.y = element_blank())+
    coord_flip()  
}

TornadoAll <-function(Strategies,Parms,Outcomes){
  opt<-which.max(colMeans(Outcomes))
  # calculate min and max vectors of the parameters (e.g., lower 2.5% and 97.5%)
  X <- as.matrix(Parms)
  y <- as.matrix(Outcomes[,opt])
  Y <- as.matrix(Outcomes)
  ymean <- mean(y)
  n <- nrow(Parms)
  nParams <- ncol(Parms)
  #paramNames <- Names[seq(8,7+nParams)]
  paramNames <- colnames(Parms)
  Parms.sorted <- apply(Parms,2,sort,decreasing=F)#Sort in increasing order each column of Parms
  lb <- 2.5
  ub <- 97.5 
  Xmean <- rep(1,nParams) %*% t(colMeans(X))
  XMin <- Xmean
  XMax <- Xmean
  paramMin <- as.vector(Parms.sorted[round(lb*n/100),])
  paramMax <- as.vector(Parms.sorted[round(ub*n/100),])
  
  diag(XMin) <- paramMin
  diag(XMax) <- paramMax
  
  XMin <- cbind(1, XMin)
  XMax <- cbind(1, XMax)
  
  X <- cbind(1,X)
  B <- solve(t(X) %*% X) %*% t(X) %*% y
  #install.packages("matrixStats")
  library(matrixStats)
  bigBeta <- solve(t(X) %*% X) %*% t(X) %*% Y
  yMin <- rowMaxs(XMin %*% bigBeta - ymean)
  yMax <- rowMaxs(XMax %*% bigBeta - ymean)
  ySize <- abs(yMax - yMin) 
  
  rankY<- order(ySize)
  xmin <- min(c(yMin, yMax)) + ymean
  xmax <- max(c(yMin, yMax)) + ymean
  
  paramNames2 <- paste(paramNames, "[", round(paramMin,2), ",", round(paramMax,2), "]")
  
  strategyNames<-Strategies
  strategyColors <- c("red","darkgreen","blue")
  
  ## Polygon graphs:
  nRect <- 0
  x1Rect <- NULL
  x2Rect <- NULL
  ylevel <- NULL
  colRect <- NULL
  
  for (p in 1:nParams){
    xMean <- colMeans(X)
    xStart = paramMin[rankY[p]]
    xEnd = paramMax[rankY[p]]
    xStep = (xEnd-xStart)/1000
    for (x in seq(xStart,xEnd, by = xStep)){
      #for each point determine which one is the optimal strategy
      xMean[rankY[p] + 1] <- x 
      yOutcomes <- xMean %*% bigBeta
      yOptOutcomes <- max(yOutcomes)
      yOpt <- which.max(yOutcomes)
      if (x == xStart){
        yOptOld <- yOpt
        y1 <- yOptOutcomes
      }
      #if yOpt changes, then plot a rectangle for that region
      if (yOpt != yOptOld | x == xEnd){
        nRect <- nRect + 1
        x1Rect[nRect] <- y1
        x2Rect[nRect] <- yOptOutcomes
        ylevel[nRect] <- p
        colRect[nRect] <- strategyColors[yOptOld]
        yOptOld <- yOpt
        y1 <- yOptOutcomes
      }
    }
  }
  
  txtsize <- 12
  d=data.frame(x1=x2Rect, x2=x1Rect, y1=ylevel-0.4, y2=ylevel+0.4, t=colRect, r = ylevel)
  ggplot(d, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = t)) +
    ggtitle("Torando Diagram") + 
    xlab("Expected NHB") +
    ylab("Parameters") + 
    geom_rect()+
    theme_bw() + 
    scale_y_continuous(limits = c(0.5, nParams + 0.5),breaks=seq(1:8), labels=paramNames2[rankY]) +
    #scale_y_discrete(breaks=seq(1:8), labels=paramNames2[rankY]) + 
    scale_fill_discrete(name="Optimal\nStrategy",
                        #breaks=c("ctrl", "trt1", "trt2"),
                        labels=strategyNames,
                        l=50) + 
    geom_vline(xintercept=ymean, linetype="dotted") + 
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}
MNL.SA <- function(parm,range,Strategies,Parms,Outcomes){
  #Extract parameter column number in Parms matrix
  x<-which(colnames(Parms)==parm)
  #Calculate the preferred strategy (i.e., optimal startegy) in terms of its cost-efectiveness 
  #for each of the simulations.
  Optimal <- data.frame(max.col(Outcomes)); names(Optimal)<-"Strategy"
  
  MNL <- data.frame(Optimal,Parms)
  MNL.vglm = vglm(Strategy ~ ., data = MNL, family=multinomial)
  
  y = seq(2.5,97.5,length=400) #vector to define 400 samples between the 2.5th and 97.5th percentiles
  j = round(y*(length(Parms[,x])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
  vector<-sort(Parms[j,x])
  #vector<-seq(from=range[1],to=range[2],length.out=201)
  
  #Generate matrix to use for prediction 
  Sim.fit <- matrix(rep(colMeans(Parms)), 
                    nrow = length(vector),
                    ncol = ncol(Parms), byrow = T)
  Sim.fit[, x] <- vector
  MNL.fit<-data.frame(Sim.fit) #Transform to data frame, the format required for predict
  colnames(MNL.fit)<-colnames(Parms) #Name data frame's columns with parameters' names
  
  #Predict Outcomes using MMMR Metamodel fit
  plotdata = data.frame(predict(MNL.vglm, newdata = MNL.fit, type = "response"))
  
  colnames(plotdata) <- Strategies
  plotdata = stack(plotdata, select=Strategies)
  plotdata = cbind(MNL.fit, plotdata) 
  
  plotdata$parm<-plotdata[,parm];
  
  txtsize<-12 #Text size for the graphs
  ggplot(data = plotdata, aes(x = parm, y = values, color = ind)) +
    geom_point(size = 2) + #maybe shape=3;18;21;124; Other shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
    geom_line() +
    ggtitle("Multinomial sensitivity analysis") + 
    xlab(parm) +
    ylab("Probability of Strategy Being Optimal") +
    scale_colour_hue("Strategy", l=50) +
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
  
  #   MNL.fit <- expand.grid(parm=vector)
  #   Sim.MNL <- array(0,dim=c(dim(MNL.fit)[1],length(Parms)))
  #   for (i in 1:(length(Parms))){
  #     if (i==x){
  #       Sim.MNL[,i] <-MNL.fit[,1]
  #     } else {
  #       Sim.MNL[,i] <-mean(Parms[,i])
  #     }
  #   }
  #   
  #   MNL.fit = data.frame(Sim.MNL)
  #   colnames(MNL.fit)<-colnames(Parms)
  #   plotdata = data.frame(predict(MNL.vglm, newdata = MNL.fit, type = "response"))
  #   colnames(plotdata) <- Strategies
  #   plotdata = stack(plotdata, select=Strategies)
  #   plotdata = cbind(MNL.fit, plotdata) 
  #   
  #   txtsize<-12 #Text size for the graphs
  #   ggplot(data = plotdata, aes(x = muDieCancer, y = values, color = ind)) +
  #     geom_point(size = 2) + #maybe shape=3;18;21;124; Other shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  #     geom_line() +
  #     ggtitle("Multinomial sensitivity analysis") + 
  #     xlab("muDieCancer") +
  #     ylab("Probability of Strategy Being Optimal") +
  #     scale_colour_hue("Strategy", l=50) +
  #     theme_bw() +
  #     theme(legend.position="bottom",legend.title=element_text(size = txtsize),
  #           legend.key = element_rect(colour = "black"),
  #           legend.text = element_text(size = txtsize),
  #           title = element_text(face="bold", size=15),
  #           axis.title.x = element_text(face="bold", size=txtsize),
  #           axis.title.y = element_text(face="bold", size=txtsize),
  #           axis.text.y = element_text(size=txtsize),
  #           axis.text.x = element_text(size=txtsize))
}
PlaneCE<-function(Strategies,Outcomes){
  library(reshape2)
  ndep<-length(Strategies)*2 #Determine number of outcomes for all starteges, i.e., cost and effectiveness
  ind_c<-seq(1,(ndep-1),by=2) #Index to extract the costs from matrix Outcomes
  ind_e<-seq(2,ndep,by=2) #Index to extract the effectiveness from matrix Outcomes
  Cost<-melt(Outcomes[,ind_c],variable.name = "Strategy")
  levels(Cost$Strategy)<-Strategies
  Eff<-melt(Outcomes[,ind_e],variable.name="Strategy")
  levels(Eff$Strategy)<-Strategies
  CE<-cbind(Cost,Eff[,2])
  colnames(CE)<-c("Strategy","Cost","Effectiveness")
  
  #Dataframe with means of strategies. Source: http://stackoverflow.com/questions/18729724/ggplot2-scatter-plot-with-overlay-of-means-and-bidirectional-sd-bars
  Means <- ddply(CE,.(Strategy),summarise,
                 N = length(Cost),
                 Cost.mean = mean(Cost),
                 Eff.mean = mean(Effectiveness))
  
  #Define ggplot object
  txtsize<-12
  ggplot(Means, aes(x = Eff.mean, y = Cost.mean, color=Strategy)) + 
    geom_point(size=4, aes(shape=Strategy)) +
    ggtitle("Cost-Effectiveness Plane") +
    scale_colour_discrete(l=50) +  # Use a slightly darker palette than normal
    scale_y_continuous(labels = dollar)+
    scale_x_continuous(breaks=number_ticks(6), labels=comma)+
    xlab("Effectiveness")+
    ylab("Cost")+
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
}
#------------#
#### CEAF ####
#------------#
CEAF <- function(wtp, strategies, Outcomes){
  library(reshape2)
  # Outcomes must be ordered in a way that for each strategy the cost must appear first then the effectiveness
  # Create scalar with number of simulations
  n.sim <- nrow(Outcomes)
  # Create scalar with number of strategies
  n.strategies <- length(strategies)
  # Vector to index costs
  costInd <- seq(1, 2*n.strategies, by = 2) 
  # Vector to index effectiveness
  effInd  <- seq(2, 2*n.strategies, by = 2) 
  # Matrix to store NHB for each strategy
  NMB <- array(0, dim = c(n.sim, n.strategies))
  colnames(NMB) <- strategies
  # Matrix to store indicator of CEAC
  cea  <- array(0, dim = c(length(wtp), n.strategies))
  # Vector to store indicator of strategy at CEAF
  ceaf.ind <- array(0, dim = c(length(wtp), 1))
  
  for(l in 1:length(wtp)){
    NMB <-  wtp[l]*Outcomes[, effInd] - Outcomes[, costInd] # Effectiveness minus Costs, with vector indexing
    # Calculate point of CEAF, i.e., the strategy with the highest expected NMB
    ceaf.ind[l, 1] <- which.max(colMeans(NMB))
    # Calculate points in CEAC, i.e, the probability that each strategy is cost-effective
    Max.NMB <- max.col(NMB)
    opt <- table(Max.NMB)
    cea[l, as.numeric(names(opt))] <- opt/n.sim
  }
  ceaf <- cea[cbind(1:length(wtp), ceaf.ind)]
  
  cea <- data.frame(cbind(wtp, cea, ceaf))
  colnames(cea) <- c("WTP", strategies, "Frontier")
  
  ceac <- melt(cea, id.vars = "WTP") 
  colnames(ceac)[2] <- "Strategy"
  return(list(ceac = ceac,
              ceaf.ind = ceaf.ind,
              ceaf = ceaf))
}

plotCEAF <- function(ceaf,
                     title = "Cost-Effectiveness Acceptability Curves and Frontier", 
                     txtsize = 12,
                     currency = "$"){
  library(ggplot2)
  strats <- 1:(length(unique(ceaf$Strategy))-1)
  point.shapes <- c(strats+14, 0) # Shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  colors <- c(gg_color_hue(length(strats)), "#696969")
  point.size <- c(rep(2, length(strats)), 4) # Trick consists on firts define size as aes then manually change it
  
  ggplot(data = ceaf, aes(x = WTP/1000, y = value)) +
    geom_point(aes(shape = Strategy, color = Strategy, size = Strategy)) +
    geom_line(aes(color = Strategy)) +
    ggtitle(title) + 
    # scale_colour_hue(l=50, values=colors) +
    scale_x_continuous(breaks=number_ticks(20))+
    scale_shape_manual(values=point.shapes) +
    #scale_shape(solid=TRUE) +
    scale_color_manual(values=colors) + 
    scale_size_manual(values = point.size) +
    #scale_alpha_manual(values=c(rep(0, length(strats)), 0.5)) + 
    xlab(paste("Willingness to pay (Thousand ", currency,"/QALY)", sep = "")) +
    ylab("Pr Cost-Effective") +
    theme_bw() +
    theme(legend.title=element_text(size = txtsize), #legend.position="right",
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          # legend.background = element_rect(fill=alpha(0.4)),
          title = element_text(face="bold", size=14),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize),
          legend.position = c(0.5, 0.5))
}

#------------#
#### CEAF ####
#------------#
ceaf <- function(v.wtp, strategies, m.e, m.c, currency = "$"){
  library(reshape2)
  library(ggplot2)
  library(scales)
  n.sim <- nrow(m.e)
  n.str  <- ncol(m.e)
  # Matrix to store indicator of CEAC
  m.cea  <- array(0, dim = c(length(v.wtp), n.str))
  # Vector to store indicator of strategy at CEAF
  v.ceaf <- array(0, dim = c(length(v.wtp), 1))
  
  for(l in 1:length(v.wtp)){
    m.nmb <-  v.wtp[l]*m.e - m.c # Effectiveness minus Costs
    # Calculate point of CEAF, i.e., the strategy with the highest expected NMB
    v.ceaf[l, 1] <- which.max(colMeans(m.nmb))
    # Calculate points in CEAC, i.e, the probability that each strategy is cost-effective
    max.nmb <- max.col(m.nmb)
    opt <- table(max.nmb)
    m.cea[l, as.numeric(names(opt))] <- opt/n.sim
  }
  m.ceaf <- m.cea[cbind(1:length(v.wtp), v.ceaf)]
  
  df.cea <- data.frame(cbind(v.wtp, m.cea, m.ceaf))
  colnames(df.cea) <- c("WTP", strategies, "Frontier")
  
  df.ceac <- melt(df.cea, 
                  id.vars = "WTP", 
                  variable.name = "Strategy") 
  
  ## Plot CEAC & CEAF
  # Format to plot frontier
  strats <- 1:(length(unique(df.ceac$Strategy))-1)
  point.shapes <- c(strats+14, 0) # Shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  colors <- c(gg_color_hue(length(strats)), "#696969")
  point.size <- c(rep(2, length(strats)), 4) # Trick consists on firts define size as aes then manually change it
  # Plot CEAC & CEAF
  print(
  ggplot(data = df.ceac, aes(x = WTP/1000, y = value)) +
    geom_point(aes(shape = Strategy, color = Strategy, size = Strategy)) +
    geom_line(aes(color = Strategy)) +
    ggtitle("Cost-Effectiveness Acceptability Curves and Frontier") + 
    # scale_colour_hue(l=50, values=colors) +
    scale_x_continuous(breaks=number_ticks(20))+
    scale_shape_manual(values=point.shapes) +
    #scale_shape(solid=TRUE) +
    scale_color_manual(values=colors) + 
    scale_size_manual(values = point.size) +
    #scale_alpha_manual(values=c(rep(0, length(strats)), 0.5)) + 
    xlab(paste("Willingness to pay (Thousand ", currency,"/QALY)", sep = "")) +
    ylab("Pr Cost-Effective") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  )
  
}


#----------------------------#
#### Expected Loss Curves ####
#----------------------------#
elc <- function(v.wtp, strategies, m.e, m.c){
  library(reshape2)
  library(ggplot2)
  library(scales)
  n.sim <- nrow(m.e)
  n.str  <- ncol(m.e)
  m.exp.loss <- matrix(0, nrow = length(v.wtp), ncol = n.str)
  for(l in 1:length(v.wtp)){
    m.nmb <- m.e*v.wtp[l] - m.c # Effectiveness minus Costs, with vector indexing
    max.str <- max.col(m.nmb)
    m.loss <- m.nmb - m.nmb[cbind(1:n.sim, max.str)]
    m.exp.loss[l, ] <- colMeans(m.loss)
  }
  # Optimal strategy based on lowest expected loss
  optimal.str <- max.col(m.exp.loss)
  # Expected loss of optimal strategy
  optimal.el <- m.exp.loss[cbind(1:length(v.wtp), optimal.str)]
  # Format expected loss for plotting
  df.exp.loss <- data.frame(cbind(v.wtp, m.exp.loss, optimal.el))
  colnames(df.exp.loss) <- c("WTP", strategies, "Frontier & EVPI")
  df.exp.loss.plot <- melt(df.exp.loss, 
                           id.vars = "WTP", 
                           variable.name = "Strategy")
  ## Plot expected losses
  # Format to plot frontier
  strats <- 1:(length(unique(df.exp.loss.plot$Strategy))-1)
  point.shapes <- c(strats+14, 0) # Shapes: http://sape.inf.usi.ch/quick-reference/ggplot2/shape
  colors <- c(gg_color_hue(length(strats)), "#696969")
  point.size <- c(rep(2, length(strats)), 4) # Trick consists on firts define size as aes then manually change it
  # Plot ELC
  print(
  ggplot(data = df.exp.loss.plot, aes(x = WTP/1000, y = -value)) +
    geom_point(aes(shape = Strategy, color = Strategy, size = Strategy)) +
    geom_line(aes(color = Strategy)) +
    ggtitle("Expected Loss Curves") + 
    #scale_colour_hue(l=50, values=colors) +
    scale_x_continuous(breaks=number_ticks(20))+
    scale_y_continuous(labels = comma, breaks = number_ticks(8)) +
    xlab("Willingness to Pay (Thousand $/QALY)") +
    ylab("Expected loss ($)") +
    scale_shape_manual(values=point.shapes) +
    scale_color_manual(values=colors) + 
    scale_size_manual(values = point.size) +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  )
}
OneWaySA <- function(strategies, y, x, 
                     parm, range = NULL,
                     poly.order = 2,
                     txtsize = 12){
  # Extract parameter column number in x matrix
  par.col <- which(colnames(x)==parm)
  dep <- length(strategies) #Number of dependent variables, i.e., strategies outcomes
  indep <- ncol(x) #Number of independent variables, i.e., parameters
  Sim <- data.frame(y,x)
  #Determine range of of the parameer to be plotted
  if (is.null(range)){ #If user does not define a range
    #Default range given by the domain of the parameter's sample
    #vector to define 400 samples between the 2.5th and 97.5th percentiles
    percentiles = seq(2.5, 97.5, length = 400) 
    j = round(percentiles*(length(x[,par.col])/100)) #indexing vector;j=round(y*n/100) where n is the size of vector of interest
    vector<-sort(x[j, par.col]) 
  }
  else{ #If user defines a range  
    vector <- seq(range[1],range[2],length.out=400)
  }
  
  #Generate a formula by pasting column names for both dependent and independent variables. Imposes a 1 level interaction
  f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep], collapse=','), 
                        ') ~ (','poly(', parm,',', poly.order,',raw=TRUE)+',
                        paste(colnames(x)[-par.col], collapse='+'),')'))
  #Run Multiple Multivariate Regression (MMR) Metamodel
  Oway.mlm = lm(f, data=Sim)
  
  # Create data frame with all combinations between both parameters of interest
  OWSA <- data.frame(parm = vector)
  
  #Generate matrix to use for prediction 
  Sim.fit <- matrix(rep(colMeans(x)),
                    nrow = length(vector), 
                    ncol = ncol(x), 
                    byrow = T)
  Sim.fit[, par.col] <- OWSA[, 1]
  # Transform to data frame, the format required for predict
  Sim.fit <- data.frame(Sim.fit)
  # Name data frame's columns with parameters' names
  colnames(Sim.fit) <- colnames(x) #Name data frame's columns with parameters' names
  
  # Predict Outcomes using MMMR Metamodel fit
  Sim.OW = data.frame(predict(Oway.mlm, newdata = Sim.fit))
  colnames(Sim.OW) <- strategies #Name the predicted outcomes columns with strategies names
  
  #Reshape dataframe for ggplot
  Sim.OW = stack(Sim.OW, select=strategies) #
  Sim.OW = cbind(Sim.fit, Sim.OW) #Append parameter's dataframe to predicted outcomes dataframe
  
  #A simple trick I use to define my variables in my functions environment
  #Borrowed from http://stackoverflow.com/questions/5106782/use-of-ggplot-within-another-function-in-r
  Sim.OW$parm<-Sim.OW[,parm];
  
  owsa <- ggplot(data = Sim.OW, aes(x = parm, y = values, color = ind)) +
    geom_line() +
    ggtitle("One-way sensitivity analysis") + #\n Net Health Benefit
    xlab(parm) +
    ylab("E[Outcome]") +
    scale_colour_hue("Strategy", l=50) +
    scale_x_continuous(breaks=number_ticks(6)) + #Adjust for number of ticks in x axis
    scale_y_continuous(breaks=number_ticks(6)) +
    theme_bw() +
    theme(legend.position="bottom",legend.title=element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=txtsize+2),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
  return(owsa)
}
TwoWaySA <- function(strategies = NULL, y, x,  
                     parm1, parm2, 
                     range1 = NULL, range2 = NULL,
                     poly.order = 2,
                     maximize = TRUE,
                     txtsize = 12){
  # Load dependencies
  require(ggplot2)
  # Create scalar with number of strategies (i.e. number of columns of 
  # `y`)
  n.strategies <- ncol(y)
  # If the name of the strategies is not provided, generate a generic vector
  if (is.null(strategies)){
    strategies <- paste(rep("Strategy_", n.strategies), seq(1, n.strategies), " ", sep = "")
  }
  # Extract parameter column number in Parms matrix
  x1  <- which(colnames(x)==parm1)
  x2  <- which(colnames(x)==parm2)
  # Number of dependent variables, i.e., strategies
  dep <- n.strategies
  # Number of independent variables, i.e., parameters
  indep <- ncol(x) 
  # Generate data frame with outputs and inputs of the model
  Sim <- data.frame(y, x)
  
  # Determine range of of the parameer to be plotted base in user's inputs
  if (is.null(range1) & is.null(range2)){ # If user doesn't define either range
    range1 <- range(x[,x1])
    range2 <- range(x[,x2])
  }
  else if (is.null(range2)){ # If user only defines first range
    range2 <- range(x[,x2])   
  }
  else if (is.null(range1)){ # If user only defines second range
    range1 <- range(x[,x1])   
  }
  # Create vectors with values to use to evaluate TWSA
  vector1 <- seq(from = range1[1], 
                 to = range1[2],
                 length.out = 301)
  vector2 <- seq(from = range2[1], 
                 to = range2[2],
                 length.out = 301)
  
  # Generate a formula by pasting column names for both dependent and independent variables
  # with interactions
  f <- as.formula(paste('cbind(',paste(colnames(Sim)[1:dep],collapse=','), ') ~ (',
                        'poly(',parm1, ',', poly.order,')*',
                        'poly(',parm2, ',', poly.order,')+',
                        paste(colnames(x)[c(-x1,-x2)], collapse='+'),')'))
  # Run Multiple Multivariate Regression (MMR) Metamodel
  Tway.mlm <- lm(f, data=Sim)
  
  # Create data frame with all combinations between both parameters of interest
  TWSA <- expand.grid(parm1 = vector1, 
                      parm2 = vector2)
  
  #Generate matrix to use for prediction 
  Sim.fit <- matrix(rep(colMeans(x)), 
                    nrow = nrow(TWSA),
                    ncol = ncol(x), byrow = T)
  Sim.fit[, x1] <- TWSA[, 1]
  Sim.fit[, x2] <- TWSA[, 2]
  # Transform to data frame, the format required for predict
  Sim.fit <- data.frame(Sim.fit) 
  # Name data frame's columns with parameters' names
  colnames(Sim.fit) <- colnames(x)
  
  # Predict Outcomes using MMMR Metamodel fit
  Sim.TW <- data.frame(predict(Tway.mlm, newdata = Sim.fit))
  # Find optimal strategy in terms of maximum expected outcome
  if (maximize){
    Optimal <- max.col(Sim.TW)
  } else { # Find optimal strategy in terms of minimum expected outcome
    Optimal <- min.col(Sim.TW)
  }
  
  # Add a variable with optimal startegy as factor
  TWSA$Strategy <- factor(Optimal, labels = strategies)
  
  twsa <- ggplot(TWSA, aes(x = parm1, y = parm2))+ 
    geom_tile(aes(fill = Strategy)) +
    theme_bw() +
    ggtitle(expression(atop("Two-way sensitivity analysis", 
                            atop("Net Health Benefit")))) + 
    scale_fill_discrete("Strategy: ", l=50) +
    scale_x_continuous(breaks = number_ticks(6)) +
    scale_y_continuous(breaks = number_ticks(6)) +
    xlab(parm1)+
    ylab(parm2)+
    theme(legend.position = "bottom", 
          legend.title = element_text(size = txtsize),
          legend.key = element_rect(colour = "black"),
          legend.text = element_text(size = txtsize),
          title = element_text(face="bold", size=15),
          axis.title.x = element_text(face="bold", size=txtsize),
          axis.title.y = element_text(face="bold", size=txtsize),
          axis.text.y = element_text(size=txtsize),
          axis.text.x = element_text(size=txtsize))
  return(twsa)
}
##############################
#### Formatting Functions ####
##############################
betaPar <- function(m, s)  # extract the  parameters of a beta distribution from mean and st. deviation
{
  a <- m * ((m * (1 - m) / s ^ 2) - 1)
  b <- (1 - m) * ((m * (1 - m) / s ^ 2) - 1)
  list(a = a, b = b)
}

Srange <- function(low, high)  # estimate the standard error from the upper and lower 95% confidence interval
{
  s = (high - low) / (2 * qnorm(0.975))
  list(s = s)
}

gammaPar <- function(mu, sigma) {   
  # Extract the parameters of a gamma distribution from mean and st. deviation 
  # mu: mean  
  # sigma: standard deviation 
  shape <- mu ^ 2 / sigma ^ 2
  scale <- sigma ^ 2 / mu
  list(shape = shape, scale = scale)
}

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

number_ticks <- function(n) {function(limits) pretty(limits, n)} #Function for number of ticks in ggplot

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}