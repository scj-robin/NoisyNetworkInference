rm(list=ls())
library(NoisySBM)
 





  
  ####################################################################### 

#----------------------------------------------------------------------------------
#---------------------------- Param SCENARIO 1 : communauté et scores très contrastés. 
#---------------------------------------------------------------------------------
nbNodes  <- 80 
blockProp <- c(1/3,1/2,1/6)
nbBlocks <- length(blockProp)
connectParam <- matrix(0.1,nbBlocks,nbBlocks)
diag(connectParam) <-  c(0.8,0.5,0.3)
connectParam <- 0.5*(connectParam + t(connectParam))
density = 0.3
connectParam <- connectParam /  c(blockProp %*% connectParam %*% blockProp)  *  density
emissionParam <- list()
nbScores <- 4
emissionParam$noEdgeParam <- list(mean = c(0.1,0.5,0,1));
emissionParam$noEdgeParam$var <- matrix(0.1,nrow = nbScores,ncol = nbScores)
diag(emissionParam$noEdgeParam$var ) <- rep(1,nbScores)
emissionParam$edgeParam <- list( mean = (1:nbScores))
emissionParam$edgeParam$var <-  diag(0.1,nrow = nbScores,ncol = nbScores)
diag(emissionParam$noEdgeParam$var ) <- rep(1,nbScores)
paramSim <- list(density = density, nbNodes  = nbNodes, emissionParam = emissionParam, connectParam = connectParam, blockProp = blockProp)
save(paramSim, file ='EtudeSimu1/dataSim/scenario1/parSim.Rdata')

#----------------------------------------------------------------------------------
#---------------------------- Param SCENARIO 2 : nested et scores très contrastés. 
#---------------------------------------------------------------------------------
nbBlocks   <- length(paramSim$blockProp)
connectParamNested <- matrix(0.1,nbBlocks,nbBlocks)
connectParamNested[1,] <-  0.7
connectParamNested[,1] <-  0.6
connectParamNested[2,2] <-  0.4
connectParamNested[2,3] <-  0.3
connectParamNested <- 0.5 * (connectParamNested+ t(connectParamNested))
connectParamNested <- connectParamNested /  c(paramSim$blockProp %*% connectParamNested %*% paramSim$blockProp)  *  density
paramSim$connectParam <- connectParamNested
save(paramSim, file ='EtudeSimu1/dataSim/scenario2/parSim.Rdata')


#----------------------------------------------------------------------------------
#---------------------------- Param SCENARIO 3 : communauté et scores peu contrastés. 
#---------------------------------------------------------------------------------
load(file ='EtudeSimu1/dataSim/scenario1/parSim.Rdata')
emissionParam <- paramSim$emissionParam
emissionParam$edgeParam$mean =  c(1,2,0.9,2)

paramSim$emissionParam <- emissionParam
save(paramSim, file ='EtudeSimu1/dataSim/scenario3/parSim.Rdata')

#----------------------------------------------------------------------------------
#---------------------------- Param SCENARIO 4 : communauté et scores peu contrastés. 
#---------------------------------------------------------------------------------
load(file ='EtudeSimu1/dataSim/scenario3/parSim.Rdata')
paramSim$connectParam <- connectParamNested
save(paramSim, file ='EtudeSimu1/dataSim/scenario4/parSim.Rdata')



################################################################""
#---------------------------- SIMULATION
################################################################""
rm(list=ls())
library(NoisySBM)
 
fSimul <- function(indSim,scenario){
  set.seed(100)
  load(file = paste('EtudeSimu1/dataSim/scenario',scenario,'/parSim.Rdata',sep = ''))
  for (s in indSim){
    dataSim <- rNoisySBM(paramSim$nbNodes,directed = FALSE, paramSim$blockProp,paramSim$connectParam,paramSim$emissionParam,seed = NULL)
    saveFile <- paste('EtudeSimu1/dataSim/scenario',scenario,'/data',s,'.Rdata',sep = '')
    save(dataSim,paramSim, file = saveFile)
  }
}
for (scenario in 1:4){fSimul(1:50,scenario)}


################################################################""
#---------------------------- Estimation
################################################################""
#------ Estim 
rm(list = ls())
library(NoisySBM)
library(mclust)
# ###################################################################### 
fEstim <- function(indSim,scenario){
  for (s in indSim) {
    dataFile <- paste('EtudeSimu1/dataSim/scenario',scenario,'/data',s,'.Rdata',sep='')
    load(dataFile)
    scoreList <- dataSim$noisyNetworks
    rm(paramSim)
    rm(dataSim)
    resEstim <- estimateNoisySBM(scoreList,directed = FALSE, estimOptions= list(verbosity = 0), monitoring =  list(lowerBound = TRUE))
    resultFile <- paste('EtudeSimu1/resSim/scenario',scenario,'/res',s,'.Rdata',sep='')
    save(resEstim, file = resultFile)
  }
}

for (scenario in 1:4){fEstim(1:50,scenario)}

