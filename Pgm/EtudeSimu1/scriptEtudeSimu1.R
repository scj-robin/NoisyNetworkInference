library(NoisySBM)

#----------------------------------------------
#---------------------------- SCENARIO 1. 
#-----------------------------------------------
nbNodes  <- 100
directed <- FALSE
blockProp <- c(1/3,1/2,1/6)
nbBlocks   <- length(blockProp)
connectParam <- matrix(0.1,nbBlocks,nbBlocks)
diag(connectParam) <-  c(0.8,0.5,0.3)
connectParam <- 0.5*(connectParam + t(connectParam))
emissionParam <- list()
nbScores <- 4
emissionParam$noEdgeParam <- list(mean=c(0.1,0.5,0,1));
emissionParam$noEdgeParam$var <- matrix(0.1,nrow = nbScores,ncol = nbScores)
diag(emissionParam$noEdgeParam$var ) <- rep(1,nbScores)
emissionParam$edgeParam <- list( mean = 1:nbScores)
emissionParam$edgeParam$var <-  diag(0.1,nrow = nbScores,ncol = nbScores)
diag(emissionParam$noEdgeParam$var ) <- rep(1,nbScores)
paramSim <- list(emissionParam = emissionParam, connectParam = connectParam, blockProp = blockProp)

#------ Simu 
for (s in 1:100) {
  dataSim <- rNoisySBM(nbNodes,directed = FALSE, blockProp,connectParam,emissionParam,seed = NULL)
  saveFile <- paste('EtudeSimu1/dataSim/scenario1/data',s,'.Rdata',sep = '')
  save(dataSim,paramSim, file = saveFile)
}


#------ Estim 
rm(list = ls())
library(NoisySBM)
for (s in 1:50) {
  dataFile <- paste('EtudeSimu1/dataSim/scenario1/data',s,'.Rdata',sep='')
  load(dataFile)
  scoreList <- dataSim$noisyNetworks
  rm(paramSim)
  rm(dataSim)
  resEstim <- estimateNoisySBM(scoreList,directed = FALSE, estimOptions= list(verbosity = 0), monitoring =  list(lowerBound = TRUE))
  resultFile <- paste('EtudeSimu1/resSim/scenario1/res',s,'.Rdata',sep='')
  save(resEstim, file = resultFile)
}




#----------------------------------------------
#---------------------------- SCENARIO 2. 
#-----------------------------------------------
nbNodes  <- 60
directed <- FALSE
blockProp <- c(1/3,1/2,1/6)
nbBlocks   <- length(blockProp)
connectParam <- matrix(0.1,nbBlocks,nbBlocks)
connectParam[1,] <-  0.7
connectParam[,1] <-  0.6
connectParam[2,2] <-  0.4
connectParam[2,3] <-  0.3


connectParam <- 0.5 * (connectParam + t(connectParam))
emissionParam <- list()
nbScores <- 4
emissionParam$noEdgeParam <- list(mean = c(0.1,0.5,0,1));
emissionParam$noEdgeParam$var <- matrix(0.1,nrow = nbScores,ncol = nbScores)
diag(emissionParam$noEdgeParam$var ) <- rep(1,nbScores)
emissionParam$edgeParam <- list( mean = 1:nbScores)
emissionParam$edgeParam$var <-  diag(0.1,nrow = nbScores,ncol = nbScores)
diag(emissionParam$noEdgeParam$var ) <- rep(1,nbScores)
paramSim <- list(emissionParam = emissionParam, connectParam = connectParam, blockProp = blockProp)

#------ Simu 
for (s in 1:50) {
  dataSim <- rNoisySBM(nbNodes,directed = FALSE, blockProp,connectParam,emissionParam,seed = NULL)
  saveFile <- paste('EtudeSimu1/dataSim/scenario2/data',s,'.Rdata',sep = '')
  save(dataSim,paramSim, file = saveFile)
}


#------ Estim 
rm(list = ls())
library(NoisySBM)
for (s in 1:50) {
  dataFile <- paste('EtudeSimu1/dataSim/scenario2/data',s,'.Rdata',sep='')
  load(dataFile)
  scoreList <- dataSim$noisyNetworks
  rm(paramSim)
  rm(dataSim)
  resEstim <- estimateNoisySBM(scoreList,directed = FALSE, estimOptions= list(verbosity = 0), monitoring =  list(lowerBound = TRUE))
  resultFile <- paste('EtudeSimu1/resSim/scenario2/res',s,'.Rdata',sep = '')
  save(resEstim, file = resultFile)
}

#----------------------------------------------
#---------------------------- SCENARIO 3. 
#-----------------------------------------------
nbNodes  <- 100
directed <- FALSE
blockProp <- c(1/3,1/2,1/6)
nbBlocks   <- length(blockProp)
connectParam <- matrix(0.1,nbBlocks,nbBlocks)
diag(connectParam) <-  c(0.8,0.5,0.3)
connectParam <- 0.5*(connectParam + t(connectParam))
emissionParam <- list()
nbScores <- 4
emissionParam$noEdgeParam <- list(mean=c(0.1,0.5,0,1));
emissionParam$noEdgeParam$var <- matrix(0.1,nrow = nbScores,ncol = nbScores)
diag(emissionParam$noEdgeParam$var ) <- rep(1,nbScores)
emissionParam$edgeParam <- list( mean = c(1,2,0.9,2))
emissionParam$edgeParam$var <-  diag(0.1,nrow = nbScores,ncol = nbScores)
diag(emissionParam$noEdgeParam$var ) <- rep(1,nbScores)
paramSim <- list(emissionParam = emissionParam, connectParam = connectParam, blockProp = blockProp)

#------ Simu 
for (s in 1:50) {
  dataSim <- rNoisySBM(nbNodes,directed = FALSE, blockProp,connectParam,emissionParam,seed = NULL)
  saveFile <- paste('EtudeSimu1/dataSim/scenario3/data',s,'.Rdata',sep = '')
  save(dataSim,paramSim, file = saveFile)
}


#------ Estim 
rm(list = ls())
library(NoisySBM)
for (s in 1:50) {
  dataFile <- paste('EtudeSimu1/dataSim/scenario3/data',s,'.Rdata',sep='')
  load(dataFile)
  scoreList <- dataSim$noisyNetworks
  rm(paramSim)
  rm(dataSim)
  resEstim <- estimateNoisySBM(scoreList,directed = FALSE, estimOptions= list(verbosity = 0), monitoring =  list(lowerBound = TRUE))
  resultFile <- paste('EtudeSimu1/resSim/scenario3/res',s,'.Rdata',sep='')
  save(resEstim, file = resultFile)
}



