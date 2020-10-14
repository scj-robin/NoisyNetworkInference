# Simulates data, fit network inference, get scores

rm(list=ls()); par(mfrow=c(1, 1), mex=.5, pch=20)
seed <- 2; set.seed(seed)
library(gtools); library(mvtnorm); library(sna); library(rggm)
library(blockmodels)
library(MASS); library(mclust); library(PLNmodels)
library(huge); library(saturnin); library(GeneNet); library(BGGM)
source('../../NoisySBM/R/tools.R')
source('Functions/funcSimuls.R')

# Dims
nbNodes <- 80; nbBlocks <- 3; directed <- FALSE; nbObs <- 100; nbReplicates <- 5
nbNodes <- 80; nbBlocks <- 3; directed <- FALSE; nbObs <- 250; nbReplicates <- 2
netDensity <- sqrt(nbNodes)/nbNodes;
simFile <- paste0('../Sim/SimulScoreSR-n', nbNodes, '-nbObs', nbObs, '-K', nbBlocks, '-rep', nbReplicates, '-seed', seed, '.Rdata')

# Parms
blockProp <- (1:nbBlocks); blockProp <- blockProp/sum(blockProp)
connectParam <- (nbBlocks:1)^2%o%(nbBlocks:1)^2
connectParam <- .2*connectParam + .8*diag(nbBlocks)
connectParam <- connectParam * (netDensity / (blockProp%*%connectParam%*%blockProp)[1, 1])
connectParam <- vect2Mat(mat2Vect(connectParam, symmetric=!directed, diag=TRUE), symmetric=!directed, diag=TRUE)

# Network simulation
netSim <- rSimSBM(nbNodes=nbNodes, directed=directed, blockProp=blockProp, connectParam=connectParam)
while(min(colMeans(netSim$network))==0){
   netSim <- rSimSBM(nbNodes=nbNodes, directed=directed, blockProp=blockProp, connectParam=connectParam)
}
sna::gplot(netSim$network, gmode='graph', vertex.col=1+netSim$membership%*%(1:nbBlocks), edge.col=8)
networkVec <- mat2Vect(netSim$network, symmetric=TRUE)
covStruc <- makeVariance(network=netSim$network)
save(netSim, file=simFile)

# Data simulation: Gaussian
dataList <- list()
for(rep in 1:nbReplicates){
   dataList[[rep]] <- rmvnorm(nbObs, sigma=covStruc$sigma)
}
save(netSim, dataList, file=simFile)

# Scores
scoreList <- list()
for(rep in 1:nbReplicates){
   cat(rep, '')
   scoreList[[rep]] <- list()
   scoreList[[rep]]$huge <- getScoreHuge(dataList[[rep]])
   scoreList[[rep]]$saturnin <- saturnin::edge.prob(20/nbObs*saturnin::lweights_gaussian(dataList[[rep]]), log=TRUE)
   # scoreList[[rep]]$gnPcor <- GeneNet::ggm.estimate.pcor(dataList[[rep]])^2
   # scoreList[[rep]]$bGGM <- (BGGM::estimate(dataList[[rep]])$parcors_mat)^2
   # scoreList[[rep]]$gnPval <- vect2Mat(1-GeneNet::network.test.edges(GeneNet::ggm.estimate.pcor(dataList[[rep]]), verbose=FALSE)$pval, symmetric=TRUE)
}
save(netSim, dataList, scoreList, file=simFile)

# Transformations
transList <- list()
for(rep in 1:nbReplicates){ # rep <- 1
   cat(rep, '')
   transList[[rep]] <- list()
   transList[[rep]]$raw <- scoreList[[rep]]
   transList[[rep]]$boxCox <- lapply(scoreList[[rep]], function(score){boxCoxTransform(score, plotit=TRUE)})
   # transList[[rep]]$gmmBC <- lapply(scoreList[[rep]], function(score){boxCoxGmmTransform(score, plotit=TRUE)})
   # transList[[rep]]$adhoc <- list(huge=log(scoreList[[rep]]$huge), saturnin=qnorm(scoreList[[rep]]$saturnin),
   # gnPcor=atanh(GeneNet::ggm.estimate.pcor(data[[rep]]))^2)
   # transList[[rep]]$adhoc <- list(huge=log(scoreList[[rep]]$huge), saturnin=qnorm(scoreList[[rep]]$saturnin),
   # gnPcor=atanh(GeneNet::ggm.estimate.pcor(data))^2, gnPval=qnorm(scoreList[[rep]]$gnPval))
   transList[[rep]]$adhoc <- list(huge=log(scoreList[[rep]]$huge), saturnin=qnorm(scoreList[[rep]]$saturnin))
   # transList[[rep]]$adhoc <- list(plnNet=log(scoreList[[rep]]$plnNet), emTree=qnorm(scoreList[[rep]]$emTree),
                          # plnPcor=atanh(cov2cor(solve(pln$model_par$Sigma)))^2)
}
save(netSim, dataList, scoreList, transList, file=simFile)
# load(file=simFile)
