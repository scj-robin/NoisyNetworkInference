# Attempts to fit np-noisySBM network inference to simulated scores.

rm(list=ls()); par(mfrow=c(1, 1), mex=.5, pch=20)
library(mvtnorm); library(sna); library(rggm); library(ks)
library(blockmodels)
library(MASS); library(mclust);library(pbmcapply)
# library(NoisySBM)
source('Functions/funcSimuls.R')
source('../../NoisySBM/R/funcVEM.R')
source('../../NoisySBM/R/tools.R')
source('../../NoisySBM/R/estimateNoisySBM.R')
source('../../NoisySBM/R/initInferenceNoisySBM.R')

################################################################################
# Data
nbNodes <- 80; nbBlocks <- 3; directed <- FALSE; nbObs <- 500; seed <- 2
load(paste0('../Sim/SimulScoreSR-n', nbNodes, '-K', nbBlocks, '-seed', seed, '.Rdata'))
networkVec <- mat2Vect(netSim$network, symmetric=!directed)

################################################################################
# Parms for np-inference
trans <- 2
scoreList <- transList[[trans]];
scoreMat <- scoreList2scoreMat(scoreList, symmetric=TRUE)
nbScores <- ncol(scoreMat); nbEdges <- nrow(scoreMat)

# Distribution
par(mfrow=c(2, nbScores))
for(s in 1:length(scoreList)){
  plot(density(scoreMat[, s]), main='')
  lines(density(scoreMat[which(networkVec==0), s]), col=4)
  lines(density(scoreMat[which(networkVec==1), s]), col=2)
}
plot(as.data.frame(scoreMat), col=1+networkVec)

# Kernel width + symmetric Gram matrix
kerSigma <- Hpi(scoreMat)
gram <- sapply(1:nbEdges, function(ij){dmvnorm(scoreMat, mean=scoreMat[ij, ], sigma=kerSigma)})
gramMarg <- list()
for(s in 1:nbScores){
  gramMarg[[s]] <- sapply(1:nbEdges, function(ij){dnorm(scoreMat[, s], mean=scoreMat[ij, s], sd=sqrt(kerSigma[s, s]))})
  }

################################################################################
# Fit np - Noisy SBM

estimateNoisySBM(transList[[trans]], directed=FALSE, nparm=TRUE)

# Init noisySBM
init <- initInferenceNoisySBM(transList[[trans]], directed=directed)

kerSigma <- Hpi(scoreList2scoreMat(transList[[trans]], symmetric=!directed))
gram <- sapply(1:nrow(scoreMat), function(ij){dmvnorm(scoreMat, mean=scoreMat[ij, ], sigma=kerSigma)})
vem <- VEMNoisySBM(scoreMat=scoreList2scoreMat(transList[[trans]], symmetric=!directed), directed=directed, 
                   qDistInit=list(psi = init$psi, tau=init$tau[[3]], eta=init$eta[[3]]), nparm=TRUE)
# np-Phi
prop <- colMeans(init$psi)
phi <- gram%*%init$psi %*% diag(prop)
colScale <- ceiling(10*(phi-min(phi)) / (max(phi)-min(phi)))

# Fit
par(mfrow=c(2, 2), mex=.5)
plot(scoreMat, col=colScale[, 1]); plot(scoreMat, col=colScale[, 2])
for(s in 1:length(scoreList)){
  plot(density(scoreMat[, s]), main='')
  points(scoreMat[, s], colMeans(gramMarg[[s]]), col=1, pch='.')
  for(g in c(0, 1)){
    dens <- density(scoreMat[which(networkVec==g), s])
    lines(dens$x, prop[g+1]*dens$y, col=2*(1+g))
    points(scoreMat[, s], gramMarg[[s]]%*%init$psi[, 1+g] / nbEdges, col=2*(1+g), pch='.')
  }
}

