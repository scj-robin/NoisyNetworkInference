# Attempts to fit np-noisySBM network inference to simulated scores.

rm(list=ls()); par(mfrow=c(1, 1), mex=.5, pch=20)
library(mvtnorm); library(sna); library(rggm); library(ks)
library(blockmodels)
library(MASS); library(mclust);library(pbmcapply)
# library(NoisySBM)
source('Functions/funcSimuls.R')

# source('D:/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/NOISY_SBM_Package/Package/NoisySBM/R/funcVEM.R')
# source('D:/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/NOISY_SBM_Package/Package/NoisySBM/R/tools.R')
# source('D:/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/NOISY_SBM_Package/Package/NoisySBM/R/estimateNoisySBM.R')
# source('D:/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/NOISY_SBM_Package/Package/NoisySBM/R/initInferenceNoisySBM.R')

source('../../NoisySBM/R/funcVEM.R')
source('../../NoisySBM/R/tools.R')
# source('../../NoisySBM/R/estimateNoisySBM.R')
# source('../../NoisySBM/R/initInferenceNoisySBM.R')
source('../../NoisySBM/R/estimateScoreSBM.R')
source('../../NoisySBM/R/initInferenceScoreSBM.R')

################################################################################
# Data
nbNodes <- 80; nbBlocks <- 3; directed <- FALSE; nbObs <- 100; nbReplicates <- 5; seed <- 2
nbNodes <- 80; nbBlocks <- 3; directed <- FALSE; nbObs <- 250; nbReplicates <- 2
simFile <- paste0('../Sim/SimulScoreSR-n', nbNodes, '-nbObs', nbObs, '-K', nbBlocks, '-rep', nbReplicates, '-seed', seed, '.Rdata')
load(file=simFile)
networkVec <- mat2Vect(netSim$network, symmetric=!directed)

################################################################################
# Parms for np-inference
trans <- 3
scoreRepList <- list(); scoreRepNum <- 0
for(rep in 1:nbReplicates){ # rep <- 1
   for(sc in 1:length(transList[[1]][[trans]])){ 
      scoreRepNum <- scoreRepNum +1
      scoreRepList[[scoreRepNum]] <- transList[[rep]][[trans]][[sc]]
   }
}
length(scoreRepList)
scoreRepMat <- scoreList2scoreMat(scoreRepList, symmetric=TRUE)
nbRepScores <- ncol(scoreRepMat); nbEdges <- nrow(scoreRepMat)

# Distribution
par(mfrow=c(ceiling(sqrt(nbRepScores)), round(sqrt(nbRepScores))))
for(s in 1:length(scoreRepList)){
  plot(density(scoreRepMat[, s]), main='')
  lines(density(scoreRepMat[which(networkVec==0), s]), col=4)
  lines(density(scoreRepMat[which(networkVec==1), s]), col=2)
}
# plot(as.data.frame(scoreMat), col=1+networkVec)

# Kernel width + symmetric Gram matrix
# kerSigma <- Hpi(scoreRepMat)
kerSigma <- diag(apply(scoreRepMat, 2, hpi))
gram <- sapply(1:nbEdges, function(ij){dmvnorm(scoreRepMat, mean=scoreRepMat[ij, ], sigma=kerSigma)})
gramMarg <- list()
for(s in 1:nbRepScores){
  gramMarg[[s]] <- sapply(1:nbEdges, function(ij){dnorm(scoreRepMat[, s], mean=scoreRepMat[ij, s], sd=sqrt(kerSigma[s, s]))})
  }

################################################################################
# Fit Noisy SBM

init <- initInferenceScoreSBM(scoreRepList, directed=directed)
initTau <- init$tau[[nbBlocks]]
fit <- estimateScoreSBM(scoreRepList, directed=FALSE, nparm=FALSE)
fitTau <- fit[[which(sapply(1:length(fit), function(dim){ncol(fit[[dim]]$qDist$tau)})==nbBlocks)]]$qDist$tau
npFit <- estimateScoreSBM(scoreRepList, directed=FALSE, nparm=TRUE, kerSigma=kerSigma)
npFitTau <- npFit[[which(sapply(1:length(fit), function(dim){ncol(npFit[[dim]]$qDist$tau)})==nbBlocks)]]$qDist$tau
round(t(netSim$membership)%*%initTau)
round(t(netSim$membership)%*%fitTau)
round(t(netSim$membership)%*%npFitTau)

# np-Phi
prop <- colMeans(init$psi)
phi <- gram%*%init$psi %*% diag(prop)
colScale <- ceiling(10*(phi-min(phi)) / (max(phi)-min(phi)))

# Fit
par(mfrow=c(2, 2), mex=.5)
plot(scoreRepMat, col=colScale[, 1]); plot(scoreRepMat, col=colScale[, 2])
for(s in 1:length(scoreRepList)){
  plot(density(scoreRepMat[, s]), main='')
  points(scoreRepMat[, s], colMeans(gramMarg[[s]]), col=1, pch='.')
  for(g in c(0, 1)){
    dens <- density(scoreRepMat[which(networkVec==g), s])
    lines(dens$x, prop[g+1]*dens$y, col=2*(1+g))
    points(scoreRepMat[, s], gramMarg[[s]]%*%init$psi[, 1+g] / nbEdges, col=2*(1+g), pch='.')
  }
}

