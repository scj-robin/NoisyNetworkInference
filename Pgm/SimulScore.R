# Simulates data, fit network inference, get scores

rm(list=ls()); par(mfrow=c(1, 1), pch=20)
library(gtools); library(sna); library(MASS)
library(mvtnorm); library(mclust); library(PLNmodels)
library(huge); library(saturnin); library(GeneNet)
source('../../NoisySBM/R/tools.R')
source('Functions/funcSimuls.R')

# Dims
nbNodes <- 50; nbBlocks <- 3; directed <- TRUE; nbObs <- 100
netDensity <- sqrt(nbNodes)/nbNodes; 

# Parms
# blockProp <- rdirichlet(1, rep(1, nbBlocks))
# connectParam <- matrix(rbeta(nbBlocks^2, 1, 1), nbBlocks, nbBlocks)
blockProp <- (1:nbBlocks); blockProp <- blockProp/sum(blockProp) 
connectParam <- (nbBlocks:1)%o%(nbBlocks:1)
connectParam <- connectParam * (netDensity / (blockProp%*%connectParam%*%blockProp)[1, 1])
connectParam <- vect2Mat(mat2Vect(connectParam, symmetric=!directed, diag=TRUE), symmetric=!directed, diag=TRUE)

# Network simulation
netSim <- rSimSBM(nbNodes=nbNodes, directed=directed, blockProp=blockProp, connectParam=connectParam)
sna::gplot(netSim$network, gmode='graph', vertex.col=netSim$membership%*%(1:nbBlocks))
networkVec <- mat2Vect(netSim$network, symmetric=TRUE)
covStruc <- makeVariance(network=netSim$network)

# Data simulation: Gaussian
data <- rmvnorm(nbObs, sigma=covStruc$sigma)

# Gaussian scores
scoreList <- list()
scoreList$huge <- getScoreHuge(data)
scoreList$saturnin <- saturnin::edge.prob(saturnin::lweights_gaussian(data))
scoreList$gnPcor <- GeneNet::ggm.estimate.pcor(data)^2
scoreList$gnPval <- vect2Mat(1-GeneNet::network.test.edges(GeneNet::ggm.estimate.pcor(data), verbose=FALSE)$pval, symmetric=TRUE)
scoreAdHocList <- list(huge=log(scoreList$huge), saturnin=qnorm(scoreList$saturnin),
                       gnPcor=atanh(GeneNet::ggm.estimate.pcor(data))^2, gnPval=qnorm(scoreList$gnPval))

# # Data simulation: PLN
# Z <- rmvnorm(nbObs, sigma=covStruc$sigma)
# mu=rep(5, nbNodes)
# dataPLN <- matrix(rpois(nbObs*nbNodes, exp(rep(1, nbObs)%o%mu + Z)), nbObs, nbNodes)
# 
# # PLN scores
# pln <- PLN(dataPLN ~ 1)
# scoreList <- list()
# scoreList$plnNet <- getScorePLNnet(dataPLN, pln)
# scoreList$emTree <- EMtree::EMtree(PLNobject=pln)$edges_prob
# scoreList$plnPcor <- cov2cor(solve(pln$model_par$Sigma))^2
# scoreAdHocList <- list(plnNet=log(scoreList$plnNet), emTree=qnorm(scoreList$emTree),
#                        plnPcor=atanh(cov2cor(solve(pln$model_par$Sigma)))^2)

# Box-Cox + GMM transformations 
par(mfcol=c(length(scoreList), 2))
scoreBCList <- lapply(scoreList, function(score){boxCoxTransform(score, plotit=TRUE)})
scoreGMMList <- lapply(scoreList, function(score){boxCoxGmmTransform(score, plotit=TRUE)})

# Histograms
par(mfrow=c(length(scoreList), 4))
sapply(1:length(scoreList), function(s){
   hist(scoreList[[s]], breaks=nbNodes, xlab=names(scoreList)[s], ylab='', main='raw')
   hist(scoreBCList[[s]], breaks=nbNodes, xlab=names(scoreList)[s], ylab='', main='BC')
   hist(scoreGMMList[[s]], breaks=nbNodes, xlab=names(scoreList)[s], ylab='', main='GMM')
   hist(scoreAdHocList[[s]], breaks=nbNodes, xlab=names(scoreList)[s], ylab='', main='ad-hoc')
})

# Box-plots
par(mfrow=c(length(scoreList), 4))
sapply(1:length(scoreList), function(s){
   boxplot(mat2Vect(scoreList[[s]], symmetric=TRUE) ~ networkVec, xlab='', ylab=names(scoreList)[s], main='raw')
   boxplot(mat2Vect(scoreBCList[[s]], symmetric=TRUE) ~ networkVec, xlab='', ylab=names(scoreList)[s], main='BC')
   boxplot(mat2Vect(scoreGMMList[[s]], symmetric=TRUE) ~ networkVec, xlab='', ylab=names(scoreList)[s], main='GMM')
   boxplot(mat2Vect(scoreAdHocList[[s]], symmetric=TRUE) ~ networkVec, xlab='', ylab=names(scoreList)[s], main='ad_hoc')
})

# PCA
par(mfrow=c(2, 2))
scoreMat <- scoreList2scoreMat(scoreList, symmetric=TRUE)
scoreBCMat <- scoreList2scoreMat(scoreBCList, symmetric=TRUE)
scoreGMMMat <- scoreList2scoreMat(scoreGMMList, symmetric=TRUE)
scoreAdHocMat <- scoreList2scoreMat(scoreAdHocList, symmetric=TRUE)

# Compare transformation: PCA
pcaList <- list(raw=prcomp(scoreMat)$x, 
                BC=prcomp(scoreBCMat)$x, 
                GMM=prcomp(scoreGMMMat)$x,
                AdHoc=prcomp(scoreAdHocMat)$x)
par(mfrow=c(2, 2))
sapply(1:length(pcaList), function(m){
   plot(pcaList[[m]][, 1:2], pch=20, col=1+networkVec, main=names(pcaList)[m])
})

# Compare transformation: GMM
gmmList <- lapply(pcaList, function(pca){Mclust(scale(pca), G=2)})

unlist(lapply(gmmList, function(gmm){gmm$loglik}))
lapply(gmmList, function(gmm){adjustedRandIndex(gmm$classification, networkVec)})

       