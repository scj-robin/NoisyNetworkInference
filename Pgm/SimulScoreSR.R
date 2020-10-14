# Simulates data, fit network inference, get scores

# devtools::install_github("jchiquet/rggm", build_vignettes=TRUE)
# dev.off(); 
rm(list=ls()); par(mfrow=c(1, 1), mex=.5, pch=20)
seed <- 2; set.seed(seed)
library(gtools); library(mvtnorm); library(sna); library(rggm)
library(blockmodels)
library(MASS); library(mclust); library(PLNmodels)
library(huge); library(saturnin); library(GeneNet); library(BGGM)
source('../../NoisySBM/R/tools.R')
source('Functions/funcSimuls.R')

# Dims
nbNodes <- 80; nbBlocks <- 3; directed <- FALSE; nbObs <- 500
netDensity <- sqrt(nbNodes)/nbNodes;

# Parms
# blockProp <- rdirichlet(1, rep(1, nbBlocks))
# connectParam <- matrix(rbeta(nbBlocks^2, 1, 1), nbBlocks, nbBlocks)
blockProp <- (1:nbBlocks); blockProp <- blockProp/sum(blockProp)
connectParam <- (nbBlocks:1)^2%o%(nbBlocks:1)^2
connectParam <- .2*connectParam + .8*diag(nbBlocks)
connectParam <- connectParam * (netDensity / (blockProp%*%connectParam%*%blockProp)[1, 1])
connectParam <- vect2Mat(mat2Vect(connectParam, symmetric=!directed, diag=TRUE), symmetric=!directed, diag=TRUE)

# Network simulation
netSim <- rSimSBM(nbNodes=nbNodes, directed=directed, blockProp=blockProp, connectParam=connectParam)
# save(netSim, file=paste0('../Sim/SimulScoreSR-n', nbNodes, '-K', nbBlocks, '-seed', seed, '.Rdata'))
load(paste0('../Sim/SimulScoreSR-n', nbNodes, '-K', nbBlocks, '-seed', seed, '.Rdata'))
while(min(colMeans(netSim$network))==0){
   netSim <- rSimSBM(nbNodes=nbNodes, directed=directed, blockProp=blockProp, connectParam=connectParam)
   }
sna::gplot(netSim$network, gmode='graph', vertex.col=1+netSim$membership%*%(1:nbBlocks), edge.col=8)
networkVec <- mat2Vect(netSim$network, symmetric=TRUE)
covStruc <- makeVariance(network=netSim$network)

# Blockmodel on edges
LBM <- BM_bernoulli(membership_type='SBM_sym', adj=netSim$network)
LBM$estimate()
Klbm <- which.max(LBM$ICL)
Zlbm <- LBM$memberships[[Klbm]]$Z
print(round(t(netSim$membership)%*%Zlbm))

# Data simulation: Gaussian
data <- rmvnorm(nbObs, sigma=covStruc$sigma)

# Scores
scoreList <- list()
scoreList$huge <- getScoreHuge(data)
scoreList$saturnin <- saturnin::edge.prob(20/nbObs*saturnin::lweights_gaussian(data), log=TRUE)
# scoreList$gnPcor <- GeneNet::ggm.estimate.pcor(data)^2
# scoreList$bGGM <- (BGGM::estimate(data)$parcors_mat)^2
# scoreList$gnPval <- vect2Mat(1-GeneNet::network.test.edges(GeneNet::ggm.estimate.pcor(data), verbose=FALSE)$pval, symmetric=TRUE)

# # Data simulation: PLN
# Z <- rmvnorm(nbObs, sigma=covStruc$sigma)
# mu=rep(5, nbNodes)
# dataPLN <- matrix(rpois(nbObs*nbNodes, exp(rep(1, nbObs)%o%mu + Z)), nbObs, nbNodes)

# # PLN scores
# pln <- PLN(dataPLN ~ 1)
# scoreList <- list()
# scoreList$plnNet <- getScorePLNnet(dataPLN, pln)
# scoreList$emTree <- EMtree::EMtree(PLNobject=pln)$edges_prob
# scoreList$plnPcor <- cov2cor(solve(pln$model_par$Sigma))^2

# Score matrices
par(mfrow=c(2, 2))
lapply(scoreList, function(score){
  image(1:nbNodes, 1:nbNodes, score[order(netSim$cluster), order(netSim$cluster)])
  abline(v=.5+cumsum(table(netSim$cluster)), h=.5+cumsum(table(netSim$cluster)))
})
save(netSim, scoreList, file=paste0('../Sim/SimulScoreSR-n', nbNodes, '-nbOds', nbObs, '-K', nbBlocks, '-seed', seed, '.Rdata'))

scoreMat <- scoreList2scoreMat(scoreList, symmetric=!directed)

# Transformations
par(mfcol=c(length(scoreList), 2))
transList <- list()
# transList$raw <- scoreList
# transList$boxCox <- lapply(scoreList, function(score){boxCoxTransform(score, plotit=TRUE)})
# transList$gmmBC <- lapply(scoreList, function(score){boxCoxGmmTransform(score, plotit=TRUE)})
# transList$adhoc <- list(huge=log(scoreList$huge), saturnin=qnorm(scoreList$saturnin),
                        # gnPcor=atanh(GeneNet::ggm.estimate.pcor(data))^2)
# transList$adhoc <- list(huge=log(scoreList$huge), saturnin=qnorm(scoreList$saturnin),
                        # gnPcor=atanh(GeneNet::ggm.estimate.pcor(data))^2, gnPval=qnorm(scoreList$gnPval))
# transList$adhoc <- list(huge=log(scoreList$huge), saturnin=qnorm(scoreList$saturnin))
# transList$adhoc <- list(plnNet=log(scoreList$plnNet), emTree=qnorm(scoreList$emTree),
#                        plnPcor=atanh(cov2cor(solve(pln$model_par$Sigma)))^2)
transList$unif <- lapply(scoreList, function(score){inverseTransform(score)})
transList$gauss <- lapply(scoreList, function(score){qnorm(inverseTransform(score))})
# scoreMat <- scoreList2scoreMat(scoreList, symmetric=!directed)
# transList$gmm <- list()
# for(s in 1:ncol(scoreMat)){
#    scoreVec <- scoreMat[, s]
#    GMM <- Mclust(scoreVec, G=2)
#    qTauVec <- qnorm(GMM$z[, which.max(GMM$parameters$mean)])
#    transList$gmm[[s]] <- vect2Mat(qTauVec, symmetric=!directed)
# } 

# Remove diags
par(mfrow=c(length(transList), length(scoreList)))
sapply(1:length(transList), function(t){sapply(1:length(transList[[t]]), function(s){
   diag(transList[[t]][[s]]) <- NA
   image(1:nbNodes, 1:nbNodes, transList[[t]][[s]][order(netSim$cluster), order(netSim$cluster)], 
         main=paste(names(transList)[t], names(scoreList)[s]))
   abline(v=.5+cumsum(table(netSim$cluster)), h=.5+cumsum(table(netSim$cluster)))
   })})
save(netSim, scoreList, transList, file=paste0('../Sim/SimulScoreSR-n', nbNodes, '-K', nbBlocks, '-seed', seed, '.Rdata'))

# Histograms
par(mfrow=c(length(scoreList), length(transList)))
sapply(1:length(scoreList), function(s){
  sapply(1:length(transList), function(t){
    hist(transList[[t]][[s]], breaks=nbNodes, xlab='', ylab='',
         main=paste(names(transList)[t], names(transList[[t]])[s]))
  })
})

# Box-plots
par(mfrow=c(length(scoreList), length(transList)))
for(s in 1:length(scoreList)){for(t in 1:length(transList)){
  cat(s, t, '\n')
  boxplot(mat2Vect(transList[[t]][[s]], symmetric=TRUE) ~ networkVec,
          xlab='', ylab='', main=paste(names(transList)[t], names(transList[[t]])[s]))
}}

# Densities
par(mfrow=c(length(scoreList), length(transList)))
for(s in 1:length(scoreList)){for(t in 1:length(transList)){
   cat(s, t, '\n')
   plot(density(mat2Vect(transList[[t]][[s]], symmetric=TRUE)[which(networkVec==0)]))
   lines(density(mat2Vect(transList[[t]][[s]], symmetric=TRUE)[which(networkVec==1)]), col=2)
}}

# Matrix form
transMatList <- lapply(transList, function(trans){scoreList2scoreMat(trans, symmetric=TRUE)})
lapply(transMatList, function(mat){plot(as.data.frame(mat), col=1+networkVec)})

# PCA
pcaList <- lapply(transMatList, function(mat){prcomp(mat)$x})
par(mfrow=c(2, 2))
sapply(1:length(pcaList), function(m){
   plot(pcaList[[m]][, 1:2], pch=20, col=1+networkVec, main=names(pcaList)[m])
})
# transList$gaussPC <- list(PC1=vect2Mat(pcaList$gauss[, 1], symmetric=TRUE),
#                           PC2=vect2Mat(pcaList$gauss[, 2], symmetric=TRUE))
# transMatList$gaussPC <- scoreList2scoreMat(transList$gaussPC, symmetric=TRUE)

# LDA / Manova
ldaList <- lapply(transMatList, function(mat){lda(networkVec ~ mat)})
par(mfrow=c(2, 2))
sapply(1:length(transMatList), function(t){
   boxplot(transMatList[[t]]%*%ldaList[[t]]$scaling ~ networkVec, main=names(transList)[t])
})
lapply(transMatList, function(mat){summary(manova(lm(mat ~ networkVec)))})

# # MclustDA
# mdaList <- lapply(transMatList, function(mat){MclustDA(mat, networkVec)})
# lapply(mdaList, function(mda){plot(mda, what='classification')})

# Compare transformation: GMM
gmmList <- lapply(transMatList, function(mat){Mclust(mat, G=2)})
lapply(gmmList, function(gmm){plot(gmm, what='classification')})
lapply(gmmList, function(gmm){print(round(t(gmm$z) %*% cbind(1-networkVec, networkVec)))})

################################################################################
# Fit Noisy SBM
# library(NoisySBM)
library(pbmcapply)
source('../../NoisySBM/R/funcVEM.R')
source('../../NoisySBM/R/initInferenceNoisySBM.R')
# source('../../NoisySBM/R/tools.R')
source('../../NoisySBM/R/estimateNoisySBM.R')

# Fit noisySBM
par(mfrow=c(2, 2)); fitList <- list()
for(t in 1:length(transList)){ # t <- 5
  nSBM <- estimateNoisySBM(scoreList=transList[[t]], directed = FALSE,
                           estimOptions = list(explorFactor=2.5, etaTol=1e-4, tauTol=1e-4, maxIterVE=1), 
                           monitoring = list())
  KList <- unlist(lapply(nSBM, function(fit){ncol(fit$qDist$tau)}))
  J <- unlist(lapply(nSBM, function(fit){max(fit$lowerBound)}))
  ICL <- unlist(lapply(nSBM, function(fit){fit$ICL}))
  pen <- unlist(lapply(nSBM, function(fit){fit$pen}))
  Jpen <- J+pen
  fitList[[t]] <- list(nSBM=nSBM, J=J, ICL=ICL, pen=pen, Jpen=Jpen, KList=KList)
}
save(netSim, scoreList, transList, fitList, file=paste0('../Sim/SimulScoreSR-n', nbNodes, '-K', nbBlocks, '-seed', seed, '.Rdata'))

# Results
par(mfrow=c(length(transList), 1))
for(t in 1:length(transList)){ # t <- 5
  fit <- fitList[[t]]
  # plot(fit$KList, fit$J, ylim=c(min(c(fit$ICL, fit$Jpen)), max(fit$J)), main=names(transList)[t])
  # points(fit$KList, fit$ICL, col=2); points(fit$KList, fit$Jpen, col='4');
  # abline(v=fit$KList[which.max(fit$Jpen)], col=4, lwd=2)
  optK <- which.max(fit$Jpen); # optK <- which(fit$KList==nbBlocks)
  plot(fit$nSBM[[optK]]$lowerBound, col=c(1, 2))
  print(round(t(netSim$membership)%*%fit$nSBM[[optK]]$qDist$tau))
}
print(round(t(netSim$membership)%*%Zlbm))
