# Attempts to fit np-noisySBM network inference to simulated scores.

rm(list=ls()); par(mfrow=c(1, 1), mex=.5, pch=20)
library(mvtnorm); library(sna); library(rggm); library(ks)
library(blockmodels)
library(MASS); library(mclust);library(pbmcapply)
library(gtools)
# library(NoisySBM)
source('Functions/funcSimuls.R')

source('D:/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/NOISY_SBM_Package/Package/NoisySBM/R/funcVEM.R')
source('D:/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/NOISY_SBM_Package/Package/NoisySBM/R/tools.R')
source('D:/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/NOISY_SBM_Package/Package/NoisySBM/R/estimateNoisySBM.R')
source('D:/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/NOISY_SBM_Package/Package/NoisySBM/R/initInferenceNoisySBM.R')

source('../../NoisySBM/R/funcVEM.R')
source('../../NoisySBM/R/tools.R')
source('../../NoisySBM/R/estimateNoisySBM.R')
source('../../NoisySBM/R/initInferenceNoisySBM.R')

################################################################################
# Data
nbNodes <- 80; nbBlocks <- 3; directed <- FALSE; nbObs <- 500; seed <- 2



# Dims
blockProp <- c(1/3, 1/3, 1/3)  # group proportions
nbBlock   <- length(blockProp) # number of blocks
connectParam <- diag(.4, nbBlock) + 0.01 # connectivity matrix: affiliation network
mySBM <- rSBM(nbNodes, connectParam, blockProp)
library(igraph)
netSim <- as_adjacency_matrix(mySBM)


clusterTrue <- vertex_attr(mySBM, "memberships")

# Network simulation

sna::gplot(as.matrix(netSim), gmode='graph', vertex.col = vertex_attr(mySBM, "memberships"))
networkVec <- mat2Vect(as.matrix(netSim), symmetric = TRUE)
#covStruc <- makeVariance(network=netSim)

# Data simulation: Gaussian
Omega <- graph2prec(mySBM, cond_var = rep(1, nbNodes), neg_prop = 0.5)
Sigma <- solve(Omega)
means <- rep(0, ncol(Sigma))
data <- rmgaussian(nbObs, means, Sigma)

# Gaussian scores
scoreList <- list()
scoreList$huge <- getScoreHuge(data)
scoreList$saturnin <- saturnin::edge.prob(saturnin::lweights_gaussian(data))
scoreList$gnPcor <- GeneNet::ggm.estimate.pcor(data)^2
scoreList$gnPval <- vect2Mat(1 - GeneNet::network.test.edges(GeneNet::ggm.estimate.pcor(data), verbose=FALSE)$pval, symmetric=TRUE)

r <- max(scoreList$saturnin[scoreList$saturnin < 0.9999])
scoreList$saturnin[scoreList$saturnin > 0.9999] = r
scoreAdHocList <- list(huge= log(scoreList$huge), saturnin=logit(scoreList$saturnin),
                       gnPcor=atanh(GeneNet::ggm.estimate.pcor(data))^2, gnPval=qnorm(scoreList$gnPval))



################################################################################
# Parms for np-inference
trans <- 2
scoreMat <- scoreList2scoreMat(scoreAdHocList, symmetric=TRUE)
nbScores <- ncol(scoreMat); nbEdges <- nrow(scoreMat)

# Distribution
par(mfrow=c(2,nbScores/2))
for(s in 1:length(scoreList)){
  plot(density(scoreMat[, s]), main = '', col="black")
  lines(density(scoreMat[which(networkVec==0), s]), col="blue")
  lines(density(scoreMat[which(networkVec==1), s]), col="red")
}

# Kernel width + symmetric Gram matrix
nbScores = 3
kerSigma <- Hpi(scoreMat[,1:nbScores])
gram <- sapply(1:nbEdges, function(ij){dmvnorm(scoreMat[, 1:nbScores], mean=scoreMat[ij, 1:nbScores], sigma = kerSigma)})
gramMarg <- list()
for(s in 1:nbScores){
  gramMarg[[s]] <- sapply(1:nbEdges, function(ij){dnorm(scoreMat[, s], mean=scoreMat[ij, s], sd=sqrt(kerSigma[s, s]))})
  }

################################################################################
# Fit np - Noisy SBM
scoreListPar<- list(scoreList[[1]], scoreList[[2]], scoreList[[3]])
estimateNoisySBM(scoreListPar, directed=FALSE, nparm = TRUE)

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

