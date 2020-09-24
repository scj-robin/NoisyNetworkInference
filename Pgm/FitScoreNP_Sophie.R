# Attempts to fit np-ScoreSBM network inference to simulated scores.

rm(list = ls()); 
par(mfrow = c(1, 1), mex = .5, pch = 20)
library(mvtnorm); library(sna); library(rggm); library(ks)
library(blockmodels)
library(MASS); library(mclust);library(pbmcapply)
library(igraph)
library(gtools)
library(sbm)
#library(ScoreSBM)

################################################################################
#------------------------------- FUNCTIONS
################################################################################
setwd("~/WORK/RECHERCHE/TRAVAUX_DE_RECHERCHE/StéphaneRobin/NoisySBM/NoisyNetworkInference/Pgm")
source('Functions/funcSimuls.R')
adresse_package <- '/home/sophie/WORK/RECHERCHE/PACKAGES_R/ScoreSBM/ScoreSBM/'
source(paste(adresse_package, 'R/funcVEM.R',sep = ''))
source(paste(adresse_package, 'R/tools.R',sep = ''))
source(paste(adresse_package, 'R/estimateScoreSBM.R',sep = ''))
source(paste(adresse_package, 'R/initInferenceScoreSBM.R',sep = ''))

# source('../../ScoreSBM/R/funcVEM.R')
# source('../../ScoreSBM/R/tools.R')
# source('../../ScoreSBM/R/estimateScoreSBM.R')
# source('../../ScoreSBM/R/initInferenceScoreSBM.R')
#source('../../ScoreSBM/R/funcVEM.R')
#source('../../ScoreSBM/R/tools.R')
#source('../../ScoreSBM/R/estimateScoreSBM.R')
#source('../../ScoreSBM/R/initInferenceScoreSBM.R')

################################################################################
# Data
################################################################################
# Param sim
nbNodes <- 60; directed <- FALSE; nbObs <- 1000; seed <- 2

blockProp <- c(1/3, 1/3, 1/3)  # group proportions
nbBlocks   <- length(blockProp) # number of blocks
connectParam <- diag(.4, nbBlocks) + 0.05 # connectivity matrix: affiliation network

# Network simulation
mySBM <- rSBM(nbNodes, connectParam, blockProp)
netSim <- as_adjacency_matrix(mySBM)
clusterTrue <- vertex_attr(mySBM, "memberships")
sna::gplot(as.matrix(netSim), gmode='graph', vertex.col = vertex_attr(mySBM, "memberships"))
networkVec <- mat2Vect(as.matrix(netSim), symmetric = TRUE)
#covStruc <- makeVariance(network=netSim)

# Data simulation: Gaussian
Omega <- graph2prec(mySBM, cond_var = rep(1, nbNodes), neg_prop = 0.5)
Sigma <- solve(Omega)
means <- rep(0, ncol(Sigma))
data <- rmvnorm(nbObs, means, Sigma)
 
O <- BM_bernoulli('SBM_sym',as.matrix(netSim))
O$estimate()
I <- which.max(O$ICL)
clusterOracle <- apply(O$memberships[[I]]$Z,1,which.max)
table(clusterOracle, clusterTrue)
################################################################################
#  Score on edges
################################################################################
scoreList <- list()
scoreList$huge <- getScoreHuge(data)
scoreList$saturnin <- saturnin::edge.prob(saturnin::lweights_gaussian(data))
hist(scoreList$saturnin)
scoreList$gnPcor <- GeneNet::ggm.estimate.pcor(data)^2
scoreList$gnPval <- vect2Mat(1 - GeneNet::network.test.edges(GeneNet::ggm.estimate.pcor(data), verbose = FALSE)$pval, symmetric = TRUE)




r <- max(scoreList$saturnin[scoreList$saturnin < 0.99999])
scoreList$saturnin[scoreList$saturnin > 0.99999] = r
r <- min(scoreList$saturnin[scoreList$saturnin > 0.00001])
scoreList$saturnin[scoreList$saturnin < 0.00001] = r

scoreTransfo <- list(
  huge = log(scoreList$huge), 
  saturnin = logit(scoreList$saturnin),
  gnPcor = atanh(GeneNet::ggm.estimate.pcor(data))^2, 
  gnPval = qnorm(1 - scoreList$gnPval)
)

# PLOTS des scores
scoreMat <- scoreList2scoreMat(scoreTransfo, symmetric = TRUE)
nbScores <- ncol(scoreMat); nbEdges <- nrow(scoreMat)
par(mfrow = c(2,1))
for (s in 1:2){
  hist(scoreMat[, s], main = '', col = "black")
  hist(scoreMat[which(networkVec == 0), s], col = "blue",add=TRUE)
  hist(scoreMat[which(networkVec == 1), s], col = "red",add=TRUE)
}

################################################################################


# Kernel width + symmetric Gram matrix
whichScores <- 1:1
myScoreList = lapply(whichScores,function(l){scoreTransfo[[l]]})
myScoreMat <- scoreList2scoreMat(myScoreList, symmetric = TRUE)
#myScoreMat[,2] <- myScoreMat[,2]+rnorm(nbEdges,0,10^(-1))

plot(as.data.frame(myScoreMat),col=networkVec + 1)
resMclust <- Mclust(myScoreMat,G = 2)
#plot(resMclust)
# 
kerSigma <- Hpi(myScoreMat)
kerSigma <- 0.5 * (kerSigma + t(kerSigma))
gram <- sapply(1:nbEdges, function(ij){dmvnorm(myScoreMat, mean = myScoreMat[ij, ], sigma = kerSigma)})
gramMarg <- list()
for (s in 1:ncol(myScoreMat)) {
  gramMarg[[s]] <- sapply(1:nbEdges, function(ij){dnorm(myScoreMat[, s], mean = myScoreMat[ij, s], sd = sqrt(kerSigma[s, s]))})
  }

# 
kerSigmaDiag <- diag(sapply(1:ncol(myScoreMat),function(s){hpi(myScoreMat[,s])}))
gramDiag <- sapply(1:nbEdges, function(ij){dmvnorm(myScoreMat, mean = myScoreMat[ij, ], sigma = kerSigmaDiag)})
if (ncol(myScoreMat) == 1){
  gramDiag <- sapply(1:nbEdges, function(ij){dnorm(myScoreMat, mean = myScoreMat[ij, ], kerSigmaDiag)})
}


################################################################################
####################### Init générale ScoreSBM
init <- initInferenceScoreSBM(myScoreList, directed = directed)
init1 <- initInferenceScoreSBM(list(myScoreList[[1]]), directed = directed)
init2 <- initInferenceScoreSBM(list(myScoreList[[2]]), directed = directed)


par(mfrow=c(1,1))
plot(Mclust(myScoreMat,G=2),what = 'classification')
length(init)
table(networkVec , init$G)
mean(init$psi[,2])
mean(networkVec)
plot(init$psi[,2],networkVec)

#### on cherche une autres init
rankAverage <- rowMeans(apply(myScoreMat,2,rank))
hist(rankAverage, main = '', col = "black")
hist(rankAverage[which(networkVec == 0)], col = "blue",add=TRUE)
hist(rankAverage[which(networkVec == 1)], col = "red",add=TRUE)
mcRank <- Mclust(rankAverage, G = 2)
#plot(mcRank)
table(networkVec ,apply(mcRank$z,1,which.max)-1)
#### 
dens <- mean(networkVec)
Ginit <- (rankAverage  > (nbEdges * (1-dens)))
table(Ginit,networkVec)


# Fit model  with 3 blocs
################################################################################

#  Init à 3 blocs
qDistInit = list(psi = init$psi, tau = init$tau[[3]], eta = init$eta[[3]])
qDistInit1 = list(psi = init1$psi, tau = init1$tau[[3]], eta = init1$eta[[3]])
qDistInit2 = list(psi = init2$psi, tau = init2$tau[[3]], eta = init2$eta[[3]])


clusterInit <- apply(qDistInit$tau,1,which.max)
table(clusterInit,clusterTrue)
table(clusterInit,clusterOracle)


#----------------- EstimNparam
vemNparm <- VEMScoreSBM(myScoreMat , directed = directed, qDistInit = qDistInit , nparm = TRUE, gram = gramDiag,monitoring = list(lowerBound = TRUE))
plot(vemNparm$lowerBound,type = 'l')


vemNparm1 <- VEMScoreSBM(myScoreMat , directed = directed, qDistInit = qDistInit1 , nparm = TRUE, gram = gram,monitoring = list(lowerBound = TRUE))
plot(vemNparm1$lowerBound,type = 'l')


boxplot(vemNparm$qDist$psi[,2] ~ networkVec)

vemNparm2 <- VEMScoreSBM(myScoreMat , directed = directed, qDistInit = qDistInit2 , nparm = TRUE, gram = gram,monitoring = list(lowerBound = TRUE))
plot(vemNparm2$lowerBound,type = 'l')


max(vemNparm$lowerBound)
max(vemNparm1$lowerBound)
max(vemNparm2$lowerBound)


GestimNparm <- apply(vemNparm$qDist$psi,1, which.max) - 1
table(apply(init$psi,1, which.max) - 1,networkVec)

table(GestimNparm,networkVec)

clusterEstimNparm <- apply(vemNparm$qDist$tau,1,which.max)
clusterEstimNparm1 <- apply(vemNparm1$qDist$tau,1,which.max)
clusterEstimNparm2 <- apply(vemNparm2$qDist$tau,1,which.max)

table(clusterEstimNparm1,clusterEstimNparm)
table(clusterEstimNparm1,clusterEstimNparm2)


table(clusterEstimNparm,clusterTrue)
table(clusterEstimNparm,clusterOracle)
table(clusterInit, clusterEstimNparm)

#----------------- EstimNparamDiag
vemNparmDiag <- VEMScoreSBM(myScoreMat , directed = directed, qDistInit = qDistInit , nparm = TRUE, gram = gramDiag,monitoring = list(lowerBound = TRUE))
plot(vemNparmDiag$lowerBound,type = 'l')

clusterEstimNparmDiag <- apply(vemNparmDiag$qDist$tau,1,which.max)
table(clusterEstimNparmDiag,clusterTrue)
table(clusterEstimNparmDiag,clusterOracle)
table(clusterInit, clusterEstimNparmDiag)
table(clusterEstimNparm, clusterEstimNparmDiag)
mclust::adjustedRandIndex(clusterInit,clusterTrue)
mclust::adjustedRandIndex(clusterOracle,clusterEstimNparm)
mclust::adjustedRandIndex(clusterOracle,clusterEstimNparmDiag)

#-----------------  Estim Gaussienne

vemGauss <- VEMScoreSBM(myScoreMat , directed = directed, qDistInit = qDistInit ,estimOptions = list(verbosity = 0),monitoring = list(lowerBound = TRUE))
clusterEstimGauss <- apply(vemGauss$qDist$tau,1,which.max)
table(clusterEstimGauss,clusterTrue)
table(clusterEstimGauss,clusterOracle)


################################################################################
# Fit model  with alls blocs
################################################################################

#----------------- Estim Non param
estimNparm <- estimateScoreSBM(myScoreList, directed = directed, nparm = TRUE,kerSigma = kerSigma, estimOptions = list(verbosity = 2),monitoring = list(lowerBound = TRUE))
bestNbBlocksNparm <- estimNparm[[1]]$nbBlocks
print(bestNbBlocksNparm)
#---------------------- Estim Non Param
estimParm <- estimateScoreSBM(myScoreList, directed = directed, nparm = FALSE,monitoring = list(lowerBound = TRUE))
bestNbBlocksParm <- estimParm[[1]]$nbBlocks
print(bestNbBlocksParm)
