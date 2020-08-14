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
source('Functions/funcSimuls.R')
adresse_package <- '~/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/ScoreSBM/'
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
nbNodes <- 60; directed <- FALSE; nbObs <- 500; seed <- 2

blockProp <- c(1/3, 1/3, 1/3)  # group proportions
nbBlocks   <- length(blockProp) # number of blocks
connectParam <- diag(.4, nbBlocks) + 0.1 # connectivity matrix: affiliation network

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
data <- rmgaussian(nbObs, means, Sigma)
 
O <- BM_bernoulli('SBM_sym',as.matrix(netSim))
O$estimate()
I <- which.max(O$ICL)
clusterOracle <- apply(O$memberships[[I]]$Z,1,which.max)

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
par(mfrow = c(2,nbScores/2))
for (s in 1:nbScores){
  plot(density(scoreMat[, s]), main = '', col = "black")
  lines(density(scoreMat[which(networkVec == 0), s]), col = "blue")
  lines(density(scoreMat[which(networkVec == 1), s]), col = "red")
}

################################################################################


# Kernel width + symmetric Gram matrix
whichScores <- 1:4
myScoreList = lapply(whichScores,function(l){scoreTransfo[[l]]})
myScoreMat <- scoreList2scoreMat(myScoreList, symmetric = TRUE)

kerSigma <- Hpi(myScoreMat)
kerSigma <- 0.5 * (kerSigma + t(kerSigma))
gram <- sapply(1:nbEdges, function(ij){dmvnorm(myScoreMat, mean = myScoreMat[ij, ], sigma = kerSigma)})
gramMarg <- list()
for (s in 1:ncol(myScoreMat)) {
  gramMarg[[s]] <- sapply(1:nbEdges, function(ij){dnorm(myScoreMat[, s], mean = myScoreMat[ij, s], sd = sqrt(kerSigma[s, s]))})
  }

################################################################################
####################### Init générale ScoreSBM
init <- initInferenceScoreSBM(myScoreList, directed = directed)
length(init)

# Fit model  with 3 blocs
################################################################################

#  Init à 3 blocs
qDistInit = list(psi = init$psi, tau = init$tau[[3]], eta = init$eta[[3]])
clusterInit <- apply(qDistInit$tau,1,which.max)
table(clusterInit,clusterTrue)
table(clusterInit,clusterOracle)


#----------------- EstimNparam
vemNparm <- VEMScoreSBM(myScoreMat , directed = directed, qDistInit = qDistInit , nparm = TRUE, gram = gram,monitoring = list(lowerBound = TRUE))
plot(vemNparm$lowerBound,type = 'l')

clusterEstimNparm <- apply(vemNparm$qDist$tau,1,which.max)
table(clusterEstimNparm,clusterTrue)
table(clusterEstimNparm,clusterOracle)


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
print(bestNbBlocksNparm)
