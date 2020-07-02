# Attempts to fit np-noisySBM network inference to simulated scores.

rm(list=ls()); par(mfrow=c(1, 1), mex=.5, pch=20)
library(mvtnorm); library(sna); library(rggm); library(ks)
library(blockmodels)
library(MASS); library(mclust);library(pbmcapply)
library(gtools)
# library(NoisySBM)
source('Functions/funcSimuls.R')

#adresse_package <- '~/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/Noisy_SBM_Package/NoisySBM/'
adresse_package <- 'D:/WORK_ALL/RECHERCHE/PACKAGES_R/SCORESBM/NOISY_SBM_Package/Package/NoisySBM/'

source(paste(adresse_package, 'R/funcVEM.R',sep = ''))
source(paste(adresse_package, 'R/tools.R',sep = ''))
source(paste(adresse_package, 'R/estimateNoisySBM.R',sep = ''))
source(paste(adresse_package, 'R/initInferenceNoisySBM.R',sep = ''))

#source('../../NoisySBM/R/funcVEM.R')
#source('../../NoisySBM/R/tools.R')
#source('../../NoisySBM/R/estimateNoisySBM.R')
#source('../../NoisySBM/R/initInferenceNoisySBM.R')

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
                       gnPcor=atanh(GeneNet::ggm.estimate.pcor(data))^2, gnPval=qnorm(1-scoreList$gnPval))



################################################################################
# Parms for np-inference

scoreList <- scoreAdHocList 

scoreMat <- scoreList2scoreMat(scoreAdHocList, symmetric=TRUE)
nbScores <- ncol(scoreMat); nbEdges <- nrow(scoreMat)

# Distribution
par(mfrow = c(2,nbScores/2))
for(s in 1:length(scoreList)){
  plot(density(scoreMat[, s]), main = '', col="black")
  lines(density(scoreMat[which(networkVec==0), s]), col="blue")
  lines(density(scoreMat[which(networkVec==1), s]), col="red")
}

# Kernel width + symmetric Gram matrix
nbScores = 4

kerSigma <- Hpi(scoreMat[,1:nbScores])
gram <- sapply(1:nbEdges, function(ij){dmvnorm(scoreMat[, 1:nbScores], mean=scoreMat[ij, 1:nbScores], sigma = kerSigma)})
gramMarg <- list()
for (s in 1:nbScores){
  gramMarg[[s]] <- sapply(1:nbEdges, function(ij){dnorm(scoreMat[, s], mean=scoreMat[ij, s], sd=sqrt(kerSigma[s, s]))})
  }

################################################################################
# Fit model  with 3 blocs
################################################################################

#----------  Init générale noisySBM
init <- initInferenceNoisySBM(scoreList, directed = directed)



#------------------  Init à 3 blocs
qDistInit = list(psi = init$psi, tau = init$tau[[3]], eta = init$eta[[3]])
clusterInit <- apply(qDistInit$tau,1,which.max)
table(clusterInit,clusterTrue)

#----------------- EstimNparam
vemNparm <- VEMNoisySBM(scoreMat , directed = directed, qDistInit = qDistInit , nparm = TRUE,gram=gram,monitoring = list(lowerBound = TRUE))
clusterEstimNparm <- apply(vemNparm$qDist$tau,1,which.max)
table(clusterEstimNparm,clusterTrue)


#-----------------  EstimNparam
vemGauss <- VEMNoisySBM(scoreMat , directed = directed, qDistInit = qDistInit ,estimOptions=list(verbosity =2),monitoring = list(lowerBound = TRUE))
clusterEstimGauss <- apply(vemGauss$qDist$tau,1,which.max)
table(clusterEstimGauss,clusterTrue)


################################################################################
# Fit model  with alls blocs
################################################################################

#----------------- EstimNparam
estimNparm <- estimateNoisySBM(scoreList, directed = directed, nparm = TRUE,kerSigma = kerSigma, estimOptions = list(verbosity = 2),monitoring = list(lowerBound = TRUE))

#---------------------- EstimPan
estimParm <- estimateNoisySBM(scoreList, directed = directed, nparm = TRUE,monitoring = list(lowerBound = TRUE))

