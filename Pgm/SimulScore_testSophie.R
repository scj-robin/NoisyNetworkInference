# Simulates data, fit network inference, get scores

rm(list=ls()); par(mfrow=c(1, 1), pch=20)

library(gtools); library(sna); library(MASS)
library(mvtnorm); library(mclust); library(PLNmodels)
library(huge); library(saturnin); library(GeneNet)
library(NoisySBM)
library(ggplot2)
library(rggm)

source('Functions/funcSimuls.R')


set.seed(13)


# Dims
nNodes  <- 60
blockProp <- c(1/3, 1/3, 1/3)  # group proportions
nbBlock   <- length(blockProp) # number of blocks
connectParam <- diag(.4, nbBlock) + 0.01 # connectivity matrix: affiliation network
mySBM <- rSBM(nNodes, connectParam, blockProp)
library(igraph)
netSim <- as_adjacency_matrix(mySBM)


clusterTrue <- vertex_attr(mySBM, "memberships")

# Network simulation

sna::gplot(as.matrix(netSim), gmode='graph', vertex.col = vertex_attr(mySBM, "memberships"))
networkVec <- mat2Vect(as.matrix(netSim), symmetric = TRUE)
#covStruc <- makeVariance(network=netSim)

# Data simulation: Gaussian
Omega <- graph2prec(mySBM, cond_var = rep(1, nNodes), neg_prop = 0.5)
Sigma <- solve(Omega)
n <- 300
means <- rep(0, ncol(Sigma))
data <- rmgaussian(n, means, Sigma)

# Gaussian scores
scoreList <- list()
scoreList$huge <- getScoreHuge(data)
scoreList$saturnin <- saturnin::edge.prob(saturnin::lweights_gaussian(data))
scoreList$gnPcor <- GeneNet::ggm.estimate.pcor(data)^2
scoreList$gnPval <- vect2Mat(1 - GeneNet::network.test.edges(GeneNet::ggm.estimate.pcor(data), verbose=FALSE)$pval, symmetric=TRUE)

r <- max(scoreList$saturnin[scoreList$saturnin < 0.9999])
scoreList$saturnin[scoreList$saturnin > 0.9999] = r
scoreAdHocList <- list(huge=log(scoreList$huge), saturnin=logit(scoreList$saturnin),
                       gnPcor=atanh(GeneNet::ggm.estimate.pcor(data))^2, gnPval=qnorm(scoreList$gnPval))


# PLOT
nbScores <- length(scoreList) 
SS <- data.frame(network = as.factor(rep(networkVec ,nbScores)))
SS$namesScores <- as.factor(rep(c('huge','saturnin','gnPcor','gnPval'),each = length(networkVec)))
U <-  sapply(1:nbScores, function(p){mat2Vect(scoreAdHocList[[p]],symmetric = TRUE)})
SS$Scores <- c(U)
ggplot(SS,aes(x = network,y = Scores,color = namesScores)) + geom_boxplot() + facet_wrap(~namesScores,scales = "free")
ggplot(SS,aes(x = Scores,group = network,color = network)) + geom_density() + facet_wrap(~namesScores,scales = "free")




# #--------------------
library(cluster)
df = U
df <- U[,-4]

# Hierarchical clustering using Complete Linkage
hc1 <- hclust(dist(df, method = "euclidean") )
gpHC <- cutree(hc1, k = 2) - 1
M1 <- mean(df[gpHC==1,])
M0 <- mean(df[gpHC==0,])
if (M1 < M0) {gpHC = 1- gpHC }
table(gpHC,networkVec)


res_Mclust <- Mclust(df, G = 2)
gpMCLUST <- res_Mclust$classification - 1
mu  <- res_Mclust$parameters$mean
test_G <- rowMeans(t(mu)) #identify G = 0  and G =1
if (test_G[1] > test_G[2]) {
  gpMCLUST = 1 - gpMCLUST
}
table(networkVec)
table(gpHC)
table(gpMCLUST)

table(gpMCLUST,networkVec)
table(gpMCLUST,gpHC)

sum(gpMCLUST == networkVec)
sum(gpHC == networkVec)

#---------------------------------
library(blockmodels)
resOracle <- BM_bernoulli("SBM",as.matrix(netSim))
resOracle$estimate()
M <- which.max(resOracle$ICL)
clusterOracle <- apply(resOracle$memberships[[M]]$Z,1,which.max) 
############################################################################################

scoreAdHocList <- list(huge=log(scoreList$huge), saturnin=qnorm(scoreList$saturnin),
                                          gnPcor=atanh(GeneNet::ggm.estimate.pcor(data))^2)

resInit <- initInferenceNoisySBM(scoreAdHocList,directed = FALSE)

clusterInit3 <- apply(resInit$tau[[3]],1,which.max)
table(clusterInit3,clusterTrue)
table(clusterOracle,clusterTrue)


res_Estim <- estimateNoisySBM(scoreAdHocList, directed= FALSE,estimOptions = list(nbBlocksRange = c(1,10),exploreFactor = 12))
M <- which(vapply(res_Estim,function(res){res$nbBlocks},1) == 3)  
plot(res_Estim[[M]]$lowerBound,type = 'l')
clusterEstim <- apply(res_Estim[[M]]$qDist$tau,1,which.max)
table(clusterEstim,clusterTrue)



