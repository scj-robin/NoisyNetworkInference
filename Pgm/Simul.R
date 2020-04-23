# Simulates data, fit network inference, get scores

library(gtools); library(sna); library(mvtnorm); library(huge); library(saturnin)
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

# Data simulation
covStruc <- makeVariance(network=netSim$network)
data <- rmvnorm(nbObs, sigma=covStruc$sigma)
plot(covStruc$sigma, cov(data), pch=20); abline(0, 1)

# Scores
scoreHuge <- getHugeScore(data); 
scoreSaturnin <- getSaturninScore(data); 

# Transformation

