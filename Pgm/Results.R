# dataSimation of edge scores
rm(list=ls())

library(mclust); library(ROCR)
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')

# Dirs
simDir = '../Simul/'

# Dims
p = 30; n = 100; K = 3; g = 2; gLassoMethod = 'mb'
simName = paste0('dataSimVEM-n', n, '-p', p, '-K', K, '-g', g, '-', gLassoMethod)

# Resultats
load(file=paste0(simDir, simName, '.Rdata'))
nbSim = length(simRes)
ARI = matrix(0, nbSim, 5); colnames(ARI) = c('Oracle', 'Glasso', 'logGlasso', 'Tree', 'logitTree')
AUC = matrix(0, nbSim, 6); colnames(AUC) = c('Glasso', 'vemGlasso', 'vemLogGlasso', 'Tree', 'vemTree', 'vemLogitTree')
for(sim in 1:nbSim){
   # ARI
   Ztrue = as.vector(simRes[[sim]]$dataSim$Z%*%(1:K))
   ARI[sim, ] = c(adjustedRandIndex(Ztrue, apply(simRes[[sim]]$sbmG$memberships[[K]]$Z, 1, which.max)),
                  adjustedRandIndex(Ztrue, apply(simRes[[sim]]$vemGlasso$tau, 1, which.max)), 
                  adjustedRandIndex(Ztrue, apply(simRes[[sim]]$vemLogGlasso$tau, 1, which.max)), 
                  adjustedRandIndex(Ztrue, apply(simRes[[sim]]$vemTree$tau, 1, which.max)), 
                  adjustedRandIndex(Ztrue, apply(simRes[[sim]]$vemLogitTree$tau, 1, which.max)))
   # AUC
   Gvec = mat_vect_low(simRes[[sim]]$dataSim$G)
   AUC[sim, ] = c(performance(prediction(mat_vect_low(simRes[[sim]]$scoreGlasso), Gvec), "auc")@y.values[1][[1]],
                  performance(prediction(simRes[[sim]]$vemGlasso$Psi1, Gvec), "auc")@y.values[1][[1]],
                  performance(prediction(simRes[[sim]]$vemLogGlasso$Psi1, Gvec), "auc")@y.values[1][[1]],
                  performance(prediction(mat_vect_low(simRes[[sim]]$scoreTree), Gvec), "auc")@y.values[1][[1]], 
                  performance(prediction(simRes[[sim]]$vemTree$Psi1, Gvec), "auc")@y.values[1][[1]],
                  performance(prediction(simRes[[sim]]$vemLogitTree$Psi1, Gvec), "auc")@y.values[1][[1]])
}
par(mfrow=c(2, 1), pch=20, mex=.6)
boxplot(ARI)
boxplot(AUC)

