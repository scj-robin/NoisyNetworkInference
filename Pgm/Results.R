# dataSimation of edge scores
rm(list=ls())

library(mclust); library(ROCR)
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')

# Dirs
simDir = '../Simul/'

# Dims
p = 30; n = 100; Ktrue = 3; g = 2; nbSim = 50
simParms = paste0('simVEM-p', p, '-n', n, '-Ktrue', Ktrue, '-g', g, '/')

# Resultats
K = Ktrue
ARI = matrix(0, nbSim, 4); colnames(ARI) = c('Oracle', 'MB', 'Glasso', 'Tree')
AUC = matrix(0, nbSim, 6); colnames(AUC) = c('MB', 'vemMB', 'Glasso', 'vemGlasso', 'Tree', 'vemTree')
lastSim = 0
for(sim in 1:nbSim){
   # sim = 1
   simLabel = paste0('simul-K', K, '-s', sim); print(paste0(simParms, simLabel))
   if(file.exists(paste0(simDir, simParms, simLabel, '.Rdata'))){
      load(file=paste0(simDir, simParms, simLabel, '.Rdata'))
      if(p>=n){simRes$scoreGlasso = simRes$scoreMB; simRes$vemGlasso = simRes$vemMB; simRes$vemLogGlasso = simRes$vemLogMB}
      # ARI
      Ztrue = as.vector(simRes$dataSim$Z%*%(1:K))
      ARI[sim, ] = c(
         adjustedRandIndex(Ztrue, apply(simRes$sbmG$memberships[[K]]$Z, 1, which.max)),
         adjustedRandIndex(Ztrue, apply(simRes$vemMB$tau, 1, which.max)), 
         # adjustedRandIndex(Ztrue, apply(simRes$vemLogMB$tau, 1, which.max)), 
         adjustedRandIndex(Ztrue, apply(simRes$vemGlasso$tau, 1, which.max)), 
         # adjustedRandIndex(Ztrue, apply(simRes$vemLogGlasso$tau, 1, which.max)), 
         adjustedRandIndex(Ztrue, apply(simRes$vemTree$tau, 1, which.max))
         # adjustedRandIndex(Ztrue, apply(simRes$vemLogitTree$tau, 1, which.max))
         )
      # AUC
      Gvec = mat_vect_low(simRes$dataSim$G)
      AUC[sim, ] = c(
         performance(prediction(mat_vect_low(simRes$scoreMB), Gvec), "auc")@y.values[1][[1]],
         performance(prediction(simRes$vemMB$Psi1, Gvec), "auc")@y.values[1][[1]],
         # performance(prediction(simRes$vemLogMB$Psi1, Gvec), "auc")@y.values[1][[1]],
         performance(prediction(mat_vect_low(simRes$scoreGlasso), Gvec), "auc")@y.values[1][[1]],
         performance(prediction(simRes$vemGlasso$Psi1, Gvec), "auc")@y.values[1][[1]],
         # performance(prediction(simRes$vemLogGlasso$Psi1, Gvec), "auc")@y.values[1][[1]],
         performance(prediction(mat_vect_low(simRes$scoreTree), Gvec), "auc")@y.values[1][[1]], 
         performance(prediction(simRes$vemTree$Psi1, Gvec), "auc")@y.values[1][[1]]
         # performance(prediction(simRes$vemLogitTree$Psi1, Gvec), "auc")@y.values[1][[1]]
         )
      lastSim = sim
   }
}
ARI = ARI[1:lastSim, ]; AUC = AUC[1:lastSim, ]; 
ARIcol = rep(0, ncol(ARI)); AUCcol = rep(0, ncol(AUC)); 
ARIcol[1] = AUCcol[c(1, 3, 5)] = 4
ARIborder = rep(1, ncol(ARI)); AUCborder = rep(1, ncol(AUC)); 
if(p>=n){ARIborder[3:4] = AUCborder[3:4] = 8}
par(mfrow=c(2, 1), pch=20, mex=.6)
boxplot(ARI, main='ARI', ylim=c(0, 1), col=ARIcol, border=ARIborder)
boxplot(AUC, main='AUC', ylim=c(.5, 1), col=AUCcol, border=AUCborder)
