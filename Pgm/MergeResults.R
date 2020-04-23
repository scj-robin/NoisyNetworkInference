# Gathering simulation results
rm(list=ls())

library(mclust); library(ROCR)
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')

# Dirs
simDir = '../Simul/'

# Dims
Ktrue = 3; g = 2; nbSim = 50; signCor = 'pos'
pList = c(20, 30, 50, 80); nList = c(20, 50, 100, 200); 
Ktrue = 3; g = 2; nbSim = 50; 
pList = c(20, 30, 50, 80); betaList = c(20, 50, 100, 200)/2; 

Parms = AUCvemAll = AUCoracleAll = ARIvemAll = ARIoracleAll = iterVemAll = c()
# Merge results
for(p in pList){
   # p = 20
   for(beta in betaList){
      # n = 20
      # simParms = paste0('simVEM-p', p, '-n', n, '-Ktrue', Ktrue, '-g', g, '-C', signCor, '/')
      simParms = paste0('simWithScores-p', p, '-beta', beta, '-Ktrue', Ktrue, '-g', g, '/')
      cat('\n', simParms, ':')
      
      # Results
      K = Ktrue
      ARIvem = matrix(0, nbSim, 7); colnames(ARIvem) = c('vemMB', 'vemGlasso', 'vemTree', 'vemLogitTree', 'vemMBTree', 
                                                      'vemnpLogitTree', 'vemnpMBTree')
      ARIoracle = matrix(0, nbSim, 4); colnames(ARIoracle) = c('sbmG', 'sbmMB', 'sbmGlasso', 'sbmTree')
      AUCvem = matrix(0, nbSim, 7); colnames(AUCvem) = c('vemMB', 'vemGlasso', 'vemTree', 'vemLogitTree', 'vemMBTree', 
                                                      'vemnpLogitTree', 'vemnpMBTree')
      AUCoracle = matrix(0, nbSim, 3); colnames(AUCoracle) = c('MB', 'Glasso', 'Tree')
      iterVem = matrix(0, nbSim, 7); colnames(iterVem) = c('vemMB', 'vemGlasso', 'vemTree', 'vemLogitTree', 'vemMBTree', 
                                                         'vemnpLogitTree', 'vemnpMBTree')
      lastSim = 0
      for(sim in 1:nbSim){
         # sim = 1
         simLabel = paste0('simul-K', K, '-s', sim); cat('', sim)
         if(file.exists(paste0(simDir, simParms, simLabel, '.Rdata'))){
            load(file=paste0(simDir, simParms, simLabel, '.Rdata'))
            if(p>=n){simRes$scoreGlasso = simRes$scoreMB; simRes$vemGlasso = simRes$vemMB; 
            simRes$sbmGlasso = simRes$sbmMB}
            if(typeof(simRes$sbmGlasso)=='character'){
               simRes$sbmGlasso$memberships = list(); 
               tauTmp = matrix(0, p, K); tauTmp[, 1] = 1
               simRes$sbmGlasso$memberships[[K]] = list(Z = tauTmp)
            }
            # Check class order
            if(simRes$vemMB$phi0[1] > simRes$vemMB$phi1[1]){cat('Pb MB: mu0 > mu1\n')}
            if(simRes$vemGlasso$phi0[1] > simRes$vemGlasso$phi1[1]){cat('Pb Glasso: mu0 > mu1\n')}
            if(simRes$vemTree$phi0[1] > simRes$vemTree$phi1[1]){cat('Pb Tree: mu0 > mu1\n')}
            if(mean(simRes$vemMBTree$phi0$mu) > mean(simRes$vemMBTree$phi1$mu)){cat('Pb MB-Tree: mu0 > mu1\n')}
            
            # ARI
            Ztrue = as.vector(simRes$dataSim$Z%*%(1:K))
            ARIvem[sim, ] = c(
               adjustedRandIndex(Ztrue, apply(simRes$vemMB$tau, 1, which.max)), 
               adjustedRandIndex(Ztrue, apply(simRes$vemGlasso$tau, 1, which.max)), 
               adjustedRandIndex(Ztrue, apply(simRes$vemTree$tau, 1, which.max)),
               adjustedRandIndex(Ztrue, apply(simRes$vemLogitTree$tau, 1, which.max)),
               adjustedRandIndex(Ztrue, apply(simRes$vemMBTree$tau, 1, which.max)), 
               adjustedRandIndex(Ztrue, apply(simRes$vemnpLogitTree$tau, 1, which.max)), 
               adjustedRandIndex(Ztrue, apply(simRes$vemnpMBTree$tau, 1, which.max)) 
            )
            ARIoracle[sim, ] = c(
               adjustedRandIndex(Ztrue, apply(simRes$sbmG$memberships[[K]]$Z, 1, which.max)),
               adjustedRandIndex(Ztrue, apply(simRes$sbmMB$memberships[[K]]$Z, 1, which.max)), 
               adjustedRandIndex(Ztrue, apply(simRes$sbmGlasso$memberships[[K]]$Z, 1, which.max)), 
               adjustedRandIndex(Ztrue, apply(simRes$sbmTree$memberships[[K]]$Z, 1, which.max))
            )
            # AUC
            Gvec = mat_vect_low(simRes$dataSim$G)
            AUCvem[sim, ] = c(
               performance(prediction(simRes$vemMB$Psi1, Gvec), "auc")@y.values[1][[1]],
               performance(prediction(simRes$vemGlasso$Psi1, Gvec), "auc")@y.values[1][[1]],
               performance(prediction(simRes$vemTree$Psi1, Gvec), "auc")@y.values[1][[1]],
               performance(prediction(simRes$vemLogitTree$Psi1, Gvec), "auc")@y.values[1][[1]],
               performance(prediction(simRes$vemMBTree$Psi1, Gvec), "auc")@y.values[1][[1]],
               performance(prediction(simRes$vemnpLogitTree$Psi1, Gvec), "auc")@y.values[1][[1]],
               performance(prediction(simRes$vemnpMBTree$Psi1, Gvec), "auc")@y.values[1][[1]]
            )
            AUCoracle[sim, ] = c(
               performance(prediction(mat_vect_low(simRes$scoreMB), Gvec), "auc")@y.values[1][[1]],
               performance(prediction(mat_vect_low(simRes$scoreGlasso), Gvec), "auc")@y.values[1][[1]],
               performance(prediction(mat_vect_low(simRes$scoreTree), Gvec), "auc")@y.values[1][[1]]
            )
            # ARI
            iterVem[sim, ] = round(c(
               length(simRes$vemMB$borne_inf), 
               length(simRes$vemGlasso$borne_inf), 
               length(simRes$vemTree$borne_inf),
               length(simRes$vemLogitTree$borne_inf),
               length(simRes$vemMBTree$borne_inf), 
               length(simRes$vemnpLogitTree$borne_inf), 
               length(simRes$vemnpMBTree$borne_inf) 
            )/3)
            lastSim = sim
         }
      }
      if(lastSim>1){
         ARIvem = ARIvem[1:lastSim, ]; AUCvem = AUCvem[1:lastSim, ]; 
         ARIoracle = ARIoracle[1:lastSim, ]; AUCoracle = AUCoracle[1:lastSim, ];
         iterVem = iterVem[1:lastSim, ]; 
         Parms = rbind(Parms, cbind(rep(p, lastSim), rep(n, lastSim), (1:lastSim)))
         ARIvemAll = rbind(ARIvemAll, ARIvem); AUCvemAll = rbind(AUCvemAll, AUCvem); 
         ARIoracleAll = rbind(ARIoracleAll, ARIoracle); AUCoracleAll = rbind(AUCoracleAll, AUCoracle);
         iterVemAll = rbind(iterVemAll, iterVem); 
      }
   }
}

# Merging
colnames(Parms) = names(Parms) = c('p', 'n', 'sim')
colnames(AUCvemAll) = paste0('AUC.', colnames(AUCvemAll))
colnames(ARIvemAll) = paste0('ARI.', colnames(ARIvemAll))
colnames(AUCoracleAll) = paste0('AUC.', colnames(AUCoracleAll))
colnames(ARIoracleAll) = paste0('ARI.', colnames(ARIoracleAll))
colnames(iterVemAll) = paste0('iter.', colnames(iterVemAll))
Res = as.data.frame(cbind(Parms, ARIvemAll, AUCvemAll, ARIoracleAll, AUCoracleAll, iterVemAll))

# Export
write.csv(Res, file=paste0(simDir, 'simVEM-Ktrue', Ktrue, '-g', g, '-C', signCor, '-K', K, '.csv'), 
          row.names=FALSE)

