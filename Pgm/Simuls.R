# Simulation of edge scores
rm(list=ls()); #par(pch=20)

library(sna); library(mvtnorm); library(huge); library(blockmodels); library(EMtree); library(mclust); 
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')

# Dirs
simDir = '../Simul/'

# Dims : p = nb nodes, n = nb replicates
# pList = c(20, 30, 50, 80); nList = c(20, 50, 100)
pList = c(20, 30, 50, 80); nList = c(20, 50, 100, 200)
Ktrue = 3; g = 2; simNb = 50

for(p in pList){
   # p = 30
   for(n in nList){
      # n = 20
      simParms = paste0('simVEM-p', p, '-n', n, '-Ktrue', Ktrue, '-g', g, '/')
      
      # Parms
      pi = (1:Ktrue); pi = pi / sum(pi)
      rho = log(p)/p
      gamma = (pi^g)[Ktrue:1]%o%(pi^g)[Ktrue:1]; gamma = rho*gamma/(t(pi)%*%gamma%*%pi)[1, 1]
      gamma[which(gamma>1)] = 1
      parms = list(pi=pi, gamma=gamma, rho=rho)
      
      # Boucle de simulation
      K = Ktrue
      for(sim in 1:simNb){
         # sim = 1
         # Deja fait ?
         simLabel = paste0('simul-K', K, '-s', sim); print(paste0(simParms, simLabel))
         if(file.exists(paste0(simDir, simParms, simLabel, '.Rdata'))){
            load(paste0(simDir, simParms, simLabel, '.Rdata')); 
         }else{simRes = list(parms=parms)}
         # Simul data
         if(!is.element("dataSim", names(simRes))){
            dataSim = SimulZGY(pi, gamma, n, p)
            simRes$dataSim = dataSim
         }
         # Oracle : SBM sur G
         if(!is.element("sbmG", names(simRes))){
            sbmG = BM_bernoulli('SBM_sym', simRes$dataSim$G, plotting='', verbosity=0, ncores=1); 
            sbmG$estimate()
            simRes$sbmG = sbmG
         }
         # MB
         if(!is.element("vemLogMB", names(simRes))){
            scoreMB = fitHuge(simRes$dataSim$Y, method='mb')
            cat(simLabel, 'MB :'); vemMB = VEM(scoreMB, K); 
            cat('\n', simLabel, 'LogMB :'); vemLogMB = VEM(log(1+scoreMB), K); 
            simRes$scoreMB = scoreMB
            simRes$vemMB = vemMB
            simRes$vemLogMB = vemLogMB
         }
         # Tree
         if(!is.element("vemTree", names(simRes))){
            scoreTree = fitEMtree(simRes$dataSim$Y)
            cat('\n', simLabel, 'Tree:'); vemTree = VEM(scoreTree, K)
            simRes$scoreTree = scoreTree
            simRes$vemTree = vemTree
         }
         # # Tiger
         # # if(!is.element("vemLogTiger", names(simRes))){
         #    scoreTiger = fitHuge(simRes$dataSim$Y, method='tiger')
         #    cat(simLabel, 'Tiger :'); vemTiger = VEM(scoreTiger, K)
         #    cat('\n', simLabel, 'LogTiger :'); vemLogTiger = VEM(log(1+scoreTiger), K)
         #    simRes$scoreTiger = scoreTiger
         #    simRes$vemTiger = vemTiger
         #    simRes$vemLogTiger = vemLogTiger
         # # }
         # Glasso
         if(p < n){
            if(!is.element("vemLogGlasso", names(simRes))){
               scoreGlasso = fitHuge(simRes$dataSim$Y, method='glasso')
               cat('\n', simLabel, 'Glasso :'); vemGlasso = VEM(scoreGlasso, K)
               cat('\n', simLabel, 'LogGlasso :'); vemLogGlasso = VEM(log(1+scoreGlasso), K)
               simRes$scoreGlasso = scoreGlasso
               simRes$vemGlasso = vemGlasso
               simRes$vemLogGlasso = vemLogGlasso
            }
         }
         # Export   
         save(simRes, file=paste0(simDir, simParms, simLabel, '.Rdata'))
         simRes = c()
         cat('\n')
      }
   }
}

