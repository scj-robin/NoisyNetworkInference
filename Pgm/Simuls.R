# Simulation of edge scores + VEM inference
rm(list=ls()); #par(pch=20)

library(sna); library(mvtnorm); library(huge); library(blockmodels); library(EMtree); library(mclust); 
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')
source('Functions/VEMmultiFunctions.R')

# Dirs
simDir = '../Simul/'

# Dims : p = nb nodes, n = nb replicates
pList = c(20, 30, 50, 80); nList = c(20, 50, 100, 200)
# pList = c(20, 30, 50, 80); nList = c(20, 50, 100)
pList = c(20, 30, 50, 80); nList = c(20, 50, 100, 200)
Ktrue = 3; g = 2; simNb = 50

for(p in pList){
   # p = 20
   for(n in nList){
      # n = 50
      simParms = paste0('simVEM-p', p, '-n', n, '-Ktrue', Ktrue, '-g', g, '/')
      
      # Parms
      pi = (1:Ktrue); pi = pi / sum(pi)
      rho = 1.5*log(p)/p
      gamma = (pi^g)[Ktrue:1]%o%(pi^g)[Ktrue:1]; gamma = rho*gamma/(t(pi)%*%gamma%*%pi)[1, 1]
      gamma[which(gamma>1)] = 1
      parms = list(pi=pi, gamma=gamma, rho=rho)
      
      # Boucle de simulation
      K = Ktrue
      for(sim in 1:simNb){
         # sim = 1
         # Deja fait ?
         cat('Parms ')
         simLabel = paste0('simul-K', K, '-s', sim); print(paste0(simParms, simLabel))
         if(file.exists(paste0(simDir, simParms, simLabel, '.Rdata'))){
            load(paste0(simDir, simParms, simLabel, '.Rdata')); 
         }else{simRes = list(parms=parms)}
         # Simul data
         cat('data ')
         if(!is.element("dataSim", names(simRes))){
            dataSim = SimulZGY(pi, gamma, n, p)
            # if(p <= 50){dataSim = SimulZGY(pi, gamma, n, p)}
            # if(p > 50){dataSim = SimulZGY(pi, gamma, n, p, connected=FALSE)}
            simRes$dataSim = dataSim
         }
         # Oracle : SBM sur G
         cat('oracle ')
         if(!is.element("sbmG", names(simRes))){
            sbmG = BM_bernoulli('SBM_sym', simRes$dataSim$G, plotting='', verbosity=0, ncores=1); 
            sbmG$estimate()
            simRes$sbmG = sbmG
         }
         # MB
         cat('MB ')
         if(!is.element("vemMB", names(simRes))){
            scoreMB = fitHuge(simRes$dataSim$Y, method='mb')
            cat(simLabel, 'MB :'); vemMB = VEM(scoreMB, K); 
            # cat('\n', simLabel, 'LogMB :'); vemLogMB = VEM(log(1+scoreMB), K); 
            simRes$scoreMB = scoreMB
            simRes$vemMB = vemMB
            # simRes$vemLogMB = vemLogMB
         }
         # Tree
         cat('Tree ')
         if(!is.element("vemTree", names(simRes))){
            scoreTree = fitEMtree(simRes$dataSim$Y)
            cat('\n', simLabel, 'Tree:'); vemTree = VEM(scoreTree, K)
            simRes$scoreTree = scoreTree
            simRes$vemTree = vemTree
         }
         cat('MB-Tree ')
         if(!is.element("vemMBTree", names(simRes))){
            Smat = cbind(mat_vect_low(simRes$scoreMB), mat_vect_low(simRes$scoreTree))
            cat('\n', simLabel, 'MB-Tree:'); vemMBTree = VEMmulti(Smat, K)
            simRes$vemMBTree = vemMBTree
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
            cat('Glasso ')
            if(!is.element("vemGlasso", names(simRes))){
               scoreGlasso = fitHuge(simRes$dataSim$Y, method='glasso')
               cat('\n', simLabel, 'Glasso :'); vemGlasso = VEM(scoreGlasso, K)
               # cat('\n', simLabel, 'LogGlasso :'); vemLogGlasso = VEM(log(1+scoreGlasso), K)
               simRes$scoreGlasso = scoreGlasso
               simRes$vemGlasso = vemGlasso
               # simRes$vemLogGlasso = vemLogGlasso
            }
            cat('MB-Glasso-Tree ')
            if(!is.element("vemMBGlassoTree", names(simRes))){
               Smat = cbind(mat_vect_low(simRes$scoreMB), mat_vect_low(simRes$scoreGlasso), 
                            mat_vect_low(simRes$scoreGlasso))
               cat('\n', simLabel, 'MB-Glasso-Tree:'); vemMBGlassoTree = VEMmulti(Smat, K)
               simRes$vemMBGlassoTree = vemMBGlassoTree
            }
         }
         cat('sbmMB ')
         if(!is.element("sbmMB", names(simRes))){
            # sbmMB : SBM on MB-inferred network
            selMB = huge.select(huge(simRes$dataSim$Y, method='mb'))
            graphMB = as.matrix(huge(simRes$dataSim$Y, method='mb', lambda=selMB$opt.lambda)$path[[1]])
            sbmMB = BM_bernoulli('SBM_sym', graphMB, plotting='', verbosity=0, ncores=1); 
            sbmMB$estimate()
            simRes$sbmMB = sbmMB
         }
         if(p < n){
            cat('sbmGlasso ')
            if(!is.element("sbmGlasso", names(simRes))){
               # sbmGlasso : SBM on glasso-inferred network
               selGlasso = huge.select(huge(simRes$dataSim$Y, method='glasso'))
               graphGlasso = 1*(selGlasso$opt.icov!=0)
               if(sum(graphGlasso)>p){
                  sbmGlasso = BM_bernoulli('SBM_sym', graphGlasso, plotting='', verbosity=0, ncores=1); 
                  sbmGlasso$estimate()
                  simRes$sbmGlasso = sbmGlasso
               }else{simRes$sbmGlasso = 'failed'}
            }
         }
         cat('sbmTree ')
         if(!is.element("sbmTree", names(simRes))){
            # sbmTree : SBM on Tree-inferred network
            sbmTree = BM_bernoulli('SBM_sym', 1*(simRes$scoreTree > 2/p), plotting='', verbosity=0, ncores=1); 
            sbmTree$estimate()
            simRes$sbmTree = sbmTree
         }
         cat('sbmScore ')
         if(!is.element("sbmScore", names(simRes))){
            # sbmScore : SBM directly on Gaussian scores
            sbmScore = BM_gaussian('SBM_sym', simRes$scoreTree, plotting='', verbosity=0, ncores=1); 
            sbmScore$estimate()
            simRes$sbmScore = sbmScore
         }
         # Export   
         save(simRes, file=paste0(simDir, simParms, simLabel, '.Rdata'))
         simRes = c()
         cat('\n')
      }
   }
}

