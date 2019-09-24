# Simulation of edge scores
rm(list=ls()); par(pch=20)

library(sna); library(mvtnorm); library(huge); library(blockmodels); library(EMtree); library(mclust); 
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')

# Dirs
simDir = '../Simul/'

# Dims : p = nb nodes, n = nb replicates
# p = 20, 30, 50, 80, n = 20, 50, 100
p = 30; n = 50; K = 3; g = 2; simNb = 50
gLassoMethod = 'mb'
simName = paste0('dataSimVEM-n', n, '-p', p, '-K', K, '-g', g, '-', gLassoMethod)

# Parms
pi = (1:K); pi = pi / sum(pi)
rho = log(p)/p
gamma = (pi^g)[K:1]%o%(pi^g)[K:1]; gamma = rho*gamma/(t(pi)%*%gamma%*%pi)[1, 1]
gamma[which(gamma>1)] = 1
parms = list(pi=pi, gamma=gamma, rho=rho)

# Boucle de simulation
if(file.exists(paste0(simDir, simName, '.Rdata'))){
   load(paste0(simDir, simName, '.Rdata')); simInit = length(simRes)+1
}else{
   simRes = list(); simInit = 1
}
for(sim in simInit:simNb){
   simLabel = paste0('n', n, '-p', p, '-s', sim)
   # Simul data
   dataSim = SimulZGY(pi, gamma, n, p)
   # Oracle : SBM sur G
   sbmG = BM_bernoulli('SBM_sym', dataSim$G, plotting='', verbosity=0, ncores=1); 
   sbmG$estimate()
   # Scores
   scoreGlasso = fitGlasso(dataSim$Y, method=gLassoMethod)
   scoreTree = fitEMtree(dataSim$Y)
   # VEM
   par(mfrow=c(2, 2), cex=.6)
   cat(simLabel, 'Glasso :'); vemGlasso = VEM(scoreGlasso, K); #plot(vemGlasso$borne_inf)
   cat('\n', simLabel, 'LogGlasso :'); vemLogGlasso = VEM(log(1+scoreGlasso), K); #plot(vemLogGlasso$borne_inf[-1])
   cat('\n', simLabel, 'Tree:'); vemTree = VEM(scoreTree, K); #plot(vemTree$borne_inf[-1])
   # logitTree supprime : resultats trop mauvais
   # cat('\nLogitTree ', sim, ':'); vemLogitTree = VEM(qlogis(scoreTree), K); ; plot(vemLogitTree$borne_inf[-1])
   vemLogitTree = c()
   # Export   
   simRes[[sim]] = list(parms=parms, dataSim=dataSim, sbmG=sbmG, 
                        scoreGlasso=scoreGlasso, scoreTree=scoreTree, 
                        vemGlasso=vemGlasso, vemTree=vemTree, 
                        vemLogGlasso=vemLogGlasso, vemLogitTree=vemLogitTree)
   save(simRes, file=paste0(simDir, simName, '.Rdata'))
}
