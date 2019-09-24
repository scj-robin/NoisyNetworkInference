# dataSimation of edge scores
rm(list=ls())

library(sna); library(mvtnorm); library(huge); library(blockmodels); library(EMtree); library(mclust); 
# library(MASS)
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')

# Dirs
simDir = '../Simul/'

# Dims
p = 50; n = 100; K = 3; g = 2; simNb = 100
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
   # Simul data
   dataSim = SimulZGY(pi, gamma, n, p)
   # Oracle : SBM sur G
   sbmG = BM_bernoulli('SBM_sym', dataSim$G, plotting='', verbosity=0); sbmG$estimate()
   # Scores
   scoreGlasso = fitGlasso(dataSim$Y, method=gLassoMethod)
   scoreTree = fitEMtree(dataSim$Y)
   # VEM
   cat('Glasso ', sim, ':'); vemGlasso = VEM(scoreGlasso, K)
   cat('\nLogGlasso ', sim, ':'); vemLogGlasso = VEM(log(1+scoreGlasso), K)
   cat('\nTree ', sim, ':'); vemTree = VEM(scoreTree, K)
   cat('\nLogitTree ', sim, ':'); vemLogitTree = VEM(qlogis(scoreTree), K)
   # Export   
   simRes[[sim]] = list(parms=parms, dataSim=dataSim, sbmG=sbmG, 
                        scoreGlasso=scoreGlasso, scoreTree=scoreTree, 
                        vemGlasso=vemGlasso, vemTree=vemTree, 
                        vemLogGlasso=vemLogGlasso, vemLogitTree=vemLogitTree)
   save(simRes, file=paste0(simDir, simName, '.Rdata'))
}
