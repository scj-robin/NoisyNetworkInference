# Simulation of multiple edge scores + VEM inference
rm(list=ls()); #par(pch=20)

library(sna); library(mvtnorm); library(huge); library(blockmodels); library(EMtree); library(mclust); 
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')
source('Functions/VEMmultiFunctions.R')

# Dirs
simDir = '../Simul/'

# Dims : p = nb nodes, n = nb replicates
p = 30; n = 50; Ktrue = 3; g = 2

# # Parms
# pi = (1:Ktrue); pi = pi / sum(pi)
# rho = 1.5*log(p)/p
# gamma = (pi^g)[Ktrue:1]%o%(pi^g)[Ktrue:1]; gamma = rho*gamma/(t(pi)%*%gamma%*%pi)[1, 1]
# gamma[which(gamma>1)] = 1
# parms = list(pi=pi, gamma=gamma, rho=rho)
# 
# # Data & scores
# dataSim = SimulZGY(pi, gamma, n, p)
# scoreMB = fitHuge(dataSim$Y, method='mb')
# scoreTree = fitEMtree(dataSim$Y)
# scoreGlasso = fitHuge(dataSim$Y, method='glasso')
# save(parms, dataSim, scoreMB, scoreTree, scoreGlasso, file='SimMulti.Rdata')
load('SimMulti.Rdata')

load('../Simul//simVEM-p20-n20-Ktrue3-g2/simul-K3-s1.Rdata')
scoreMB = simRes$scoreMB; scoreGlasso = simRes$scoreGlasso; scoreTree = simRes$scoreTree

# Scores
plot(as.data.frame(cbind(log(mat_vect_low(scoreMB)), log(mat_vect_low(scoreGlasso)), 
                   qlogis(mat_vect_low(scoreTree)))), col=1+mat_vect_low(dataSim$G), pch=20)
plot(as.data.frame(cbind((mat_vect_low(scoreMB)), (mat_vect_low(scoreGlasso)),
                   (mat_vect_low(scoreTree)))), col=1+mat_vect_low(dataSim$G), pch=20)
Smat = cbind(mat_vect_low(scoreMB), mat_vect_low(scoreTree))
# Smat = cbind(mat_vect_low(scoreMB), mat_vect_low(scoreGlasso), mat_vect_low(scoreTree))

# Multi VEM
K = Ktrue
vemMulti = VEMmulti(Smat, K)
