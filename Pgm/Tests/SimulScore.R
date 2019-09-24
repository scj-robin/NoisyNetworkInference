# Simulation of edge scores
rm(list=ls()); par(pch=20)

library(sna); library(mvtnorm); library(huge); library(blockmodels); library(EMtree); library(mclust); library(saturnin)
# library(MASS)
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')

# Dims
p = 30; n = 100; K = 3
gLassoMethod = 'mb'

# Parms
pi = (1:K); pi = pi / sum(pi)
rho = log(p)/p
gamma = (pi^2)[K:1]%o%(pi^2)[K:1]; gamma = rho*gamma/(t(pi)%*%gamma%*%pi)[1, 1]

# Simul data
simul = SimulZGY(pi, gamma, n, p)
Zvec = simul$Z%*%(1:K); Gvec = mat_vect_low(simul$G)
# gplot(G, gmode='graph', vertex.col=Zvec); plot(Sigma, cov(Y)); abline(0, 1)

# Oracle : SBM sur G
sbm = BM_bernoulli('SBM_sym', simul$G); sbm$estimate()
tauOracle = sbm$memberships[[K]]$Z; Zoracle = apply(tauOracle, 1, which.max)

# Glasso
scoreGlasso = fitGlasso(simul$Y)
# boxplot(mat_vect_low(log(1+S)) ~ Gvec)
par(mfrow=c(2, 2))
hist(mat_vect_low(scoreGlasso), breaks=p); plot(density(mat_vect_low(scoreGlasso)))
hist(log(mat_vect_low(scoreGlasso)), breaks=p); plot(density(log(mat_vect_low(scoreGlasso))))
boxplot(mat_vect_low(scoreGlasso) ~ Gvec); boxplot(log(mat_vect_low(scoreGlasso)) ~ Gvec); 

# EMtree
scoreTree = fitEMtree(simul$Y)
par(mfrow=c(2, 2))
hist(mat_vect_low(scoreTree), breaks=p); plot(density(mat_vect_low(scoreTree)))
hist(qlogis(mat_vect_low(scoreTree)), breaks=p); plot(density(qlogis(mat_vect_low(scoreTree))))
boxplot(mat_vect_low(scoreTree) ~ Gvec); boxplot(qlogis(mat_vect_low(scoreTree)) ~ Gvec); 

# VEM
vemGlasso = VEM(scoreGlasso, K)
# boxplot(res$Psi1 ~ Gvec)
Zinit = apply(res$tau_init, 1, which.max); Zhat = apply(res$tau, 1, which.max)
adjustedRandIndex(Zvec, Zoracle); adjustedRandIndex(Zvec, Zinit); adjustedRandIndex(Zvec, Zhat); 
