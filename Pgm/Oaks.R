# Analyse des données de Barents
rm(list=ls()); par(pch=20)

library(sna); library(igraph)
library(blockmodels); library(EMtree); library(mclust); 
library(huge); library(PLNmodels); library(MRFcov)
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')

# Dirs
dataDir = '../Data/'

# Parms
Kmax = 4

# Donnees
dataName = 'oaks'
load(paste0(dataDir, dataName, '.Rdata'))
n = nrow(Data$count); p = ncol(Data$count)

# # MRFcov
# X = glm(Data$count[, 1] ~ as.factor(Data$covariates$tree) + Data$covariates$distTOground + Data$covariates$distTOtrunk, family='poisson', x=TRUE)$x[, -1]
# d = ncol(X)
# YX = cbind(Data$count, X)
# mrf = MRFcov(YX, symmetrise='min', prep_covariates=TRUE, n_nodes=p, n_covariates=d, family='poisson')
# mrf = MRFcov(Data$count, symmetrise='min', prep_covariates=TRUE, n_nodes=p, n_covariates=0, family='poisson')
# gplot(mrf$graph, gmode='graph')
# mrf$direct_coefs
# bmrf = bootstrap_MRF(Data$count, n_bootstrap=1e1, symmetrise='min', n_nodes=p, family='poisson')

# PLN
# pln = PLN(Data$count ~ as.factor(Data$covariates$tree) + Data$covariates$distTOground + 
             # Data$covariates$distTOtrunk)
# save(pln, file=paste0(dataDir, dataName, '-pln.Rdata'))
load(paste0(dataDir, dataName, '-pln.Rdata'))
Sigma = pln$model_par$Sigma
Omega = solve(Sigma)
Beta = t(pln$model_par$Theta); d = nrow(Beta)

# PLN net
# plnnet = PLNnetwork(Data$count ~ as.factor(Data$covariates$tree) + Data$covariates$distTOground + 
#                        Data$covariates$distTOtrunk)
# save(plnnet, file=paste0(dataDir, dataName, '-plnnet.Rdata'))
load(paste0(dataDir, dataName, '-plnnet.Rdata'))
plnnet$plot()
plnnetBest = plnnet$getBestModel()
Omega = plnnetBest$model_par$Omega
Sigma = plnnetBest$model_par$Sigma
Beta = t(plnnetBest$model_par$Theta); d = nrow(Beta)

# Scores MB
# scoreMB = fitHuge(Sigma, method='mb'); hist(scoreMB, breaks=p)
# vemMB = list(); for(K in 1:Kmax){vemMB[[K]] = VEM(S=scoreMB, K)}
# save(scoreMB, vemMB, file=paste0(dataDir, dataName, '-vemMB.Rdata'))
load(paste0(dataDir, dataName, '-vemMB.Rdata'))

# Scores Tree
# weightTree = -n/2*log(1-cov2cor(Sigma)^2)
# scoreTree = EdgeProba(weightTree); hist(scoreTree, breaks=p)
# vemTree = list(); for(K in 1:Kmax){vemTree[[K]] = VEM(S=scoreTree, K)}
# save(scoreTree, vemTree, file=paste0(dataDir, dataName, '-vemTree.Rdata'))
load(paste0(dataDir, dataName, '-vemTree.Rdata'))

# Choose K
lbMB = unlist(lapply(vemMB, function(v){v$borne_inf[length(v$borne_inf)]}))
lbTree = unlist(lapply(vemTree, function(v){v$borne_inf[length(v$borne_inf)]}))
plot(lbMB, type='b'); plot(lbTree, type='b')
Kopt = which.max(lbTree); abline(v=Kopt, col=2)

# Network
Psi = vect_mat_low(vemTree[[Kopt]]$Psi1); hist(Psi, breaks=p)
Ghat = 1*(Psi > .5)

# Node classif
Tau = vemTree[[Kopt]]$tau; Zhat = apply(Tau, 1, which.max); groupSize = table(Zhat)
vemTree[[Kopt]]$Pi_hat; vemTree[[Kopt]]$Gamma_hat

# Network plot
gplot(Ghat, gmode='graph', vertex.col=Zhat)
plot(graph_from_adjacency_matrix(Ghat, mode='undirected'))

# Covariance
nodeOrder = order(Tau%*%(1:Kopt))
image(1:p, 1:p, cov2cor(Sigma[nodeOrder, nodeOrder])); abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
image(1:p, 1:p, cov2cor(Omega[nodeOrder, nodeOrder])); abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
image(2:d, 1:p, Beta[-1, nodeOrder]); abline(h=cumsum(groupSize)+.5)
