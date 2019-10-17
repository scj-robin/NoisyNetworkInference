# Analyse des données de Barents & oaks
rm(list=ls()); par(mfrow=c(1, 1), pch=20)

library(sna); library(igraph)
library(blockmodels); library(EMtree); library(mclust); 
library(huge); library(PLNmodels); library(gplots)
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/ResultFunctions.R')
source('Functions/VEMFunctions.R')
source('Functions/VEMmultiFunctions.R')

# Dirs
dataDir = '../Data/'
dataName = 'BarentsFish'; Kmax = 6
dataName = 'oaks'; Kmax = 12
# dataName = 'oaksFungi'; Kmax = 6
lwd=3; seed = 1; set.seed(seed); alpha = 1/3

###############################################################################
# Exemples
###############################################################################
load(paste0(dataDir, dataName, '.Rdata'))
n = nrow(Data$count); p = ncol(Data$count)
circlePos = 2*(1:p)*acos(-1)/p; circlePos = cbind(cos(circlePos), sin(circlePos))

# Donnees
Y = Data$count; O = log(Data$offset);  speciesName = colnames(Y)
if(dataName=='BarentsFish'){
   O = 0*O; 
   X = scale(Data$covariate)
   speciesName = as.vector(sapply(speciesName, function(s){gsub('_', '.', s)}))
   }
if((dataName=='oaks') | (dataName=='oaksFungi')){X = glm(Y[, 1] ~ as.factor(Data$covariates$tree) + 
                                Data$covariates$distTOground + Data$covariates$distTOtrunk, 
                             family='poisson', x=TRUE)$x[, -1]
                     X[, 3:4] = scale(X[, 3:4])
                     speciesName = as.vector(sapply(speciesName, function(s){gsub('_', '', s)}))
                     eaNum = Data$numEA
                     }
d = ncol(X)
speciesSel = (6:9)
Fvec(speciesName[speciesSel])
Ftab(Y[1:5, speciesSel])

###############################################################################
# # MRFcov
###############################################################################
# YX = cbind(Y, X)
# mrf = MRFcov(YX, symmetrise='min', prep_covariates=TRUE, n_nodes=p, n_covariates=d, family='poisson')
# mrf = MRFcov(Data$count, symmetrise='min', prep_covariates=TRUE, n_nodes=p, n_covariates=0, family='poisson')
# gplot(mrf$graph, gmode='graph')
# mrf$direct_coefs
# bmrf = bootstrap_MRF(Data$count, n_bootstrap=1e1, symmetrise='min', n_nodes=p, family='poisson')

###############################################################################
# PLN models
###############################################################################
# PLN : no covariate
if(!file.exists(paste0(dataDir, dataName, '-plnNone.Rdata'))){
   plnNone = PLN(Y ~ 1)
   save(plnNone, file=paste0(dataDir, dataName, '-plnNone.Rdata'))
}else{load(paste0(dataDir, dataName, '-plnNone.Rdata'))}
SigmaNone = plnNone$model_par$Sigma
BetaNone = t(plnNone$model_par$Theta)

# PLN : all covariate
if(!file.exists(paste0(dataDir, dataName, '-plnAll.Rdata'))){
   plnAll = PLN(Y ~ X)
   save(plnAll, file=paste0(dataDir, dataName, '-plnAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-plnAll.Rdata'))}
SigmaAll = plnAll$model_par$Sigma
BetaAll = t(plnAll$model_par$Theta)

# PLN net : no covariates
if(!file.exists(paste0(dataDir, dataName, '-plnNetNone.Rdata'))){
   plnNetNone = PLNnetwork(Y ~ 1 + offset(O))
   logLambda = log(plnNetNone$penalties); lMin = min(logLambda); lMax = max(logLambda)
   lambda = exp(seq(lMax, 2*lMin-lMax, length.out=5*p))
   plnNetNone = PLNnetwork(Y ~ 1 + offset(O), penalties=lambda)
   save(plnNetNone, file=paste0(dataDir, dataName, '-plnNetNone.Rdata'))
   }else{load(paste0(dataDir, dataName, '-plnNetNone.Rdata'))}
# plnNetNone$plot()
plnNetBestNone = plnNetNone$getBestModel()
SigmaNetNone = plnNetBestNone$model_par$Sigma
OmegaNetNone = plnNetBestNone$model_par$Omega
BetaNetNone = t(plnNetBestNone$model_par$Theta)

# PLN net : all covariates
if(!file.exists(paste0(dataDir, dataName, '-plnNetAll.Rdata'))){
   plnNetAll = PLNnetwork(Y ~ X + offset(O))
   logLambda = log(plnNetAll$penalties); lMin = min(logLambda); lMax = max(logLambda)
   lambda = exp(seq(lMax, 2*lMin-lMax, length.out=5*p))
   plnNetAll = PLNnetwork(Y ~ X + offset(O), penalties=lambda)
   save(plnNetAll, file=paste0(dataDir, dataName, '-plnNetAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-plnNetAll.Rdata'))}
# plnNetAll$plot()
plnNetBestAll = plnNetAll$getBestModel()
SigmaNetAll = plnNetBestAll$model_par$Sigma
OmegaNetAll = plnNetBestAll$model_par$Omega
BetaNetAll = t(plnNetBestAll$model_par$Theta)

###############################################################################
# Bad practice : 
###############################################################################
# Ghat all
GhatAll = 1*(plnNetAll$getBestModel()$model_par$Omega !=0)
pdf(paste0(dataDir, dataName, '-GhatAll.pdf')); par(mex=.01)
nodePos = gplot(GhatAll, gmode='graph', label=speciesName, label.pos=5, vertex.col=0, vertex.cex=2)
dev.off()
# Ghat none
GhatNone = 1*(plnNetNone$getBestModel()$model_par$Omega !=0)
pdf(paste0(dataDir, dataName, '-GhatNone.pdf')); par(mex=.01)
gplot(GhatNone, gmode='graph', label=speciesName, label.pos=5, vertex.col=0, vertex.cex=2, coord=nodePos)
dev.off()
# SBM none
if(!file.exists(paste0(dataDir, dataName, '-SBMNone.Rdata'))){
   SBMNone = BM_bernoulli('SBM_sym', GhatNone)
   SBMNone$estimate(); 
   save(SBMNone, file=paste0(dataDir, dataName, '-SBMNone.Rdata'))
}else{load(paste0(dataDir, dataName, '-SBMNone.Rdata'))}
KNone = which.max(SBMNone$ICL); ZhatNone = apply(SBMNone$memberships[[KNone]]$Z, 1, which.max)
# SBM all
pdf(paste0(dataDir, dataName, '-GhatNoneSBM.pdf')); par(mex=.01)
gplot(GhatNone, gmode='graph', label=speciesName, label.pos=5, vertex.col=2*ZhatNone, vertex.cex=2, coord=nodePos)
dev.off()
if(!file.exists(paste0(dataDir, dataName, '-SBMAll.Rdata'))){
   SBMAll = BM_bernoulli('SBM_sym', GhatAll)
   SBMAll$estimate(); 
   save(SBMAll, file=paste0(dataDir, dataName, '-SBMAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-SBMAll.Rdata'))}
KAll = which.max(SBMAll$ICL); ZhatAll = apply(SBMAll$memberships[[KAll]]$Z, 1, which.max)
pdf(paste0(dataDir, dataName, '-GhatAllSBM.pdf')); par(mex=.01)
gplot(GhatAll, gmode='graph', label=speciesName, label.pos=5, vertex.col=2*ZhatAll, vertex.cex=2, coord=nodePos)
dev.off()

###############################################################################
# Scores
# Scores Glasso none
if(!file.exists(paste0(dataDir, dataName, '-vemGlassoNone.Rdata'))){
   scoreGlassoNone = fitHuge(SigmaNone, method='glasso'); 
   vemGlassoNone = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemGlassoNone[[K]] = VEM(S=scoreGlassoNone, K, explorFact=3)}
   save(scoreGlassoNone, vemGlassoNone, file=paste0(dataDir, dataName, '-vemGlassoNone.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemGlassoNone.Rdata'))}

# Scores Tree none
if(!file.exists(paste0(dataDir, dataName, '-vemTreeNone.Rdata'))){
   scoreTreeNone = EdgeProba(-n/2*log(1-cov2cor(SigmaNone)^2)); 
   vemTreeNone = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemTreeNone[[K]] = VEM(S=scoreTreeNone, K, explorFact=3)}
   save(scoreTreeNone, vemTreeNone, file=paste0(dataDir, dataName, '-vemTreeNone.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemTreeNone.Rdata'))}

# Scores Glasso all
if(!file.exists(paste0(dataDir, dataName, '-vemGlassoAll.Rdata'))){
   scoreGlassoAll = fitHuge(SigmaAll, method='glasso'); 
   vemGlassoAll = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemGlassoAll[[K]] = VEM(S=scoreGlassoAll, K, explorFact=3)}
   save(scoreGlassoAll, vemGlassoAll, file=paste0(dataDir, dataName, '-vemGlassoAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemGlassoAll.Rdata'))}

# Scores Tree all
if(!file.exists(paste0(dataDir, dataName, '-vemTreeAll.Rdata'))){
   scoreTreeAll = EdgeProba(-n/2*log(1-cov2cor(SigmaAll)^2)); 
   vemTreeAll = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemTreeAll[[K]] = VEM(S=scoreTreeAll, K, explorFact=3)}
   save(scoreTreeAll, vemTreeAll, file=paste0(dataDir, dataName, '-vemTreeAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemTreeAll.Rdata'))}

# Scores Glasso-Tree all
if(!file.exists(paste0(dataDir, dataName, '-vemGlassoTreeAll.Rdata'))){
   scoreMat = cbind(mat_vect_low(scoreGlassoAll), mat_vect_low(scoreTreeAll))
   vemGlassoTreeAll = list(); 
   for(K in 1:Kmax){cat('\n', K, ':'); vemGlassoTreeAll[[K]] = VEMmulti(S=scoreMat, K, explorFact=3)}
   save(vemGlassoTreeAll, file=paste0(dataDir, dataName, '-vemGlassoTreeAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemGlassoTreeAll.Rdata'))}

###############################################################################   
# Model selection
par(mfrow=c(2, 2))
critGlassoNone = ModelSelect(vemGlassoNone)
critGlassoAll = ModelSelect(vemGlassoAll); 
critTreeNone = ModelSelect(vemTreeNone)
critTreeAll = ModelSelect(vemTreeAll); 
pdf(paste0(dataDir, dataName, '-critTreeAll.pdf')); par(mex=.6, pch=20, cex=2, lwd=2)
plot(critTreeAll$Klist, critTreeAll$lb, type='b', 
     ylim=c(min(critTreeAll$iclZG), max(critTreeAll$lb)), xlab='', ylab=''); 
lines(critTreeAll$Klist, critTreeAll$bic, type='b', col=4); 
lines(critTreeAll$Klist, critTreeAll$iclZ, type='b', col=3); 
lines(critTreeAll$Klist, critTreeAll$iclG, type='b', col=5); 
lines(critTreeAll$Klist, critTreeAll$iclZG, type='b', col=2); 
abline(v = critTreeAll$K, col=2)
dev.off()

###############################################################################   
# Effect of the covariates
betaAll = plnAll$model_par$Theta[, -1]
plnAllClus = heatmap.2(t(betaAll), dendrogram="col", key=F, Rowv=F, trace='none')
speciesBetaOrder = plnAllClus$colInd
image(vemGlassoNone[[critGlassoNone$K]]$tau[speciesBetaOrder, ])

# Comparing classification (model = all)
round(t(vemGlassoAll[[critGlassoAll$K]]$tau) %*% vemTreeAll[[critTreeAll$K]]$tau, 1)

###############################################################################   
# Results from model = all, method = Tree
# Parms (model = all, method = Tree)
Sigma = SigmaAll; SigmaNet = SigmaNetAll; OmegaNet = OmegaNetAll 
K = critTreeAll$K
pi = vemTreeAll[[K]]$Pi_hat;  gamma = vemTreeAll[[K]]$Gamma_hat; tau = vemTreeAll[[K]]$tau; 
piOrder = order(pi); pi = pi[piOrder]; gamma = gamma[piOrder, piOrder]; tau=tau[, piOrder]
Fvec(round(100*pi, 1)); Ftab(round(100*gamma, 1))

# Species classification
Z = apply(tau, 1, which.max); groupSize = table(Z)
speciesTauOrder = order(tau%*%(1:K)); 
speciesTauRank = rank(tau%*%(1:K)); 
plot((tau%*%(1:K))[speciesTauOrder], xlab='', ylab='', pch=20); 
abline(h=1:K, lty=2) # abline(v=cumsum(groupSize)+.5, lty=2)
pdf(paste0(dataDir, dataName, '-nodeClustTreeAll.pdf')); 
plot(tau[speciesTauOrder, 1], type='b', col=1, pch=20, ylim=c(0, 1), xlab='', ylab='')
sapply(2:K, function(q){points(tau[speciesTauOrder, q], type='b', col=q, pch=20)})
dev.off()
classifTable = data.frame(speciesName = speciesName); classifTable$cluster = Z
save(classifTable, file=paste0(dataDir, dataName, '-nodeClustTreeAll.Rdata'))

# Network 
par(mfrow=c(1, 1))
psi = vemTreeAll[[K]]$Psi1
pdf(paste0(dataDir, dataName, '-histPsiTreeAll.pdf')); 
hist(psi, breaks=p, main='', xlab='' ,ylab='', cex.axis=2)
dev.off()
G = vect_mat_low(1*(psi > .5))
pdf(paste0(dataDir, dataName, '-graphTreeAll.pdf')); par(mex=.01)
gplot(G, gmode='graph', label=speciesName, label.pos=5, vertex.col=1+Z, vertex.border=ceiling((1+Z)/2), vertex.cex=2, edge.col=8)
dev.off()
?pdf(paste0(dataDir, dataName, '-circleGraphTreeAll.pdf')); par(mex=.01)
gplot(G[speciesTauOrder,speciesTauOrder], gmode='graph', label=speciesName, label.pos=5, 
      vertex.col=1+Z[speciesTauOrder], vertex.cex=2, coord=circlePos)
dev.off()
if(substr(dataName, 1, 3)=='oak'){
   pdf(paste0(dataDir, dataName, '-circleGraphTreeAll-Ea.pdf')); par(mex=.01)
   vertexCol = rep(0, p); vertexCol[eaNum] = 2
   edgeCol = 8*G; edgeCol[eaNum, ] = edgeCol[, eaNum] = G[eaNum, ]
   gplot(G[speciesTauOrder,speciesTauOrder], gmode='graph', label.pos=5, vertex.border=1+Z[speciesTauOrder], 
         vertex.col=vertexCol[speciesTauOrder], vertex.cex=2, 
         edge.col=edgeCol[speciesTauOrder,speciesTauOrder], coord=circlePos)
   dev.off()
}

# Covariance
Omega = solve(Sigma); corSigma = cov2cor(Sigma); corOmega = -cov2cor(Omega); diag(corOmega) = 1
corSigmaNet = cov2cor(SigmaNet); corOmegaNet = -cov2cor(OmegaNet); diag(corOmegaNet) = 1
pdf(paste0(dataDir, dataName, '-corrTreeAll.pdf')); par(mex=.3)
image(1:p, 1:p, abs(corSigma[speciesTauOrder, speciesTauOrder])^(1/3), xlab='', ylab='', zlim=c(0, 1)); 
abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
if(substr(dataName, 1, 3)=='oak'){abline(v=speciesTauRank[eaNum]+c(-.5, .5), h=speciesTauRank[eaNum]+c(-.5, .5), col=4)}
dev.off()
pdf(paste0(dataDir, dataName, '-partCorrTreeAll.pdf')); par(mex=.3)
image(1:p, 1:p, abs(corOmega[speciesTauOrder, speciesTauOrder])^(1/3), xlab='', ylab='', zlim=c(0, 1)); 
abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
if(substr(dataName, 1, 3)=='oak'){abline(v=speciesTauRank[eaNum]+c(-.5, .5), h=speciesTauRank[eaNum]+c(-.5, .5), col=4)}
dev.off()
pdf(paste0(dataDir, dataName, '-corrNetTreeAll.pdf')); par(mex=.3)
image(1:p, 1:p, abs(corSigmaNet[speciesTauOrder, speciesTauOrder])^(1/3), xlab='', ylab='', zlim=c(0, 1)); 
abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
if(substr(dataName, 1, 3)=='oak'){abline(v=speciesTauRank[eaNum]+c(-.5, .5), h=speciesTauRank[eaNum]+c(-.5, .5), col=4)}
dev.off()
pdf(paste0(dataDir, dataName, '-partCorrNetTreeAll.pdf')); par(mex=.3)
image(1:p, 1:p, abs(corOmegaNet[speciesTauOrder, speciesTauOrder])^(1/3), xlab='', ylab='', zlim=c(0, 1)); 
abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
if(substr(dataName, 1, 3)=='oak'){abline(v=speciesTauRank[eaNum]+c(-.5, .5), h=speciesTauRank[eaNum]+c(-.5, .5), col=4)}
dev.off()
# image(1:d, 1:p, t(betaAll[speciesTauOrder, ]), xlab='', ylab=''); abline(h=cumsum(groupSize)+.5)

# Score distribution
phi0 = vemTreeAll[[K]]$phi0; phi1 = vemTreeAll[[K]]$phi1 
p1 = (pi%*%gamma%*%(pi))[1, 1]; p0 = 1 - p1
s = seq(min(scoreTreeAll), max(scoreTreeAll), length.out=500)
f0 = dnorm(s, mean=phi0[1], sd=sqrt(phi0[2]))
f1 = dnorm(s, mean=phi1[1], sd=sqrt(phi1[2]))
g = p0*f0 + p1*f1
pdf(paste0(dataDir, dataName, '-scoreDistTreeAll.pdf')); par(mex=.5)
H = hist(mat_vect_low(scoreTreeAll), breaks=p, main='', xlab='', ylab='', lwd=2)
lines(s, sum(H$counts)*mean(diff(H$breaks))*p0*f0, lwd=4, col=4)
lines(s, sum(H$counts)*mean(diff(H$breaks))*p1*f1, lwd=4, col=2)
lines(s, sum(H$counts)*mean(diff(H$breaks))*g, lwd=4, col=1, lty=2)
dev.off()




# scoreGlassoNone = fitHuge(SigmaNone, method='glasso'); 
# vemGlassoNone = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemGlassoNone[[K]] = VEM(S=scoreGlassoNone, K, explorFact=3)}
# save(scoreGlassoNone, vemGlassoNone, file=paste0(dataDir, dataName, '-vemGlassoNone.Rdata'))
# 
# 
# scoreTreeNone = EdgeProba(-n/2*log(1-cov2cor(SigmaNone)^2)); 
# vemTreeNone = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemTreeNone[[K]] = VEM(S=scoreTreeNone, K, explorFact=3)}
# save(scoreTreeNone, vemTreeNone, file=paste0(dataDir, dataName, '-vemTreeNone.Rdata'))
# 
# scoreGlassoAll = fitHuge(SigmaAll, method='glasso'); 
# vemGlassoAll = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemGlassoAll[[K]] = VEM(S=scoreGlassoAll, K, explorFact=3)}
# save(scoreGlassoAll, vemGlassoAll, file=paste0(dataDir, dataName, '-vemGlassoAll.Rdata'))
