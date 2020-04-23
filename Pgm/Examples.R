# Analyse des données de Barents & oaks
rm(list=ls()); par(mfrow=c(1, 1), pch=20)

library(sna); library(igraph)
library(blockmodels); library(EMtree); library(mclust); 
library(huge); library(PLNmodels); library(gplots); library(KernSmooth)
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/ResultFunctions.R')
source('Functions/VEMFunctions.R')
source('Functions/VEMmultiFunctions.R')
source('Functions/VEMnparmFunctions.R')

# Dirs
dataDir = '../Data/'
dataName = 'BarentsFish'; Kmax = 6; explorFact=6
dataName = 'oaks'; Kmax = 12; explorFact=3
# dataName = 'oaksFungi'; Kmax = 6
lwd=3; seed = 1; set.seed(seed); alpha = 1/3
plotPdf = FALSE

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
if(plotPdf){pdf(paste0(dataDir, dataName, '-GhatAll.pdf'))}; par(mex=.01)
nodePos = gplot(GhatAll, gmode='graph', label=speciesName, label.pos=5, vertex.col=0, vertex.cex=2)
if(plotPdf){dev.off()}
# Ghat none
GhatNone = 1*(plnNetNone$getBestModel()$model_par$Omega !=0)
if(plotPdf){pdf(paste0(dataDir, dataName, '-GhatNone.pdf'))}; par(mex=.01)
gplot(GhatNone, gmode='graph', label=speciesName, label.pos=5, vertex.col=0, vertex.cex=2, coord=nodePos)
if(plotPdf){dev.off()}
# SBM none
if(!file.exists(paste0(dataDir, dataName, '-SBMNone.Rdata'))){
   SBMNone = BM_bernoulli('SBM_sym', GhatNone)
   SBMNone$estimate(); 
   save(SBMNone, file=paste0(dataDir, dataName, '-SBMNone.Rdata'))
}else{load(paste0(dataDir, dataName, '-SBMNone.Rdata'))}
KNone = which.max(SBMNone$ICL); ZhatNone = apply(SBMNone$memberships[[KNone]]$Z, 1, which.max)
# SBM all
if(plotPdf){pdf(paste0(dataDir, dataName, '-GhatNoneSBM.pdf'))}; par(mex=.01)
gplot(GhatNone, gmode='graph', label=speciesName, label.pos=5, vertex.col=2*ZhatNone, vertex.cex=2, coord=nodePos)
if(plotPdf){dev.off()}
if(!file.exists(paste0(dataDir, dataName, '-SBMAll.Rdata'))){
   SBMAll = BM_bernoulli('SBM_sym', GhatAll)
   SBMAll$estimate(); 
   save(SBMAll, file=paste0(dataDir, dataName, '-SBMAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-SBMAll.Rdata'))}
KAll = which.max(SBMAll$ICL); ZhatAll = apply(SBMAll$memberships[[KAll]]$Z, 1, which.max)
if(plotPdf){pdf(paste0(dataDir, dataName, '-GhatAllSBM.pdf'))}; par(mex=.01)
gplot(GhatAll, gmode='graph', label=speciesName, label.pos=5, vertex.col=2*ZhatAll, vertex.cex=2, coord=nodePos)
if(plotPdf){dev.off()}

###############################################################################
# Scores
# Scores Glasso none
if(!file.exists(paste0(dataDir, dataName, '-vemGlassoNone.Rdata'))){
   scoreGlassoNone = fitHuge(SigmaNone, method='glasso'); 
   vemGlassoNone = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemGlassoNone[[K]] = VEM(S=scoreGlassoNone, K, explorFact= explorFact)}
   save(scoreGlassoNone, vemGlassoNone, file=paste0(dataDir, dataName, '-vemGlassoNone.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemGlassoNone.Rdata'))}

# Scores Tree none
if(!file.exists(paste0(dataDir, dataName, '-vemTreeNone.Rdata'))){
   scoreTreeNone = EdgeProba(-n/2*log(1-cov2cor(SigmaNone)^2)); 
   vemTreeNone = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemTreeNone[[K]] = VEM(S=scoreTreeNone, K, explorFact= explorFact)}
   save(scoreTreeNone, vemTreeNone, file=paste0(dataDir, dataName, '-vemTreeNone.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemTreeNone.Rdata'))}

# Scores MB all
if(!file.exists(paste0(dataDir, dataName, '-vemMBAll.Rdata'))){
   scoreMBAll = fitHuge(SigmaAll, method='mb'); 
   vemMBAll = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemMBAll[[K]] = VEM(S=scoreMBAll, K, explorFact= explorFact)}
   save(scoreMBAll, vemMBAll, file=paste0(dataDir, dataName, '-vemMBAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemMBAll.Rdata'))}

# Scores Glasso all
if(!file.exists(paste0(dataDir, dataName, '-vemGlassoAll.Rdata'))){
   scoreGlassoAll = fitHuge(SigmaAll, method='glasso'); 
   vemGlassoAll = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemGlassoAll[[K]] = VEM(S=scoreGlassoAll, K, explorFact= explorFact)}
   save(scoreGlassoAll, vemGlassoAll, file=paste0(dataDir, dataName, '-vemGlassoAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemGlassoAll.Rdata'))}

# Scores Tree all
if(!file.exists(paste0(dataDir, dataName, '-vemTreeAll.Rdata'))){
   scoreTreeAll = EdgeProba(-n/2*log(1-cov2cor(SigmaAll)^2)); 
   vemTreeAll = list(); for(K in 1:Kmax){cat('\n', K, ':'); vemTreeAll[[K]] = VEM(S=scoreTreeAll, K, explorFact= explorFact)}
   save(scoreTreeAll, vemTreeAll, file=paste0(dataDir, dataName, '-vemTreeAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemTreeAll.Rdata'))}

# Scores logitTree all : nparm
if(!file.exists(paste0(dataDir, dataName, '-vemnpLogitTreeAll.Rdata'))){
   scoreLogitTreeAll = qlogis(scoreTreeAll)
   vemnpLogitTreeAll = list(); Svec = mat_vect_low(scoreLogitTreeAll)
   par(mfrow=c(ceiling(Kmax/2), 4));
   for(K in 1:Kmax){
      cat('\n', K, ':');
      vemnpLogitTreeAll[[K]] = VEMnparm(S=scoreLogitTreeAll, K, niter=200, explorFact=6)
      F_Hist(Svec, vemnpLogitTreeAll[[K]]);
      hist(vemnpLogitTreeAll[[K]]$Psi1, breaks=p, 
           main=round(vemnpLogitTreeAll[[K]]$borne_inf[length(vemnpLogitTreeAll[[K]]$borne_inf)], 2))
      }
   save(scoreLogitTreeAll, vemnpLogitTreeAll, file=paste0(dataDir, dataName, '-vemnpLogitTreeAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemnpLogitTreeAll.Rdata'))}

# Scores Glasso-Tree all
if(!file.exists(paste0(dataDir, dataName, '-vemGlassoTreeAll.Rdata'))){
   scoreArray = array(dim=c(p, p, 2))
   scoreArray[, , 1] = scoreGlassoAll; scoreArray[, , 2] = scoreTreeAll
   vemGlassoTreeAll = list(); 
   for(K in 1:Kmax){cat('\n', K, ':'); vemGlassoTreeAll[[K]] = VEMmulti(S=scoreArray, K, explorFact= explorFact)}
   save(vemGlassoTreeAll, file=paste0(dataDir, dataName, '-vemGlassoTreeAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemGlassoTreeAll.Rdata'))}

# Scores MB-Tree all
if(!file.exists(paste0(dataDir, dataName, '-vemMBTreeAll.Rdata'))){
   scoreArray[, , 1] = scoreMBAll; scoreArray[, , 2] = scoreTreeAll
   vemGlassoTreeAll = list(); 
   for(K in 1:Kmax){cat('\n', K, ':'); vemGlassoTreeAll[[K]] = VEMmulti(S=scoreArray, K, explorFact= explorFact)}
   save(vemGlassoTreeAll, file=paste0(dataDir, dataName, '-vemMBTreeAll.Rdata'))
}else{load(paste0(dataDir, dataName, '-vemMBTreeAll.Rdata'))}

###############################################################################   
# Model selection
par(mfrow=c(3, 2))
critGlassoNone = ModelSelect(vemGlassoNone)
critGlassoAll = ModelSelect(vemGlassoAll); 
critTreeNone = ModelSelect(vemTreeNone)
critTreeAll = ModelSelect(vemTreeAll); 
critnpLogitTreeAll = ModelSelect(vemnpLogitTreeAll); 

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
modelName = 'TreeAll'; scoreModel = scoreTreeAll; vemModel = vemTreeAll; critModel = critTreeAll
# modelName = 'npLogitTreeAll'; scoreModel = scoreLogitTreeAll; vemModel = vemnpLogitTreeAll; critModel = critnpLogitTreeAll

# Criteria
if(plotPdf){pdf(paste0(dataDir, dataName, '-crit', modelName, '.pdf'))}; par(mex=.6, pch=20, cex=2, lwd=2)
plot(critModel$Klist, critModel$lb, type='b', 
     ylim=c(min(critModel$iclZG), max(critModel$lb)), xlab='', ylab=''); 
lines(critModel$Klist, critModel$bic, type='b', col=4); 
lines(critModel$Klist, critModel$iclZ, type='b', col=3); 
lines(critModel$Klist, critModel$iclG, type='b', col=5); 
lines(critModel$Klist, critModel$iclZG, type='b', col=2); 
abline(v = critModel$K, col=2)
if(plotPdf){dev.off()}

# Parms (model = all, method = Tree)
Sigma = SigmaAll; SigmaNet = SigmaNetAll; OmegaNet = OmegaNetAll 
K = critModel$K
pi = vemModel[[K]]$Pi_hat;  gamma = vemModel[[K]]$Gamma_hat; tau = vemModel[[K]]$tau; 
piOrder = order(pi); pi = pi[piOrder]; gamma = gamma[piOrder, piOrder]; tau=tau[, piOrder]
Fvec(round(100*pi, 1)); Ftab(round(100*gamma, 1))

# Species classification
Z = apply(tau, 1, which.max); groupSize = table(Z)
speciesTauOrder = order(tau%*%(1:K)); 
speciesTauRank = rank(tau%*%(1:K)); 
plot((tau%*%(1:K))[speciesTauOrder], xlab='', ylab='', pch=20); 
abline(h=1:K, lty=2) # abline(v=cumsum(groupSize)+.5, lty=2)
if(plotPdf){pdf(paste0(dataDir, dataName, '-nodeClust', modelName, '.pdf'))}; 
plot(tau[speciesTauOrder, 1], type='b', col=1, pch=20, ylim=c(0, 1), xlab='', ylab='')
sapply(2:K, function(q){points(tau[speciesTauOrder, q], type='b', col=q, pch=20)})
if(plotPdf){dev.off()}
classifTable = data.frame(speciesName = speciesName); classifTable$cluster = Z
save(classifTable, file=paste0(dataDir, dataName, '-nodeClust', modelName, '.Rdata'))

# Network 
par(mfrow=c(1, 1))
psi = vemModel[[K]]$Psi1
if(plotPdf){pdf(paste0(dataDir, dataName, '-histPsi', modelName, '.pdf'))}; 
hist(psi, breaks=p, main='', xlab='' ,ylab='', cex.axis=2)
if(plotPdf){dev.off()}
G = vect_mat_low(1*(psi > .5))
if(plotPdf){pdf(paste0(dataDir, dataName, '-graph', modelName, '.pdf'))}; par(mex=.01)
gplot(G, gmode='graph', label=speciesName, label.pos=5, vertex.col=1+Z, vertex.border=ceiling((1+Z)/2), vertex.cex=2, edge.col=8)
if(plotPdf){dev.off()}
if(plotPdf){pdf(paste0(dataDir, dataName, '-circleGraph', modelName, '.pdf'))}; par(mex=.01)
gplot(G[speciesTauOrder,speciesTauOrder], gmode='graph', label=speciesName, label.pos=5, 
      vertex.col=1+Z[speciesTauOrder], vertex.cex=2, coord=circlePos)
if(plotPdf){dev.off()}
if(substr(dataName, 1, 3)=='oak'){
   if(plotPdf){pdf(paste0(dataDir, dataName, '-circleGraph', modelName, '-Ea.pdf'))}; par(mex=.01)
   vertexCol = rep(0, p); vertexCol[eaNum] = 2
   edgeCol = 8*G; edgeCol[eaNum, ] = edgeCol[, eaNum] = G[eaNum, ]
   gplot(G[speciesTauOrder,speciesTauOrder], gmode='graph', label.pos=5, vertex.border=1+Z[speciesTauOrder], 
         vertex.col=vertexCol[speciesTauOrder], vertex.cex=2, 
         edge.col=edgeCol[speciesTauOrder,speciesTauOrder], coord=circlePos)
   if(plotPdf){dev.off()}
}

# Covariance
Omega = solve(Sigma); corSigma = cov2cor(Sigma); corOmega = -cov2cor(Omega); diag(corOmega) = 1
corSigmaNet = cov2cor(SigmaNet); corOmegaNet = -cov2cor(OmegaNet); diag(corOmegaNet) = 1
if(plotPdf){pdf(paste0(dataDir, dataName, '-corr', modelName, '.pdf'))}; par(mex=.3)
image(1:p, 1:p, abs(corSigma[speciesTauOrder, speciesTauOrder])^(1/3), xlab='', ylab='', zlim=c(0, 1)); 
abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
if(substr(dataName, 1, 3)=='oak'){abline(v=speciesTauRank[eaNum]+c(-.5, .5), h=speciesTauRank[eaNum]+c(-.5, .5), col=4)}
if(plotPdf){dev.off()}
if(plotPdf){pdf(paste0(dataDir, dataName, '-partCorr', modelName, '.pdf'))}; par(mex=.3)
image(1:p, 1:p, abs(corOmega[speciesTauOrder, speciesTauOrder])^(1/3), xlab='', ylab='', zlim=c(0, 1)); 
abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
if(substr(dataName, 1, 3)=='oak'){abline(v=speciesTauRank[eaNum]+c(-.5, .5), h=speciesTauRank[eaNum]+c(-.5, .5), col=4)}
if(plotPdf){dev.off()}
if(plotPdf){pdf(paste0(dataDir, dataName, '-corrNet', modelName, '.pdf'))}; par(mex=.3)
image(1:p, 1:p, abs(corSigmaNet[speciesTauOrder, speciesTauOrder])^(1/3), xlab='', ylab='', zlim=c(0, 1)); 
abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
if(substr(dataName, 1, 3)=='oak'){abline(v=speciesTauRank[eaNum]+c(-.5, .5), h=speciesTauRank[eaNum]+c(-.5, .5), col=4)}
if(plotPdf){dev.off()}
if(plotPdf){pdf(paste0(dataDir, dataName, '-partCorrNet', modelName, '.pdf'))}; par(mex=.3)
image(1:p, 1:p, abs(corOmegaNet[speciesTauOrder, speciesTauOrder])^(1/3), xlab='', ylab='', zlim=c(0, 1)); 
abline(v=cumsum(groupSize)+.5, h=cumsum(groupSize)+.5)
if(substr(dataName, 1, 3)=='oak'){abline(v=speciesTauRank[eaNum]+c(-.5, .5), h=speciesTauRank[eaNum]+c(-.5, .5), col=4)}
if(plotPdf){dev.off()}
# image(1:d, 1:p, t(betaAll[speciesTauOrder, ]), xlab='', ylab=''); abline(h=cumsum(groupSize)+.5)

# Score distribution (Gaussian)
s = mat_vect_low(scoreModel)
if(substr(modelName, 1, 2)=='np'){
   p0 = mean(vemModel[[K]]$Psi0); p1 = mean(vemModel[[K]]$Psi1)
   f0 = vemModel[[K]]$fu[, 1]; f1 = vemModel[[K]]$fu[, 2]
   g = p0*f0 + p1*f1
}else{
   phi0 = vemModel[[K]]$phi0; phi1 = vemModel[[K]]$phi1 
   p1 = (pi%*%gamma%*%(pi))[1, 1]; p0 = 1 - p1
   f0 = dnorm(s, mean=phi0[1], sd=sqrt(phi0[2]))
   f1 = dnorm(s, mean=phi1[1], sd=sqrt(phi1[2]))
   g = p0*f0 + p1*f1
}
if(plotPdf){pdf(paste0(dataDir, dataName, '-scoreDist', modelName, '.pdf'))}; par(mex=.5)
H = hist(mat_vect_low(scoreModel), breaks=p, main='', xlab='', ylab='', lwd=2)
points(s, sum(H$counts)*mean(diff(H$breaks))*p0*f0, lwd=4, col=4, pch=20)
points(s, sum(H$counts)*mean(diff(H$breaks))*p1*f1, lwd=4, col=2, pch=20)
points(s, sum(H$counts)*mean(diff(H$breaks))*g, lwd=4, col=1, lty=2, pch=20)
if(plotPdf){dev.off()}
