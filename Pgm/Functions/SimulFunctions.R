############################################################
############ Fonctions de bases à utiliser dans R ####################
############ pour les simulations ####################
############################################################

##########################################################################
# Simul Z, G, et Y selon un modèle gaussien
SimulZGY <- function(pi, gamma, n, p, connected=TRUE){
   # Blocks
   Z = t(rmultinom(p, 1, pi)); 
   # Graph
   G = vect_mat_low(mat_vect_low(matrix(rbinom(p^2, 1, Z%*%gamma%*%t(Z)), p, p), diag=T), diag=T); 
   if(connected){
      while(is.connected(G, connected='weak')==FALSE){
         G = vect_mat_low(mat_vect_low(matrix(rbinom(p^2, 1, Z%*%gamma%*%t(Z)), p, p), diag=T), diag=T)
      }
   }
   diag(G) = 0
   # Sigma
   lambda = 1; Omega = lambda*diag(colSums(G)) - G
   while(min(eigen(Omega)$values) < 1e-10){lambda = 1.1*lambda; Omega = lambda*diag(colSums(G)) - G}
   Sigma = solve(Omega)
   # Observed
   Y = rmvnorm(n, sigma=Sigma)
   return(list(Z=Z, G=G, Y=Y, Omega=Omega, Sigma=Sigma))
}

##########################################################################
# Recupere les scores du lasso
fitHuge <- function(Y, method='mb'){
   # penalites
   resHuge = huge(Y, method=method, verbose=FALSE); lMin = min(resHuge$lambda); lMax = max(resHuge$lambda)
   while(min(resHuge$sparsity)>0){lMax = 2*lMax; resHuge = huge(Y, lambda=lMax, method=method, verbose=FALSE)}
   while(max(resHuge$sparsity)<1){lMin = lMin/2; resHuge = huge(Y, lambda=lMin, method=method, verbose=FALSE)}
   lambda = exp(seq(log(lMax), log(lMin), length.out=p^2))
   # fit "glasso"
   resHuge = huge(Y, lambda=lambda, method=method)
   # scores   
   S = matrix(0, p, p)
   sapply(1:length(resHuge$lambda), function(l){
      S[which((as.matrix(resHuge$path[[l]])==1) & (S==0))] <<- resHuge$lambda[l] 
   })
   return(S)
}

##########################################################################
# Recupere les scores de EMTree
fitEMtree <- function(Y){
   weightTree = -n/2*log(1-cor(Y)^2)
   scoreTree = EdgeProba(weightTree)
   return(scoreTree)
}
