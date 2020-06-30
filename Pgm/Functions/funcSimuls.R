# Functions for noisy network inference simulations

###############################################################################
# Network simulation
rSimSBM <- function(nbNodes, directed, blockProp, connectParam){
   # Simulates a network according to SBM
   nbBlocks <- length(blockProp)
   membership <- t(rmultinom(nbNodes, 1, blockProp))
   edgeProb <- membership %*% connectParam %*% t(membership)
   network <- matrix(rbinom(nbNodes^2, 1, edgeProb), nbNodes, nbNodes)
   network <- vect2Mat(mat2Vect(network, symmetric=!directed, diag=FALSE), symmetric=!directed, diag=FALSE)
   return(list(membership=membership, network=network))
}

makeVariance <- function(network, signed=FALSE, omegaTol=1e-10){
   # Builds precision and variance matrices given a network
   nbNodes <- nrow(network)
   # Add signs
   if(signed){network <- network * matrix((1-2*rbinom(nbNodes^2, 1, .5)), nbNodes, nbNodes)}
   network <- vect2Mat(mat2Vect(network, symmetric=TRUE, diag=FALSE), symmetric=TRUE, diag=FALSE)
   # Creat positive definite Omega faithfull to the network
   omega <- network; diagOmega <- diag(rowSums(network)); lambda <- 1.1
   omega <- lambda*diagOmega - network
   while(min(Re(eigen(omega)$values)) < omegaTol){
      diagOmega <- lambda*diagOmega 
      omega <- lambda*diagOmega - network
   }
   return(list(sigma=solve(omega), omega=omega))
}

###############################################################################
# Network inference
library('huge')
getScoreHuge <- function(data){
   nbNodes <- ncol(data)
   # Choses lambda penalization grid
   fit <- huge(x=data); lambda <- fit$lambda; logRange <- range(log(fit$lambda))
   while(sum(as.matrix(fit$path[[1]])) > 0){
      logRange[2] <- logRange[2] + diff(logRange)
      lambda <- exp(seq(logRange[1], logRange[2], length.out=length(lambda)))
      fit <- huge(x=data, lambda=lambda)
   }
   while(sum(as.matrix(fit$path[[length(fit$path)]])) < nbNodes*(nbNodes-1)){
      logRange[1] <- logRange[1] - diff(logRange)
      lambda <- exp(seq(logRange[2], logRange[1], length.out=length(lambda)))
      fit <- huge(x=data, lambda=lambda)
   }
   # Final fit
   lambda <- exp(seq(logRange[2], logRange[1], length.out=min(nbNodes^2, 1e3)))
   fit <- huge(x=data, lambda=lambda)
   # Get scores
   score <- matrix(NA, nbNodes, nbNodes)
   for(s in 1:length(fit$path)){
      score[which(is.na(score) & (as.matrix(fit$path[[s]]) !=0))] <- lambda[s]
   }
   return(score)
}

###############################################################################
# Score tranformation
boxCoxFunction <- function(x, lambda){
   geoMean <- exp(mean(log(x)))
   if(lambda>0){
      y <- (x^lambda-1) / (lambda * geoMean^(lambda-1)) 
   }else{
      y <- geoMean * log(x)
   }
}

boxCoxTransform <- function(score, plotit=FALSE){
   scoreVec <- mat2Vect(score, symmetric=TRUE)
   fitBC <- boxcox(scoreVec ~ 1, plotit=plotit); 
   scoreBC <- boxCoxFunction(scoreVec, lambda=fitBC$x[which.max(fitBC$y)])
   return(vect2Mat(scoreBC, symmetric=TRUE))
}

boxCoxGMMTransform <- function(score, nbLambda=20, plotit=FALSE){
   scoreVec <- mat2Vect(score, symmetric=TRUE)
   fitBC <- boxcox(scoreVec ~ 1, plotit=plotit); 
   lambda <- seq(min(fitBC$x), max(fitBC$x), length.out=nbLambda)
   cat('lambda =')
   logLikGMM <- sapply(1:length(lambda), function(l){
      cat(lambda[l], '')
      Mclust(boxCoxFunction(scoreVec, lambda[l]), verbose=FALSE)$loglik
   })
   cat('\n')
   if(plotit){plot(lambda, logLikGMM, type='b'); abline(v=lambda[which.max(logLikGMM)])}
   scoreGMM <- boxCoxFunction(scoreVec, lambda=lambda[which.max(logLikGMM)])
   return(vect2Mat(scoreGMM, symmetric=TRUE))
}
   