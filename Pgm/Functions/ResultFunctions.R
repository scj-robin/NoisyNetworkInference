############################################################
############ Fonctions de bases à utiliser dans R ####################
############ pour analyser les resultats de simulation ####################
############################################################

xLogX <- function(x){x * log(x + 1*(x==0))}

ModelSelect <- function(vemRes, plot=TRUE){
   Kmax = length(vemRes); p = length(vemRes[[1]]$tau)
   Klist = (1:Kmax); 
   pen = .5*((Klist-1)*log(p) + Klist*(Klist+1)/2*log(p*(p+1)/2))
   entZ = sapply(1:Kmax, function(K){-sum(xLogX(vemRes[[K]]$tau))})
   entG = sapply(1:Kmax, function(K){-sum(xLogX(vemRes[[K]]$Psi0)) -sum(xLogX(vemRes[[K]]$Psi1))})
   lb = unlist(lapply(vemRes, function(v){v$borne_inf[length(v$borne_inf)]}))
   bic = lb - pen; iclZ = bic - entZ; iclG = bic - entG; iclZG = bic - entZ - entG; 
   if(plot){
      plot(lb, type='b', pch=20, ylim=c(min(iclZG), max(lb))); 
      lines(bic, type='b', col=4, pch=20); 
      lines(iclZ, type='b', col=3, pch=20); 
      lines(iclG, type='b', col=5, pch=20); 
      lines(iclZG, type='b', col=2, pch=20); 
      abline(v = which.max(iclZG), col=2)
   }
   return(list(Klist=Klist, lb=lb, entZ=entZ, entG=entG, bic=bic, 
               iclZ=iclZ, iclG=iclG, iclZG=iclZG, K=which.max(iclZG)))
}

# Latex display for a vector
Fvec <- function(a, linebreak=T){
   sapply(1:(length(a)-1), function(i){cat(a[i], '& ')})
   if(linebreak){cat(a[length(a)], '\\\\ \n')}else{cat(a[length(a)], '\n')}
}

# Latex display for a table
Ftab <- function(A){
   sapply(1:(nrow(A)-1), function(i){Fvec(A[i, ])})
   Fvec(A[nrow(A), ], linebreak=F)
}

