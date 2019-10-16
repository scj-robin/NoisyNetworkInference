# Plotting simulation results
rm(list=ls())

# Dirs & parms
simDir = '../Simul/'
K = Ktrue = 3; g = 2
simName = paste0('simVEM-Ktrue', Ktrue, '-g', g, '-K', K)
plot = FALSE

# Results
Res = read.table(paste0(simDir, simName, '.csv'), sep=',', h=TRUE)
Parms = Res[, 1:3]; Parms$p = as.factor(Parms$p); Parms$n = as.factor(Parms$n)
ARItab = Res[, 4:7]; names(ARItab) = gsub('ARI.', '', names(ARItab))
AUCtab = Res[, 8:13]; names(AUCtab) = gsub('AUC.', '', names(AUCtab))

# Vectorization
ARIvec = c(); AUCvec = c(); 
invisible(sapply(1:ncol(ARItab), function(c){
   Method = rep(names(ARItab)[c], nrow(ARItab)); ARI = ARItab[, c]
   ARIvec <<- rbind(ARIvec, cbind(Parms, Method, ARI))
   }))
invisible(sapply(1:ncol(AUCtab), function(c){
   Method = rep(names(AUCtab)[c], nrow(AUCtab)); AUC = AUCtab[, c]
   AUCvec <<- rbind(AUCvec, cbind(Parms, Method, AUC))
}))

# Plots
lwd=2; cex.axis=2; mex=.4
ARIlim = c(min(ARIvec$ARI[which(ARIvec$Method!='Oracle')]), max(ARIvec$ARI[which(ARIvec$Method!='Oracle')]))
AUClim = c(min(AUCvec$AUC[which(substr(AUCvec$Method, 1, 3) == 'vem')]), 1)
pList = levels(Parms$p)
par(mfrow=c(2, length(pList)))
sapply(pList, function(p){
   # ARI
   selARI = which((ARIvec$p==p) & (ARIvec$Method!='Oracle'))
   colMethod = rep(2:4, each=nlevels(ARIvec$n))
   if(plot){pdf(paste0(simDir, simName, '-ARI-p', p, '.pdf')); par(mex=mex)}
   boxplot(ARI ~ n + Method, data=ARIvec, subset=selARI, drop=TRUE, border=colMethod,
           xlab='', ylab='', lwd=lwd, cex.axis=cex.axis, ylim=ARIlim); 
   abline(v = .5+nlevels(ARIvec$n)*(1:2))
   if(plot){dev.off()}
   # AUC
   selAUC = which((AUCvec$p==p) & (substr(AUCvec$Method, 1, 3) == 'vem'))
   colMethod = rep(2:4, each=nlevels(AUCvec$n))
   if(plot){pdf(paste0(simDir, simName, '-AUC-p', p, '.pdf')); par(mex=mex)}
   boxplot(AUC ~ n + Method, data=AUCvec, subset=selAUC, drop=TRUE, border=colMethod, 
           xlab='', ylab='', lwd=lwd, cex.axis=cex.axis, ylim=AUClim); 
   abline(v = .5+nlevels(AUCvec$n)*(1:2))
   if(plot){dev.off()}
})
