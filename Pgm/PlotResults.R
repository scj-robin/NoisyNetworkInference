# Plotting simulation results
rm(list=ls())
source('Functions/ResultFunctions.R')

# Dirs & parms
simDir = '../Simul/'
K = Ktrue = 3; g = 2; signCor = 'neg'
simType = 'WithScores'
# simType = 'WithData'
simName = paste0('simVEM-Ktrue', Ktrue, '-g', g, '-C', signCor, '-K', K)
simName = paste0('sim', simType, '-Ktrue', Ktrue, '-g', g, '-K', K)
figDir = '../Figs/'
plotPdf = FALSE

# Results
Res = read.table(paste0(simDir, simName, '.csv'), sep=',', h=TRUE)
Parms = Res[, 1:3]; Parms$p = as.factor(Parms$p); 
if(simType=='WithScores'){
   Parms$beta = as.factor(Parms$beta)
   ARIvemTab = Res[, 4:7]; AUCvemTab = Res[, 8:11]; 
   ARIoracleTab = Res[, 12:14]; AUCoracleTab = Res[, 15:16]
   iterVemTab = Res[, 17:20]; 
}else{
   Parms$n = as.factor(Parms$n)
   ARIvemTab = Res[, 4:10]; AUCvemTab = Res[, 11:17]; 
   ARIoracleTab = Res[, 18:21]; AUCoracleTab = Res[, 22:24]; 
   iterVemTab = Res[, 25:31]; 
}
names(ARIvemTab) = gsub('ARI.vem', '', names(ARIvemTab)); names(AUCvemTab) = gsub('AUC.vem', '', names(AUCvemTab))
names(ARIoracleTab) = gsub('ARI.', '', names(ARIoracleTab)); names(AUCoracleTab) = gsub('AUC.', '', names(AUCoracleTab))
names(iterVemTab) = gsub('iter.vem', '', names(iterVemTab)); 
names(ARIvemTab); names(AUCvemTab); names(ARIoracleTab); names(AUCoracleTab); names(iterVemTab)

# Method selection
vemSel = (1:ncol(ARIvemTab)); oracleARIsel = (1:ncol(ARIoracleTab)); oracleAUCsel = (1:ncol(AUCoracleTab)); 
# vemSel = c(1, 3, 5, 7); oracleARIsel = c(1, 2, 4); # oracleAUCsel = vemSel
ARIvemTab = ARIvemTab[, vemSel]; AUCvemTab = AUCvemTab[, vemSel]; iterVemTab = iterVemTab[, vemSel]; 
ARIoracleTab = ARIoracleTab[, oracleARIsel]; AUCoracleTab = AUCoracleTab[, oracleAUCsel]; 
names(ARIvemTab); names(AUCvemTab); names(ARIoracleTab); names(AUCoracleTab); names(iterVemTab)

# Vectorization
ARIvemVec = F_Tab2Vec(Parms, ARIvemTab); ARIoracleVec = F_Tab2Vec(Parms, ARIoracleTab)
AUCvemVec = F_Tab2Vec(Parms, AUCvemTab); AUCoracleVec = F_Tab2Vec(Parms, AUCoracleTab)
iterVemVec = F_Tab2Vec(Parms, iterVemTab)

# Plots 
pList = levels(Parms$p); 
lwd=2; cex.axis=2; mex=.4; par(mfrow=c(2, length(pList)))
sapply(pList, function(p){
   F_BoxPlotBeta(ARIvemVec, which(ARIvemVec$p==p), title=paste0('ARI vem:', p))
})
sapply(pList, function(p){
   F_BoxPlotBeta(ARIoracleVec, which(ARIoracleVec$p==p), title=paste0('ARI oracle:', p), borderShift=0)
})
sapply(pList, function(p){
   F_BoxPlotBeta(AUCvemVec, which(AUCvemVec$p==p), title=paste0('AUC vem:', p))
})
sapply(pList, function(p){
   F_BoxPlotBeta(AUCoracleVec, which(AUCoracleVec$p==p), title=paste0('AUC oracle:', p))
})
sapply(pList, function(p){
   F_BoxPlotBeta(iterVemVec, which(iterVemVec$p==p), title=paste0('Iter vem:', p), perfLim = c(0, 100))
})

# Export pdf
if(plotPdf){
   sapply(pList, function(p){
      pdf(paste0(figDir, simName, '-ARIvem-p', p, '.pdf'))
      F_BoxPlotBeta(ARIvemVec, which(ARIvemVec$p==p))
      dev.off()
   })
   sapply(pList, function(p){
      pdf(paste0(figDir, simName, '-ARIoracle-p', p, '.pdf'))
      F_BoxPlotBeta(ARIoracleVec, which(ARIoracleVec$p==p), borderShift=0)
      dev.off()
   })
   sapply(pList, function(p){
      pdf(paste0(figDir, simName, '-AUCvem-p', p, '.pdf'))
      F_BoxPlotBeta(AUCvemVec, which(AUCvemVec$p==p))
      dev.off()
   })
   sapply(pList, function(p){
      pdf(paste0(figDir, simName, '-AUCoracle-p', p, '.pdf'))
      F_BoxPlotBeta(AUCoracleVec, which(AUCoracleVec$p==p))
      dev.off()
   })
}

