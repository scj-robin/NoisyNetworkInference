# dataSimation of edge scores
rm(list=ls())

library(mclust); library(ROCR)
source('Functions/BasicFunctions.R')
source('Functions/SimulFunctions.R')
source('Functions/VEMFunctions.R')

# Dirs
simDir = '../Simul/'

# Parms
Ktrue = 3; g = 2
load(paste0(simDir, 'simVEM-Ktrue', Ktrue, '-g', g, '-K', K, '.Rdata'))
