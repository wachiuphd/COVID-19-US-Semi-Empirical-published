source("setup_data.R")

chainnum<-3

niter<-5000
runname <- file.path(mcmcpath,"chain5000")
set.seed(1414214)
source(file.path(functionfolder,"Run_chain_fast.R"))

niter<-5000
restartname <- file.path(mcmcpath,"chain5000")
runname <- file.path(mcmcpath,"chain2x5000")
set.seed(1414214+10)
source(file.path(functionfolder,"Restart_chain_fast.R"))

niter<-10000
restartname <- file.path(mcmcpath,"chain2x5000")
runname <- file.path(mcmcpath,"chain4x5000")
set.seed(1414214+20)
source(file.path(functionfolder,"Restart_chain_fast.R"))

