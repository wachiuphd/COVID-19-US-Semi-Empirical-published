source("setup_data.R")

chainnum<-1

niter<-5000
runname <- file.path(mcmcpath,"chain5000")
set.seed(314159)
source(file.path(functionfolder,"Run_chain_fast.R"))

niter<-5000
restartname <- file.path(mcmcpath,"chain5000")
runname <- file.path(mcmcpath,"chain2x5000")
set.seed(314159+10)
source(file.path(functionfolder,"Restart_chain_fast.R"))

niter<-10000
restartname <- file.path(mcmcpath,"chain2x5000")
runname <- file.path(mcmcpath,"chain4x5000")
set.seed(314159+20)
source(file.path(functionfolder,"Restart_chain_fast.R"))

