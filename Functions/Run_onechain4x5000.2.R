source("setup_data.R")

chainnum<-2

niter<-5000
runname <- file.path(mcmcpath,"chain5000")
set.seed(2718282)
source(file.path(functionfolder,"Run_chain_fast.R"))

niter<-5000
restartname <- file.path(mcmcpath,"chain5000")
runname <- file.path(mcmcpath,"chain2x5000")
set.seed(2718282+10)
source(file.path(functionfolder,"Restart_chain_fast.R"))

niter<-10000
restartname <- file.path(mcmcpath,"chain2x5000")
runname <- file.path(mcmcpath,"chain4x5000")
set.seed(2718282+20)
source(file.path(functionfolder,"Restart_chain_fast.R"))

