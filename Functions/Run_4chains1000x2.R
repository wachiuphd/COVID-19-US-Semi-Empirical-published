source("setup_data.R")

niter<-1000
runname <- file.path(mcmcpath,"chain1000")

set.seed(314159)
chainnum<-1
source(file.path(functionfolder,"Run_chain_fast.R"))

set.seed(2718282)
chainnum<-2
source(file.path(functionfolder,"Run_chain_fast.R"))

set.seed(1414214)
chainnum<-3
source(file.path(functionfolder,"Run_chain_fast.R"))

set.seed(1732051)
chainnum<-4
source(file.path(functionfolder,"Run_chain_fast.R"))

## load into mcmc.list
if (Fixed_n) {
  allparms.list <-mcmc.list(
    list(mcmc(fread(paste0(runname,".1.csv"))[,c(4:7,60:110)]),
         mcmc(fread(paste0(runname,".2.csv"))[,c(4:7,60:110)]),
         mcmc(fread(paste0(runname,".3.csv"))[,c(4:7,60:110)]),
         mcmc(fread(paste0(runname,".4.csv"))[,c(4:7,60:110)]))
  )
} else {
  allparms.list <-mcmc.list(
    list(mcmc(fread(paste0(runname,".1.csv"))[,c(2:7,9:110)]),
         mcmc(fread(paste0(runname,".2.csv"))[,c(2:7,9:110)]),
         mcmc(fread(paste0(runname,".3.csv"))[,c(2:7,9:110)]),
         mcmc(fread(paste0(runname,".4.csv"))[,c(2:7,9:110)]))
  )
}
## Gelman 
print(gr.diag<-gelman.diag(allparms.list))

restartname <- file.path(mcmcpath,"chain1000")
runname <- file.path(mcmcpath,"chain2x1000")

set.seed(314159+10)
chainnum<-1
source(file.path(functionfolder,"Restart_chain_fast.R"))

set.seed(2718282+10)
chainnum<-2
source(file.path(functionfolder,"Restart_chain_fast.R"))

set.seed(1414214+10)
chainnum<-3
source(file.path(functionfolder,"Restart_chain_fast.R"))

set.seed(1732051+10)
chainnum<-4
source(file.path(functionfolder,"Restart_chain_fast.R"))

## load into mcmc.list
if (Fixed_n) {
  allparms.list <-mcmc.list(
    list(mcmc(fread(paste0(runname,".1.csv"))[,c(4:7,60:110)][-(1:200),]),
         mcmc(fread(paste0(runname,".2.csv"))[,c(4:7,60:110)][-(1:200),]),
         mcmc(fread(paste0(runname,".3.csv"))[,c(4:7,60:110)][-(1:200),]),
         mcmc(fread(paste0(runname,".4.csv"))[,c(4:7,60:110)][-(1:200),]))
  )
} else {
  allparms.list <-mcmc.list(
    list(mcmc(fread(paste0(runname,".1.csv"))[,c(2:7,9:110)][-(1:200),]),
         mcmc(fread(paste0(runname,".2.csv"))[,c(2:7,9:110)][-(1:200),]),
         mcmc(fread(paste0(runname,".3.csv"))[,c(2:7,9:110)][-(1:200),]),
         mcmc(fread(paste0(runname,".4.csv"))[,c(2:7,9:110)][-(1:200),]))
  )
}
## Gelman 
print(gr.diag<-gelman.diag(allparms.list,autoburnin = FALSE))

source(file.path(functionfolder,"summarize_mcmc_results.R"))

