source(file.path("International","setup_intl_data.R"))

mcmcpath <- "MCMC.n0.5"
sampparms<-fread(file.path(mcmcpath,"sampled_parameters.csv"))
sampparms<-as.data.frame(sampparms)
nsamp<-200

#### n=0.5

preddat.samp<-data.frame()
for (j.samp in 1:nsamp) {
  if ((j.samp %% 10)==1) print(j.samp)
  Tinf.vec.samp<-sampparms$Tinf[j.samp]
  ncoef.vec.samp<-inv.logit(sampparms$mu_ncoef[j.samp])
  Roffset.vec.samp<-100*inv.logit(sampparms$mu_Roffset[j.samp])
  preddat.tmp <- get_all_post_pred(ncoef.vec.samp,
                                   Roffset.vec.samp,
                                   Tinf.vec.samp)[,c("t","numDate","location","id","I.pct","SP.pct","dI.pct","Itot.pct")]
  preddat.tmp$j.samp <- j.samp
  preddat.samp <- rbind(preddat.samp,preddat.tmp)
}
preddat.samp$date <- as.Date(preddat.samp$t,origin=datezero)
predquant.I.pct <- aggregate(I.pct ~ date + location,
                             data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.I.pct<-do.call(data.frame,predquant.I.pct)
predquant.Itot.pct <- aggregate(Itot.pct ~ date + location,
                             data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.Itot.pct<-do.call(data.frame,predquant.Itot.pct)
predquant.SP.pct <- aggregate(SP.pct ~ date + location,
                              data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.SP.pct<-do.call(data.frame,predquant.SP.pct)
predquant.dI.pct <- aggregate(dI.pct ~ date + location,
                              data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.dI.pct<-do.call(data.frame,predquant.dI.pct)

predquant<-merge.data.frame(merge.data.frame(merge.data.frame(predquant.I.pct,predquant.SP.pct),predquant.dI.pct),
predquant.Itot.pct)
fwrite(predquant,file.path(intlpath,"Intl_n0.5_Posterior_I_SP_dI_Quantiles.csv"))


mcmcpath <- "MCMC"
sampparms<-fread(file.path(mcmcpath,"sampled_parameters.csv"))
sampparms<-as.data.frame(sampparms)
nsamp<-200

### Fixed effects no offset

preddat.samp<-data.frame()
for (j.samp in 1:nsamp) {
  if ((j.samp %% 10)==1) print(j.samp)
  Tinf.vec.samp<-sampparms$Tinf[j.samp]
  ncoef.vec.samp<-inv.logit(sampparms$mu_ncoef[j.samp])
  Roffset.vec.samp<-0
  preddat.tmp <- get_all_post_pred(ncoef.vec.samp,
                                   Roffset.vec.samp,
                                   Tinf.vec.samp)[,c("t","numDate","location","id","I.pct","SP.pct","dI.pct","Itot.pct")]
  preddat.tmp$j.samp <- j.samp
  preddat.samp <- rbind(preddat.samp,preddat.tmp)
}
preddat.samp$date <- as.Date(preddat.samp$t,origin=datezero)
predquant.I.pct <- aggregate(I.pct ~ date + location,
                             data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.I.pct<-do.call(data.frame,predquant.I.pct)
predquant.Itot.pct <- aggregate(Itot.pct ~ date + location,
                                data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.Itot.pct<-do.call(data.frame,predquant.Itot.pct)
predquant.SP.pct <- aggregate(SP.pct ~ date + location,
                              data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.SP.pct<-do.call(data.frame,predquant.SP.pct)
predquant.dI.pct <- aggregate(dI.pct ~ date + location,
                              data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.dI.pct<-do.call(data.frame,predquant.dI.pct)

predquant<-merge.data.frame(merge.data.frame(merge.data.frame(predquant.I.pct,predquant.SP.pct),predquant.dI.pct),
                            predquant.Itot.pct)
fwrite(predquant,file.path(intlpath,
                           "Intl_FixedEffects_NoOffset_Posterior_I_SP_dI_Quantiles.csv"))

#### Fixed effects

preddat.samp<-data.frame()
for (j.samp in 1:nsamp) {
  if ((j.samp %% 10)==1) print(j.samp)
  Tinf.vec.samp<-sampparms$Tinf[j.samp]
  ncoef.vec.samp<-inv.logit(sampparms$mu_ncoef[j.samp])
  Roffset.vec.samp<-100*inv.logit(sampparms$mu_Roffset[j.samp])
  preddat.tmp <- get_all_post_pred(ncoef.vec.samp,
                                   Roffset.vec.samp,
                                   Tinf.vec.samp)[,c("t","numDate","location","id","I.pct","SP.pct","dI.pct","Itot.pct")]
  preddat.tmp$j.samp <- j.samp
  preddat.samp <- rbind(preddat.samp,preddat.tmp)
}
preddat.samp$date <- as.Date(preddat.samp$t,origin=datezero)
predquant.I.pct <- aggregate(I.pct ~ date + location,
                             data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.I.pct<-do.call(data.frame,predquant.I.pct)
predquant.Itot.pct <- aggregate(Itot.pct ~ date + location,
                                data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.Itot.pct<-do.call(data.frame,predquant.Itot.pct)
predquant.SP.pct <- aggregate(SP.pct ~ date + location,
                              data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.SP.pct<-do.call(data.frame,predquant.SP.pct)
predquant.dI.pct <- aggregate(dI.pct ~ date + location,
                              data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.dI.pct<-do.call(data.frame,predquant.dI.pct)

predquant<-merge.data.frame(merge.data.frame(merge.data.frame(predquant.I.pct,predquant.SP.pct),predquant.dI.pct),
                            predquant.Itot.pct)
fwrite(predquant,file.path(intlpath,"Intl_FixedEffects_Posterior_I_SP_dI_Quantiles.csv"))

### Random effects
randomstate <- data.frame(
  ncoef=inv.logit(rnorm(nsamp,m=sampparms$mu_ncoef,sd=exp(sampparms$logSigma_ncoef))),
  Roffset=100*inv.logit(rnorm(nsamp,m=sampparms$mu_Roffset,sd=exp(sampparms$logSigma_Roffset))),
  Tinf=sampparms$Tinf,
  sig2err=sampparms$sig2err
)

# randomstate2 <- randomstate
# idrand <- sample.int(n.id, nsamp, replace = TRUE)
# for (j in 1:nsamp) {
#   randomstate2$ncoef[j] <- sampparms[j,paste0("ncoef.",idrand[j])]
#   randomstate2$Roffset[j] <- sampparms[j,paste0("Roffset.",idrand[j])]
# }

preddat.samp<-data.frame()
for (j.samp in 1:nsamp) {
  if ((j.samp %% 10)==1) print(j.samp)
  Tinf.vec.samp<-randomstate$Tinf[j.samp]
  ncoef.vec.samp<-randomstate$ncoef[j.samp]
  Roffset.vec.samp<-randomstate$Roffset[j.samp]
  preddat.tmp <- get_all_post_pred(ncoef.vec.samp,
                                   Roffset.vec.samp,
                                   Tinf.vec.samp)[,c("t","numDate","location","id","I.pct","SP.pct","dI.pct","Itot.pct")]
  preddat.tmp$j.samp <- j.samp
  preddat.samp <- rbind(preddat.samp,preddat.tmp)
}
preddat.samp$date <- as.Date(preddat.samp$t,origin=datezero)
predquant.I.pct <- aggregate(I.pct ~ date + location,
                             data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.I.pct<-do.call(data.frame,predquant.I.pct)
predquant.Itot.pct <- aggregate(Itot.pct ~ date + location,
                                data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.Itot.pct<-do.call(data.frame,predquant.Itot.pct)
predquant.SP.pct <- aggregate(SP.pct ~ date + location,
                              data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.SP.pct<-do.call(data.frame,predquant.SP.pct)
predquant.dI.pct <- aggregate(dI.pct ~ date + location,
                              data=preddat.samp,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))
predquant.dI.pct<-do.call(data.frame,predquant.dI.pct)

predquant<-merge.data.frame(merge.data.frame(merge.data.frame(predquant.I.pct,predquant.SP.pct),predquant.dI.pct),
                            predquant.Itot.pct)
fwrite(predquant,file.path(intlpath,"Intl_Posterior_I_SP_dI_Quantiles.csv"))