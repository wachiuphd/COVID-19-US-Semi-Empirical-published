source("setup_data.R")

if (Fixed_n) {
  runname <- file.path(mcmcpath,"chain2x1000")
  fitparms.list <-mcmc.list(
    list(mcmc(fread(paste0(runname,".1.csv"))[,c(4:7,60:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".2.csv"))[,c(4:7,60:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".3.csv"))[,c(4:7,60:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".4.csv"))[,c(4:7,60:110)][-(1:100),],start=101))
  )
  scaledparms.list <-mcmc.list(
    list(mcmc(fread(paste0(runname,".1.csv"))[,-c(2:7,9:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".2.csv"))[,-c(2:7,9:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".3.csv"))[,-c(2:7,9:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".4.csv"))[,-c(2:7,9:110)][-(1:100),],start=101))
  )
  allparms <-rbind(rbind(rbind(fread(paste0(runname,".1.csv"))[-(1:100),],
                               fread(paste0(runname,".2.csv"))[-(1:100),]),
                         fread(paste0(runname,".3.csv"))[-(1:100),]),
                   fread(paste0(runname,".4.csv"))[-(1:100),])
} else {
  runname <- file.path(mcmcpath,"chain4x5000")
  fitparms.list <-mcmc.list(
    list(mcmc(fread(paste0(runname,".1.csv"))[,c(2:7,9:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".2.csv"))[,c(2:7,9:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".3.csv"))[,c(2:7,9:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".4.csv"))[,c(2:7,9:110)][-(1:100),],start=101))
  )
  scaledparms.list <-mcmc.list(
    list(mcmc(fread(paste0(runname,".1.csv"))[,-c(2:7,9:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".2.csv"))[,-c(2:7,9:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".3.csv"))[,-c(2:7,9:110)][-(1:100),],start=101),
         mcmc(fread(paste0(runname,".4.csv"))[,-c(2:7,9:110)][-(1:100),],start=101))
  )
  allparms <-rbind(rbind(rbind(fread(paste0(runname,".1.csv"))[-(1:100),],
                               fread(paste0(runname,".2.csv"))[-(1:100),]),
                         fread(paste0(runname,".3.csv"))[-(1:100),]),
                   fread(paste0(runname,".4.csv"))[-(1:100),])
}
capture.output(summary(fitparms.list),
               file=file.path(mcmcpath,"FitParameter.Summary.txt"))
pdf(paste0(runname,".fitparms.coda.pdf"),height=8,width=12)
plot(fitparms.list)
dev.off()
capture.output(print(gr.fitparms.diag<-gelman.diag(fitparms.list,autoburnin = FALSE)),
               file=file.path(mcmcpath,"FitParameter.psrf.txt"))
capture.output(summary(scaledparms.list),file=file.path(mcmcpath,"Parameter.Summary.txt"))
pdf(paste0(runname,".scaledparms.coda.pdf"),height=8,width=12)
plot(scaledparms.list)
dev.off()
capture.output(print(gr.scaledparms.diag<-gelman.diag(scaledparms.list,autoburnin = FALSE,multivariate = FALSE)),
               file=file.path(mcmcpath,"Parameter.psrf.txt"))

allparms <- as.data.frame(allparms)
allparms$ncoef.0 <- inv.logit(allparms$mu_ncoef)
allparms$Roffset.0 <- inv.logit(allparms$mu_Roffset)*100
set.seed(123456)
nsamp<-2000
sampparms <- allparms[sample(nrow(allparms),size=nsamp),]
fwrite(sampparms,file.path(mcmcpath,"sampled_parameters.csv"))

### WAIC
lnlikesamps <- as.matrix(allparms[,grep("lnlike.",names(allparms))])
llpd <- sum(log(apply(exp(lnlikesamps),2,mean)))
pwaic <- sum(apply(lnlikesamps,2,var))
waic <- data.frame(mcmcpath=mcmcpath,
                   llpd=llpd,pWAIC=pwaic,WAIC=-2*(llpd - pwaic))
print(waic)

### update posterior predictions
source(file.path(functionfolder,"update_posterior_predictions_fast.R"))
predquant<-fread(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
predquant<-as.data.frame(predquant)
predquant$date <- as.Date(predquant$date)
predquant <- subset(predquant,date < today)

fitdat<-merge.data.frame(serodat,predquant)
fitresid<-log(fitdat$Data.50.)-log(fitdat$SP.pct.50.)
fitmad <- mad(fitresid)
sigerr.med <- sqrt(median(allparms$sig2err))
res<-lm(log(Data.50.)~log(SP.pct.50.),data=fitdat)
ressum<-summary(res)
adjr2 <- signif(ressum$adj.r.squared,2)

fitdat.val<-merge.data.frame(cdcseroval,predquant)
fitdat.val<- subset(fitdat.val,Data.50.>0)
fitresid.val<-log(fitdat.val$Data.50.)-log(fitdat.val$SP.pct.50.)
fitmad.val <- mad(fitresid.val)
sigerr.med.val <- sd(fitresid.val)
res.val<-lm(log(Data.50.)~log(SP.pct.50.),data=fitdat.val)
ressum.val<-summary(res.val)
adjr2.val <- signif(ressum.val$adj.r.squared,2)

### Extended SEIR predictions
predfile<-file.path(datafolder,"SEIR-ModelPredictions.csv")
seirpred<-fread(predfile)
seirpred$model <- "Extended-SEIR"
seirpred$date<-as.Date(seirpred$date)
### Imperial predictions
imperialpredfile<-file.path(datafolder,"ImperialPredictions.csv")
imperialpred<-fread(imperialpredfile)
imperialpred$model <- "Imperial"
imperialpred$date<-as.Date(imperialpred$date)
imperialpred <- subset(imperialpred,N_I_tot.50.>0)

### comparisons
predquant.1 <- subset(predquant,date<=max(imperialpred$date,seirpred$date) &
                        Itot.pct.50. > 0)
seirscat.df <- merge.data.frame(seirpred[,c("state","date","N_I_tot.2.5.",
                                            "N_I_tot.50.","N_I_tot.97.5.","N_I_tot.Random")],
                                predquant.1[,c("state","date","Itot.pct.50.")])
seirresid <- log(seirscat.df$N_I_tot.50.)-log(seirscat.df$Itot.pct.50.)
seirmad <- mad(seirresid)
seirerr <- sd(seirresid)
res.seir<-lm(log(N_I_tot.50.)~log(Itot.pct.50.),data=seirscat.df)
ressum.seir <-summary(res.seir)
adjr2.seir <- signif(ressum.seir$adj.r.squared,2)

imperialscat.df <- merge.data.frame(imperialpred[,c("state","date","N_I_tot.2.5.","N_I_tot.50.","N_I_tot.97.5.","N_I_tot.Random")],
                                    predquant.1[,c("state","date","Itot.pct.50.")])
imperialresid <- log(imperialscat.df$N_I_tot.50.)-log(imperialscat.df$Itot.pct.50.)
imperialmad <- mad(imperialresid)
imperialerr <- sd(imperialresid)
res.imperial<-lm(log(N_I_tot.50.)~log(Itot.pct.50.),data=imperialscat.df)
ressum.imperial <-summary(res.imperial)
adjr2.imperial <- signif(ressum.imperial$adj.r.squared,2)

fitstats <- cbind(waic,data.frame(
  fitmad=fitmad,
  sigerr.med=sigerr.med,
  adjr2=adjr2,
  fitmad.val = fitmad.val,
  adjr2.val=adjr2.val,
  seirmad=seirmad,
  seirerr=seirerr,
  adjr2.seir=adjr2.seir,
  imperialmad=imperialmad,
  imperialerr=imperialerr,
  adjr2.imperial=adjr2.imperial
))

print(fitstats)
fwrite(fitstats,file.path(mcmcpath,"FitStats.csv"))

