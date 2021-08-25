library(data.table)
library(tidyr)
library(boot)
library(coda)
datezero <- "2019-12-31"
maxdatadate <- "2020-12-31"
figfolder<-"Figures"
tabfolder<-"Tables"
mcmcpath <- "MCMC"
intlpath <- "International"
runname <- file.path(mcmcpath,"chain4x5000")

### Seroprevalence
intl_serofile<-file.path(intlpath,"InternationalSeroprevalenceData.csv")
intl_serodat<-fread(intl_serofile)
intl_serodat$Output_Var <- "Seroprevalence"
intl_serodat$date <- as.Date((as.numeric(as.Date(intl_serodat$mindate))+
                                     as.numeric(as.Date(intl_serodat$maxdate)))/2 -
                                    as.numeric(as.Date(datezero)),
                                  origin=datezero)
locationsvec <- unique(intl_serodat$location)
# If zero, make CI symmetric
iszero<-intl_serodat$Data.2.5.==0
intl_serodat$Data.2.5.[iszero] <- intl_serodat$Data.50.[iszero]/
  (intl_serodat$Data.97.5.[iszero]/intl_serodat$Data.50.[iszero])

## Generate data files
source(file.path(intlpath,"generate_intl_data_functions.R"))
generate_intl_testpositivity(locationsvec = locationsvec)
## Load data files
### Test positivity - England
englfile<-file.path(intlpath,"EnglandTestPositivityData.csv")
engldat<-fread(englfile)
engldat$date <- as.Date(engldat$date)
engldat<-subset(engldat,date <= maxdatadate)
engldat$Cases_14_pct <- 100*engldat$Cases_14/engldat$population

### Test positivity - owid
owidfile<-file.path(intlpath,"OWIDTestPositivityData.csv")
owiddat<-fread(owidfile)
owiddat$date <- as.Date(owiddat$date)
owiddat<-subset(owiddat,date <= maxdatadate)
owiddat$Cases_14_pct <- 100*owiddat$Cases_14/owiddat$population

### Combined
combnames <- intersect(names(engldat),names(owiddat))
intl_dat <- rbind(as.data.frame(engldat)[,combnames],
                       as.data.frame(owiddat)[,combnames])


# ## Global data frames
pred.lag <- 7
n.id <- length(unique(intl_serodat$location))
# ### COVID testing data
testdat<-intl_dat
testdat$id <- as.numeric(factor(testdat$location))
testdat$date <- testdat$date-pred.lag
testdat$t <- testdat$numDate-pred.lag
testdat <- testdat[order(testdat$id,testdat$t),]
testdat$I.pct<-0
testdat$SP.pct <- 0
### Serological data
intl_serodat <- subset(intl_serodat, !is.na(intl_serodat$Data.50.))
intl_serodat.df <- data.frame(id=as.numeric(factor(intl_serodat$location)),
                                   t=floor(as.numeric(as.Date(intl_serodat$date))-
                                             as.numeric(as.Date(datezero))),
                                   y=log(intl_serodat$Data.50.),
                                   y.se=log(intl_serodat$Data.97.5./intl_serodat$Data.50.)/(qnorm(0.975)))
intl_serodat.df <- intl_serodat.df[order(intl_serodat.df$id,intl_serodat.df$t,intl_serodat.df$y),]

# get_one_post_pred <- function(idnow,ncoef,Roffset,Tinf,sig2err) {
#   # Only select current id
#   indx <- testdat$id==idnow
#   dattmp <- testdat[indx,]
#   ## Calculate predictions
#   I.pct <- (dattmp$PositivePct_14)^(1-ncoef)*(dattmp$Cases_14_pct)^ncoef
#   I.pct[is.na(I.pct)]<-0
#   SP <- Roffset+c(0,cumsum(I.pct[-length(I.pct)]))/Tinf
#   dattmp$I.pct<-I.pct
#   dattmp$SP.pct<-SP
#   return(dattmp)
# }
# 
# get_one_fitpred <- function(idnow,ncoef,Roffset,Tinf,sig2err) {
#   
#   dattmp <- get_one_pred(idnow,ncoef,Roffset,Tinf,sig2err)
#   ## Merge
#   fitdat<-merge.data.frame(subset(serodat.df,id==idnow),
#                            dattmp[,c("date","t","state","id","fips","I.pct","SP.pct")])
#   # Reorder to match original data
#   fitdat <- fitdat[order(fitdat$id,fitdat$t,fitdat$y),]
#   fitdat$ypred <- log(fitdat$SP.pct)
#   return(fitdat)
# }
# 
# get_one_likelihood <- function(idnow,ncoef,Roffset,Tinf,sig2err) {
#   fitdat <- get_one_fitpred(idnow,ncoef,Roffset,Tinf,sig2err)
#   fitdat$lnlike <- dnorm(fitdat$y,mean=fitdat$ypred,sd=sqrt(fitdat$y.se^2+sig2err),log=TRUE)
#   lnlike <- sum(fitdat$lnlike)
#   return(lnlike)
# }

get_all_post_pred <- function(ncoef.vec=0.5,
                              Roffset.vec=0,
                              Tinf.vec=14,
                              sig2err=0.1,
                              addCumulPos = 1,
                              Trec.vec = Tinf.vec) {
  ## Calculate predictions
  dattmp.all <- testdat
  n.id <- max(dattmp.all$id)
  if (length(ncoef.vec)==1) ncoef.vec <- rep(ncoef.vec,n.id)
  if (length(Roffset.vec)==1) Roffset.vec <- rep(Roffset.vec,n.id)
  if (length(Tinf.vec)==1) Tinf.vec <- rep(Tinf.vec,n.id)
  if (length(Trec.vec)==1) Trec.vec <- rep(Trec.vec,n.id)
  for (idnow in unique(dattmp.all$id)) {
    ncoef <- ncoef.vec[idnow]
    Roffset <- Roffset.vec[idnow]
    Tinf <- Tinf.vec[idnow]
    Trec <- Trec.vec[idnow]
    indx <- dattmp.all$id==idnow
    dattmp <- dattmp.all[indx,]
    I.pct <- (dattmp$PositivePct_14)^(1-ncoef)*(dattmp$Cases_14_pct)^ncoef
    I.pct[is.na(I.pct)]<-0
    ## Correction for diagnosed cases
    I.pct.minus.Cases <- I.pct - dattmp$Cases_14_pct
    I.pct.minus.Cases[is.na(I.pct.minus.Cases)] <- 0
    I.pct.minus.Cases[(I.pct.minus.Cases < 0)] <- 0
    ## Cumulative diagnosed 
    dTcumul <- pred.lag+ceiling(Tinf)
    CumulPosPct.1 <- c(rep(0,dTcumul),dattmp$CumulPosPct[1:(nrow(dattmp)-dTcumul)])
    # SP <- Roffset+c(0,cumsum(I.pct[-length(I.pct)]))/Tinf+
    #   addCumulPos*CumulPosPct.1
    gamma.vec <- I.pct/I.pct.minus.Cases
    Teff <- gamma.vec*(1-gamma.vec^Tinf)/(1-gamma.vec)
    Teff[is.na(Teff)] <- Tinf
    Teff[Teff == 0] <- Tinf
    I.pct.Teff <- I.pct/Teff
    I.pct.Teff[is.na(I.pct.Teff)]<-0
    # ### Seroprevalence total
    SP <- Roffset+c(0,cumsum(I.pct.Teff[-length(I.pct.Teff)]))+
      addCumulPos*CumulPosPct.1
    dattmp.all[indx,"I.pct"]<-I.pct
    dattmp.all[indx,"SP.pct"]<-SP
    I.pct.tplus1 <- c(I.pct[-1],NA)
    dI.pct <- I.pct.tplus1-I.pct.minus.Cases*(1-1/Tinf)
    dattmp.all[indx,"dI.pct"]<-dI.pct
    # - add +1 so 'current day' not double counted
    CumulPosPct.0 <- c(rep(0,pred.lag+1),dattmp$CumulPosPct[1:(nrow(dattmp)-(pred.lag+1))])
    dTrec <- pred.lag+ceiling(Trec)
    CumulPosPct.2 <- c(rep(0,dTrec),dattmp$CumulPosPct[1:(nrow(dattmp)-dTrec)])
    ID.pct <- CumulPosPct.0 - CumulPosPct.2
    dattmp.all[indx,"ID.pct"]<-ID.pct
    Itot.pct <- I.pct + ID.pct
    dattmp.all[indx,"Itot.pct"]<-Itot.pct
  }
  return(dattmp.all)
}

# get_all_fitpred <- function(ncoef.vec=rep(0.4,length(statesvec)),
#                             Roffset.vec=rep(0,length(statesvec)),
#                             Tinf.vec=rep(10,length(statesvec)),
#                             sig2err=0.1) {
#   dattmp.all <- get_all_pred(ncoef.vec,Roffset.vec,Tinf.vec,sig2err)
#   ## Merge
#   fitdat<-merge.data.frame(serodat.df,
#                            dattmp.all[,c("date","t","state","id","fips","I.pct","SP.pct")])
#   # Reorder to match original data
#   fitdat <- fitdat[order(fitdat$id,fitdat$t,fitdat$y),]
#   fitdat$ypred <- log(fitdat$SP.pct)
#   return(fitdat)
# }
# 
# get_all_likelihood <- function(ncoef.vec=rep(0.4,length(statesvec)),
#                                Roffset.vec=rep(0,length(statesvec)),
#                                Tinf.vec=rep(10,length(statesvec)),
#                                sig2err=0.1) {
#   fitdat <- get_all_fitpred(ncoef.vec,Roffset.vec,Tinf.vec,sig2err)
#   fitdat$lnlike <- dnorm(fitdat$y,mean=fitdat$ypred,sd=sqrt(fitdat$y.se^2+sig2err),log=TRUE)
#   lnlike <- sum(fitdat$lnlike)
#   return(lnlike)
# }


