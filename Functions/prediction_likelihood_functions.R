get_one_pred <- function(idnow,ncoef,Roffset,Tinf,sig2err,
                         addCumulPos = 1,Trec=10) {
  # Only select current id
  indx <- testdat$id==idnow
  dattmp <- testdat[indx,]
  ## Calculate predictions
  I.pct <- (dattmp$PositivePct_tau)^(1-ncoef)*(dattmp$Cases_tau_pct)^ncoef
  I.pct[is.na(I.pct)]<-0
  ## Correction for diagnosed cases
  I.pct.minus.Cases <- I.pct - dattmp$Cases_tau_pct
  I.pct.minus.Cases[is.na(I.pct.minus.Cases)] <- 0
  I.pct.minus.Cases[(I.pct.minus.Cases < 0)] <- 0
  ## Cumulative diagnosed 
  dTcumul <- pred.lag+ceiling(Tinf)
  CumulPosPct.1 <- c(rep(0,dTcumul),dattmp$CumulPosPct[1:(nrow(dattmp)-dTcumul)])
  ## SP calculation
  # SP <- Roffset+c(0,cumsum(I.pct.minus.Cases[-length(I.pct)]))/Tinf+
  #   addCumulPos*CumulPosPct.1
  ### [Alternative]
  gamma.vec <- I.pct/I.pct.minus.Cases
  Teff <- gamma.vec*(1-gamma.vec^Tinf)/(1-gamma.vec)
  Teff[is.na(Teff)] <- Tinf
  Teff[Teff == 0] <- Tinf
  I.pct.Teff <- I.pct/Teff
  I.pct.Teff[is.na(I.pct.Teff)]<-0
  # ### Seroprevalence total
  SP <- Roffset+c(0,cumsum(I.pct.Teff[-length(I.pct.Teff)]))+
    addCumulPos*CumulPosPct.1
  I.pct.tplus1 <- c(I.pct[-1],NA)
  dI.pct <- I.pct.tplus1-(I.pct.minus.Cases-I.pct.Teff)
  # add +1 so 'current day' not double counted
  CumulPosPct.0 <- c(rep(0,pred.lag+1),dattmp$CumulPosPct[1:(nrow(dattmp)-(pred.lag+1))])
  dTrec <- pred.lag+ceiling(Trec)
  CumulPosPct.2 <- c(rep(0,dTrec),dattmp$CumulPosPct[1:(nrow(dattmp)-dTrec)])
  ID.pct <- CumulPosPct.0 - CumulPosPct.2
  Itot.pct <- I.pct+ID.pct
  dattmp$I.pct<-I.pct
  dattmp$SP.pct<-SP
  dattmp$dI.pct<-dI.pct
  dattmp$ID.pct<-ID.pct
  dattmp$Itot.pct<-Itot.pct
  return(dattmp)
}

get_one_fitpred <- function(idnow,ncoef,Roffset,Tinf,sig2err,addCumulPos = 1,Trec=10) {
  
  dattmp <- get_one_pred(idnow,ncoef,Roffset,Tinf,sig2err,addCumulPos = addCumulPos,Trec=Trec)
  ## Merge
  fitdat<-merge.data.frame(subset(serodat.df,id==idnow),
                           dattmp[,c("date","t","state","id","fips","I.pct","SP.pct")])
  # Reorder to match original data
  fitdat <- fitdat[order(fitdat$id,fitdat$t,fitdat$y),]
  fitdat$ypred <- log(fitdat$SP.pct)
  return(fitdat)
}

get_one_likelihood <- function(idnow,ncoef,Roffset,Tinf,sig2err,addCumulPos = 1,Trec=10) {
  fitdat <- get_one_fitpred(idnow,ncoef,Roffset,Tinf,sig2err,addCumulPos = addCumulPos,Trec=Trec)
  fitdat$lnlike <- dnorm(fitdat$y,mean=fitdat$ypred,sd=sqrt(fitdat$y.se^2+sig2err),log=TRUE)
  lnlike <- sum(fitdat$lnlike)
  return(lnlike)
}

get_all_pred <- function(ncoef.vec=rep(0.4,length(statesvec)),
                         Roffset.vec=rep(0,length(statesvec)),
                         Tinf.vec=rep(10,length(statesvec)),
                         sig2err=0.1,
                         addCumulPos = 1,
                         Trec.vec = rep(10,length(statesvec))) {
  ## Calculate predictions
  dattmp.all <- testdat
  for (idnow in unique(dattmp.all$id)) {
    ncoef <- ncoef.vec[idnow]
    Roffset <- Roffset.vec[idnow]
    Tinf <- Tinf.vec[idnow]
    Trec <- Trec.vec[idnow]
    indx <- dattmp.all$id==idnow
    dattmp <- dattmp.all[indx,]
    I.pct <- (dattmp$PositivePct_tau)^(1-ncoef)*(dattmp$Cases_tau_pct)^ncoef
    I.pct[is.na(I.pct)]<-0
    ## Correction for diagnosed cases
    I.pct.minus.Cases <- I.pct - dattmp$Cases_tau_pct
    I.pct.minus.Cases[is.na(I.pct.minus.Cases)] <- 0
    I.pct.minus.Cases[(I.pct.minus.Cases < 0)] <- 0
    ## Cumulative diagnosed 
    dTcumul <- pred.lag+ceiling(Tinf) 
    CumulPosPct.1 <- c(rep(0,dTcumul),dattmp$CumulPosPct[1:(nrow(dattmp)-dTcumul)])
    # ## SP calculation
    # SP <- Roffset+c(0,cumsum(I.pct.minus.Cases[-length(I.pct)]))/Tinf+
    #   addCumulPos*CumulPosPct.1
    ### [Alternative]
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
    dI.pct <- I.pct.tplus1-(I.pct.minus.Cases-I.pct.Teff)
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

get_all_pred_fast <- function(ncoef.vec=rep(0.4,n.id),
                         Roffset.vec=rep(0,n.id),
                         Tinf.vec=10,
                         sig2err=0.1,
                         addCumulPos = 1,
                         Trec.vec=10) {
  ## Parameters
  ndates<-nrow(Cases_tau.pct.mat)
  ncoef.mat <- matrix(rep(ncoef.vec,ndates),ncol=n.id,byrow=T)
  Roffset.mat <- matrix(rep(Roffset.vec,ndates),ncol=n.id,byrow=T)
  ## Predictions
  I.pct.mat <- (PositivePct_tau.mat^(1-ncoef.mat))*(Cases_tau.pct.mat)^ncoef.mat
  I.pct.mat[is.na(I.pct.mat)]<-0
  ## Correction for diagnosed cases
  I.pct.minus.Cases.mat <- I.pct.mat - Cases_tau.pct.mat
  I.pct.minus.Cases.mat[is.na(I.pct.minus.Cases.mat)] <- 0
  I.pct.minus.Cases.mat[(I.pct.minus.Cases.mat < 0)] <- 0
  ### Offset by one, redo dates
  I.pct.minus.Cases.mat.m1 <- rbind(rep(0,n.id),I.pct.minus.Cases.mat[-ndates,])
  rownames(I.pct.minus.Cases.mat.m1)<-rownames(I.pct.minus.Cases.mat)
  ### Cumulative - offset by Tinf = Seroprevalence diagnosed
  dTcumul <- pred.lag+ceiling(mean(Tinf.vec)) # All Tinf the same 
  CumulPosPct.mat.m1 <- rbind(matrix(data=rep(0,n.id*dTcumul),
                                     nrow = dTcumul,ncol=n.id),
                              CumulPosPct.mat[-((ndates-dTcumul+1):ndates),]) 
  rownames(CumulPosPct.mat.m1) <- rownames(I.pct.mat)
  ### Seroprevalence total
  # SP.pct.mat <- Roffset.mat + apply(I.pct.minus.Cases.mat.m1,2,cumsum)/Tinf + 
  #   addCumulPos*CumulPosPct.mat.m1
  ### [Alternative]
  gamma.mat <- I.pct.mat/I.pct.minus.Cases.mat
  Teff.mat <- gamma.mat*(1-gamma.mat^mean(Tinf.vec))/(1-gamma.mat)
  Teff.mat[is.na(Teff.mat)] <- mean(Tinf.vec)
  Teff.mat[Teff.mat == 0] <- mean(Tinf.vec)
  I.pct.Teff.mat <- I.pct.mat/Teff.mat
  # Offset by 1, redo dates
  I.pct.Teff.mat.m1 <- rbind(rep(0,n.id),I.pct.Teff.mat[-ndates,])
  rownames(I.pct.Teff.mat.m1)<-rownames(I.pct.mat)
  I.pct.Teff.mat.m1[is.na(I.pct.Teff.mat.m1)] <- 0
  # ### Seroprevalence total
  SP.pct.mat <- Roffset.mat + apply(I.pct.Teff.mat.m1,2,cumsum) +
    addCumulPos*CumulPosPct.mat.m1
  ## Offset by one, redo dates
  I.pct.mat.p1 <- rbind(I.pct.mat[-1,],rep(NA,n.id))
  rownames(I.pct.mat.p1)<-rownames(I.pct.mat)
  ### Daily new undiagnosed cases 
  dI.pct.mat <- I.pct.mat.p1-(I.pct.minus.Cases.mat-I.pct.Teff.mat)
  ### Diagnosed prevalence = Cumulative Diagnosed - Recovered
  ### Cumulative Diagnosed - up to t - (pred.lag+1)
  # - add +1 so 'current day' not double counted
  CumulPosPct.mat.m0 <- rbind(matrix(data=rep(0,n.id*(pred.lag+1)),
                                     nrow = (pred.lag+1),ncol=n.id),
                              CumulPosPct.mat[-((ndates-(pred.lag+1)+1):ndates),]) 
  rownames(CumulPosPct.mat.m0) <- rownames(I.pct.mat)
  ### Cumulative Recovered - up to t - (pred.lag + Trec)
  dTrec <- pred.lag+ceiling(mean(Trec.vec))
  CumulPosPct.mat.m2 <- rbind(matrix(data=rep(0,n.id*dTrec),
                                     nrow = dTrec,ncol=n.id),
                              CumulPosPct.mat[-((ndates-dTrec+1):ndates),]) 
  rownames(CumulPosPct.mat.m2) <- rownames(I.pct.mat)
  ID.pct.mat <- CumulPosPct.mat.m0 - CumulPosPct.mat.m2
  rownames(ID.pct.mat) <- rownames(I.pct.mat)
  ### Total prevalence
  Itot.pct.mat <- I.pct.mat.p1 + ID.pct.mat
  ## Reshape
  # I.pct.df <- as.data.frame(as.table(I.pct.mat),stringsAsFactors = FALSE)
  # SP.pct.df <- as.data.frame(as.table(SP.pct.mat),stringsAsFactors = FALSE)
  # dI.pct.df <- as.data.frame(as.table(dI.pct.mat),stringsAsFactors = FALSE)
  # dattmp.all <- cbind(I.pct.df,SP.pct.df[["Freq"]],dI.pct.df[["Freq"]])
  # names(dattmp.all)<-c("date","id","I.pct","SP.pct","dI.pct")
  # return(dattmp.all)
  return(list(I.pct.mat,SP.pct.mat,dI.pct.mat,ID.pct.mat,Itot.pct.mat))
}


get_all_fitpred_fast <- function(ncoef.vec=rep(0.4,length(statesvec)),
                            Roffset.vec=rep(0,length(statesvec)),
                            Tinf.vec=rep(10,length(statesvec)),
                            sig2err=0.1,addCumulPos = 1,
                            Trec.vec=rep(10,length(statesvec))) {
  dattmp.list <- get_all_pred_fast(ncoef.vec,Roffset.vec,Tinf.vec,sig2err,addCumulPos = addCumulPos,Trec.vec=Trec.vec)
  I_fast <- as.data.frame.table(dattmp.list[[1]], responseName = "I.pct")
  names(I_fast)[1:2]<-c("date","id")
  SP_fast <- as.data.frame.table(dattmp.list[[2]], responseName = "SP.pct")
  names(SP_fast)[1:2]<-c("date","id")
  dI_fast <- as.data.frame.table(dattmp.list[[3]], responseName = "dI.pct")
  names(dI_fast)[1:2]<-c("date","id")
  dattmp.all <- inner_join(inner_join(I_fast,SP_fast,by=c("date","id")),
                           dI_fast,by=c("date","id"))
  dattmp.all$date <- as.Date(dattmp.all$date)
  dattmp.all$id <- as.numeric(dattmp.all$id)
  dattmp.all <- inner_join(testdat[,1:(ncol(testdat)-2)],dattmp.all,
                           by=c("date","id"))
  ## Merge
  fitdat<-inner_join(serodat.df,
                     dattmp.all[,c("date","t","state","id","fips","I.pct","SP.pct")],
                     by=c("id","t"))
  # Reorder to match original data
  fitdat <- fitdat[order(fitdat$id,fitdat$t,fitdat$y),]
  fitdat$ypred <- log(fitdat$SP.pct)
  return(fitdat)
}

get_all_fitpred <- function(ncoef.vec=rep(0.4,length(statesvec)),
                            Roffset.vec=rep(0,length(statesvec)),
                            Tinf.vec=rep(10,length(statesvec)),
                            sig2err=0.1,addCumulPos = 1,Trec.vec=rep(10,length(statesvec))) {
  dattmp.all <- get_all_pred(ncoef.vec,Roffset.vec,Tinf.vec,sig2err,addCumulPos = addCumulPos,Trec.vec=Trec.vec)
  ## Merge
  fitdat<-merge.data.frame(serodat.df,
                           dattmp.all[,c("date","t","state","id","fips","I.pct","SP.pct")])
  # Reorder to match original data
  fitdat <- fitdat[order(fitdat$id,fitdat$t,fitdat$y),]
  fitdat$ypred <- log(fitdat$SP.pct)
  return(fitdat)
}


get_all_likelihood_fast <- function(ncoef.vec=rep(0.4,length(statesvec)),
                               Roffset.vec=rep(0,length(statesvec)),
                               Tinf.vec=rep(10,length(statesvec)),
                               sig2err=0.1,addCumulPos = 1,
                               Trec.vec=rep(10,length(statesvec))) {
  fitdat <- get_all_fitpred_fast(ncoef.vec,Roffset.vec,Tinf.vec,sig2err,addCumulPos = addCumulPos,Trec.vec=Trec.vec)
  fitdat$lnlike <- dnorm(fitdat$y,mean=fitdat$ypred,sd=sqrt(fitdat$y.se^2+sig2err),log=TRUE)
  lnlike <- sum(fitdat$lnlike)
  return(lnlike)
}

get_all_likelihood <- function(ncoef.vec=rep(0.4,length(statesvec)),
                               Roffset.vec=rep(0,length(statesvec)),
                               Tinf.vec=rep(10,length(statesvec)),
                               sig2err=0.1,addCumulPos = 1,Trec.vec=rep(10,length(statesvec))) {
  fitdat <- get_all_fitpred(ncoef.vec,Roffset.vec,Tinf.vec,sig2err,addCumulPos = addCumulPos,Trec.vec=Trec.vec)
  fitdat$lnlike <- dnorm(fitdat$y,mean=fitdat$ypred,sd=sqrt(fitdat$y.se^2+sig2err),log=TRUE)
  lnlike <- sum(fitdat$lnlike)
  return(lnlike)
}
