# tuning sd for component jumps
#   eta.tuning 
# hyperparameters: 
#   mu_ncoef - logit of m_ncoef has flat prior
#   logSigma_ncoef - log of sd_ncoef has flat prior
#   mu_Roffset - logit of m_Roffset/100 has flat prior
#   logSigma_Roffset - log of sd_ncoef has flat prior
# global parameters: 
#   Tinf - has normal prior (m=14, sd=3.5)
#   logsig2err - log of sig^2 error has flat prior 
# population model: 
#   logit(ncoef) ~ N(mu_ncoef,exp(logSigma_ncoef))
#   logit(Roffset/100) ~ N(mu_Roffset,exp(logSigma_Roffset))
# Algorithm
#   (1) for mu_ncoef, logSigma_ncoef, mu_Roffset, logSigma_Roffset, (likelihoods don't change)
#       (i) add q ~ rnorm(0,eta.tuning)
#       (ii) calculate population prior'
#       (iii) accept with probability min(1, (prior' / prior))
#       (iv) repeat 10 times for better mixing
#   (2) Tinf, logsig2err (no population model, but affects prior and likelihood)
#       (i) add q ~ rnorm(0,eta.tuning) 
#       (ii) calculate prior' and overall likelihood'
#       (iii) accept with probability min(1, (prior' * likelihood') / (prior * likelihood))
#   (3) for each state, separately for ncoeff and Roffset
#       (i) add q ~ rnorm(0,eta.tuning)
#       (ii) calculate population prior'
#       (iii) calculate individual likelihood'
#       (iv) accept with probability min(1, (prior' * likelihood') / (prior * likelihood))
#   (4) repeat

eta.tuning <- 0.1
eta.tuning.hyper <- eta.tuning*c(10,4,1,4) 
eta.tuning.global <- eta.tuning*c(10,4)
eta.tuning.indiv <- eta.tuning*c(rep(5,n.id),rep(20,n.id))

# Default values
m_ncoef <- 0.5
Sigma_ncoef <- 0.1
m_Roffset <- 0.3/100
Sigma_Roffset <- 0.1
m_Tinf <- 14
sd_Tinf <- 3.5
sig2err <- 0.5^2

# Check if fixed n
if (Fixed_n) {
  m_ncoef <- n_default # Fix value of ncoef
  Sigma_ncoef <- 0 # No population variation
  eta.tuning.hyper[1:2] <- 0 # Do not jump population
  eta.tuning.indiv[1:n.id] <- 0 # Do not jump individual
}

hyperparms.0 <- c(logit(m_ncoef),log(Sigma_ncoef),logit(m_Roffset),log(Sigma_Roffset))
n.hyperparms <- length(hyperparms.0)
globalparms.0 <- c(m_Tinf,log(sig2err))
n.globalparms <- length(globalparms.0)
indivparms.0 <- c(rep(hyperparms.0[1],n.id),
                  rep(hyperparms.0[3],n.id))
n.indivparms <- length(indivparms.0)

# restart conditions
allparms<-data.frame(fread(paste0(restartname,".",chainnum,".csv"),colClasses="numeric"))
lastparms<-allparms[nrow(allparms),]
iter0 <- lastparms$iter

hyperparms.current <- as.numeric(lastparms[,2:5])
globalparms.current <- as.numeric(lastparms[,6:7])
indivparms.current <- as.numeric(lastparms[,9:110])

# initial priors
lnprior.current <- dnorm(globalparms.current[1],m=m_Tinf,sd=sd_Tinf,log=TRUE)
if (Fixed_n) {
  lnpopprior.current <- 
    sum(dnorm(indivparms.current[n.id+(1:n.id)],
              m=hyperparms.current[3],sd=exp(hyperparms.current[4]),log=TRUE))
} else {
  lnpopprior.current <- 
    sum(dnorm(indivparms.current[1:n.id],
              m=hyperparms.current[1],sd=exp(hyperparms.current[2]),log=TRUE))+
    sum(dnorm(indivparms.current[n.id+(1:n.id)],
              m=hyperparms.current[3],sd=exp(hyperparms.current[4]),log=TRUE))
}

# scaling
ncoef.vec.current <- inv.logit(indivparms.current[1:n.id])
Roffset.vec.current <- 100*inv.logit(indivparms.current[n.id+(1:n.id)])
Tinf.vec.current <- rep(globalparms.current[1],n.id)
sig2err.current <- exp(globalparms.current[2])

# initial likelihoods
lnlike.current <- get_all_likelihood_fast(ncoef.vec.current,
                                     Roffset.vec.current,
                                     Tinf.vec.current,
                                     sig2err.current,
                                     addCumulPos = addCumulPos,
                                     Trec.vec=Trec.vec.global)
lnlike.indiv.current <- numeric(n.id)
for (j in 1:n.id) lnlike.indiv.current[j] <- 
  get_one_likelihood(j,ncoef.vec.current[j],
                     Roffset.vec.current[j],
                     Tinf.vec.current[j],
                     sig2err.current,
                     addCumulPos = addCumulPos,
                     Trec=Trec.vec.global[j])

lnpost.current <- lnprior.current+lnpopprior.current+lnlike.current
# # all parameters
# allparms <- c(0,
#               hyperparms.current,
#               globalparms.current,
#               sig2err.current,
#               indivparms.current,
#               ncoef.vec.current,
#               Roffset.vec.current,
#               lnlike.indiv.current,
#               lnprior.current,
#               lnpopprior.current,
#               lnlike.current,
#               lnpost.current)
# names(allparms) <- c("iter",
#                      "mu_ncoef", "logSigma_ncoef", "mu_Roffset","logSigma_Roffset",
#                      "Tinf", "logsig2err","sig2err",
#                      paste("logit.ncoef",1:n.id,sep="."),
#                      paste("logit.Roffset.100",1:n.id,sep="."),
#                      paste("ncoef",1:n.id,sep="."),
#                      paste("Roffset",1:n.id,sep="."),
#                      paste("lnlike",1:n.id,sep="."),
#                      "lnprior","lnpopprior","lnlike","lnpost")
# 
# allparms <- as.data.frame(t(allparms))

accept.hyper<-0
total.hyper<-0
accept.global<-0
total.global<-0
accept.indiv<-0
total.indiv<-0
for (iter in (iter0+(1:niter))) {
  if ((iter %% 10)==0) print(iter)
  for (subiter in 1:10) { # do 10 steps here
    for (j in 1:n.hyperparms) { # only affects population prior
      if (eta.tuning.hyper[j] > 0) { # don't run if no jumps
        total.hyper<-total.hyper+1
        hyperparms.proposed <- hyperparms.current
        globalparms.proposed <- globalparms.current
        indivparms.proposed <- indivparms.current
        q <- rnorm(1,m=0,sd=eta.tuning.hyper[j])
        hyperparms.proposed[j] <- hyperparms.proposed[j] + q 
        if (Fixed_n) {
          lnpopprior.proposed <-
            sum(dnorm(indivparms.proposed[n.id+(1:n.id)],
                      m=hyperparms.proposed[3],sd=exp(hyperparms.proposed[4]),log=TRUE))
        } else {
          lnpopprior.proposed <-
            sum(dnorm(indivparms.proposed[1:n.id],
                      m=hyperparms.proposed[1],sd=exp(hyperparms.proposed[2]),log=TRUE))+
            sum(dnorm(indivparms.proposed[n.id+(1:n.id)],
                      m=hyperparms.proposed[3],sd=exp(hyperparms.proposed[4]),log=TRUE))
        }          
        lndiff <- lnpopprior.proposed - lnpopprior.current
        if (!is.nan(lndiff)) {
          if (log(runif(1)) < lndiff) {
            hyperparms.current <- hyperparms.proposed
            lnpopprior.current <- lnpopprior.proposed
            accept.hyper<-accept.hyper+1
          }
        }
        # cat("hyper",paste(c(hyperparms.current,lnpopprior.current)),"\n")
      }
    }
  }
  for (j in 1:n.globalparms) { # only affects prior and likelihood
    total.global<-total.global+1
    hyperparms.proposed <- hyperparms.current
    globalparms.proposed <- globalparms.current
    indivparms.proposed <- indivparms.current
    q <- rnorm(1,m=0,sd=eta.tuning.global[j])
    globalparms.proposed[j] <- globalparms.proposed[j] + q 
    # if (j==2) print(paste(globalparms.current[j],globalparms.proposed[j]))
    ncoef.vec.proposed <- inv.logit(indivparms.proposed[1:n.id])
    Roffset.vec.proposed <- 100*inv.logit(indivparms.proposed[n.id+(1:n.id)])
    Tinf.vec.proposed <- rep(globalparms.proposed[1],n.id)
    if (min(Tinf.vec.proposed)<0) Tinf.vec.proposed <- Tinf.vec.current # truncate at 0
    sig2err.proposed <- exp(globalparms.proposed[2])
    # proposed and current prior for Tinf
    lnprior.proposed <- dnorm(globalparms.proposed[1],m=m_Tinf,sd=sd_Tinf,log=TRUE)
    # proposed likelihoods
    lnlike.proposed <- get_all_likelihood_fast(ncoef.vec.proposed,
                                          Roffset.vec.proposed,
                                          Tinf.vec.proposed,
                                          sig2err.proposed,
                                          addCumulPos = addCumulPos,
                                          Trec.vec=Trec.vec.global)
    lndiff <- lnprior.proposed - lnprior.current + lnlike.proposed - lnlike.current
    if (!is.nan(lndiff)) {
      if (log(runif(1)) < lndiff) {
        globalparms.current <- globalparms.proposed
        lnprior.current <- lnprior.proposed
        lnlike.current <- lnlike.proposed
        Tinf.vec.current <- Tinf.vec.proposed
        sig2err.current <- sig2err.proposed
        for (idnow in 1:n.id) lnlike.indiv.current[idnow] <-
          get_one_likelihood(idnow,ncoef.vec.proposed[idnow],
                             Roffset.vec.proposed[idnow],
                             Tinf.vec.proposed[idnow],
                             sig2err.proposed,
                             addCumulPos = addCumulPos,
                             Trec=Trec.vec.global[idnow])
        accept.global<-accept.global+1
      }
    }
    # cat("global",paste(c(globalparms.current,lnlike.current)),"\n")
  }
  for (j in 1:n.indivparms) { # affects population prior and likelihood
    if (eta.tuning.indiv[j] > 0) { # Don't run if no jumps
      total.indiv<-total.indiv+1
      if (j > n.id) {idnow <- (j-n.id)} else {idnow <- j}
      hyperparms.proposed <- hyperparms.current
      globalparms.proposed <- globalparms.current
      indivparms.proposed <- indivparms.current
      q <- rnorm(1,m=0,sd=eta.tuning.indiv[j])
      indivparms.proposed[j] <- indivparms.proposed[j] + q
      ncoef.vec.proposed <- inv.logit(indivparms.proposed[1:n.id])
      Roffset.vec.proposed <- 100*inv.logit(indivparms.proposed[n.id+(1:n.id)])
      Tinf.vec.proposed <- rep(globalparms.proposed[1],n.id)
      if (min(Tinf.vec.proposed)<0) Tinf.vec.proposed <- Tinf.vec.current # truncate at 0
      sig2err.proposed <- exp(globalparms.proposed[2])
      if (Fixed_n) {
        lnpopprior.proposed <-
          sum(dnorm(indivparms.proposed[n.id+(1:n.id)],
                    m=hyperparms.proposed[3],sd=exp(hyperparms.proposed[4]),log=TRUE))
      } else {
        lnpopprior.proposed <-
          sum(dnorm(indivparms.proposed[1:n.id],
                    m=hyperparms.proposed[1],sd=exp(hyperparms.proposed[2]),log=TRUE))+
          sum(dnorm(indivparms.proposed[n.id+(1:n.id)],
                    m=hyperparms.proposed[3],sd=exp(hyperparms.proposed[4]),log=TRUE))
      }
      lndiff <- lnpopprior.proposed - lnpopprior.current
      lnlike.indiv.proposed <- lnlike.indiv.current
      lnlike.indiv.proposed[idnow] <-
        get_one_likelihood(idnow,ncoef.vec.proposed[idnow],
                           Roffset.vec.proposed[idnow],
                           Tinf.vec.proposed[idnow],
                           sig2err.proposed,
                           addCumulPos = addCumulPos,
                           Trec=Trec.vec.global[idnow])
      lndiff <- lndiff + lnlike.indiv.proposed[idnow] - lnlike.indiv.current[idnow]
      if (!is.nan(lndiff)) {
        if (log(runif(1)) < lndiff) {
          indivparms.current <- indivparms.proposed
          lnpopprior.current <- lnpopprior.proposed
          lnlike.indiv.current <- lnlike.indiv.proposed
          lnlike.current <- sum(lnlike.indiv.proposed)
          ncoef.vec.current <- ncoef.vec.proposed
          Roffset.vec.current <- Roffset.vec.proposed
          accept.indiv<-accept.indiv+1
        }
      }
      #cat("indiv",paste(c(lnpopprior.current,sum(lnlike.indiv.current))),"\n")
    }
  }
  # print(lnlike.current-get_all_likelihood_fast(ncoef.vec.current,
  #                           Roffset.vec.current,
  #                           Tinf.vec.current,
  #                           sig2err.current,
  #                           addCumulPos = addCumulPos,
  #                           Trec.vec=Trec.vec.global))
  lnpost.current <- lnprior.current+lnpopprior.current+lnlike.current
  newparms <- c(iter,
                hyperparms.current, 
                globalparms.current, 
                sig2err.current,
                indivparms.current,
                ncoef.vec.current, 
                Roffset.vec.current,
                lnlike.indiv.current,
                lnprior.current,
                lnpopprior.current,
                lnlike.current,
                lnpost.current)
  allparms <- rbind(allparms,newparms)
}
if (Fixed_n) {
  pairs(allparms[-(1:10),c(4:7,264:267)])
} else {
  pairs(allparms[-(1:10),c(2:7,264:267)])
}
plot((allparms$lnpost)[-(1:10)])
print(paste("Hyper acceptance rate:",accept.hyper/total.hyper))
print(paste("Global acceptance rate:",accept.global/total.global))
print(paste("Indiv acceptance rate:",accept.indiv/total.indiv))
accept.parms <- numeric(1+n.hyperparms+n.globalparms+n.indivparms)
names(accept.parms)<-names(allparms)[2:(2+n.hyperparms+n.globalparms+n.indivparms)]
for (j in 2:(2+n.hyperparms+n.globalparms+n.indivparms)) accept.parms[j-1]<-length(unique(allparms[[j]]))/(niter+1)
plot(accept.parms)

fwrite(allparms,paste0(runname,".",chainnum,".csv"))
# if (Fixed_n) {
#   allparms.coda <- mcmc(allparms[-(1:10),-(1:3)],start=11)
# } else {
#   allparms.coda <- mcmc(allparms[-(1:10),-1],start=11)
# }
# pdf(paste0(runname,".coda.",chainnum,".pdf"),height=8,width=12)
# plot(allparms.coda)
# dev.off()
