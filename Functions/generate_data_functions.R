library(tidyverse)
library(jsonlite)
library(rvest)
library(data.table)
library(RCurl)
datezero <- "2019-12-31"
get_meantestdata <- function(datezero="2019-12-31",
                             mindate="2020-03-19",
                             nmean=7,align="right",
                             cleandat=TRUE,
                             datafolder="Data") {
  datfile <- file.path(datafolder,"daily.csv")
  if (file.exists(datfile)) {
    dat <- read.csv(datfile,as.is=TRUE)
  } else {
    download <- getURL("https://api.covidtracking.com/v1/states/daily.csv")
    dat <- read.csv(text = download,as.is=TRUE)
    write.csv(dat,datfile,row.names = FALSE)
  }
  dat$date <- as.Date(as.character(dat$date),tryFormats = c("%Y%m%d"))
  dat$numDate <- as.numeric(dat$date)-as.numeric(as.Date(datezero))
  meandat <- data.frame()
  for (statenow in unique(dat$state)) {
    tmpdat <- subset(dat,state==statenow)
    # if (cleandat) {
    #   indx <- tmpdat$positiveIncrease<0 |
    #     tmpdat$totalTestResultsIncrease<=0
    #   tmpdat$positiveIncrease[indx] <- NA
    #   tmpdat$deathIncrease[indx] <- NA
    #   tmpdat$totalTestResultsIncrease[indx] <- NA
    # }
    tmpdat <- tmpdat[order(tmpdat$numDate),]
    tmpdat$meanpositiveIncrease <- 
      frollmean(tmpdat$positiveIncrease,nmean,align=align,na.rm=TRUE)
    meanpositiveIncrease <- tmpdat$meanpositiveIncrease 
    meanpositiveIncrease[is.na(meanpositiveIncrease)]<-0
    tmpdat$meanpositiveCumul <- cumsum(meanpositiveIncrease)
    tmpdat$meandeathIncrease <- 
      frollmean(tmpdat$deathIncrease,nmean,align=align,na.rm=TRUE)
    tmpdat$meantotalTestResultsIncrease <- 
      frollmean(tmpdat$totalTestResultsIncrease,nmean,align=align,na.rm=TRUE)
    tmpdat$meantotalTestResultsIncrease[tmpdat$meantotalTestResultsIncrease<=0]<-NA
    # # Using posNeg instead
    # tmpdat$posNegIncrease <- c(tmpdat$posNeg[1],diff(tmpdat$posNeg))
    # tmpdat$meanposNegIncrease <- 
    #   frollmean(tmpdat$posNegIncrease,nmean,align=align,na.rm=TRUE)
    # tmpdat$meanposNegIncrease[tmpdat$meanposNegIncrease<=0]<-NA
    meandat <- rbind(meandat,tmpdat)
  }
  dat<-meandat
  dat <- subset(dat,date >= as.Date(mindate))
  dat.df <- melt(as.data.table(dat[,c(
    "date","numDate","state","fips",
    "meanpositiveIncrease",
    "meandeathIncrease",
    "meantotalTestResultsIncrease",
    "meanpositiveCumul")]),
    id.vars=1:4)
  dat.df <- dat.df[order(dat.df$variable,dat.df$numDate),]
  return(dat.df)
}

generate_testpositivity <- function(datafolder="Data",forceupdate=FALSE,
                                    tau=14,mcmcpath="MCMC") {
  datfile<-file.path(mcmcpath,"TestPositivityData.csv")
  if (file.exists(datfile) & !forceupdate) {
    alldat<-fread(datfile)
  } else {
    alldat_tau.df <- get_meantestdata(align="right",nmean = tau)
    alldat_tau <- spread(as.data.table(alldat_tau.df), variable, value)
    names(alldat_tau)[c(5,6,7)] <- c("Cases_tau","Deaths_tau","TotalTests_tau")
    alldat_tau$PositivePct_tau <- 100*alldat_tau$Cases_tau/(alldat_tau$TotalTests_tau)
    alldat_tau$TestsPer1000_tau <- (alldat_tau$TotalTests_tau)*1000/
      fips_table[as.character(alldat_tau$state),"pop"]  
    alldat<-alldat_tau
    fwrite(alldat,file=datfile)
  }
}


generate_imperial <- function(datafolder="Data") {
  ## Imperial college model data
  imperialorigfile<-file.path(datafolder,"Imperial-data-model-estimates.csv")
  imperialpredfile<-file.path(datafolder,"ImperialPredictions.csv")
  set.seed(exp(2))
  imperialdat <- fread(imperialorigfile)
  imperialdat$date <- as.Date(imperialdat$date)
  imperialdat$numDate <- as.numeric(imperialdat$date)-as.numeric(as.Date(datezero))
  popscale <- 100000 / fips_table[imperialdat$state,"pop"]
  imperialdat$reported_cases <- imperialdat$reported_cases * popscale
  imperialdat$predicted_infections_mean <- imperialdat$predicted_infections_mean * popscale
  imperialdat$predicted_infections_lower_CI_95 <- imperialdat$predicted_infections_lower_CI_95 * popscale
  imperialdat$predicted_infections_higher_CI_95_cumulative <- imperialdat$predicted_infections_higher_CI_95_cumulative * popscale
  mu <- log(sqrt(imperialdat$predicted_infections_lower_CI_95*imperialdat$predicted_infections_higher_CI_95_cumulative))
  sig <- log((imperialdat$predicted_infections_higher_CI_95_cumulative/imperialdat$predicted_infections_lower_CI_95))/(2*qnorm(0.975))
  imperialdat$predicted_infections_random <- exp(mu + sig*rnorm(nrow(imperialdat)))
  imperialdat$reported_deaths <- imperialdat$reported_deaths * popscale
  imperialdat$estimated_deaths_mean <- imperialdat$estimated_deaths_mean * popscale
  imperialdat$estimated_deaths_lower_CI_95 <- imperialdat$estimated_deaths_lower_CI_95 * popscale
  imperialdat$estimated_deaths_higher_CI_95 <- imperialdat$estimated_deaths_higher_CI_95 * popscale
  mu <- log(sqrt(imperialdat$estimated_deaths_lower_CI_95*imperialdat$estimated_deaths_higher_CI_95))
  sig <- log((imperialdat$estimated_deaths_higher_CI_95/imperialdat$estimated_deaths_lower_CI_95))/(2*qnorm(0.975))
  imperialdat$estimated_deaths_random <- exp(mu + sig*rnorm(nrow(imperialdat)))
  imperialdat$mean_infectious_individuals <- imperialdat$mean_infectious_individuals * popscale
  imperialdat$infectious_individuals_lower_CI_95 <- imperialdat$infectious_individuals_lower_CI_95 * popscale
  imperialdat$infectious_individuals_higher_CI_95 <- imperialdat$infectious_individuals_higher_CI_95 * popscale
  mu <- log(sqrt(imperialdat$infectious_individuals_lower_CI_95*imperialdat$infectious_individuals_higher_CI_95))
  sig <- log((imperialdat$infectious_individuals_higher_CI_95/imperialdat$infectious_individuals_lower_CI_95))/(2*qnorm(0.975))
  imperialdat$infectious_individuals_random <- exp(mu + sig*rnorm(nrow(imperialdat)))
  
  mu.Rt <- log(sqrt(imperialdat$`time_varying_reproduction_number_R(t)_lower_CI_95`*imperialdat$`time_varying_reproduction_number_R(t)_higher_CI_95`))
  sig.Rt <- log((imperialdat$`time_varying_reproduction_number_R(t)_higher_CI_95`/imperialdat$`time_varying_reproduction_number_R(t)_lower_CI_95`))/(2*qnorm(0.975))
  imperialdat$`time_varying_reproduction_number_R(t)_random` <- exp(mu.Rt + sig.Rt*rnorm(nrow(imperialdat)))
  
  setnames(imperialdat,
           old=c("reported_cases",
                 "predicted_infections_mean",
                 "predicted_infections_lower_CI_95",
                 "predicted_infections_higher_CI_95_cumulative",
                 "predicted_infections_random",
                 "reported_deaths",
                 "estimated_deaths_mean",
                 "estimated_deaths_lower_CI_95",
                 "estimated_deaths_higher_CI_95",
                 "estimated_deaths_random",
                 "mean_time_varying_reproduction_number_R(t)",
                 "time_varying_reproduction_number_R(t)_lower_CI_95",
                 "time_varying_reproduction_number_R(t)_higher_CI_95",
                 "time_varying_reproduction_number_R(t)_random",
                 "mean_infectious_individuals",
                 "infectious_individuals_lower_CI_95",
                 "infectious_individuals_higher_CI_95",
                 "infectious_individuals_random"),
           new=c("Cases",
                 "N_dtCumInfected.50.",
                 "N_dtCumInfected.2.5.",
                 "N_dtCumInfected.97.5.",
                 "N_dtCumInfected.Random",
                 "Deaths",
                 "D_pos.50.",
                 "D_pos.2.5.",
                 "D_pos.97.5.",
                 "D_pos.Random",
                 "Rt.50.",
                 "Rt.2.5.",
                 "Rt.97.5.",
                 "Rt.Random",
                 "N_I_tot.50.",
                 "N_I_tot.2.5.",
                 "N_I_tot.97.5.",
                 "N_I_tot.Random")
  )
  
  ## Cumulative infected (per 100000) calculated by cumulative sum of daily infected
  imperialdat$CumInfected.2.5. <- imperialdat$N_dtCumInfected.2.5.
  imperialdat$CumInfected.50. <- imperialdat$N_dtCumInfected.50.
  imperialdat$CumInfected.97.5. <- imperialdat$N_dtCumInfected.97.5.
  for (statenow in unique(imperialdat$state)) {
    indx <- imperialdat$state==statenow
    imperialdat$CumInfected.2.5.[indx]<-cumsum(imperialdat$CumInfected.2.5.[indx])
    imperialdat$CumInfected.50.[indx]<-cumsum(imperialdat$CumInfected.50.[indx])
    imperialdat$CumInfected.97.5.[indx]<-cumsum(imperialdat$CumInfected.97.5.[indx])
  }
  mu.cumI <- log(sqrt(imperialdat$CumInfected.2.5.*imperialdat$CumInfected.97.5.))
  sig.cumI <- log((imperialdat$CumInfected.97.5./imperialdat$CumInfected.2.5.))/(2*qnorm(0.975))
  imperialdat$CumInfected.Random <- exp(mu.cumI + sig.cumI*rnorm(nrow(imperialdat)))
  
  ## Susceptible (per 100000)
  imperialdat$S.2.5. <- 100000 - imperialdat$CumInfected.97.5.
  imperialdat$S.50. <- 100000 - imperialdat$CumInfected.50.
  imperialdat$S.97.5. <- 100000 - imperialdat$CumInfected.2.5.
  mu.S <- log(sqrt(imperialdat$S.2.5.*imperialdat$S.97.5.))
  sig.S <- log((imperialdat$S.97.5./imperialdat$S.2.5.))/(2*qnorm(0.975))
  imperialdat$S.Random <- exp(mu.S + sig.S*rnorm(nrow(imperialdat)))
  
  # Assume S and Reff independent - convert S back to per capita
  mu.Reff <- mu.Rt+mu.S-log(100000)
  sig.Reff <- sqrt(sig.Rt+sig.S)
  imperialdat$Refft.2.5. <- exp(mu.Reff+qnorm(0.025)*sig.Reff)
  imperialdat$Refft.50. <- exp(mu.Reff)
  imperialdat$Refft.97.5. <- exp(mu.Reff+qnorm(0.975)*sig.Reff)
  imperialdat$Refft.Random <- exp(mu.Reff + sig.Reff*rnorm(nrow(imperialdat)))
  
  fwrite(imperialdat,imperialpredfile)
}

