library(tidyverse)
library(jsonlite)
library(rvest)
library(data.table)
library(RCurl)

get_england_meantestdata <- function(datezero="2019-12-31",
                                     mindate="2020-03-01",
                                     nmean=14,align="right",
                                     include_negtest=TRUE,
                                     cleandat=FALSE) {
  englandfile <- file.path("International","EnglandTestingData.csv")
  if (!file.exists(englandfile)) {
    download <- getURL("https://api.coronavirus.data.gov.uk/v2/data?areaType=nation&areaCode=E92000001&metric=cumPillarOneTwoTestsByPublishDate&metric=newCasesByPublishDate&metric=newPillarOneTwoTestsByPublishDate&metric=cumCasesByPublishDate&format=csv")
    dat <- read.csv(text = download,as.is=TRUE)
    write.csv(dat,file=englandfile,row.names = FALSE)
  }
  dat <- read.csv(englandfile,as.is=TRUE)
  dat$date <- as.Date(as.character(dat$date))#,tryFormats = c("%Y%m%d"))
  dat$numDate <- as.numeric(dat$date)-as.numeric(as.Date(datezero))
  meandat <- dat
  if (cleandat) {
    indx <- meandat$new_cases<0 |
      meandat$new_tests<0
    meandat$new_cases[indx] <- NA
    meandat$new_tests[indx] <- NA
  }
  meandat <- meandat[order(meandat$numDate),]
  meandat$meannew_cases <- 
    frollmean(meandat$newCasesByPublishDate,nmean,align=align,na.rm=TRUE)
  meandat$meannew_cases[is.na(meandat$meannew_cases)]<-0
  meandat$meannew_casesCumul <- cumsum(meandat$meannew_cases)
  meandat$meannew_tests <- 
    frollmean(meandat$newPillarOneTwoTestsByPublishDate,nmean,align=align,na.rm=TRUE)
  dat<-meandat
  dat <- subset(dat,date >= as.Date(mindate))
  dat$population <-56287000
  if (include_negtest) {
    dat.df<-pivot_longer(dat[,c(
      "date","numDate","areaName","areaCode","population",
      "newCasesByPublishDate",
      "newPillarOneTwoTestsByPublishDate",
      "meannew_cases",
      "meannew_tests","meannew_casesCumul")],
      cols = 6:10,names_to="variable")
  } else {
    dat.df<-pivot_longer(dat[,c(
      "date","numDate","areaName","areaCode","population",
      "newCasesByPublishDate",
      "meannew_cases",
      "meannew_casesCumul")],
      cols = 6:8,names_to="variable")
  }
  dat.df <- dat.df[order(dat.df$variable,dat.df$numDate),]
  return(dat.df)
}

get_owid_meantestdata <- function(datezero="2019-12-31",
                                  mindate="2020-03-01",
                                  nmean=14,align="right",
                                  include_negtest=TRUE,
                                  cleandat=FALSE) {
  download <- getURL("https://covid.ourworldindata.org/data/owid-covid-data.csv")
  dat <- read.csv(text = download,as.is=TRUE)
  dat$date <- as.Date(as.character(dat$date))#,tryFormats = c("%Y%m%d"))
  dat$numDate <- as.numeric(dat$date)-as.numeric(as.Date(datezero))
  meandat <- data.frame()
  for (locationnow in unique(dat$location)) {
    tmpdat <- subset(dat,location==locationnow)
    if (cleandat) {
      indx <- tmpdat$new_cases<0 |
        tmpdat$new_tests<0
      tmpdat$new_cases[indx] <- NA
      tmpdat$new_deaths[indx] <- NA
      tmpdat$new_tests[indx] <- NA
    }
    tmpdat <- tmpdat[order(tmpdat$numDate),]
    tmpdat$meannew_cases <- 
      frollmean(tmpdat$new_cases,nmean,align=align,na.rm=TRUE)
    tmpdat$meannew_cases[is.na(tmpdat$meannew_cases)]<-0
    tmpdat$meannew_casesCumul <- cumsum(tmpdat$meannew_cases)
    tmpdat$meannew_deaths <- 
      frollmean(tmpdat$new_deaths,nmean,align=align,na.rm=TRUE)
    tmpdat$meannew_tests <- 
      frollmean(tmpdat$new_tests,nmean,align=align,na.rm=TRUE)
    meandat <- rbind(meandat,tmpdat)
  }
  dat<-meandat
  dat <- subset(dat,date >= as.Date(mindate))
  if (include_negtest) {
    # dat.df <- melt(as.data.table(dat[,c(
    #   "date","numDate","location","iso_code",
    #   "new_cases",
    #   "new_deaths",
    #   "new_tests",
    #   "meannew_cases",
    #   "meannew_deaths",
    #   "meannew_tests")]),
    #   id.vars=1:4)
    dat.df<-pivot_longer(dat[,c(
      "date","numDate","location","iso_code","population",
      "new_cases",
      "new_deaths",
      "new_tests",
      "meannew_cases",
      "meannew_deaths",
      "meannew_tests",
      "meannew_casesCumul")],
      cols = 6:12,names_to="variable")
  } else {
    # dat.df <- melt(as.data.table(dat[,c(
    #   "date","numDate","location","iso_code","population",
    #   "new_cases",
    #   "new_deaths",
    #   "meannew_cases",
    #   "meannew_deaths")]),
    #   id.vars=1:5)
    dat.df<-pivot_longer(dat[,c(
      "date","numDate","location","iso_code","population",
      "new_cases",
      "new_deaths",
      "meannew_cases",
      "meannew_deaths",
      "meannew_casesCumul")],
      cols = 6:10,names_to="variable")
  }
  dat.df <- dat.df[order(dat.df$variable,dat.df$numDate),]
  return(dat.df)
}


generate_intl_testpositivity <- function(datafolder="International",forceupdate=FALSE,
                                         locationsvec = 
                                           c("Denmark","Greece","Hungary","Luxembourg","Iceland","Slovenia")
) {
  
  datfile<-file.path(datafolder,"EnglandTestPositivityData.csv")
  if (file.exists(datfile) & !forceupdate) {
    englanddat<-fread(datfile)
  } else {
    englanddat.df <- get_england_meantestdata(include_negtest=TRUE,align="right",nmean = 14)
    englanddat <- spread(as.data.table(englanddat.df), variable, value)
    englanddat <- rename(englanddat,
                         location = areaName,
                         new_cases = newCasesByPublishDate,
                         new_tests = newPillarOneTwoTestsByPublishDate,
                         Cases_14 = meannew_cases,
                         Tests_14 = meannew_tests)
    englanddat$PositivePct_14 <- 100*englanddat$Cases_14/englanddat$Tests_14
    englanddat$TestsPer1000_14 <- englanddat$Tests_14*1000/englanddat$population
    englanddat$Cases_14_pop <- englanddat$Cases_14/englanddat$population*100000
    englanddat$Cases_pop <- englanddat$new_cases/englanddat$population*100000
    englanddat$TestsPer1000 <- englanddat$new_tests*1000/englanddat$population
    englanddat$PositivePct <- englanddat$Cases_14_pop/englanddat$TestsPer1000
    englanddat$PositivePct[is.infinite(englanddat$PositivePct)] <- NA
    englanddat$PositivePct[englanddat$PositivePct>100] <- NA
    englanddat$CumulPosPct <- 100*englanddat$meannew_casesCumul/englanddat$population
    
    englanddat <- subset(englanddat,!is.na(new_tests))
    fwrite(englanddat,file=datfile)
  }
  
  datfile<-file.path(datafolder,"OWIDTestPositivityData.csv")
  if (file.exists(datfile) & !forceupdate) {
    owiddat<-fread(datfile)
  } else {
    owiddat.df <- get_owid_meantestdata(include_negtest=TRUE,align="right",nmean = 14)
    owiddat <- spread(as.data.table(owiddat.df), variable, value)
    owiddat <- rename(owiddat,
                      Cases_14 = meannew_cases,
                      Deaths_14 = meannew_deaths,
                      Tests_14 = meannew_tests)
    owiddat$PositivePct_14 <- 100*owiddat$Cases_14/owiddat$Tests_14
    owiddat$TestsPer1000_14 <- owiddat$Tests_14*1000/owiddat$population
    owiddat$Cases_14_pop <- owiddat$Cases_14/owiddat$population*100000
    owiddat$Cases_pop <- owiddat$new_cases/owiddat$population*100000
    owiddat$TestsPer1000 <- owiddat$new_tests*1000/owiddat$population
    owiddat$PositivePct <- owiddat$Cases_14_pop/owiddat$TestsPer1000
    owiddat$PositivePct[is.infinite(owiddat$PositivePct)] <- NA
    owiddat$PositivePct[owiddat$PositivePct>100] <- NA
    owiddat$CumulPosPct <- 100*owiddat$meannew_casesCumul/owiddat$population
    owiddat <- subset(owiddat,(location %in% locationsvec) & !is.na(new_tests))
    fwrite(owiddat,file=datfile)
  }
}



