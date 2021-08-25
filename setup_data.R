library(data.table)
library(tidyr)
library(boot)
library(coda)
library(googledrive)
datezero <- "2019-12-31"
datafolder<-"Data"
functionfolder<-"Functions"
figfolder<-"Figures"
tabfolder<-"Tables"
fips_table <- read.csv(file.path(datafolder,"FIPS_TABLE.csv"),colClasses=c(
  rep("character",4),rep("numeric",2)
))
rownames(fips_table)<-as.character(fips_table$Alpha.code)
statesvec <- fips_table$Alpha.code[2:52]

## Generate data files
source(file.path(functionfolder,"generate_data_functions.R"))
source(file.path(functionfolder,"prediction_likelihood_functions.R"))

# Google drive path
googledrivepath<-"~/COVID-19-US-Semi-Empirical"
# App path
apppath<-"COVID-19-Prevalence-and-Seroprevalence"

## Model options
addCumulPos <- 1
Trec.global <- 10 
Trec.vec.global <- rep(Trec.global,length(statesvec))
tau <- 14 # 14 primary, range 7 - 28
pred.lag <- ceiling(tau/2)
Fixed_n <- FALSE # FALSE primary, simpler model TRUE
n_default <- 0.5 # For fixed, also include 0.4 and 0.6
mcmcpath <- "MCMC"
if (Trec.global != 10) mcmcpath <- paste0(mcmcpath,".Trec",Trec.global)
if (tau != 14) mcmcpath <- paste0(mcmcpath,".tau",tau)
if (pred.lag != ceiling(tau/2)) mcmcpath <- paste0(mcmcpath,".lag",pred.lag)
if (Fixed_n) mcmcpath <- paste0(mcmcpath,".n",n_default)

if (!dir.exists(mcmcpath)) dir.create(mcmcpath)
if (!dir.exists(file.path(apppath,mcmcpath))) dir.create(file.path(apppath,mcmcpath))

today<- as.Date("2021-01-05") + tau 

## Load data files
### Test positivity
datfile<-file.path(mcmcpath,"TestPositivityData.csv")
if (file.exists(datfile)) {
  alldat<-fread(datfile)
  alldat$date <- as.Date(alldat$date)
  maxdatadate <- max(alldat$date)
} else {
  generate_testpositivity(tau=tau,mcmcpath=mcmcpath)
  alldat<-fread(datfile)
  alldat$date <- as.Date(alldat$date)
  maxdatadate <- max(alldat$date)
  fwrite(alldat,file.path(apppath,datfile))
  # drive_put(datfile,path=file.path(googledrivepath,datfile))
}

maxdatadate <- as.Date(today)-1
alldat<-subset(alldat,date <= maxdatadate & state %in% statesvec)
alldat$Cases_tau_pct <- 100*alldat$Cases_tau/fips_table[alldat$state,"pop"]

### Seroprevalence
serofile<-file.path(datafolder,"SeroprevalenceData.csv")
serodat<-fread(serofile)
serodat <- subset(serodat, state %in% statesvec)
serodat$Output_Var <- "Seroprevalence"
serodat$date <- as.Date((as.numeric(as.Date(serodat$mindate))+
                           as.numeric(as.Date(serodat$maxdate)))/2 -
                          as.numeric(as.Date(datezero)),
                        origin=datezero)

## Global data frames
n.id <- length(unique(serodat$state))
### COVID testing data - date of testdat includes lag
testdat<-alldat
testdat$id <- as.numeric(factor(testdat$state,levels=statesvec))
testdat$date <- testdat$date-pred.lag
testdat$t <- testdat$numDate-pred.lag
testdat$CumulPosPct <- 100*testdat$meanpositiveCumul/fips_table[testdat$state,"pop"]
testdat <- testdat[order(testdat$id,testdat$t),]
testdat$I.pct<-0
testdat$SP.pct <- 0

### COVID Testing data in matrix-like form
Cases_tau.pct.df <- pivot_wider(testdat[,c("date","id","Cases_tau_pct")],id_cols=1,names_from=2,values_from=3)
Cases_tau.pct.mat<-as.matrix(Cases_tau.pct.df[,-1])
rownames(Cases_tau.pct.mat)<-as.character(Cases_tau.pct.df$date)
PositivePct_tau.df <- pivot_wider(testdat[,c("date","id","PositivePct_tau")],id_cols=1,names_from=2,values_from=3)
PositivePct_tau.mat<-as.matrix(PositivePct_tau.df[,-1])
rownames(PositivePct_tau.mat)<-as.character(PositivePct_tau.df$date)
# For cumulative - use original date
CumulPosPct.df <- pivot_wider(testdat[,c("date","id","CumulPosPct")],id_cols=1,names_from=2,values_from=3)
CumulPosPct.df$date <- CumulPosPct.df$date + pred.lag
CumulPosPct.mat<-as.matrix(CumulPosPct.df[,-1])
rownames(CumulPosPct.mat)<-as.character(CumulPosPct.df$date)

### Serological data
serodat <- subset(serodat, !is.na(serodat$Data.50.))
serodat.df <- data.frame(id=as.numeric(factor(serodat$state,levels=statesvec)),
                         t=floor(as.numeric(as.Date(serodat$date))-
                                   as.numeric(as.Date(datezero))),
                         y=log(serodat$Data.50.),
                         y.se=log(serodat$Data.97.5./serodat$Data.50.)/(qnorm(0.975)))
serodat.df <- serodat.df[order(serodat.df$id,serodat.df$t,serodat.df$y),]

### CDC Seropositivity API for validation data

cdcserofile <- file.path(datafolder,"CDC_Sero_Validation.csv")
if (file.exists(cdcserofile)) {
  cdcseroval <- fread(cdcserofile)
} else {
  #cdcserorawfile <- getURL("https://data.cdc.gov/resource/d2tw-32xv.csv")
  #cdcseroval  <- read.csv(text = cdcserorawfile,as.is=TRUE)
  #write.csv(cdcseroval,file.path(datafolder,"d2tw-32xv.csv"),row.names = FALSE)
  cdcseroval <- fread(file.path(datafolder,"d2tw-32xv.csv"))
  daterangemat <- str_split(cdcseroval$date_range_of_specimen," - ",simplify = TRUE)
  cdcseroval$maxdate <- as.character(as.Date(trimws(daterangemat[,2]),tryFormats = "%b %d, %Y"))
  cdcseroval$mindate <- as.character(as.Date(trimws(paste(daterangemat[,1],
                                                          format(as.Date(cdcseroval$maxdate),", %Y"))),tryFormats = "%b %d , %Y"))
  cdcseroval$date <- as.Date((as.numeric(as.Date(cdcseroval$mindate))+
                                as.numeric(as.Date(cdcseroval$maxdate)))/2 -
                               as.numeric(as.Date(datezero)),
                             origin=datezero)
  cdcseroval<-rename(cdcseroval,
                     state=site,
                     Data.50.=rate_cumulative_prevalence,
                     Data.2.5.=lower_ci_cumulative_prevalence,
                     Data.97.5.=upper_ci_cumulative_prevalence)
  cdcseroval <- subset(cdcseroval, !is.na(Data.50.) & Data.50. <= 100 &
                         state %in% statesvec & date <= today & round > 7 &
                         round < 12)
  cdcseroval$Source <- "Validation CDC (2021)"
  write.csv(cdcseroval,file.path(datafolder,"CDC_Sero_Validation.csv"),row.names = FALSE)
}
