library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggstance)
library(lubridate)
source("setup_data.R")

for (mcmcpath in c("MCMC","MCMC.n0.5")) {
  predquant <- as.data.frame(fread(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv")))
  predquant$date <- as.Date(predquant$date)
  uspredquant <- as.data.frame(fread(file.path(mcmcpath,"US_Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv")))
  uspredquant$date <- as.Date(uspredquant$date)
  
  alldat$pop <- fips_table[alldat$state,"pop"]
  alldat$posPct <- 100*alldat$meanpositiveCumul/alldat$pop
  alldat$Cases_tau_pct.weighted <- alldat$Cases_tau_pct*alldat$pop
  alldat$TestsPer1000_tau.weighted <- alldat$TestsPer1000_tau*alldat$pop
  # 
  maxdate <- as.Date("2020-12-31") # max(predquant$date)
  predquant <- subset(predquant,date<=maxdate)
  uspredquant <- subset(uspredquant,date<=maxdate)
  alldat <- subset(alldat,date<=(maxdate+7))
  
  usdat <- aggregate(. ~ date,data=alldat[,c("date","meanpositiveCumul","Cases_tau_pct.weighted",
                                             "TestsPer1000_tau.weighted")],sum)
  usdat<-rename(usdat,Cases_tau_pct=Cases_tau_pct.weighted)
  usdat<-rename(usdat,TestsPer1000_tau=TestsPer1000_tau.weighted)
  usdat$state<-"US"
  usdat$posPct<- 100*usdat$meanpositiveCumul/sum(fips_table[statesvec,"pop"])
  usdat$PositivePct_tau <- 1000*usdat$Cases_tau_pct/usdat$TestsPer1000_tau
  usdat$Cases_tau_pct<- usdat$Cases_tau_pct/sum(fips_table[statesvec,"pop"])
  
  uspredquant$origdate <- uspredquant$date + pred.lag
  usdat$origdate <- usdat$date
  usbiasquant <- merge.data.frame(usdat[,c("origdate","state","TestsPer1000_tau",
                                          "PositivePct_tau","Cases_tau_pct")],
                                uspredquant[,-1])
  usbiasquant$b_pos.2.5. <- usbiasquant$PositivePct_tau/usbiasquant$I.pct.97.5.
  usbiasquant$b_pos.50. <- usbiasquant$PositivePct_tau/usbiasquant$I.pct.50.
  usbiasquant$b_pos.97.5. <- usbiasquant$PositivePct_tau/usbiasquant$I.pct.2.5.
  usbiasquant$month <- as.factor(month(usbiasquant$origdate,label=TRUE))
  print(aggregate(b_pos.50. ~ month,data=usbiasquant,median))
  
  uspredquant$serodate <- uspredquant$date
  usdat$serodate <- usdat$date + 14
  
  uscumbiasquant <- merge.data.frame(usdat[,c("serodate","state","meanpositiveCumul")],
                                   uspredquant[,-1])
  uscumbiasquant$CumCasesPct <- uscumbiasquant$meanpositiveCumul/sum(fips_table[statesvec,"pop"])*100
  
  uscumbiasquant$b_cumcase.2.5. <- uscumbiasquant$CumCasesPct/uscumbiasquant$SP.pct.97.5.
  uscumbiasquant$b_cumcase.50. <- uscumbiasquant$CumCasesPct/uscumbiasquant$SP.pct.50.
  uscumbiasquant$b_cumcase.97.5. <- uscumbiasquant$CumCasesPct/uscumbiasquant$SP.pct.2.5.
  
  uscumbiasquant$month <- as.factor(month(uscumbiasquant$serodate,label=TRUE))
  
  print(aggregate(b_cumcase.50.~month,data=uscumbiasquant,median))
  
  ### Prevalence
  predquant$origdate <- predquant$date + pred.lag
  alldat$origdate <- alldat$date
  
  biasquant <- merge.data.frame(alldat[,c("origdate","state","TestsPer1000_tau",
                                          "PositivePct_tau","Cases_tau_pct")],
                                predquant[,-1])
  biasquant <- subset(biasquant,(origdate>="2020-04-01") & (origdate <= (maxdate-7)))
  biasquant$month <- as.factor(month(biasquant$origdate,label=TRUE))
  
  bdf <- data.frame(b=(1/c(1,10,100)),
                    TestingBias=factor(paste(c(1,10,100)),
                                       levels=paste(c(1,10,100))),
                    CumCaseBias=factor(paste(c(1,0.1,0.01)),
                                       levels=paste(c(1,0.1,0.01)))) 
  
  pconcept<-ggplot(bdf)+
    geom_errorbar(data=biasquant,aes(x=PositivePct_tau,ymin=I.pct.2.5.,ymax=I.pct.97.5.,color=month),alpha=0.25)+
    geom_point(data=biasquant,aes(x=PositivePct_tau,y=I.pct.50.,color=month),size=0.2,alpha=0.25)+
    geom_label(data=subset(biasquant,origdate==max(biasquant$origdate)),
               aes(x=PositivePct_tau,y=I.pct.50.,label=state),alpha=0.5,size=2,label.size=0)+
    geom_abline(aes(intercept=b,slope=1,linetype=TestingBias),
                color="grey")+
    scale_x_log10(limits=c(0.01,120),
                  breaks=10^seq(-2,2),
                  labels=10^seq(-2,2))+
    scale_y_log10(limits=c(0.001,12),
                  breaks=10^seq(-3,1),
                  labels=10^seq(-3,1))+
    # geom_boxploth(aes(y=1/1000,PositivePct_tau,
    #                   color="U.S. State Positivity Rates"),
    #               data=subset(alldat,TestsPer1000_tau>0),width=0.1,
    #               outlier.size=0.05)+
    annotation_logticks(short=unit(0.05, "cm"),
                        mid = unit(0.1, "cm"),
                        long = unit(0.15, "cm"))+
    scale_linetype_discrete(name="Positivity\nRate Bias (b)")+
    scale_color_viridis_d(begin=1,end=0.5,name=" Date")+
    theme_bw()+
    guides(linetype = guide_legend(order=1),
           color = guide_legend(order=2))+
    # theme(legend.position = "bottom",
    #                  legend.title = element_text(size = 10),
    #                  legend.text = element_text(size = 10),
    #                  legend.box="vertical")+
    xlab("Test Positivity % (2 wk average)")+
    ylab("Estimated Prevalence of Undiagnosed Infection %")
  
  uspconcept<-ggplot(bdf)+
    geom_errorbar(data=usbiasquant,aes(x=PositivePct_tau,ymin=I.pct.2.5.,ymax=I.pct.97.5.,color=month),alpha=0.25)+
    geom_point(data=usbiasquant,aes(x=PositivePct_tau,y=I.pct.50.,color=month),alpha=0.25)+
    geom_label(data=subset(usbiasquant,origdate==max(biasquant$origdate)),
               aes(x=PositivePct_tau,y=I.pct.50.,label=state),alpha=0.5)+
    geom_abline(aes(intercept=b,slope=1,linetype=TestingBias),
                color="grey")+
    scale_x_log10(limits=c(0.01,120),
                  breaks=10^seq(-2,2),
                  labels=10^seq(-2,2))+
    scale_y_log10(limits=c(0.001,12),
                  breaks=10^seq(-3,1),
                  labels=10^seq(-3,1))+
    # geom_boxploth(aes(y=1/1000,PositivePct_tau,
    #                   color="U.S. State Positivity Rates"),
    #               data=subset(alldat,TestsPer1000_tau>0),width=0.1,
    #               outlier.size=0.05)+
    annotation_logticks(short=unit(0.05, "cm"),
                        mid = unit(0.1, "cm"),
                        long = unit(0.15, "cm"))+
    scale_linetype_discrete(name="Positivity\nRate Bias (b)")+
    scale_color_viridis_d(begin=1,end=0.5,name=" Date")+
    theme_bw()+
    guides(linetype = guide_legend(order=1),
           color = guide_legend(order=2))+
    # theme(legend.position = "bottom",
    #                  legend.title = element_text(size = 10),
    #                  legend.text = element_text(size = 10),
    #                  legend.box="vertical")+
    xlab("Test Positivity % (2 wk average)")+
    ylab("Estimated Prevalence of Undiagnosed Infection %")
  
  
  # ggplot(biasquant) + 
  #   geom_point(aes(x=PositivePct_tau,y=I.pct.50.,color=origdate))+
  #   scale_y_log10()+scale_x_log10()+
  #   scale_color_viridis_c(trans="date",begin=1,end=0.5)
  # 
  biasquant$b_pos.2.5. <- biasquant$PositivePct_tau/biasquant$I.pct.97.5.
  biasquant$b_pos.50. <- biasquant$PositivePct_tau/biasquant$I.pct.50.
  biasquant$b_pos.97.5. <- biasquant$PositivePct_tau/biasquant$I.pct.2.5.
  # 
  # biasquant$b_case.2.5. <- biasquant$Cases_tau_pct/biasquant$I.pct.97.5.
  # biasquant$b_case.50. <- biasquant$Cases_tau_pct/biasquant$I.pct.50.
  # biasquant$b_case.97.5. <- biasquant$Cases_tau_pct/biasquant$I.pct.2.5.
  # 
  
  print(aggregate(b_pos.50. ~ month,data=biasquant,median))
  
  print(aggregate(b_pos.50. ~ state,data=subset(biasquant,month=="Dec"),median))
  
  
  
  # biastrend<- ggplot(subset(biasquant,month>="2020-04-01"))+
  #   geom_boxplot(aes(x=month,y=b_pos.50.,color="Positivity (14-day average)",group=month),
  #                outlier.shape = NA)+
  #   geom_boxplot(aes(x=month,y=b_case.50.,color="Reported Cases (14-day average)",group=month),
  #                outlier.shape = NA)+
  #   scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)))+
  #   annotation_logticks(sides="l")+
  #   scale_color_viridis_d(begin=0.2,end=0.7,name="")+
  #   ggtitle("Bias in estimating infection prevalence")+
  #   geom_hline(yintercept=1)+
  #   ylab("Bias (posterior median)")+
  #   scale_x_date(minor_breaks = "month")+theme_bw()+theme(legend.position = "bottom")
  
  ### Seroprevalence
  
  predquant$serodate <- predquant$date
  alldat$serodate <- alldat$date + 14
  
  cumbiasquant <- merge.data.frame(alldat[,c("serodate","state","meanpositiveCumul")],
                                   predquant[,-1])
  cumbiasquant <- subset(cumbiasquant,(serodate>="2020-04-01") & (serodate <= maxdate))
  cumbiasquant$CumCasesPct <- cumbiasquant$meanpositiveCumul/fips_table[cumbiasquant$state,"pop"]*100
  
  cumbiasquant$b_cumcase.2.5. <- cumbiasquant$CumCasesPct/cumbiasquant$SP.pct.97.5.
  cumbiasquant$b_cumcase.50. <- cumbiasquant$CumCasesPct/cumbiasquant$SP.pct.50.
  cumbiasquant$b_cumcase.97.5. <- cumbiasquant$CumCasesPct/cumbiasquant$SP.pct.2.5.
  
  cumbiasquant$month <- as.factor(month(cumbiasquant$serodate,label=TRUE))
  
  print(aggregate(b_cumcase.50.~month,data=
                    aggregate(b_cumcase.50. ~ month+state,data=subset(cumbiasquant),median),
                  range))
  
  # cumbiastrend<- ggplot(cumbiasquant)+
  #   geom_boxplot(aes(x=month,y=b_cumcase.50.,color="Cumulative Reported Cases (14-day delay)",group=month),
  #                outlier.shape = NA)+
  #   scale_y_log10(breaks=trans_breaks("log10", function(x) 10^x),
  #                 labels = trans_format("log10", math_format(10^.x)))+
  #   annotation_logticks(sides="l")+
  #   scale_color_viridis_d(begin=0.2,end=0.7,name="")+
  #   ggtitle("Bias in estimating seroprevalence")+
  #   geom_hline(yintercept=1)+
  #   ylab("Bias (posterior median)")+
  #   scale_x_date(minor_breaks = "month")+theme_bw()+theme(legend.position = "bottom")
  
  
  pserobias<-ggplot(bdf)+
    geom_errorbar(data=cumbiasquant,aes(x=CumCasesPct,ymin=SP.pct.2.5.,ymax=SP.pct.97.5.,color=month),alpha=0.25)+
    geom_point(data=cumbiasquant,aes(x=CumCasesPct,y=SP.pct.50.,color=month),size=0.2,alpha=0.25)+
    geom_label(data=subset(cumbiasquant,serodate==max(cumbiasquant$serodate)),
               aes(x=CumCasesPct,y=SP.pct.50.,label=state),alpha=0.5,size=2,label.size=0)+
    geom_abline(aes(intercept=1/b,slope=1,linetype=CumCaseBias),
                color="grey")+
    scale_x_log10(limits=c(0.01,120),
                  breaks=10^seq(-2,2),
                  labels=10^seq(-2,2))+
    scale_y_log10(limits=c(0.01,120),
                  breaks=10^seq(-2,2),
                  labels=10^seq(-2,2))+
    # geom_boxploth(aes(y=1/1000,PositivePct_tau,
    #                   color="U.S. State Positivity Rates"),
    #               data=subset(alldat,TestsPer1000_tau>0),width=0.1,
    #               outlier.size=0.05)+
    annotation_logticks(short=unit(0.05, "cm"),
                        mid = unit(0.1, "cm"),
                        long = unit(0.15, "cm"))+
    scale_linetype_discrete(name="Cumulative\nCase Under-\nreporting")+
    scale_color_viridis_d(begin=1,end=0.5,name=" Date")+
    theme_bw()+
    guides(linetype = guide_legend(order=1),
           color = guide_legend(order=2))+
    # theme(legend.position = "bottom",
    #                  legend.title = element_text(size = 10),
    #                  legend.text = element_text(size = 10),
    #                  legend.box="vertical")+
    xlab("Cumulative Cases (2-wk past) %")+
    ylab("Estimated Seroprevalence %")
  
  uspserobias<-ggplot(bdf)+
    geom_errorbar(data=uscumbiasquant,aes(x=CumCasesPct,ymin=SP.pct.2.5.,ymax=SP.pct.97.5.,color=month),alpha=0.25)+
    geom_point(data=uscumbiasquant,aes(x=CumCasesPct,y=SP.pct.50.,color=month),alpha=0.25)+
    geom_label(data=subset(uscumbiasquant,serodate==max(uscumbiasquant$serodate)),
               aes(x=CumCasesPct,y=SP.pct.50.,label=state),alpha=0.5)+
    geom_abline(aes(intercept=1/b,slope=1,linetype=CumCaseBias),
                color="grey")+
    scale_x_log10(limits=c(0.01,120),
                  breaks=10^seq(-2,2),
                  labels=10^seq(-2,2))+
    scale_y_log10(limits=c(0.01,120),
                  breaks=10^seq(-2,2),
                  labels=10^seq(-2,2))+
    # geom_boxploth(aes(y=1/1000,PositivePct_tau,
    #                   color="U.S. State Positivity Rates"),
    #               data=subset(alldat,TestsPer1000_tau>0),width=0.1,
    #               outlier.size=0.05)+
    annotation_logticks(short=unit(0.05, "cm"),
                        mid = unit(0.1, "cm"),
                        long = unit(0.15, "cm"))+
    scale_linetype_discrete(name="Cumulative\nCase Under-\nreporting")+
    scale_color_viridis_d(begin=1,end=0.5,name=" Date")+
    theme_bw()+
    guides(linetype = guide_legend(order=1),
           color = guide_legend(order=2))+
    # theme(legend.position = "bottom",
    #                  legend.title = element_text(size = 10),
    #                  legend.text = element_text(size = 10),
    #                  legend.box="vertical")+
    xlab("Cumulative Cases (2-wk past) %")+
    ylab("Estimated Seroprevalence %")
  
  figbias<-ggarrange(pconcept,uspconcept,ncol=2,legend="right",
                     common.legend = TRUE,labels=c("A","B"))
  figserobias<-ggarrange(pserobias,uspserobias,ncol=2,legend="right",
                         common.legend = TRUE,labels=c("C","D"))
  figbiasall<-ggarrange(figbias,figserobias,ncol=1)
  print(figbiasall)
  
  if (mcmcpath == "MCMC") {
    figname.pdf <- file.path(figfolder,"SupFig_Biases.pdf")
    figname.jpg <- file.path(figfolder,"SupFig_Biases.jpeg")
  } else {
    figname.pdf <- file.path(figfolder,
                             paste0("SupFig_Biases.",mcmcpath,".pdf"))
    figname.jpg <- file.path(figfolder,
                             paste0("SupFig_Biases.",mcmcpath,".jpeg"))
  }
  
  ggsave(figname.pdf,figbiasall,height=9,width=10,scale=1)
  ggsave(figname.jpg,figbiasall,height=9,width=10,scale=1,dpi=600)
}


