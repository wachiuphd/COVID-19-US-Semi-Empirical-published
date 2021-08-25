source("setup_data.R")
library(ggstance)
library(ggpubr)
for (mcmcpath in c("MCMC","MCMC.n0.5")) {
  predquant <- as.data.frame(fread(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv")))
  predquant$date <- as.Date(predquant$date)
  uspredquant <- as.data.frame(fread(file.path(mcmcpath,"US_Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv")))
  uspredquant$date <- as.Date(uspredquant$date)
  maxdate <- as.Date("2020-12-31") # max(predquant$date)
  
  pprev<-ggplot() + 
    geom_ribbon(data=predquant,
                aes(x=date,ymin=I.pct.2.5.,ymax=I.pct.97.5.,
                    fill="Posterior 95% CrI"),alpha=0.5)+
    geom_line(data=predquant,
              aes(x=date,y=I.pct.50.,color="Posterior median"))+
    facet_wrap(~state)+
    xlab("Date")+
    ylab("Infection prevalence %")+
    scale_x_date(date_minor_breaks="1 month",limits=c(as.Date("2020-09-01"),maxdate))+
    #  scale_y_log10()+annotation_logticks(side="l")+coord_cartesian(ylim=c(0.1,100))+
    ggtitle("State-wide Infection prevalence")+
    scale_color_viridis_d(begin=0.3,name="Prevalence Predictions")+
    scale_fill_viridis_d(begin=0.5,name="")+
    guides(shape = guide_legend(order=1),
           color = guide_legend(order=2),
           fill = guide_legend(order=3))+
    theme_bw()+
    theme(legend.position="right")
  # print(pprev)
  # ggsave(file.path(figfolder,"RecentPrevalence.pdf"),pprev,height=8,width=10)
  
  
  psero<-ggplot() + 
    geom_ribbon(data=predquant,
                aes(x=date,ymin=SP.pct.2.5.,ymax=SP.pct.97.5.,
                    fill="Posterior 95% CrI"),alpha=0.5)+
    geom_line(data=predquant,
              aes(x=date,y=SP.pct.50.,color="Posterior median"))+
    facet_wrap(~state)+
    xlab("Date")+
    ylab("Cumulative % Infected")+
    scale_x_date(date_minor_breaks="1 month",limits=c(as.Date("2020-09-01"),maxdate))+
    #  scale_y_log10()+annotation_logticks(side="l")+coord_cartesian(ylim=c(0.1,100))+
    ggtitle("State-wide Seroprevalence")+
    scale_shape_discrete(name="Seroprevalence\nData Source")+
    scale_color_viridis_d(begin=0.3,name="Seroprevalence Predictions")+
    scale_fill_viridis_d(begin=0.5,name="")+
    guides(shape = guide_legend(order=1),
           color = guide_legend(order=2),
           fill = guide_legend(order=3))+
    theme_bw()+
    theme(legend.position="right")
  # print(psero)
  # ggsave(file.path(figfolder,"RecentSeroprevalence.pdf"),psero,height=8,width=10)
  
  
  currentpred <- rbind(subset(predquant,date==maxdate),
                       subset(uspredquant,date==maxdate))
  rownames(currentpred)<-currentpred$state
  currentpred <- currentpred[c(statesvec,"US"),]
  currentpred$state <- factor(currentpred$state,levels=c(statesvec,"US"))
  pprevnow<-ggplot(currentpred,aes(y=state,xmin = I.pct.2.5., xlower = I.pct.25., 
                                   xmiddle = I.pct.50., xupper = I.pct.75., xmax = I.pct.97.5.))+
    geom_boxploth(stat="identity",aes(fill=(state=="US")))+
    scale_y_discrete(limits=rev)+#xlim(0,1)+
    xlab(paste("Undiagnosed Infection Prevalence",maxdate))+theme_bw()+theme(legend.position = "none") + 
    scale_fill_viridis_d(begin=1,end=0.7,option = "magma")
  print(pprevnow)
  
  pprevtotnow<-ggplot(currentpred,aes(y=state,xmin = Itot.pct.2.5., xlower = Itot.pct.25., 
                                   xmiddle = Itot.pct.50., xupper = Itot.pct.75., xmax = Itot.pct.97.5.))+
    geom_boxploth(stat="identity",aes(fill=(state=="US")))+
    scale_y_discrete(limits=rev)+#xlim(0,1)+
    xlab(paste("Total Infection Prevalence",maxdate))+theme_bw()+theme(legend.position = "none") + 
    scale_fill_viridis_d(begin=1,end=0.7,option = "magma")
  print(pprevtotnow)
  
  
  currentpred$isus <- ifelse(currentpred$state=="US",
                             "US total\nposterior","State-level\nposterior")
  
  
  pprevbothnow<-ggplot(currentpred,aes(y=state,xmin = Itot.pct.2.5., xlower = Itot.pct.25., 
                         xmiddle = Itot.pct.50., xupper = Itot.pct.75., xmax = Itot.pct.97.5.))+
    geom_boxploth(stat="identity",aes(fill=isus,color="Total"))+
    geom_boxploth(stat="identity",aes(y=state,xmin = I.pct.2.5., xlower = I.pct.25., 
                                      xmiddle = I.pct.50., xupper = I.pct.75., xmax = I.pct.97.5.,fill=isus,color="Undiagnosed"))+
    scale_y_discrete(limits=rev)+#xlim(0,1)+
    xlab(paste("Total Infection Prevalence",maxdate))+theme_bw()+ 
    guides(color=guide_legend(order=1),fill = "none")+
    theme(legend.position=c(0.75,0.05),
          legend.background = element_rect(fill=alpha("white",0.75)))+
    scale_fill_viridis_d(name="",begin=1,end=0.7,option = "magma")+
    scale_color_viridis_d(name="",begin=0,end=0.7,option = "plasma")
  
  tmp.14<-subset(alldat,date==(currentpred$date[1]-14))
  tmp.14$CumCase.14 <- 100*tmp.14$meanpositiveCumul/fips_table[tmp.14$state,"pop"]
  tmp.7<-subset(alldat,date==(currentpred$date[1]-7))
  tmp.7$CumCase.7 <- 100*tmp.7$meanpositiveCumul/fips_table[tmp.7$state,"pop"]
  tmp.21<-subset(alldat,date==(currentpred$date[1]-21))
  tmp.21$CumCase.21 <- 100*tmp.21$meanpositiveCumul/fips_table[tmp.21$state,"pop"]
  
  
  
  
  current_cumcase <- merge.data.frame(
    merge.data.frame(tmp.14[,c("state","CumCase.14")],
                     tmp.7[,c("state","CumCase.7")]),
    tmp.21[,c("state","CumCase.21")])
  current_cumcase <- rbind(current_cumcase,
                           data.frame(state="US",
                                      CumCase.14=
                                        weighted.mean(current_cumcase$CumCase.14,
                                                      fips_table[current_cumcase$state,
                                                                 "pop"]),
                                      CumCase.7=
                                        weighted.mean(current_cumcase$CumCase.7,
                                                      fips_table[current_cumcase$state,
                                                                 "pop"]),
                                      CumCase.21=
                                        weighted.mean(current_cumcase$CumCase.21,
                                                      fips_table[current_cumcase$state,
                                                                 "pop"])
                           ))
  
  pseronow<-ggplot(currentpred)+
    geom_boxploth(stat="identity",aes(y=state,xmin = SP.pct.2.5., xlower = SP.pct.25., 
                                      xmiddle = SP.pct.50., xupper = SP.pct.75., xmax = SP.pct.97.5.,
                                      fill=isus))+
    geom_point(aes(x=CumCase.14,
                   y=state,alpha="Reported\nCases"),
               data=current_cumcase)+
    geom_errorbarh(aes(y=state,xmin = CumCase.21,
                       xmax = CumCase.7,alpha="Reported\nCases"),
                   data=current_cumcase)+
    scale_y_discrete(limits=rev)+
    xlab(paste("Seroprevalence",maxdate))+theme_bw()+
    scale_alpha_discrete(name="",range=c(0.5,1))+
    xlim(NA,30)+
    guides(alpha=guide_legend(order=1),fill = "none")+
    theme(legend.position=c(0.75,0.05),
          #legend.box.background = element_rect(colour = "grey"),
          legend.background = element_rect(fill=alpha("white",0.75)))+
    scale_fill_viridis_d(begin=1,end=0.7,option = "magma",name="")
  print(pseronow)
  
  print(ggarrange(pprevbothnow,pseronow,ncol=2,labels="AUTO"))
  if (mcmcpath == "MCMC") {
    figname.pdf <- file.path(figfolder,"SupFig_CurrentEstimates.pdf")
    figname.jpg <- file.path(figfolder,"SupFig_CurrentEstimates.jpeg")
  } else {
    figname.pdf <- file.path(figfolder,
                             paste0("SupFig_CurrentEstimates.",mcmcpath,".pdf"))
    figname.jpg <- file.path(figfolder,
                             paste0("SupFig_CurrentEstimates.",mcmcpath,".jpeg"))
  }
  
  ggsave(figname.pdf,
         plot=ggarrange(pprevbothnow,pseronow,ncol=2,labels="AUTO"),height=8,width=6)
  ggsave(figname.jpg,
         plot=ggarrange(pprevbothnow,pseronow,ncol=2,labels="AUTO"),height=8,width=6,dpi=600)
}
