library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggstance)
source("setup_data.R")

mcmcpath <- "MCMC"
predquant<-fread(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
predquant<-as.data.frame(predquant)
predquant$date <- as.Date(predquant$date)

mcmcpath.n0.5 <- "MCMC.n0.5"
predquant.n0.5<-fread(file.path(mcmcpath.n0.5,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
predquant.n0.5<-as.data.frame(predquant.n0.5)
predquant.n0.5$date <- as.Date(predquant.n0.5$date)

exdat <- subset(alldat,state %in% fips_table$Alpha.code & date <= "2021-01-15")
exdat$date <- as.Date(exdat$date)
exdat$TestsPct_tau <- exdat$TestsPer1000_tau/10

exdat.n0.5 <- merge.data.frame(exdat,predquant.n0.5)
exdat <- merge.data.frame(exdat,predquant)

explot<-ggplot(subset(exdat,date>=as.Date("2020-04-01")))+
  geom_line(aes(x=date,y=PositivePct_tau,color="Positivity %"))+
  geom_line(aes(x=date,y=Cases_tau_pct,color="Cases %"))+
  geom_line(aes(x=date,y=TestsPct_tau,color="Tests %"),linetype="dotted")+
  geom_line(aes(x=date,y=Ipct.50.,color="Prevalence-Undiagnosed"),linetype="dashed")+
  geom_line(aes(x=date,y=Itot.pct.50.,color="Prevalence-Total"),linetype="dashed")+
#  geom_line(aes(x=date,y=I_tot.gm.lag/10000,color="I_tot per 10 (GM)"),linetype="dashed")+
  scale_color_viridis_d(end=0.8)+
  scale_y_log10()+
  ylab("")+
  facet_wrap(~state,scales="free_y")
#print(explot)
#ggsave(file.path(figfolder,"ExamplePlots.pdf"),explot,height=8,width=12)
explot2<-ggplot(subset(exdat,date>=as.Date("2020-09-01")))+
  geom_point(aes(x=Cases_tau_pct,y=PositivePct_tau,color=date))+
  #geom_line(aes(x=date,y=Cases_tau_pct*1000,color="1000*Cases %"))+
#  geom_line(aes(x=date,y=TestsPct_tau,color="Tests %"),linetype="dotted")+
#  geom_line(aes(x=date,y=Itot.pct.50.,color="Prevalence"),linetype="dashed")+
  #  geom_line(aes(x=date,y=I_tot.gm.lag/10000,color="I_tot per 10 (GM)"),linetype="dashed")+
  scale_color_viridis_c(end=0.8)+
  scale_x_log10()+scale_y_log10()+
  ylab("")+
  facet_wrap(~state,scales="free_y")
# NY Times - increasing ND, HI, KS, IL, DE, SD, WY, ME, VT
# NY Times - same TN, MO, CA, IA, KY, *WI*, NE, IN, *MN*, *VA*, *MT*, *RI*, *WA*, DC, MI, WV, OR, PA, NJ, *NY*, CT
# NY Times - decreasing TX, MS, GA, NV, FL, ID, AL, AR, LA, OK, SC, NC, AZ, *AK*, UT, MD, OH, NM, CO, MA, NH
# c("WI","MN","VA","MT","RI","WA","NY","AK")

pex.list <- list()
for (statenow in c("MN","VA","WI","KY","TN")) {
  if (statenow == "MN" | statenow == "VA" | statenow == "WI") {
    dat.ex <- subset(exdat, state == statenow & date >= "2020-04-15" & date < "2020-06-15")
    dat.ex.n0.5 <- subset(exdat.n0.5, state == statenow & date >= "2020-04-15" & date < "2020-06-15")
  }
  if (statenow == "TN" | statenow == "KY") {
    dat.ex <- subset(exdat, state == statenow & date >= "2020-11-01" & date < "2021-01-01")
    dat.ex.n0.5 <- subset(exdat.n0.5, state == statenow & date >= "2020-11-01" & date < "2021-01-01")
  }
  # neff.ex <- coef(sep.res$res.re)$state[statenow,]
  # dat.ex$I_tot.re <- c(dat.ex$I_tot.re.lag[-(1:pred.lag)],rep(NA,pred.lag))
  # dat.ex$I_tot.gm <- c(dat.ex$I_tot.gm.lag[-(1:pred.lag)],rep(NA,pred.lag))
  dat.ex.df <- pivot_longer(dat.ex[,c("date","state",
                                      "Cases_tau_pct","PositivePct_tau","TestsPct_tau",
                                      "I.pct.50.","I.pct.2.5.","I.pct.97.5.",
                                      "Itot.pct.50.","Itot.pct.2.5.","Itot.pct.97.5.")],
                            cols = 3:11)
  dat.ex.df$name <- factor(dat.ex.df$name,levels = 
                             c("Cases_tau_pct","PositivePct_tau",
                               "Itot.pct.50.","I.pct.50.","TestsPct_tau",
                               "Itot.pct.2.5.","Itot.pct.97.5.",
                               "I.pct.2.5.","I.pct.97.5."))
  dat.ex.n0.5.df <- pivot_longer(dat.ex.n0.5[,c("date","state",
                                      "Cases_tau_pct","PositivePct_tau","TestsPct_tau",
                                      "I.pct.50.","I.pct.2.5.","I.pct.97.5.",
                                      "Itot.pct.50.","Itot.pct.2.5.","Itot.pct.97.5.")],
                            cols = 3:11)
  dat.ex.n0.5.df$name <- factor(dat.ex.n0.5.df$name,levels = 
                             c("Cases_tau_pct","PositivePct_tau",
                               "Itot.pct.50.","I.pct.50.","TestsPct_tau",
                               "Itot.pct.2.5.","Itot.pct.97.5.",
                               "I.pct.2.5.","I.pct.97.5."))
  labvec <- c(Cases_tau_pop="Reported Cases %",
              PositivePct_tau="Positivity %",
              Itot.pct.50.="Total\n Prevalence %",
              I.pct.50.="Undiagnosed\n Prevalence %",
              TestsPct_tau="Daily Tests %",
              Itot.pct.2.5.="Total\n Prevalence %",
              Itot.pct.97.5.="Total\n Prevalence %",
              I.pct.2.5.="Undiagnosed\n Prevalence %",
              I.pct.97.5.="Undiagnosed\n Prevalence %")
  linevec<-c(Cases_tau_pop="Data",
             PositivePct_tau="Data",
             Itot.pct.50.="Posterior",
             I.pct.50.="Posterior",
             TestsPct_tau="Data",
             Itot.pct.2.5.="Posterior",
             Itot.pct.97.5.="Posterior",
             I.pct.2.5.="Posterior",
             I.pct.97.5.="Posterior")
  dat.ex.df$label <- factor(labvec[dat.ex.df$name],levels=
                              c("Total\n Prevalence %",
                                "Undiagnosed\n Prevalence %",
                                "Reported Cases %",
                                "Positivity %",
                                "Daily Tests %"))
  dat.ex.n0.5.df$label <- factor(labvec[dat.ex.n0.5.df$name],levels=
                              c("Total\n Prevalence %",
                                "Undiagnosed\n Prevalence %",
                                "Reported Cases %",
                                "Positivity %",
                                "Daily Tests %"))
  
  dat.ex.df$colorlab <- factor(linevec[dat.ex.df$name],
                               levels=c("Posterior","Data"))
  dat.ex.n0.5.df$colorlab <- factor(linevec[dat.ex.n0.5.df$name],
                               levels=c("Posterior","Data"))
  dat.ex.df$Divergence <- ""
  dat.ex.n0.5.df$Divergence <- ""
  if (statenow != "TN" & statenow != "KY") {
    dat.ex.df$Divergence[(dat.ex.df$date >= "2020-05-01") &
                               (dat.ex.df$date <= "2020-05-31")] <- "May"
    dat.ex.n0.5.df$Divergence[(dat.ex.n0.5.df$date >= "2020-05-01") &
                           (dat.ex.n0.5.df$date <= "2020-05-31")] <- "May"
  } else {
    if (statenow == "TN") {
      dat.ex.df$Divergence[(dat.ex.df$date >= "2020-12-20") &
                           (dat.ex.df$date <= "2020-12-31")] <- "December"
      dat.ex.n0.5.df$Divergence[(dat.ex.n0.5.df$date >= "2020-12-20") &
                             (dat.ex.n0.5.df$date <= "2020-12-31")] <- "December"
    } else {
      dat.ex.df$Divergence[(dat.ex.df$date >= "2020-12-10") &
                             (dat.ex.df$date <= "2020-12-31")] <- "December"
      dat.ex.n0.5.df$Divergence[(dat.ex.n0.5.df$date >= "2020-12-10") &
                             (dat.ex.n0.5.df$date <= "2020-12-31")] <- "December"
    }
  }
  dat.ex.df$Divergence<-factor(dat.ex.df$Divergence,
                               levels=c("May","December",""))
  dat.ex.df$linetype <- NA
  dat.ex.df$linetype[dat.ex.df$Divergence != ""] <- "Divergence"
  dat.ex.df$linetype[dat.ex.df$name == "Itot.pct.2.5." | dat.ex.df$name == "Itot.pct.97.5." |
                       dat.ex.df$name == "I.pct.2.5." | dat.ex.df$name == "I.pct.97.5."] <- "95% CrI"
  
  dat.ex.n0.5.df$Divergence<-factor(dat.ex.n0.5.df$Divergence,
                               levels=c("May","December",""))
  dat.ex.n0.5.df$linetype <- "n=0.5"
  pex<-ggplot(subset(dat.ex.df,(name != "Itot.pct.2.5.") & (name != "Itot.pct.97.5.") &
                       (name != "I.pct.2.5.") & (name != "I.pct.97.5.")))+
    geom_ribbon(aes(x=date,ymin=0,ymax=value,fill=colorlab))+
    geom_line(aes(x=date,y=value,color=colorlab,linetype=linetype),
              subset(dat.ex.n0.5.df,(name != "Itot.pct.2.5.") & (name != "Itot.pct.97.5.") &
                       (name != "I.pct.2.5.") & (name != "I.pct.97.5.")),alpha=0.5)+
    geom_line(aes(x=date,y=value,color=colorlab,linetype=linetype),
              #size=1.5,
              subset(dat.ex.df,Divergence=="May" & (name != "Itot.pct.2.5.") & (name != "Itot.pct.97.5.") &
                       (name != "I.pct.2.5.") & (name != "I.pct.97.5.")))+
    geom_line(aes(x=date,y=value,color=colorlab,linetype=linetype),
              #size=1.5,
              subset(dat.ex.df,Divergence=="December" & (name != "Itot.pct.2.5.") & (name != "Itot.pct.97.5.") &
                       (name != "I.pct.2.5.") & (name != "I.pct.97.5.")))+
    geom_line(aes(x=date,y=value,color=colorlab),
              size=1,
              subset(dat.ex.df,Divergence=="May" & (name != "Itot.pct.2.5.") & (name != "Itot.pct.97.5.") &
                       (name != "I.pct.2.5.") & (name != "I.pct.97.5.")))+
    geom_line(aes(x=date,y=value,color=colorlab),
              size=1,
              subset(dat.ex.df,Divergence=="December" & (name != "Itot.pct.2.5.") & (name != "Itot.pct.97.5.") &
                       (name != "I.pct.2.5.") & (name != "I.pct.97.5.")))+
    ggtitle(statenow)+ylim(0,NA)+
    scale_fill_viridis_d(alpha=0.6,begin=0.1,end=0.6,name="")+
    scale_color_viridis_d(begin=0,end=0.6,name="")+
    scale_x_date(expand=c(0,0))+
    ylab("")+
    theme_bw()+theme(legend.position = "right")+
    facet_wrap(~label,scales="free_y",ncol=1)
    pex <- pex + geom_line(aes(x=date,y=value,color=colorlab,linetype=linetype),
                               subset(dat.ex.df,name == "Itot.pct.2.5." | 
                                        name == "I.pct.2.5."),alpha=0.5)+
      geom_line(aes(x=date,y=value,color=colorlab,linetype=linetype),
                subset(dat.ex.df,name == "Itot.pct.97.5." |
                         name == "I.pct.97.5."),alpha=0.5)+
      scale_linetype_manual(name="",values=c(3,1,2))
  # print(pex)
  pex.list[[statenow]]<-pex
}
pexstates <- ggarrange(plotlist=pex.list,nrow = 1,common.legend = TRUE,legend="right")
print(pexstates)

ggsave(file.path(figfolder,"SupFig_Divergences.pdf"),pexstates,height=6,width=8,scale=1.7)
ggsave(file.path(figfolder,"SupFig_Divergences.jpeg"),pexstates,height=6,width=8,scale=1.7,dpi=600)
