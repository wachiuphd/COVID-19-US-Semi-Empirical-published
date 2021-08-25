library(ggplot2)
library(ggpmisc)
library(ggpubr)
library(lme4)
library(ggforce)
library(ggstance)
library(ggrepel)
library(usmap)
library(lubridate)
library(scales)
source("setup_data.R")

allparms.all <- data.frame()
predquant.all <- data.frame()
mcmcpaths <- c("MCMC","MCMC.tau7","MCMC.tau28")
tauvec <- c(14,7,28)
for (i in 1:length(mcmcpaths)) {
  allparms <- fread(file.path(mcmcpaths[i],"sampled_parameters.csv"))
  allparms <- as.data.frame(allparms)
  allparms$mcmcpath <- mcmcpaths[i]
  allparms$tau <- tauvec[i]
  allparms.all <- rbind(allparms.all,allparms)
  predquant <- as.data.frame(fread(file.path(mcmcpaths[i],"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv")))
  predquant$date <- as.Date(predquant$date)
  predquant$tau <- tauvec[i]
  predquant.all <- rbind(predquant.all,predquant)
}
predquant.all$tauf <- factor(predquant.all$tau,levels=c(14,7,28))
ncoef <- allparms.all[,c(paste0("ncoef.",0:51),"tau","iter")]
ncoef.df <- pivot_longer(ncoef,1:52,names_to="state.num",values_to="value")
ncoef.df$state.num <- as.numeric(gsub("ncoef.","",ncoef.df$state.num))
ncoef.df$state <- ""
ncoef.df$state[ncoef.df$state.num>0] <-statesvec[ncoef.df$state.num[ncoef.df$state.num>0]]
ncoef.df$state[ncoef.df$state.num==0] <- "F.E."
ncoef.df$state <- factor(ncoef.df$state,levels=c("F.E.",statesvec))
pncoef<-ggplot(ncoef.df) + geom_boxplot(aes(x=value,y=factor(tau)),outlier.shape = NA) +
  facet_wrap(~state)+xlab("n")+ylab("tau (days)")+theme_bw()+
  scale_x_continuous(limits=c(0,1),breaks=c(0,0.5,1))

Roffset <- allparms.all[,c(paste0("Roffset.",0:51),"tau","iter")]
Roffset.df <- pivot_longer(Roffset,1:52,names_to="state.num",values_to="value")
Roffset.df$state.num <- as.numeric(gsub("Roffset.","",Roffset.df$state.num))
Roffset.df$state <- ""
Roffset.df$state[Roffset.df$state.num>0] <-statesvec[Roffset.df$state.num[Roffset.df$state.num>0]]
Roffset.df$state[Roffset.df$state.num==0] <- "F.E."
Roffset.df$state <- factor(Roffset.df$state,levels=c("F.E.",statesvec))
pSP0<-ggplot(Roffset.df) + geom_boxplot(aes(x=log10(value),y=factor(tau)),outlier.shape = NA) +
  facet_wrap(~state)+xlab("log10 SP0 (%)")+ylab("tau (days)")+theme_bw()+
  coord_cartesian(xlim=c(-2,2))

pboth <- ggarrange(pncoef,pSP0,ncol = 1,labels = "AUTO")

ggsave(file.path(figfolder,"SupFig_TauSensitivity_parameters.pdf"),pboth,height=9,width=7.5)
ggsave(file.path(figfolder,"SupFig_TauSensitivity_parameters.jpeg"),pboth,height=9,width=7.5,dpi=600)

pI <- ggplot(predquant.all) + 
  geom_line(aes(x=date,y=I.pct.50.,color=tauf,linetype=tauf)) + 
  scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
               date_labels = "%b")+  
  theme_bw()+theme(legend.position="bottom")+
  scale_color_viridis_d(end=0.7,name="tau")+
  scale_linetype_discrete(name="tau")+
  ggtitle("Sensitivity of Undiagnosed Infection Prevalence to Averaging Time")+
  xlab("Date")+
  ylab("Prevalence of Undiagnosed Infection %")+
  facet_wrap(~state,scales="free_y",ncol=6)
print(pI)

pItot <- ggplot(predquant.all) + 
  geom_line(aes(x=date,y=Itot.pct.50.,color=tauf,linetype=tauf)) + 
  scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
               date_labels = "%b")+  
  theme_bw()+theme(legend.position="bottom")+
  scale_color_viridis_d(end=0.7,name="tau")+
  scale_linetype_discrete(name="tau")+
  ggtitle("Sensitivity of Total Infection Prevalence to Averaging Time")+
  xlab("Date")+
  ylab("Prevalence of Total Infection %")+
  facet_wrap(~state,scales="free_y",ncol=6)
print(pItot)

pSP <- ggplot(predquant.all) + 
  geom_line(aes(x=date,y=SP.pct.50.,color=tauf,linetype=tauf)) + 
  scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
               date_labels = "%b")+  
  theme_bw()+theme(legend.position="bottom")+
  scale_color_viridis_d(end=0.7,name="tau")+
  scale_linetype_discrete(name="tau")+
  ggtitle("Sensitivity of Seroprevalence to Averaging Time")+
  xlab("Date")+
  ylab("Seroprevalence %")+
  facet_wrap(~state,scales="free_y",ncol=6)
print(pSP)

ggsave(file.path(figfolder,"SupFig_TauSensitivity_I.pdf"),pI,height=9,width=7.5,scale=1.35)
ggsave(file.path(figfolder,"SupFig_TauSensitivity_I.jpeg"),pI,height=9,width=7.5,scale=1.35,dpi=600)

ggsave(file.path(figfolder,"SupFig_TauSensitivity_Itot.pdf"),pItot,height=9,width=7.5,scale=1.35)
ggsave(file.path(figfolder,"SupFig_TauSensitivity_Itot.jpeg"),pItot,height=9,width=7.5,scale=1.35,dpi=600)

ggsave(file.path(figfolder,"SupFig_TauSensitivity_SP.pdf"),pSP,height=9,width=7.5,scale=1.35)
ggsave(file.path(figfolder,"SupFig_TauSensitivity_SP.jpeg"),pSP,height=9,width=7.5,scale=1.35,dpi=600)
