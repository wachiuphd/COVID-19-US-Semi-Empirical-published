library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggstance)
source("setup_data.R")

predquant <- fread(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
predquant<-as.data.frame(predquant)
predquant$date <- as.Date(predquant$date)


mcmcpath.n0.5 <- "MCMC.n0.5"
predquant.n0.5<-fread(file.path(mcmcpath.n0.5,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
predquant.n0.5<-as.data.frame(predquant.n0.5)
predquant.n0.5$date <- as.Date(predquant.n0.5$date)

pratio <- ggplot(predquant)+
  geom_ribbon(data=predquant,aes(x=date,ymin=IfracU.pct.2.5.,
                                 ymax=IfracU.pct.97.5.,fill="Posterior 95% CrI"),
              alpha=0.5)+
  geom_line(data=predquant,aes(x=date,y=IfracU.pct.50.,
                               color="Posterior median"),size=1)+
  geom_line(data=predquant.n0.5,aes(x=date,y=IfracU.pct.50.,
                               color="Posterior median (n=0.5)"))+
  scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
               date_labels = "%b")+  
  scale_fill_viridis_d(name="")+
  scale_color_viridis_d(name="Semi-empirical Model")+
  guides(color = guide_legend(order=1),
         fill = guide_legend(order=2))+
  theme_bw()+theme(legend.position="bottom")+
  ggtitle("Percent of prevalence undiagnosed")+
  coord_cartesian(ylim=c(0,100))+
  xlab("Date")+
  ylab("Percent of prevalence undiagnosed")+
  facet_wrap(~state,ncol=6)
print(pratio)
ggsave(file.path(figfolder,"FracPrevUndiagnosed.pdf"),pratio,height=9,width=7.5,scale=1.35)

### Extended SEIR predictions
predfile<-file.path(datafolder,"SEIR-ModelPredictions.csv")
seirpred<-fread(predfile)
seirpred$model <- "Extended-SEIR"
seirpred$date<-as.Date(seirpred$date)
### Imperial predictions
imperialpredfile<-file.path(datafolder,"ImperialPredictions.csv")
imperialpred<-fread(imperialpredfile)
imperialpred$model <- "Imperial"
imperialpred$date<-as.Date(imperialpred$date)

predquant <- subset(predquant,date<=max(imperialpred$date,seirpred$date))
predquant.n0.5 <- subset(predquant.n0.5,date<=max(imperialpred$date,seirpred$date))

## Version with SEIR & Imperial CrIs
pboth<-ggplot()+
  geom_ribbon(data=seirpred,
              aes(x=date,ymin=N_I_tot.2.5./1000,ymax=N_I_tot.97.5./1000,fill="Extended SEIR CrI"),
              alpha=0.25)+
  geom_line(data=seirpred,
            aes(x=date,y=N_I_tot.50./1000,linetype="Extended SEIR median"))+
  geom_ribbon(data=imperialpred,
              aes(x=date,ymin=N_I_tot.2.5./1000,ymax=N_I_tot.97.5./1000,fill="Imperial CrI"),
              alpha=0.25)+
  geom_line(data=imperialpred,
            aes(x=date,y=N_I_tot.50./1000,linetype="Imperial median"))+
  # geom_ribbon(data=predquant,aes(x=date,ymin=Itot.pct.2.5.,
  #                                ymax=Itot.pct.97.5.,fill="Posterior 95% CrI"),
  #             alpha=0.5)+
  geom_line(data=predquant,aes(x=date,y=Itot.pct.50.,
                               color="Posterior median"),size=1.5)+
  geom_line(data=predquant.n0.5,aes(x=date,y=Itot.pct.50.,
                               color="Posterior median (n=0.5)"),size=1)+
  scale_fill_viridis_d(begin=0.3,end=0.7,name="")+
  scale_linetype_manual(name="Epidemiologic Models",values=c(3,2,1))+
  scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
               date_labels = "%b")+  
  theme_bw()+theme(legend.position="bottom")+
  scale_color_viridis_d(begin=0.3,name="Semi-empirical Model")+
  guides(linetype = guide_legend(order=1,nrow=2),
         color = guide_legend(order=3,nrow=2),
         fill = guide_legend(order=2,nrow=3))+
  # theme(legend.title = element_text(size = 10),
  #                  legend.text = element_text(size = 10),
  #                  #legend.box="vertical",
  #                  axis.text.x = element_text(angle = 90))+
  ggtitle("State-wide Infection Prevalence")+
  xlab("Date")+
  ylab("Prevalence of Active Infection %")

pbothnolog <- pboth+facet_wrap(~state,scales="free_y",ncol=6)
print(pbothnolog)

ggsave(file.path(figfolder,"Figure3_v1_prevalence_posteriors_page.pdf"),pbothnolog,height=9,width=7.5,scale=1.35)
ggsave(file.path(figfolder,"Figure3_v1_prevalence_posteriors_page.jpeg"),pbothnolog,height=9,width=7.5,scale=1.35,dpi=600)
  
## Version with semi-empirical CrIs 
pboth<-ggplot()+
  # geom_ribbon(data=imperialpred,
  #             aes(x=date,ymin=N_I_tot.2.5./1000,ymax=N_I_tot.97.5./1000,fill="Imperial CrI"),
  #             alpha=0.5)+
  geom_line(data=imperialpred,
            aes(x=date,y=N_I_tot.50./1000,linetype="Imperial median"))+
  # geom_ribbon(data=seirpred,
  #             aes(x=date,ymin=N_I_tot.2.5./1000,ymax=N_I_tot.97.5./1000,fill="Extended SEIR CrI"),
  #             alpha=0.5)+
  geom_line(data=seirpred,
            aes(x=date,y=N_I_tot.50./1000,linetype="Extended SEIR median"))+
  geom_ribbon(data=predquant,aes(x=date,ymin=Itot.pct.2.5.,
                                 ymax=Itot.pct.97.5.,fill="Posterior 95% CrI"),
              alpha=0.5)+
  geom_line(data=predquant,aes(x=date,y=Itot.pct.50.,
                               color="Posterior median"),size=1)+
  geom_line(data=predquant.n0.5,aes(x=date,y=Itot.pct.50.,
                                    color="Posterior median (n=0.5)"))+
  scale_fill_viridis_d(begin=0.5,name="")+
  scale_linetype_manual(name="Epidemiologic Models",values=c(3,2,1))+
  scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
               date_labels = "%b")+  
  theme_bw()+theme(legend.position="bottom")+
  scale_color_viridis_d(begin=0.3,name="Semi-empirical Model")+
  guides(linetype = guide_legend(order=1,nrow=2),
         color = guide_legend(order=2,nrow=2),
         fill = guide_legend(order=3))+
  # theme(legend.title = element_text(size = 10),
  #                  legend.text = element_text(size = 10),
  #                  #legend.box="vertical",
  #                  axis.text.x = element_text(angle = 90))+
  ggtitle("State-wide Infection Prevalence")+
  xlab("Date")+
  ylab("Prevalence of Active Infection %")


# pbothnolog <- pboth+facet_wrap(~state,scales="free_y")
# print(pbothnolog)
# 
# pbothlog <- pboth+facet_wrap(~state)+scale_y_log10(breaks=10^(-3:1))+
#   coord_cartesian(ylim=c(0.001,10))
# print(pbothlog)
# ggsave(file.path(figfolder,"Figure3_prevalence_posteriors.pdf"),pbothnolog,height=8,width=10)
# ggsave(file.path(figfolder,"Figure3_prevalence_posteriors.jpeg"),pbothnolog,height=8,width=10,dpi=600)

pbothnolog <- pboth+facet_wrap(~state,scales="free_y",ncol=6)
print(pbothnolog)

ggsave(file.path(figfolder,"Figure3_v2_prevalence_posteriors_page.pdf"),pbothnolog,height=9,width=7.5,scale=1.35)
ggsave(file.path(figfolder,"Figure3_v2_prevalence_posteriors_page.jpeg"),pbothnolog,height=9,width=7.5,scale=1.35,dpi=600)

