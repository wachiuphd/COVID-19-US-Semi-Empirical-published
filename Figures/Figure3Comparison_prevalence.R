library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggstance)
source("setup_data.R")

predquant0<-fread(file.path(mcmcpath,"Posterior_I_SP_dI_Quantiles copy 3.csv"))
predquant0<-as.data.frame(predquant0)
predquant0$date <- as.Date(predquant0$date)

predquant <- fread(file.path(mcmcpath,"Posterior_I_SP_dI_Quantiles.csv"))
predquant<-as.data.frame(predquant)
predquant$date <- as.Date(predquant$date)


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

predquant0 <- subset(predquant0,date<=max(imperialpred$date,seirpred$date))
predquant <- subset(predquant,date<=max(imperialpred$date,seirpred$date))


pboth<-ggplot()+
  # geom_ribbon(data=predquant0,aes(x=date,ymin=I.pct.2.5.,
  #                                ymax=I.pct.97.5.,fill="Posterior 95% CrI (n fit)"),
  #             alpha=0.5)+
  geom_line(data=predquant0,aes(x=date,y=I.pct.50.,
                               color="Posterior median (n fit)"),size=2)+
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
  # geom_ribbon(data=predquant,aes(x=date,ymin=I.pct.2.5.,
  #                                ymax=I.pct.97.5.,fill="Posterior 95% CrI (n=0.5)"),
  #             alpha=0.5)+
  geom_line(data=predquant,aes(x=date,y=I.pct.50.,
                               color="Posterior median (n=0.5)"),size=1)+
  scale_fill_viridis_d(begin=0.5,name="")+
  scale_linetype_manual(name="Epidemiologic Models",values=c(3,2,1))+
  scale_x_date(date_minor_breaks="1 month",date_breaks = "2 months",
               date_labels = "%b")+  
  theme_bw()+theme(legend.position="bottom")+
  scale_color_viridis_d(begin=0.3,end=0.95,name="Semi-empirical Model",direction=-1)+
  guides(linetype = guide_legend(order=1,nrow=2),
         color = guide_legend(order=2,nrow=2),
         fill = guide_legend(order=3,nrow=2))+
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
ggsave(file.path(figfolder,"Figure3Comparison_prevalence_posteriors.pdf"),pbothnolog,height=8,width=10)