library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggstance)
source("setup_data.R")

predquant0<-fread(file.path(mcmcpath,"Posterior_I_SP_dI_Quantiles copy 3.csv"))
predquant0<-as.data.frame(predquant0)
predquant0$date <- as.Date(predquant0$date)
predquant0 <- subset(predquant0,date < today)

predquant<-fread(file.path(mcmcpath,"Posterior_I_SP_dI_Quantiles.csv"))
predquant<-as.data.frame(predquant)
predquant$date <- as.Date(predquant$date)
predquant <- subset(predquant,date < today)

psero<-ggplot() + 
  geom_line(data=predquant0,
            aes(x=date,y=SP.pct.50.,color="Posterior median (n fit)"),
            size=2)+
  geom_line(data=predquant,
            aes(x=date,y=SP.pct.50.,color="Posterior median (n=0.5) "),
            size=1)+
  # geom_ribbon(data=predquant,
  #             aes(x=date,ymin=SP.pct.2.5.,ymax=SP.pct.97.5.,
  #                 fill="Posterior 95% CrI"),alpha=0.5)+
  geom_point(data=serodat,aes(x=date,y=Data.50.,shape=Source))+
  geom_errorbar(data=serodat,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.))+
  facet_wrap(~state,scales="free_y",ncol=6)+
  xlab("Date")+
  ylab("Seroprevalence %")+
  scale_x_date(date_minor_breaks="1 month")+
  #scale_y_log10()+annotation_logticks(side="l")+coord_cartesian(ylim=c(0.1,100))+
  ggtitle("State-wide Seroprevalence")+
  #  scale_shape_discrete(name="Seroprevalence\nData Source")+
  scale_shape_manual(name="Seroprevalence\nData Source",
                     values=c(16,17,15,3,7,1))+
  scale_color_viridis_d(begin=0.3,end=0.95,name="Semi-empirical Model",direction=-1)+
  scale_fill_viridis_d(begin=0.5,name="")+
  guides(shape = guide_legend(order=1,nrow=3),
         color = guide_legend(order=2,nrow=2),
         fill = guide_legend(order=3))+
  theme_bw()+
  theme(legend.position="bottom")
print(psero)
pseroval <- psero+geom_point(data=cdcseroval,aes(x=date,y=Data.50.,shape=Source))+
  geom_errorbar(data=cdcseroval,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.))
print(pseroval)
ggsave(file.path(figfolder,"Figure2Comparison_Seroprevalence_posteriors_validation_page.pdf"),pseroval,height=9,width=7.5,scale=1.35)
