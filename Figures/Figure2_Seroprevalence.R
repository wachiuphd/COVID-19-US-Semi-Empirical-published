library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggstance)
source("setup_data.R")
predquant<-fread(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
predquant<-as.data.frame(predquant)
predquant$date <- as.Date(predquant$date)
predquant <- subset(predquant,date < today)

mcmcpath.n0.5 <- "MCMC.n0.5"
predquant.n0.5<-fread(file.path(mcmcpath.n0.5,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
predquant.n0.5<-as.data.frame(predquant.n0.5)
predquant.n0.5$date <- as.Date(predquant.n0.5$date)
predquant.n0.5 <- subset(predquant.n0.5,date < today)

# psero<-ggplot() + geom_point(data=serodat,aes(x=date,y=Data.50.,shape=Source))+
#   geom_errorbar(data=serodat,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.))+
#   geom_line(data=predquant,
#             aes(x=date,y=SP.pct.50.,color="Posterior median"),size=1)+
#   geom_ribbon(data=predquant,
#               aes(x=date,ymin=SP.pct.2.5.,ymax=SP.pct.97.5.,
#                   fill="Posterior 95% CrI"),alpha=0.5)+
#   geom_line(data=predquant.n0.5,
#             aes(x=date,y=SP.pct.50.,color="Posterior median (n=0.5)"))+
#   facet_wrap(~state,scales="free_y")+
#   xlab("Date")+
#   ylab("Seroprevalence %")+
#   scale_x_date(date_minor_breaks="1 month")+
#   #scale_y_log10()+annotation_logticks(side="l")+coord_cartesian(ylim=c(0.1,100))+
#   ggtitle("State-wide Seroprevalence")+
#   scale_shape_manual(name="Seroprevalence\nData Source",
#                      values=c(18,17,15,16,3,7,1))+
#   scale_color_viridis_d(begin=0.3,name="Semi-empirical Model")+
#   scale_fill_viridis_d(begin=0.5,name="")+
#   guides(shape = guide_legend(order=1,nrow=4),
#          color = guide_legend(order=2,nrow=2),
#          fill = guide_legend(order=3))+
#   theme_bw()+
#   theme(legend.position="bottom")
# print(psero)
# pseroval <- psero+geom_point(data=cdcseroval,aes(x=date,y=Data.50.,shape=Source))+
#   geom_errorbar(data=cdcseroval,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.))
# print(pseroval)
# ggsave(file.path(figfolder,"Figure2_Seroprevalence_posteriors_validation.pdf"),pseroval,height=8,width=10)
# ggsave(file.path(figfolder,"Figure2_Seroprevalence_posteriors_validation.jpeg"),pseroval,height=8,width=10,dpi=600)

psero<-ggplot() + geom_point(data=serodat,aes(x=date,y=Data.50.,shape=Source))+
  geom_errorbar(data=serodat,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.))+
  geom_line(data=predquant,
            aes(x=date,y=SP.pct.50.,color="Posterior median"),size=1)+
  geom_ribbon(data=predquant,
              aes(x=date,ymin=SP.pct.2.5.,ymax=SP.pct.97.5.,
                  fill="Posterior 95% CrI"),alpha=0.5)+
  geom_line(data=predquant.n0.5,
            aes(x=date,y=SP.pct.50.,color="Posterior median (n=0.5)"))+
  facet_wrap(~state,scales="free_y",ncol=6)+
  xlab("Date")+
  ylab("Seroprevalence %")+
  scale_x_date(date_minor_breaks="1 month")+
  #scale_y_log10()+annotation_logticks(side="l")+coord_cartesian(ylim=c(0.1,100))+
  ggtitle("State-wide Seroprevalence")+
  scale_shape_manual(name="Seroprevalence\nData Source",
                     values=c(18,17,15,16,3,7,1))+
  scale_color_viridis_d(begin=0.3,name="Semi-empirical Model")+
  scale_fill_viridis_d(begin=0.5,name="")+
  guides(shape = guide_legend(order=1,nrow=4),
         color = guide_legend(order=2,nrow=2),
         fill = guide_legend(order=3))+
  theme_bw()+
  theme(legend.position="bottom")
#print(psero)
pseroval <- psero+geom_point(data=cdcseroval,aes(x=date,y=Data.50.,shape=Source))+
  geom_errorbar(data=cdcseroval,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.))
print(pseroval)
ggsave(file.path(figfolder,"Figure2_Seroprevalence_posteriors_validation_page.pdf"),pseroval,height=9,width=7.5,scale=1.35)
ggsave(file.path(figfolder,"Figure2_Seroprevalence_posteriors_validation_page.jpeg"),pseroval,height=9,width=7.5,dpi=600,scale=1.35)
