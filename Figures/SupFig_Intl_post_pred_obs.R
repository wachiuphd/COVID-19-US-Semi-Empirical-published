library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggstance)
source(file.path("International","setup_intl_data.R"))

p<-ggplot() + geom_point(data=intl_serodat,aes(x=date,y=Data.50.))+
  geom_errorbar(data=intl_serodat,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.))+
  facet_wrap(~location)+
  ylab("Seroprevalence %")+
  xlab("Date")+
  theme_bw()+
  theme(legend.position="right")+
  guides(shape = guide_legend(order=1),
         color = guide_legend(order=2),
         fill = guide_legend(order=3))+
  scale_x_date(date_minor_breaks="1 month")+
  scale_y_log10()+annotation_logticks(side="l")+coord_cartesian(ylim=c(0.01,100))+
  facet_wrap(~location)+
  ggtitle("International Prevalence Data")
print(p)


infile <- c("Intl_FixedEffects_NoOffset_Posterior_I_SP_dI_Quantiles.csv",
            "Intl_FixedEffects_Posterior_I_SP_dI_Quantiles.csv",
            "Intl_Posterior_I_SP_dI_Quantiles.csv")
pdffile <- c("Intl_post_FE_NoOffset_nolog.pdf","Intl_post_FE_nolog.pdf","SupFig_Intl_post_RE_nolog.pdf")
jpegfile <- c("Intl_post_FE_NoOffset_nolog.jpeg","Intl_post_FE_nolog.jpeg","SupFig_Intl_post_RE_nolog.jpeg")
scatfile <- c("Intl_post_FE_NoOffset_scat.pdf","Intl_post_FE_scat.pdf","Intl_post_RE_scat.pdf")
scatjpegfile <- c("Intl_post_FE_NoOffset_scat.jpeg","Intl_post_FE_scat.jpeg","Intl_post_RE_scat.jpeg")
posttype <- c("(Fixed Effects, no offset)","Fixed Effects","Random Effects")

predquant.n0.5<-fread(file.path(intlpath,"Intl_n0.5_Posterior_I_SP_dI_Quantiles.csv"))
predquant.n0.5<-as.data.frame(predquant.n0.5)
predquant.n0.5$date <- as.Date(predquant.n0.5$date)
fitdat<-merge.data.frame(intl_serodat,predquant.n0.5)

res<-lm(log(Data.50.)~log(SP.pct.50.),data=fitdat)
ressum<-summary(res)
adjr2 <- signif(ressum$adj.r.squared,2)
pscat<-ggplot(data=fitdat,
              aes(x=SP.pct.50.,y=Data.50.,
                  label=location))+
  geom_errorbar(aes(ymin=Data.2.5.,ymax=Data.97.5.),alpha=0.7,color="grey")+
  geom_label(aes(fill=date),size=2)+
  scale_x_log10(limits=c(0.1,100))+scale_y_log10(limits=c(0.1,100))+
  annotation_logticks()+
  scale_fill_viridis_c(trans="date",begin=1,end=0.5)+
  xlab("Median posterior prediction (n=0.5)")+ylab("Observed")+
  geom_abline()+
  geom_label(aes(x=10,y=50,label=paste0("Adj R-sq=",adjr2)))+
  ggtitle(paste0("International Prevalence Median Predictions (n=0.5) vs. Data"))+
  theme_bw()
print(pscat)

ggsave(file.path(figfolder,"Intl_post_n0.5_scat.pdf"),pscat,height=7.5,width=8.5)
ggsave(file.path(figfolder,"Intl_post_n0.5_scat.jpeg"),pscat,height=7.5,width=8.5,dpi=600)

for (j in 1:3) {
  predquant<-fread(file.path(intlpath,infile[j]))
  predquant<-as.data.frame(predquant)
  predquant$date <- as.Date(predquant$date)
  
  psero<-ggplot() + 
    geom_ribbon(data=predquant,
                aes(x=date,ymin=SP.pct.2.5.,ymax=SP.pct.97.5.,
                    fill="Posterior 95% CrI"),alpha=0.5)+
    geom_point(data=intl_serodat,aes(x=date,y=Data.50.))+
    geom_errorbar(data=intl_serodat,aes(x=date,ymin=Data.2.5.,ymax=Data.97.5.))+
    #  scale_shape_discrete(name="Seroprevalence\nData Source")+
    geom_line(data=predquant,
              aes(x=date,y=SP.pct.50.,color="Posterior median"),size=1)+
    geom_line(data=predquant.n0.5,
              aes(x=date,y=SP.pct.50.,color="Posterior median (n=0.5)"))+
    facet_wrap(~location,scales="free_y")+
    ylab("Seroprevalence %")+
    xlab("Date")+
    theme_bw()+
    theme(legend.position="right")+
    scale_color_viridis_d(begin=0.3,name=
                            paste0("Seroprevalence Posteriors\n",posttype[j]))+
    scale_fill_viridis_d(begin=0.5,name="")+
    guides(shape = guide_legend(order=1),
           color = guide_legend(order=2),
           fill = guide_legend(order=3))+
    scale_x_date(date_minor_breaks="1 month")+
    #scale_y_log10()+annotation_logticks(side="l")+coord_cartesian(ylim=c(0.01,100))+
    facet_wrap(~location)+
    ggtitle(paste0("International Prevalence Predictions ",posttype[j]," vs. Data"))
  print(psero)
  
  ggsave(file.path(figfolder,pdffile[j]),psero,height=8,width=10)
  ggsave(file.path(figfolder,jpegfile[j]),psero,height=8,width=10,dpi=600)
  
  fitdat<-merge.data.frame(intl_serodat,predquant)
  
  res<-lm(log(Data.50.)~log(SP.pct.50.),data=fitdat)
  ressum<-summary(res)
  adjr2 <- signif(ressum$adj.r.squared,2)
  pscat<-ggplot(data=fitdat,
                aes(x=SP.pct.50.,y=Data.50.,
                    label=location))+
    geom_errorbar(aes(ymin=Data.2.5.,ymax=Data.97.5.),alpha=0.7,color="grey")+
    geom_label(aes(fill=date),size=2)+
    scale_x_log10(limits=c(0.1,100))+scale_y_log10(limits=c(0.1,100))+
    annotation_logticks()+
    scale_fill_viridis_c(trans="date",begin=1,end=0.5)+
    xlab("Median posterior prediction")+ylab("Observed")+
    geom_abline()+
    geom_label(aes(x=10,y=50,label=paste0("Adj R-sq=",adjr2)))+
    ggtitle(paste0("International Prevalence Median Predictions ",posttype[j]," vs. Data"))+
    theme_bw()
  print(pscat)
  
  ggsave(file.path(figfolder,scatfile[j]),pscat,height=7.5,width=8.5)
  ggsave(file.path(figfolder,scatjpegfile[j]),pscat,height=7.5,width=8.5,dpi=600)
}


