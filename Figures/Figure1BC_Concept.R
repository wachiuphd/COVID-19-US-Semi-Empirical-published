library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggstance)
library(ggforce)

datezero <- "2019-12-31"
maxdatadate <- "2020-12-31"
datafolder<-"Data"
mcmcpath <- "MCMC"
figfolder<-"Figures"
fips_table <- read.csv(file.path(datafolder,"FIPS_TABLE.csv"),colClasses=c(
  rep("character",4),rep("numeric",2)
))
rownames(fips_table)<-as.character(fips_table$Alpha.code)
statesvec <- fips_table$Alpha.code[2:52]
# datfile<-file.path(datafolder,"TestPositivityData.csv")
# alldat<-fread(datfile)
# alldat <- subset(alldat,date <= maxdatadate)


bdf <- data.frame(b=(1/c(1,10,100)),
                  TestingBias=factor(paste(c(1,10,100)),
                                     levels=paste(c(1,10,100)))) 

pconcept<-ggplot(bdf)+
  geom_abline(aes(intercept=b,slope=1,linetype=TestingBias))+
  scale_x_log10(limits=c(0.01,120),
                breaks=10^seq(-2,2),
                labels=10^seq(-2,2))+
  scale_y_log10(limits=c(0.6,12000)/1000,
                breaks=10^seq(-3,1),
                labels=paste(10^seq(-3,1)))+
  # geom_boxploth(aes(y=1/1000,PositivePct_tau,
  #                   color="U.S. State Positivity Rates"),
  #               data=subset(alldat,TestsPer1000_tau>0),width=0.1,
  #               outlier.size=0.05)+
  coord_cartesian(expand=FALSE)+
  annotation_logticks(short=unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm"))+
  scale_linetype_discrete(name="Positivity\nRate Bias (b)")+
  # scale_color_viridis_d(name="",begin=0.4)+
  xlab("Test Positivity %")+
  ylab(expression(atop(Prevalence~of~Undiagnosed,Infection~(I[U]/N)~{"%"})))+
  guides(color = guide_legend(order=1),
         fill = guide_legend(order=2))+
  theme_bw()+theme(legend.position = "bottom",
                   legend.title = element_text(size = 10),
                   legend.text = element_text(size = 10),
                   legend.box="vertical")
# print(pconcept)

x<-10^seq(-4,0,0.1)
bdf2 <- data.frame(TestingRate=x*1000,Bias=x^(-0.5),
                   BiasMin=(x)^(-0.1),BiasMax=x^(-0.9))
bdf2$BiasMax[bdf2$BiasMax>100]<-100
pconcept2 <- ggplot(bdf2)+
  geom_ribbon(aes(x=TestingRate,ymin=BiasMin,ymax=BiasMax,
                  fill="Positivity Rate Bias (b)\n(power parameter n=0.5 [0.1-0.9])"),alpha=0.3)+
  geom_line(aes(x=TestingRate,y=Bias))+#,linetype="dotted",alpha=0.3)+
  xlab("Daily Testing Rate (per 1000)")+
  ylab("Positivity Rate Bias (b)")+
  scale_y_log10(limits=c(0.6,100),expand=c(0,0),breaks=c(1,10,100))+
  annotation_logticks(side="l",short=unit(0.05, "cm"),
                      mid = unit(0.1, "cm"),
                      long = unit(0.15, "cm"))+
  # geom_boxploth(aes(y=1,TestsPer1000_tau,
  #                   color="U.S. State Testing Rates"),
  #               data=subset(alldat,TestsPer1000_tau>0),width=0.1,
  #               outlier.size=0.05)+
  # facet_zoom(xlim=c(0,50))+
  scale_fill_viridis_d(name="",end=0.7)+
  # scale_color_viridis_d(name="",begin=0.7)+
  # guides(color = guide_legend(order=1),
  #        fill = guide_legend(order=2))+
  theme_bw()+theme(legend.position = "bottom",
                   legend.title = element_text(size = 10),
                   legend.text = element_text(size = 10),
                   legend.box="vertical")
# print(pconcept2)

figconcept<-ggarrange(pconcept,pconcept2,ncol=2,labels=c("B","C"))
print(figconcept)
ggsave(file.path(figfolder,"Figure1BC_Concept.pdf"),figconcept,height=6.5,width=12,scale=0.6)
ggsave(file.path(figfolder,"Figure1BC_Concept.jpeg"),figconcept,height=6.5,width=12,scale=0.6,dpi=600)
