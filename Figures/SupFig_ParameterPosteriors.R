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

for (mcmcpath in c("MCMC","MCMC.n0.5")) {
  allparms <- fread(file.path(mcmcpath,"sampled_parameters.csv"))
  allparms <- as.data.frame(allparms)
  
  ncoef.quants <- as.data.frame(t(apply(
    allparms[,startsWith(names(allparms),"ncoef.")],
    2,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))))
  names(ncoef.quants)<-make.names(names(ncoef.quants))
  ncoef.quants$id <- as.numeric(str_remove(rownames(ncoef.quants),"ncoef."))
  ncoef.quants$state[ncoef.quants$id!=0] <- statesvec[ncoef.quants$id]
  ncoef.quants$state[ncoef.quants$id==0] <- " F.E."
  pncoef<-ggplot(ncoef.quants,aes(y=state,xmin = X2.5., xlower = X25., 
                                  xmiddle = X50., xupper = X75., xmax = X97.5.))+
    geom_boxploth(stat="identity",aes(fill=(state==" F.E.")))+
    scale_y_discrete(limits=rev)+xlim(0,1)+
    xlab("n")+theme_bw()+theme(legend.position = "none") + 
    geom_vline(xintercept=ncoef.quants["ncoef.0","X50."],linetype=2)+
    geom_label(x=0.8,y=52,label=paste0("F.E.=",signif(ncoef.quants["ncoef.0","X50."],2),
                                       " [",signif(ncoef.quants["ncoef.0","X2.5."],2),
                                       "-",signif(ncoef.quants["ncoef.0","X97.5."],2),
                                       "]"),
               size=2.5,label.size=0)+
    scale_fill_viridis_d(begin=1,end=0.7,option = "magma")
  print(pncoef)
  
  Roffset.quants <- as.data.frame(t(apply(
    allparms[,startsWith(names(allparms),"Roffset.")],
    2,quantile,prob=c(0.025,0.25,0.5,0.75,0.975))))
  names(Roffset.quants)<-make.names(names(Roffset.quants))
  Roffset.quants$id <- as.numeric(str_remove(rownames(Roffset.quants),"Roffset."))
  Roffset.quants$state[Roffset.quants$id!=0] <- statesvec[Roffset.quants$id]
  Roffset.quants$state[Roffset.quants$id==0] <- " F.E."
  pRoffset<-ggplot(Roffset.quants,aes(y=state,xmin = X2.5., xlower = X25., 
                                      xmiddle = X50., xupper = X75., xmax = X97.5.))+
    geom_boxploth(stat="identity",aes(fill=(state==" F.E.")))+
    scale_y_discrete(limits=rev)+scale_x_log10()+
    coord_cartesian(xlim=c(0.01,30),expand = 0)+
    annotation_logticks(sides="b")+
    xlab("Offset (%)")+theme_bw()+theme(legend.position = "none") + 
    geom_vline(xintercept=Roffset.quants["Roffset.0","X50."],linetype=2)+
    geom_label(x=0.8,y=52,label=paste0("F.E.=",signif(Roffset.quants["Roffset.0","X50."],2),
                                       " [",signif(Roffset.quants["Roffset.0","X2.5."],2),
                                       "-",signif(Roffset.quants["Roffset.0","X97.5."],2),
                                       "]"),
               size=2.5,label.size=0)+
    scale_fill_viridis_d(begin=1,end=0.7,option = "magma")
  print(pRoffset)
  
  set.seed(123456)
  nsamp<-500
  sampparms <- allparms[sample(nrow(allparms),size=nsamp),]
  Tinf.prior.density <- data.frame(x=seq(0,28,0.01),y=dnorm(seq(0,28,0.01),m=14,sd=3.5))
  pTinf<-ggplot()+geom_histogram(data=sampparms,
                                 aes(Tinf,after_stat(density),
                                     color="Posterior"),bins=20,fill=NA)+
    geom_line(data=Tinf.prior.density,aes(x,y,color="Prior"))+theme_bw()+
    scale_y_continuous(minor_breaks=0.01)+  
    scale_color_viridis_d(option="magma",name="",begin=0,end=0.7)+
    annotate("label",x=5,y=0.1,label=paste0("Prior:\nm=14.0\nsd=3.5"),label.size=0,hjust=0)+
    annotate("label",x=20,y=0.1,label=paste0("Posterior:\nm=",signif(mean(allparms$Tinf),3),
                                             "\nsd=",signif(sd(allparms$Tinf),3)),label.size=0,hjust=0)
  print(pTinf)
  
  if (mcmcpath == "MCMC") {
    figname.pdf <- file.path(figfolder,"SupFig_RandomEffects.pdf")
    figname.jpg <- file.path(figfolder,"SupFig_RandomEffects.jpeg")
  } else {
    figname.pdf <- file.path(figfolder,
                             paste0("SupFig_RandomEffects.",mcmcpath,".pdf"))
    figname.jpg <- file.path(figfolder,
                             paste0("SupFig_RandomEffects.",mcmcpath,".jpeg"))
  }
  if (grepl(".n",mcmcpath)) {
    ggsave(figname.pdf,plot=pRoffset,height=8,width=3)
    ggsave(figname.jpg,plot=pRoffset,height=8,width=3,dpi=600)
  } else {
    ggsave(figname.pdf,plot=ggarrange(pncoef,pRoffset,ncol=2),height=8,width=6)
    ggsave(figname.jpg,plot=ggarrange(pncoef,pRoffset,ncol=2),height=8,width=6,dpi=600)
  }
  
  if (mcmcpath == "MCMC") {
    figname.pdf <- file.path(figfolder,"Tinf.pdf")
    figname.jpg <- file.path(figfolder,"Tinf.jpeg")
  } else {
    figname.pdf <- file.path(figfolder,
                             paste0("Tinf.",mcmcpath,".pdf"))
    figname.jpg <- file.path(figfolder,
                             paste0("Tinf.",mcmcpath,".jpeg"))
  }
  ggsave(figname.pdf,plot=pTinf,height=4,width=6)
  ggsave(figname.jpg,plot=pTinf,height=4,width=6,dpi=600)
}
