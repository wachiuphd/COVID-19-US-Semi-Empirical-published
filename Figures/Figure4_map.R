library(data.table)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggstance)
library(usmap)
library(ggrepel)
source("setup_data.R")

currentdate <- as.Date("2020-12-31") # static
statexy <- read.csv(file.path(datafolder,"statexy.csv"),as.is=TRUE)
rownames(statexy) <- statexy$abbr

for (mcmcpath in c("MCMC","MCMC.n0.5")) {
  uspredquant <- fread(file.path(mcmcpath,"US_Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
  uspredquant <- as.data.frame(uspredquant)
  uspredquant$date <- as.Date(uspredquant$date)
  
  predquant <- fread(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
  predquant<-as.data.frame(predquant)
  predquant$date <- as.Date(predquant$date)
  
  uscurrentpred <- subset(uspredquant,date==currentdate)
  uspastpred <- subset(uspredquant,date==(currentdate-14))
  
  currentpred <- subset(predquant,date==currentdate)
  currentpred <- currentpred[order(currentpred$state),]
  
  pastpred <- subset(predquant,date==(currentdate-14))
  pastpred <- pastpred[order(pastpred$state),]
  
  currentpred$x <- statexy[currentpred$state,"x"]
  currentpred$y <- statexy[currentpred$state,"y"]
  pastpred$x <- statexy[pastpred$state,"x"]
  pastpred$y <- statexy[pastpred$state,"y"]
  
  ### Undiagnosed prevalence
  
  uscurrentpred$var <- uscurrentpred$I.pct.50.
  uscurrentpred$var.97.5. <- uscurrentpred$I.pct.97.5.
  uscurrentpred$var.2.5. <- uscurrentpred$I.pct.2.5.
  
  uspastpred$var <- uspastpred$I.pct.50.
  currentpred$var <- currentpred$I.pct.50.
  pastpred$var <- pastpred$I.pct.50.
  
  uscurrentpred$trend <- 100*(uscurrentpred$var/uspastpred$var-1)
  uscurrentpred$trend[uscurrentpred$trend>100]<-100
  uscurrentpred$dir<-ifelse(uscurrentpred$trend>0,"Increasing","Decreasing")
  uscurrentpred$dir <- factor(uscurrentpred$dir,levels=c("Increasing","Decreasing"))
  
  currentpred$trend <- 100*(currentpred$var/pastpred$var-1)
  currentpred$trend[currentpred$trend>100]<-100
  currentpred$dir<-ifelse(currentpred$trend>0,"Increasing","Decreasing")
  currentpred$dir <- factor(currentpred$dir,levels=c("Increasing","Decreasing"))
  
  uscurrentpred$state <- "US"
  
  mapstates<-statesvec
  
  prevmap<-plot_usmap(include=mapstates)+
    geom_point(data=currentpred,
               aes(x = x, y = y, 
                   size = var,
                   shape=dir,
                   color=trend,
                   stroke=2
               ))+
    geom_label_repel(data=currentpred,
                     aes(x=x,y=y,label=paste0(state,"\n",signif(var,2),"%")),
                     point.padding=NA,force=0.1,
                     hjust=0.5,vjust=0.5,size=3,label.size=0.1,alpha=0.8,seed=314)+
    geom_label_repel(data=currentpred,
                     aes(x=x,y=y,label=paste0(state,"\n",signif(var,2),"%")),
                     point.padding=NA,force=0.1,
                     hjust=0.5,vjust=0.5,size=3,label.size=0.1,fill=NA,seed=314)+
    # geom_text(data=currentpred,
    #           aes(x=x,y=y,label=paste0(state,"\n",signif(var,2),"%")),
    #           hjust=0.5,vjust=0.5,size=3)+
    scale_size(range=c(0.1,24),limits = c(0,ceiling(max(currentpred$var))),
               name="Semi-empirical\nPrevalence %")+
    scale_color_viridis_c(alpha=0.8,option="plasma",direction=-1,
                          limits=c(-100,100),
                          breaks=c(-100,-50,0,50,100),
                          labels=c("-100","-50","0","+50","+\u2265 100"),
                          name="2-wk Change (%)")+
    scale_shape_manual(values=c(16,1),name="2-wk Trend")+
    guides(size = guide_legend(order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 12))
           #color = guide_legend(order = 2),
    )+
    ggtitle(paste0("Undiagnosed Infection Prevalence as of ",currentdate,
                   "\nUS overall: ",round(uscurrentpred$var,2),
                   "% [95% CI: ",round(uscurrentpred$var.2.5.,2),
                   "%-",round(uscurrentpred$var.97.5.,2),
                   "%] (2-wk trend: ",uscurrentpred$dir,")"))+
    theme(legend.position = "right",
          legend.box="vertical")
  
  print(prevmap)
  
  ### Total prevalence
  
  uscurrentpred$vartot <- uscurrentpred$Itot.pct.50.
  uscurrentpred$vartot.97.5. <- uscurrentpred$Itot.pct.97.5.
  uscurrentpred$vartot.2.5. <- uscurrentpred$Itot.pct.2.5.
  
  uspastpred$vartot <- uspastpred$Itot.pct.50.
  currentpred$vartot <- currentpred$Itot.pct.50.
  pastpred$vartot <- pastpred$Itot.pct.50.
  
  uscurrentpred$trendtot <- 100*(uscurrentpred$vartot/uspastpred$vartot-1)
  uscurrentpred$trendtot[uscurrentpred$trendtot>100]<-100
  uscurrentpred$dirtot<-ifelse(uscurrentpred$trendtot>0,"Increasing","Decreasing")
  uscurrentpred$dirtot <- factor(uscurrentpred$dirtot,levels=c("Increasing","Decreasing"))
  
  currentpred$trendtot <- 100*(currentpred$vartot/pastpred$vartot-1)
  currentpred$trendtot[currentpred$trendtot>100]<-100
  currentpred$dirtot<-ifelse(currentpred$trendtot>0,"Increasing","Decreasing")
  currentpred$dirtot <- factor(currentpred$dirtot,levels=c("Increasing","Decreasing"))
  
  uscurrentpred$state <- "US"
  
  mapstates<-statesvec
  
  prevmaptot<-plot_usmap(include=mapstates)+
    geom_point(data=currentpred,
               aes(x = x, y = y, 
                   size = vartot,
                   shape=dir,
                   color=trend,
                   stroke=2
               ))+
    geom_label_repel(data=currentpred,
                     aes(x=x,y=y,label=paste0(state,"\n",signif(vartot,2),"%")),
                     point.padding=NA,force=0.1,
                     hjust=0.5,vjust=0.5,size=3,label.size=0.1,alpha=0.8,seed=314)+
    geom_label_repel(data=currentpred,
                     aes(x=x,y=y,label=paste0(state,"\n",signif(vartot,2),"%")),
                     point.padding=NA,force=0.1,
                     hjust=0.5,vjust=0.5,size=3,label.size=0.1,fill=NA,seed=314)+
    # geom_text(data=currentpred,
    #           aes(x=x,y=y,label=paste0(state,"\n",signif(vartot,2),"%")),
    #           hjust=0.5,vjust=0.5,size=3)+
    scale_size(range=c(0.1,24),limits = c(0,ceiling(max(currentpred$vartot))),
               name="Semi-empirical\nPrevalence %")+
    scale_color_viridis_c(alpha=0.8,option="plasma",direction=-1,
                          limits=c(-100,100),
                          breaks=c(-100,-50,0,50,100),
                          labels=c("-100","-50","0","+50","+\u2265 100"),
                          name="2-wk Change (%)")+
    scale_shape_manual(values=c(16,1),name="2-wk Trend")+
    guides(size = guide_legend(order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 12))
           #color = guide_legend(order = 2),
    )+
    ggtitle(paste0("Total Infection Prevalence as of ",currentdate,
                   "\nUS overall: ",round(uscurrentpred$vartot,2),
                   "% [95% CI: ",round(uscurrentpred$vartot.2.5.,2),
                   "%-",round(uscurrentpred$vartot.97.5.,2),
                   "%] (2-wk trend: ",uscurrentpred$dir,")"))+
    theme(legend.position = "right",
          legend.box="vertical")
  
  print(prevmaptot)
  
  ### Seroprevalence
  
  seromap<-plot_usmap(include=mapstates)+
    geom_point(data=currentpred,
               aes(x = x, y = y, 
                   size = SP.pct.50.,
                   stroke=1,color=SP.pct.50.
               ),shape=16)+
    geom_label_repel(data=currentpred,
                     aes(x=x,y=y,label=paste0(state,"\n",signif(SP.pct.50.,2),"%")),
                     point.padding=NA,force=0.1,
                     hjust=0.5,vjust=0.5,size=3,label.size=0,alpha=0.2,seed=314)+
    geom_label_repel(data=currentpred,
                     aes(x=x,y=y,label=paste0(state,"\n",signif(SP.pct.50.,2),"%")),
                     point.padding=NA,force=0.1,
                     hjust=0.5,vjust=0.5,size=3,label.size=0,fill=NA,seed=314)+
    scale_size(range=c(0.1,24),limits = c(0,ceiling(max(currentpred$SP.pct.50.))),
               name="Semi-empirical\nSeroprevalence %")+
    scale_color_viridis_c(alpha=0.8,option="plasma",begin=1,end=0,name="",limits=c(0,100))+
    guides(size = guide_legend(order = 1),
           shape = guide_legend(order = 2,
                                override.aes = list(size = 12))
           #color = guide_legend(order = 2),
    )+
    ggtitle(paste0("Seroprevalence as of ",currentdate,
                   "\nUS overall: ",round(uscurrentpred$SP.pct.50.,1),
                   "% [95% CI: ",round(uscurrentpred$SP.pct.2.5.,1),
                   "%-",round(uscurrentpred$SP.pct.97.5.,1),
                   "%]"))+
    theme(legend.position = "right",
          legend.box="vertical")
  
  print(seromap)
  pmaps <- ggarrange(prevmap,prevmaptot,seromap,ncol=1,nrow=3,labels="AUTO")
  print(pmaps)
  
  allcurrentnames<-c("date","state","I.pct.2.5.","I.pct.50.","I.pct.97.5.",
                     "Itot.pct.2.5.","Itot.pct.50.","Itot.pct.97.5.",
                     "SP.pct.2.5.","SP.pct.50.","SP.pct.97.5.","trend","dir",
                     "trendtot","dirtot")
  allcurrentpred <- rbind(
    currentpred[,allcurrentnames],
    uscurrentpred[,allcurrentnames]
  )
  allcurrentpred<-rename(allcurrentpred,trend.2wk.pct.change=trend)
  allcurrentpred<-rename(allcurrentpred,trend.2wk.direction=dir)
  allcurrentpred<-rename(allcurrentpred,trendtot.2wk.pct.change=trendtot)
  allcurrentpred<-rename(allcurrentpred,trendtot.2wk.direction=dirtot)
  allcurrentpred$trend.2wk.pct.change<-allcurrentpred$I.pct.50.*
    (allcurrentpred$trend.2wk.pct.change/(100+allcurrentpred$trend.2wk.pct.change))
  allcurrentpred$trendtot.2wk.pct.change<-allcurrentpred$Itot.pct.50.*
    (allcurrentpred$trendtot.2wk.pct.change/(100+allcurrentpred$trendtot.2wk.pct.change))
  
  fwrite(allcurrentpred,file.path(mcmcpath,"CurrentEstimates.csv"))
  if (mcmcpath == "MCMC") {
    figname.pdf <- file.path(figfolder,"Figure4_Maps.pdf")
    figname.jpg <- file.path(figfolder,"Figure4_Maps.jpeg")
  } else {
    figname.pdf <- file.path(figfolder,
                             paste0("Figure4_Maps.",mcmcpath,".pdf"))
    figname.jpg <- file.path(figfolder,
                             paste0("Figure4_Maps.",mcmcpath,".jpeg"))
  }
  ggsave(figname.pdf,pmaps,height=12,width=8,scale=1.5)
  ggsave(figname.jpg,pmaps,height=12,width=8,dpi=600,scale=1.5)
}

