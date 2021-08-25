source("setup_data.R")

for (mcmcpath in c("MCMC","MCMC.n0.5")) {
  
  allparms <- fread(file.path(mcmcpath,"sampled_parameters.csv"))
  allparms <- as.data.frame(allparms)
  
  predquant <- fread(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
  predquant<-as.data.frame(predquant)
  predquant$date <- as.Date(predquant$date)
  fitdat<-merge.data.frame(cdcseroval,predquant)
  fitdat<- subset(fitdat,Data.50.>0)

  sigerr.med <- sd(log(fitdat$Data.50.)-log(fitdat$SP.pct.50.))
  
  res<-lm(log(Data.50.)~log(SP.pct.50.),data=fitdat)
  ressum<-summary(res)
  adjr2 <- signif(ressum$adj.r.squared,2)
  figtitle <- "Seroprevalence Validation Data vs. Median prediction"
  xtitle <- "Median posterior prediction"
  if (mcmcpath == "MCMC.n0.5") {
    figtitle <- paste(figtitle,"(n=0.5)")
    xtitle <- paste(xtitle,"(n=0.5)")
  }
  pscatval<-ggplot(data=fitdat,
                aes(x=SP.pct.50.,y=Data.50.,
                    label=state))+
    geom_errorbar(aes(ymin=Data.2.5.,ymax=Data.97.5.),alpha=0.7,color="grey")+
    geom_label(aes(fill=date),size=2)+
    scale_x_log10(limits=c(0.1,100))+scale_y_log10(limits=c(0.1,100))+
    annotation_logticks()+
    scale_fill_viridis_c(trans="date",begin=1,end=0.5)+
    xlab(xtitle)+ylab("Observed")+
    geom_abline(intercept=1)+
    geom_abline(intercept=exp(c(-1,1)*sigerr.med),linetype="dashed")+
    geom_abline(intercept=exp(c(-1,1)*(qnorm(0.975)*sigerr.med)),linetype="dotted")+
    geom_label(aes(x=50,y=30,label=paste0("Residual SE:\n",
                                          signif(exp(sigerr.med),3),
                                          "-fold")))+
    geom_label(aes(x=50,y=50,label=paste0("Residual 2.5-97.5% CrI\nrange: ",
                                          signif(exp(2*qnorm(0.975)*sigerr.med),3),
                                          "-fold")))+
    geom_label(aes(x=10,y=50,label=paste0("Adj R-sq=",adjr2)))+
    ggtitle(figtitle)+
    theme_bw()
  print(pscatval)
  
  if (mcmcpath == "MCMC") {
    figname.pdf <- file.path(figfolder,"SupFig_Bayesian_scatter_validation.pdf")
    figname.jpg <- file.path(figfolder,"SupFig_Bayesian_scatter_validation.jpeg")
  } else {
    figname.pdf <- file.path(figfolder,
                             paste0("SupFig_Bayesian_scatter_validation.",mcmcpath,".pdf"))
    figname.jpg <- file.path(figfolder,
                             paste0("SupFig_Bayesian_scatter_validation.",mcmcpath,".jpeg"))
  }
  
  ggsave(pscatval,filename=figname.pdf,height=7.5,width=8.5)
  ggsave(pscatval,filename=figname.jpg,height=7.5,width=8.5,dpi=600)
}
