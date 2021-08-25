source("setup_data.R")

for (mcmcpath in c("MCMC","MCMC.n0.5")) {
  
  predquant <- fread(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
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
  
  seirscat.df <- merge.data.frame(seirpred[,c("state","date","N_I_tot.2.5.",
                                              "N_I_tot.50.","N_I_tot.97.5.","N_I_tot.Random")],
                                  predquant[,c("state","date","Itot.pct.50.")])
  pseirscat<-ggplot(seirscat.df)+
    geom_point(aes(x=Itot.pct.50.*1000,y=N_I_tot.50.,color="Extended SEIR"))+
    facet_wrap(~state)+scale_x_log10(limits=c(1,1e4))+scale_y_log10(limits=c(1,1e4))+geom_abline(intercept=1)+
    geom_abline(intercept=2^c(-1,1),linetype="dashed")+
    scale_color_viridis_d(begin=1,end=0.5,name="")+
    xlab("Semi-empirical (median)")+ylab("Extended SEIR (random)")+
    theme_bw()+theme(axis.text.x = element_text(angle = 90))
  print(pseirscat)
  
  imperialscat.df <- merge.data.frame(imperialpred[,c("state","date","N_I_tot.2.5.","N_I_tot.50.","N_I_tot.97.5.","N_I_tot.Random")],
                                      predquant[,c("state","date","Itot.pct.50.")])
  pimperialscat<-ggplot(imperialscat.df)+
    geom_point(aes(x=Itot.pct.50.*1000,y=N_I_tot.50.,color="Imperial"))+
    facet_wrap(~state)+scale_x_log10()+scale_y_log10()+geom_abline(intercept=1)+
    geom_abline(intercept=2^c(-1,1),linetype="dashed")+
    scale_color_viridis_d(begin=1,end=0.5,name="")+
    xlab("Semi-empirical (median)")+ylab("Imperial (random)")+
    theme_bw()+theme(axis.text.x = element_text(angle = 90))
  print(pimperialscat)
  
  pbothscat<-ggplot()+
    geom_errorbar(data=imperialscat.df,
                  aes(x=Itot.pct.50.*1000,ymin=N_I_tot.2.5.,ymax=N_I_tot.97.5.),alpha=0.5,color="grey")+
    geom_point(data=imperialscat.df,aes(x=Itot.pct.50.*1000,y=N_I_tot.50.,color="Imperial"),
               alpha=0.5,size=0.5)+
    geom_errorbar(data=seirscat.df,
                  aes(x=Itot.pct.50.*1000,ymin=N_I_tot.2.5.,ymax=N_I_tot.97.5.),alpha=0.5,color="grey")+
    geom_point(data=seirscat.df,aes(x=Itot.pct.50.*1000,y=N_I_tot.50.,color="Extended SEIR"),
               alpha=0.5,size=0.5)+
    facet_wrap(~state)+scale_x_log10(limits=c(1,1e4))+scale_y_log10(limits=c(1,1e4))+geom_abline(intercept=1)+
    #geom_abline(intercept=3^c(-1,1),linetype="dashed")+
    scale_color_viridis_d(begin=1,end=0.5,name="")+
    xlab("Semi-empirical (median)")+ylab("Epidemiologic model posterior")+
    theme_bw()+theme(axis.text.x = element_text(angle = 90))
  print(pbothscat)
  
  pbothscat_comb <- ggplot()+
    geom_point(data=imperialscat.df,aes(x=Itot.pct.50.*1000,y=N_I_tot.Random,color="Imperial"),
               alpha=0.5,size=0.5)+
    geom_point(data=seirscat.df,aes(x=Itot.pct.50.*1000,y=N_I_tot.Random,color="Extended SEIR"),
               alpha=0.5,size=0.5)+
    scale_x_log10(limits=c(1,1e4))+scale_y_log10(limits=c(1,1e4))+geom_abline(intercept=1)+
    #geom_abline(intercept=3^c(-1,1),linetype="dashed")+
    scale_color_viridis_d(begin=1,end=0.5,name="")+
    xlab("Semi-empirical (median)")+ylab("Epidemiologic model (random)")+
    theme_bw()
  
  seirres1<-lm(log10(N_I_tot.50.)~log10(Itot.pct.50.*1000),
               data=subset(seirscat.df,Itot.pct.50.>0))
  imperialres1<-lm(log10(N_I_tot.50.)~log10(Itot.pct.50.*1000),
                   data=subset(imperialscat.df,Itot.pct.50.>0 & N_I_tot.50.>0))
  seirres1sum<-summary(seirres1)
  imperialres1sum<-summary(imperialres1)
  
  seirres <- lm(log10(N_I_tot.50.)-log10(Itot.pct.50.*1000)~0,
                data=subset(seirscat.df,Itot.pct.50.>0))
  imperialres <- lm(log10(N_I_tot.50.)-log10(Itot.pct.50.*1000)~0,
                    data=subset(imperialscat.df,Itot.pct.50.>0 & N_I_tot.50.>0))
  seirressum<-summary(seirres)
  imperialressum<-summary(imperialres)
  
  pbothscat_comb_med <- ggplot()+
    geom_errorbar(data=imperialscat.df,
                  aes(x=Itot.pct.50.*1000,ymin=N_I_tot.2.5.,ymax=N_I_tot.97.5.),alpha=0.25,color="grey")+
    geom_errorbar(data=seirscat.df,
                  aes(x=Itot.pct.50.*1000,ymin=N_I_tot.2.5.,ymax=N_I_tot.97.5.),alpha=0.25,color="grey")+
    geom_point(data=imperialscat.df,aes(x=Itot.pct.50.*1000,y=N_I_tot.50.,color="Imperial"),
               alpha=0.25,size=0.5)+
    geom_point(data=seirscat.df,aes(x=Itot.pct.50.*1000,y=N_I_tot.50.,color="Extended SEIR"),
               alpha=0.25,size=0.5)+
    annotate("text",x=100,y=5000,size=3,label=paste0("Extended SEIR RSE=",signif(10^seirressum$sigma,3),"-fold",
                                                     "\nAdj R-sq=",signif(seirres1sum$adj.r.squared,2)))+
    annotate("text",x=1000,y=10,size=3,label=paste0("Imperial RSE=",signif(10^imperialressum$sigma,3),"-fold",
                                                    "\nAdj R-sq=",signif(imperialres1sum$adj.r.squared,2)))+
    scale_x_log10(limits=c(1,1e4))+scale_y_log10(limits=c(1,1e4))+geom_abline(intercept=1)+
    scale_color_viridis_d(begin=1,end=0.5,name="")+
    xlab("Semi-empirical (median)")+ylab("Epidemiologic model posterior")+
    theme_bw()
  
  print(pbothscat_comb_med+annotation_logticks())
  
  if (mcmcpath == "MCMC") {
    figname.pdf <- file.path(figfolder,"SupFig_Bayesian_prevalence_scatter_comb.pdf")
    figname.jpg <- file.path(figfolder,"SupFig_Bayesian_prevalence_scatter_comb.jpeg")
  } else {
    figname.pdf <- file.path(figfolder,
                             paste0("SupFig_Bayesian_prevalence_scatter_comb.",mcmcpath,".pdf"))
    figname.jpg <- file.path(figfolder,
                             paste0("SupFig_Bayesian_prevalence_scatter_comb.",mcmcpath,".jpeg"))
  }
  
  ggsave(figname.pdf,pbothscat_comb_med,height=4.5,width=6)
  ggsave(figname.jpg,pbothscat_comb_med,height=4.5,width=6,dpi=600)
  
  print(rse_seir_self<-sd(log10(seirscat.df$N_I_tot.Random/seirscat.df$N_I_tot.50.)))
  print(rse_imperial_self<-sd(log10(imperialscat.df$N_I_tot.Random/imperialscat.df$N_I_tot.50.),na.rm=T))
  
  seirscat.df$seir.50. <- seirscat.df$N_I_tot.50.
  imperialscat.df$imperial.50. <- imperialscat.df$N_I_tot.50.
  seirimperial <- merge.data.frame(seirscat.df[,-c(3:7)],imperialscat.df[,-c(3:7)])
  seirimperial <- subset(seirimperial,seir.50.>0 & imperial.50.>0)
  print(rse_seir_imperial <- sd(log10(seirimperial$seir.50.)-log10(seirimperial$imperial.50.),na.rm=T))
}