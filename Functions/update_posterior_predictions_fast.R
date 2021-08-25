library(googledrive)
# Google drive path
googledrivepath<-"~/COVID-19-US-Semi-Empirical"
# Update each state
sampparms<-fread(file.path(mcmcpath,"sampled_parameters.csv"))
sampparms<-as.data.frame(sampparms)
nsamp<-2000
ncoefnames<-paste0("ncoef.",1:n.id)
Roffsetnames<-paste0("Roffset.",1:n.id)
ndates<-nrow(Cases_tau.pct.mat)
I.pct.arr<-array(0,dim=c(ndates,n.id,nsamp),
                 dimnames=list(date=as.character(Cases_tau.pct.df$date),
                               id=1:n.id,
                               j.samp=1:nsamp))
SP.pct.arr<-I.pct.arr
dI.pct.arr<-I.pct.arr
ID.pct.arr<-I.pct.arr
Itot.pct.arr<-I.pct.arr
IfracU.pct.arr <- I.pct.arr
pop.weights.arr<-I.pct.arr
pop.weights.mat<-matrix(rep(fips_table[statesvec,"pop"],ndates),
                        ncol=n.id,nrow=ndates,byrow=TRUE)/sum(fips_table[statesvec,"pop"])
  
for (j.samp in 1:nsamp) {
  if ((j.samp %% 100)==1) print(j.samp)
  Tinf.samp<-sampparms[j.samp,"Tinf"]
  ncoef.vec.samp<-as.numeric(sampparms[j.samp,ncoefnames])
  Roffset.vec.samp<-as.numeric(sampparms[j.samp,Roffsetnames])
  preddat.tmp <- get_all_pred_fast(ncoef.vec.samp,
                              Roffset.vec.samp,
                              Tinf.samp,
                              addCumulPos = addCumulPos, Trec.vec=Trec.global)
  pop.weights.arr[,,j.samp]<-pop.weights.mat
  I.pct.arr[,,j.samp]<-preddat.tmp[[1]]
  SP.pct.arr[,,j.samp]<-preddat.tmp[[2]]
  dI.pct.arr[,,j.samp]<-preddat.tmp[[3]]
  ID.pct.arr[,,j.samp]<-preddat.tmp[[4]]
  Itot.pct.arr[,,j.samp]<-preddat.tmp[[5]]
  IfracU.pct.arr[,,j.samp]<- 100*preddat.tmp[[1]]/preddat.tmp[[5]]
}
predquant.I.pct.arr<-apply(I.pct.arr,1:2,quantile,
                           prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
predquant.I.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  predquant.I.pct.arr,stringsAsFactors = FALSE),id_cols = 2:3,names_from=1,values_from=4))
tmpnames<-names(predquant.I.pct)
names(predquant.I.pct)<-sub("X","I.pct.",make.names(tmpnames))
predquant.I.pct$state<-statesvec[as.numeric(predquant.I.pct$id)]

predquant.SP.pct.arr<-apply(SP.pct.arr,1:2,quantile,
                            prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
predquant.SP.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  predquant.SP.pct.arr,stringsAsFactors = FALSE),id_cols = 2:3,names_from=1,values_from=4))
tmpnames<-names(predquant.SP.pct)
names(predquant.SP.pct)<-sub("X","SP.pct.",make.names(tmpnames))
predquant.SP.pct$state<-statesvec[as.numeric(predquant.SP.pct$id)]

predquant.dI.pct.arr<-apply(dI.pct.arr,1:2,quantile,
                            prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
predquant.dI.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  predquant.dI.pct.arr,stringsAsFactors = FALSE),id_cols = 2:3,names_from=1,values_from=4))
tmpnames<-names(predquant.dI.pct)
names(predquant.dI.pct)<-sub("X","dI.pct.",make.names(tmpnames))
predquant.dI.pct$state<-statesvec[as.numeric(predquant.dI.pct$id)]

predquant.ID.pct.arr<-apply(ID.pct.arr,1:2,quantile,
                            prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
predquant.ID.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  predquant.ID.pct.arr,stringsAsFactors = FALSE),id_cols = 2:3,names_from=1,values_from=4))
tmpnames<-names(predquant.ID.pct)
names(predquant.ID.pct)<-sub("X","ID.pct.",make.names(tmpnames))
predquant.ID.pct$state<-statesvec[as.numeric(predquant.ID.pct$id)]

predquant.Itot.pct.arr<-apply(Itot.pct.arr,1:2,quantile,
                            prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
predquant.Itot.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  predquant.Itot.pct.arr,stringsAsFactors = FALSE),id_cols = 2:3,names_from=1,values_from=4))
tmpnames<-names(predquant.Itot.pct)
names(predquant.Itot.pct)<-sub("X","Itot.pct.",make.names(tmpnames))
predquant.Itot.pct$state<-statesvec[as.numeric(predquant.Itot.pct$id)]

predquant.IfracU.pct.arr<-apply(IfracU.pct.arr,1:2,quantile,
                              prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
predquant.IfracU.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  predquant.IfracU.pct.arr,stringsAsFactors = FALSE),id_cols = 2:3,names_from=1,values_from=4))
tmpnames<-names(predquant.IfracU.pct)
names(predquant.IfracU.pct)<-sub("X","IfracU.pct.",make.names(tmpnames))
predquant.IfracU.pct$state<-statesvec[as.numeric(predquant.IfracU.pct$id)]

predquant<-left_join(left_join(left_join(left_join(
  inner_join(predquant.I.pct,predquant.SP.pct),predquant.dI.pct),
  predquant.ID.pct),predquant.Itot.pct),predquant.IfracU.pct)

predquant <- predquant[,-grep("id",names(predquant))]
fwrite(predquant,file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
fwrite(predquant,file.path(apppath,mcmcpath,"Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
# drive_put(file.path(mcmcpath,"Posterior_I_SP_dI_ID_Itot_Quantiles.csv"),
#           path=file.path(googledrivepath,mcmcpath,"Posterior_I_SP_dI_ID_Itot_Quantiles.csv"))

# Update US totals

uspred.I.pct.arr<-apply(I.pct.arr*pop.weights.arr,c(1,3),sum,na.rm=T)
uspred.SP.pct.arr<-apply(SP.pct.arr*pop.weights.arr,c(1,3),sum,na.rm=T)
uspred.dI.pct.arr<-apply(dI.pct.arr*pop.weights.arr,c(1,3),sum,na.rm=T)
uspred.ID.pct.arr<-apply(ID.pct.arr*pop.weights.arr,c(1,3),sum,na.rm=T)
uspred.Itot.pct.arr<-apply(Itot.pct.arr*pop.weights.arr,c(1,3),sum,na.rm=T)
uspred.IfracU.pct.arr<-apply(IfracU.pct.arr*pop.weights.arr,c(1,3),sum,na.rm=T)

uspredquant.I.pct.arr<-apply(uspred.I.pct.arr,1,quantile,
                           prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
uspredquant.I.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  uspredquant.I.pct.arr,stringsAsFactors = FALSE),id_cols = 2,names_from=1,values_from=3))
tmpnames<-names(uspredquant.I.pct)
names(uspredquant.I.pct)<-sub("X","I.pct.",make.names(tmpnames))
uspredquant.I.pct$state<-"US"

uspredquant.SP.pct.arr<-apply(uspred.SP.pct.arr,1,quantile,
                             prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
uspredquant.SP.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  uspredquant.SP.pct.arr,stringsAsFactors = FALSE),id_cols = 2,names_from=1,values_from=3))
tmpnames<-names(uspredquant.SP.pct)
names(uspredquant.SP.pct)<-sub("X","SP.pct.",make.names(tmpnames))
uspredquant.SP.pct$state<-"US"

uspredquant.dI.pct.arr<-apply(uspred.dI.pct.arr,1,quantile,
                             prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
uspredquant.dI.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  uspredquant.dI.pct.arr,stringsAsFactors = FALSE),id_cols = 2,names_from=1,values_from=3))
tmpnames<-names(uspredquant.dI.pct)
names(uspredquant.dI.pct)<-sub("X","dI.pct.",make.names(tmpnames))
uspredquant.dI.pct$state<-"US"

uspredquant.ID.pct.arr<-apply(uspred.ID.pct.arr,1,quantile,
                              prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
uspredquant.ID.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  uspredquant.ID.pct.arr,stringsAsFactors = FALSE),id_cols = 2,names_from=1,values_from=3))
tmpnames<-names(uspredquant.ID.pct)
names(uspredquant.ID.pct)<-sub("X","ID.pct.",make.names(tmpnames))
uspredquant.ID.pct$state<-"US"

uspredquant.Itot.pct.arr<-apply(uspred.Itot.pct.arr,1,quantile,
                              prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
uspredquant.Itot.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  uspredquant.Itot.pct.arr,stringsAsFactors = FALSE),id_cols = 2,names_from=1,values_from=3))
tmpnames<-names(uspredquant.Itot.pct)
names(uspredquant.Itot.pct)<-sub("X","Itot.pct.",make.names(tmpnames))
uspredquant.Itot.pct$state<-"US"

uspredquant.IfracU.pct.arr<-apply(uspred.IfracU.pct.arr,1,quantile,
                                prob=c(0.025,0.25,0.5,0.75,0.975),na.rm=TRUE)
uspredquant.IfracU.pct<-as.data.frame(pivot_wider(as.data.frame.table(
  uspredquant.IfracU.pct.arr,stringsAsFactors = FALSE),id_cols = 2,names_from=1,values_from=3))
tmpnames<-names(uspredquant.IfracU.pct)
names(uspredquant.IfracU.pct)<-sub("X","IfracU.pct.",make.names(tmpnames))
uspredquant.IfracU.pct$state<-"US"

uspredquant<-left_join(left_join(left_join(left_join(inner_join(
  uspredquant.I.pct,uspredquant.SP.pct),uspredquant.dI.pct),
  uspredquant.ID.pct),uspredquant.Itot.pct),uspredquant.IfracU.pct)

fwrite(uspredquant,file.path(mcmcpath,"US_Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))
fwrite(uspredquant,file.path(apppath,mcmcpath,"US_Posterior_I_SP_dI_ID_Itot_IfracU_Quantiles.csv"))

# drive_put(file.path(mcmcpath,"US_Posterior_I_SP_dI_ID_Itot_Quantiles.csv"),
#           path=file.path(googledrivepath,mcmcpath,"US_Posterior_I_SP_dI_ID_Itot_Quantiles.csv"))


