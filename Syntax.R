## Loading Packages -----------------------------------------------------------------
library(lme4)
library(corrplot)
library(RColorBrewer)
library(mvtnorm)
library(ggplot2)
library(rstan)
library(rstanarm)
library(brms)
library(reshape2)
library(dplyr)
source(file.path(getwd(),"functions.R"))

## Generate Graphs -----------------------------------------------------------------
## ==========================================
## Illustration plot Fig 2.1
## ==========================================
(sce4a=setup(a=NA,b=NA,b.sd=0.03,a.sd=1,rho=0.8,p0=0.999,p1=0.9,lower=95,sl=24))
(sce4b=setup(a=NA,b=NA,b.sd=0.05,a.sd=1,rho=0.7,p0=0.999,p1=0.9,lower=95,sl=24))
# pmvnorm(lower=c(98,95),mean=c(a,m.sl),sigma=cc)
# pmvnorm(lower=c(98,95),mean=c(a,m.sl),corr=cov2cor(cc))
pd=curvy(sce=sce4b,eta.grid=seq(90,105,by=0.2),choice=1,lower2=95,seed0=201888)
pd2=curvy(sce=sce4a,eta.grid=seq(90,105,by=0.2),choice=1,lower2=95,seed0=20188,spar=0.45)
pd3=curvy(sce=sce4b,eta.grid=seq(90,105,by=0.2),choice=1,lower2=95,seed0=201888,margin=T)
pd4=curvy(sce=sce4a,eta.grid=seq(90,105,by=0.2),choice=1,lower2=95,seed0=201888,margin=T)
eta1=approx(x=pd$yy,y=pd$eta,xout=0.95)$y
eta2=approx(x=pd2$yy,y=pd2$eta,xout=0.95)$y

#windows(width=7,height=5)
png(file="Figure_1.png",width=1500, height=550,type="quartz",pointsize=25)
par(mfrow=c(1,2),font.main=6,font.lab=6,font.axis=6,cex.main=1,cex.axis=0.9,cex.lab=1,
    mgp=c(2.1,0.3,0),tcl=-0.15,mar=c(4,4,1,1))
## subplot 1
plot(pd$yy~pd$eta,type="l",lwd=2,xlim=c(92,102),ylim=c(0.7,1),col="navyblue",
     yaxt="n",xlab=bquote(""~eta[CoI]),ylab="Probability",main="Modestly Big Noise")
axis(2,las=1)
abline(h=0.95,lty=2,lwd=1,col=rgb(0,0,0,100,max=255))
abline(v=eta1,lty=2,lwd=1,col=rgb(0,0,0,100,max=255))
text(x=eta1,y=0.95,labels=bquote("("~eta[CoI]~","~q~")"),cex=1.1,col="deeppink",pos=4)
points(x=eta1,y=0.95,pch=16,cex=1.2,col="deeppink")
text(x=eta1,y=1-pnorm(eta1,sce4b$a,sqrt(sce4b$v.0)),
     labels=bquote("Pr("~Y[i0]>=eta[CoI]~")"),
     cex=0.9,col="deeppink",pos=4)
points(x=eta1,y=1-pnorm(eta1,sce4b$a,sqrt(sce4b$v.0)),
       pch=1,cex=1.2,col="deeppink")
lines(pd3$yy~pd3$eta,lty=4,lwd=2,col="navyblue")
## subplot 2
plot(pd2$yy~pd$eta,type="l",lwd=2,xlim=c(92,102),ylim=c(0.7,1),col="coral3",
     yaxt="n",xlab=bquote(""~eta[CoI]),ylab="Probability",main="Modestly Small Noise")
axis(2,las=1)
abline(h=0.95,lty=2,lwd=1,col=rgb(0,0,0,100,max=255))
abline(v=eta2,lty=2,lwd=1,col=rgb(0,0,0,100,max=255))
lines(pd4$yy~pd2$eta,lty=4,lwd=2,col="coral3")
text(x=eta2,y=0.95,labels=bquote("("~eta[CoI]~","~q~")"),cex=1.1,col="deeppink",pos=4)
points(x=eta2,y=0.95,pch=16,cex=1.2,col="deeppink")
text(x=eta2,y=1-pnorm(eta2,sce4a$a,sqrt(sce4a$v.0)),
     labels=bquote("Pr("~Y[i0]>=eta[CoI]~")"),
     cex=0.9,col="deeppink",pos=4)
points(x=eta2,y=1-pnorm(eta2,sce4a$a,sqrt(sce4a$v.0)),
       pch=1,cex=1.2,col="deeppink")
dev.off()
## ==========================================
## Simulation fo CoI
## ==========================================
sl=36
delta_sl=c(2,3)
sd_sl_wT=seq(0.3,1,by=0.1)
sd_e=seq(0.1,1.5,by=0.2)
q=c(0.95,0.99)
cases=expand.grid(sd_sl_wT,sd_e)
names(cases)=c("sd_sl_wT","sd_e")

res_perQ_3=lapply(q,function(q_ind){
  temp=cases; temp$eta=""; temp$pr0_eta=""
  for(j in 1:nrow(temp)){
    if(j%%8==0) print(paste0("Running ",j,"/",nrow(temp)))
    test_sce=setup(a=100,b=-delta_sl[2]/sl,b.sd=temp$sd_sl_wT[j]/sl,
                   a.sd=0.5,rho=0.25/(0.25+temp$sd_e[j]^2),p0=NA,p1=NA,lower=95,sl=sl)
    if(q_ind<test_sce$p1){
      pr_0_eta=NA; eta_coi="<95"  
    }else{
      test_res=curvy(sce=test_sce,eta.grid=seq(90,110,by=0.05),choice=1,lower2=95,seed0=1234,spar=0.7)
      if(max(test_res$yy)<q_ind){
        pr_0_eta=0; eta_coi=">110"
      }else{
        pid=min(which(test_res$yy>q_ind))
        pid= c(pid-1,pid+1)
        if(pid[1]<1) pid[1]=1
        if(pid[1]>length(test_res$yy)) pid[2]=length(test_res$yy)
        eta_coi=approx(x=test_res$yy[pid[1]:pid[2]],y=test_res$eta[pid[1]:pid[2]],xout=q_ind)$y
        pr_0_eta=1-pnorm(eta_coi,test_sce$a,sqrt(test_sce$v.0)) 
      }
      rm(list="test_res")
    }
    temp$eta[j]=as.character(eta_coi)
    temp$pr0_eta[j]=as.character(pr_0_eta)
    rm(list=c("eta_coi","pr_0_eta","test_sce"))
  }
  temp
})
names(res_perQ_3)=paste0("q=",q)
res_perQ_2=lapply(q,function(q_ind){
  temp=cases; temp$eta=""; temp$pr0_eta=""
  for(j in 1:nrow(temp)){
    if(j%%8==0) print(paste0("Running ",j,"/",nrow(temp)))
    test_sce=setup(a=100,b=-delta_sl[1]/sl,b.sd=temp$sd_sl_wT[j]/sl,
                   a.sd=0.5,rho=0.25/(0.25+temp$sd_e[j]^2),p0=NA,p1=NA,lower=95,sl=sl)
    if(q_ind<test_sce$p1){
      pr_0_eta=NA; eta_coi="<95"  
    }else{
      test_res=curvy(sce=test_sce,eta.grid=seq(90,110,by=0.05),choice=1,lower2=95,seed0=1234,spar=0.7)
      if(max(test_res$yy)<q_ind){
        pr_0_eta=0; eta_coi=">110"
      }else{
        pid=min(which(test_res$yy>q_ind))
        pid= c(pid-1,pid+1)
        if(pid[1]<1) pid[1]=1
        if(pid[1]>length(test_res$yy)) pid[2]=length(test_res$yy)
        eta_coi=approx(x=test_res$yy[pid[1]:pid[2]],y=test_res$eta[pid[1]:pid[2]],xout=q_ind)$y
        pr_0_eta=1-pnorm(eta_coi,test_sce$a,sqrt(test_sce$v.0)) 
      }
      rm(list="test_res")
    }
    temp$eta[j]=as.character(eta_coi)
    temp$pr0_eta[j]=as.character(pr_0_eta)
    rm(list=c("eta_coi","pr_0_eta","test_sce"))
  }
  temp
})
names(res_perQ_2)=paste0("q=",q)

res_perQ_3b=lapply(q,function(q_ind){
  temp=cases; temp$eta=""; temp$pr0_eta=""
  for(j in 1:nrow(temp)){
    if(j%%8==0) print(paste0("Running ",j,"/",nrow(temp)))
    test_sce=setup(a=98.5,b=-delta_sl[2]/sl,b.sd=temp$sd_sl_wT[j]/sl,
                   a.sd=0.5,rho=0.25/(0.25+temp$sd_e[j]^2),p0=NA,p1=NA,lower=95,sl=sl)
    if(q_ind<test_sce$p1){
      pr_0_eta=NA; eta_coi="<95"  
    }else{
      test_res=curvy(sce=test_sce,eta.grid=seq(90,110,by=0.05),choice=1,lower2=95,seed0=1234,spar=0.7)
      if(max(test_res$yy)<q_ind){
        pr_0_eta=0; eta_coi=">110"
      }else{
        pid=min(which(test_res$yy>q_ind))
        pid= c(pid-1,pid+1)
        if(pid[1]<1) pid[1]=1
        if(pid[1]>length(test_res$yy)) pid[2]=length(test_res$yy)
        eta_coi=approx(x=test_res$yy[pid[1]:pid[2]],y=test_res$eta[pid[1]:pid[2]],xout=q_ind)$y
        pr_0_eta=1-pnorm(eta_coi,test_sce$a,sqrt(test_sce$v.0)) 
      }
      rm(list="test_res")
    }
    temp$eta[j]=as.character(eta_coi)
    temp$pr0_eta[j]=as.character(pr_0_eta)
    rm(list=c("eta_coi","pr_0_eta","test_sce"))
  }
  temp
})
names(res_perQ_3b)=paste0("q=",q)
res_perQ_2b=lapply(q,function(q_ind){
  temp=cases; temp$eta=""; temp$pr0_eta=""
  for(j in 1:nrow(temp)){
    if(j%%8==0) print(paste0("Running ",j,"/",nrow(temp)))
    test_sce=setup(a=98.5,b=-delta_sl[1]/sl,b.sd=temp$sd_sl_wT[j]/sl,
                   a.sd=0.5,rho=0.25/(0.25+temp$sd_e[j]^2),p0=NA,p1=NA,lower=95,sl=sl)
    if(q_ind<test_sce$p1){
      pr_0_eta=NA; eta_coi="<95"  
    }else{
      test_res=curvy(sce=test_sce,eta.grid=seq(90,110,by=0.05),choice=1,lower2=95,seed0=1234,spar=0.7)
      if(max(test_res$yy)<q_ind){
        pr_0_eta=0; eta_coi=">110"
      }else{
        pid=min(which(test_res$yy>q_ind))
        pid= c(pid-1,pid+1)
        if(pid[1]<1) pid[1]=1
        if(pid[1]>length(test_res$yy)) pid[2]=length(test_res$yy)
        eta_coi=approx(x=test_res$yy[pid[1]:pid[2]],y=test_res$eta[pid[1]:pid[2]],xout=q_ind)$y
        pr_0_eta=1-pnorm(eta_coi,test_sce$a,sqrt(test_sce$v.0)) 
      }
      rm(list="test_res")
    }
    temp$eta[j]=as.character(eta_coi)
    temp$pr0_eta[j]=as.character(pr_0_eta)
    rm(list=c("eta_coi","pr_0_eta","test_sce"))
  }
  temp
})
names(res_perQ_2b)=paste0("q=",q)


coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)
colfunc=colorRampPalette(c("#D35099", "#E181B5"))
coul <- c(colfunc(60),coul[6:length(coul)])

# Plot for A=100 --------------------------------------------------------
par(mfrow=c(2,2),cex.main=0.8)

raw_zz=res_perQ_2$`q=0.95`
raw_zz$pr0_eta=as.numeric(raw_zz$pr0_eta)
zz =unstack(form=formula(raw_zz[,c("pr0_eta","sd_e")]),x=raw_zz)
colnames(zz)= gsub("X", "s_e=",colnames(zz)," ")
rownames(zz)= paste0("T*s_b=",sd_sl_wT," ")
pmat=is.na(zz)
zz[which(is.na(zz),arr.ind=T)]=1
corrplot(as.matrix(zz),col=coul,order="original",is.corr=F,na.label="X",
         addCoef.col = "black",tl.col="black", tl.cex=0.7,tl.srt=90,method="color",
         addgrid.col="lightgrey", number.cex=0.7,number.digits = 3,number.font =6,
         title="Modest Degradation till T (bT=2), Assurance request q=0.95",p.mat=pmat,sig.level=0.1,
         mar=c(0,0,3,0),cl.cex = 0.7,cl.pos = "n", cl.lim=c(0,1)) #tl.srt=45,

raw_zz=res_perQ_2$`q=0.99`
raw_zz$pr0_eta=as.numeric(raw_zz$pr0_eta)
zz =unstack(form=formula(raw_zz[,c("pr0_eta","sd_e")]),x=raw_zz)
colnames(zz)= gsub("X", "s_e=",colnames(zz)," ")
rownames(zz)= paste0("T*s_b=",sd_sl_wT," ")
pmat=is.na(zz)
zz[which(is.na(zz),arr.ind=T)]=1
corrplot(as.matrix(zz),col=coul,order="original",is.corr=F,na.label="X",
         addCoef.col = "black",tl.col="black", tl.cex=0.7,tl.srt=90,method="color",
         addgrid.col="lightgrey", number.cex=0.7,number.digits = 3,number.font =6,
         title="Modest Degradation till T (bT=2), Assurance request q=0.99",p.mat=pmat,sig.level=0.1,
         mar=c(0,0,3,0),cl.cex = 0.7,cl.offset = 3) #tl.srt=45,

raw_zz=res_perQ_3$`q=0.95`
raw_zz$pr0_eta=as.numeric(raw_zz$pr0_eta)
zz =unstack(form=formula(raw_zz[,c("pr0_eta","sd_e")]),x=raw_zz)
colnames(zz)= gsub("X", "s_e=",colnames(zz)," ")
rownames(zz)= paste0("T*s_b=",sd_sl_wT," ")
pmat=is.na(zz)
zz[which(is.na(zz),arr.ind=T)]=1
corrplot(as.matrix(zz),col=coul,order="original",is.corr=F,na.label="X",
         addCoef.col = "black",tl.col="black", tl.cex=0.7,tl.srt=90,method="color",
         addgrid.col="lightgrey", number.cex=0.7,number.digits = 3,number.font =6,
         title="Large Degradation till T (bT=3), Assurance request q=0.95",p.mat=pmat,sig.level=0.1,
         mar=c(0,0,3,0),cl.cex = 0.7,cl.offset = 3) #tl.srt=45,

raw_zz=res_perQ_3$`q=0.99`
raw_zz$pr0_eta=as.numeric(raw_zz$pr0_eta)
zz =unstack(form=formula(raw_zz[,c("pr0_eta","sd_e")]),x=raw_zz)
colnames(zz)= gsub("X", "s_e=",colnames(zz)," ")
rownames(zz)= paste0("T*s_b=",sd_sl_wT," ")
pmat=is.na(zz)
zz[which(is.na(zz),arr.ind=T)]=1
corrplot(as.matrix(zz),col=coul,order="original",is.corr=F,na.label="X",
         addCoef.col = "black",tl.col="black", tl.cex=0.7,tl.srt=90,method="color",
         addgrid.col="lightgrey", number.cex=0.7,number.digits = 3,number.font =6,
         title="Large Degradation till T (bT=3), Assurance request q=0.99",p.mat=pmat,sig.level=0.1,
         mar=c(0,0,3,0),cl.cex = 0.7,cl.offset = 3) #tl.srt=45,

# Plot for A=98.5 -------------------------------------------------------------------------
par(mfrow=c(2,2),cex.main=0.8)

raw_zz=res_perQ_2b$`q=0.95`
raw_zz$pr0_eta=as.numeric(raw_zz$pr0_eta)
zz =unstack(form=formula(raw_zz[,c("pr0_eta","sd_e")]),x=raw_zz)
colnames(zz)= gsub("X", "s_e=",colnames(zz)," ")
rownames(zz)= paste0("T*s_b=",sd_sl_wT," ")
pmat=is.na(zz)
zz[which(is.na(zz),arr.ind=T)]=1
corrplot(as.matrix(zz),col=coul,order="original",is.corr=F,na.label="X",
         addCoef.col = "black",tl.col="black", tl.cex=0.7,tl.srt=90,method="color",
         addgrid.col="lightgrey", number.cex=0.7,number.digits = 3,number.font =6,
         title="Modest Degradation till T (bT=2), Assurance request q=0.95",p.mat=pmat,sig.level=0.1,
         mar=c(0,0,3,0),cl.cex = 0.7,cl.offset = 3) #tl.srt=45,

raw_zz=res_perQ_2b$`q=0.99`
raw_zz$pr0_eta=as.numeric(raw_zz$pr0_eta)
zz =unstack(form=formula(raw_zz[,c("pr0_eta","sd_e")]),x=raw_zz)
colnames(zz)= gsub("X", "s_e=",colnames(zz)," ")
rownames(zz)= paste0("T*s_b=",sd_sl_wT," ")
pmat=is.na(zz)
zz[which(is.na(zz),arr.ind=T)]=1
corrplot(as.matrix(zz),col=coul,order="original",is.corr=F,na.label="X",
         addCoef.col = "black",tl.col="black", tl.cex=0.7,tl.srt=90,method="color",
         addgrid.col="lightgrey", number.cex=0.7,number.digits = 3,number.font =6,
         title="Modest Degradation till T (bT=2), Assurance request q=0.99",p.mat=pmat,sig.level=0.1,
         mar=c(0,0,3,0),cl.cex = 0.7,cl.offset = 3) #tl.srt=45,

raw_zz=res_perQ_3b$`q=0.95`
raw_zz$pr0_eta=as.numeric(raw_zz$pr0_eta)
zz =unstack(form=formula(raw_zz[,c("pr0_eta","sd_e")]),x=raw_zz)
colnames(zz)= gsub("X", "s_e=",colnames(zz)," ")
rownames(zz)= paste0("T*s_b=",sd_sl_wT," ")
pmat=is.na(zz)
zz[which(is.na(zz),arr.ind=T)]=1
corrplot(as.matrix(zz),col=coul,order="original",is.corr=F,na.label="X",
         addCoef.col = "black",tl.col="black", tl.cex=0.7,tl.srt=90,method="color",
         addgrid.col="lightgrey", number.cex=0.7,number.digits = 3,number.font =6,
         title="Large Degradation till T (bT=3), Assurance request q=0.95",p.mat=pmat,sig.level=0.1,
         mar=c(0,0,3,0),cl.pos="n",cl.lim=c(0,1)) #tl.srt=45,

raw_zz=res_perQ_3b$`q=0.99`
raw_zz$pr0_eta=as.numeric(raw_zz$pr0_eta)
zz =unstack(form=formula(raw_zz[,c("pr0_eta","sd_e")]),x=raw_zz)
colnames(zz)= gsub("X", "s_e=",colnames(zz)," ")
rownames(zz)= paste0("T*s_b=",sd_sl_wT," ")
pmat=is.na(zz)
zz[which(is.na(zz),arr.ind=T)]=1
corrplot(as.matrix(zz),col=coul,order="original",is.corr=F,na.label="X",
         addCoef.col = "black",tl.col="black", tl.cex=0.7,tl.srt=90,method="color",
         addgrid.col="lightgrey", number.cex=0.7,number.digits = 3,number.font =6,
         title="Large Degradation till T (bT=3), Assurance request q=0.99",p.mat=pmat,sig.level=0.1,
         mar=c(0,0,3,0),cl.pos="n",cl.lim=c(0,1)) #tl.srt=45,



## ==========================================
## Example Calculation
## ==========================================
(sce.u1=setup(a=NA,b=NA,b.sd=0.01,a.sd=0.5,rho=0.8,p0=0.9999,p1=0.9,lower=95,sl=36))
(sce.u2=setup(a=NA,b=NA,b.sd=0.04,a.sd=0.5,rho=0.4,p0=0.9999,p1=0.9,lower=95,sl=36))

#:-> Solution in Stan
out1=output(sce.u1,xgrid=seq(90,105,length=100))
out2=output(sce.u2,xgrid=seq(85,115,length=100))
#:-> different prior setting
out1b=output(sce.u1,xgrid=seq(90,105,length=100),use_prior = "nc-0.1-info")
out1c=output(sce.u1,xgrid=seq(90,105,length=100),use_prior = "nc-0.01")
out1d=output(sce.u1,xgrid=seq(90,105,length=100),use_prior = "nc-0.01-info")
out2b=output(sce.u2,xgrid=seq(85,115,length=100),use_prior = "nc-0.1-info")
out2c=output(sce.u2,xgrid=seq(85,115,length=100),use_prior = "nc-0.01")
out2d=output(sce.u2,xgrid=seq(85,115,length=100),use_prior = "nc-0.01-info")

# summary across different prior settings for comparison
estSummary1=lapply(c("1","1b","1c","1d"),function(yy){
  tmp = get(paste0("out",yy))
  as.matrix(extMCMC(tmp))
})
estSummary1=simplify2array(estSummary1)
dimnames(estSummary1)[[3]]=c("non-info-cauchy-0.1","info-cauchy-0.1","non-info-cauchy-0.01","info-cauchy-0.01")
for(i in 1:dim(estSummary1)[1]){
  cat(paste("\n==========================[",dimnames(estSummary1)[[1]][i],"]\n"))
  print(round(estSummary1[i,,],4))
}

estSummary2=lapply(c("2","2b","2c","2d"),function(yy){
  tmp = get(paste0("out",yy))
  as.matrix(extMCMC(tmp))
})
estSummary2=simplify2array(estSummary2)
dimnames(estSummary2)[[3]]=c("non-info-cauchy-0.1","info-cauchy-0.1","non-info-cauchy-0.01","info-cauchy-0.01")
for(i in 1:dim(estSummary2)[1]){
  cat(paste("\n==========================[",dimnames(estSummary2)[[1]][i],"]\n"))
  print(round(estSummary2[i,,],4))
}

#:-> use weakly informative prior, sd use half-cauchy with 0 location and 0.01 scale
# posterior diagnosis
mcmc_plot(out1d$mcmc,type="acf")
mcmc_plot(out1d$mcmc,type="trace")
coda::geweke.plot(as.mcmc(out1d$mcmc))
mcmc_plot(out2d$mcmc,type="acf")
mcmc_plot(out2d$mcmc,type="trace")
coda::geweke.plot(as.mcmc(out2d$mcmc))
launch_shinystan(out1d$mcmc)
launch_shinystan(out2d$mcmc)

# make plot 
leg1=c(paste("ADG =",round(out1b$adg.irl,2)),
       paste("CoT =",round(out1b$cot.irl,2)),
       paste("CoI-B1 =",round(out1b$coi.b1,2)),
       paste("CoI-B2 =",round(out1b$coi.b2,2))
)
leg2=c(paste("ADG =",round(out2b$adg.irl,2)),
       paste("CoT =",round(out2b$cot.irl,2)),
       paste("CoI-B1 =",round(out2b$coi.b1,2)),
       paste("CoI-B2 =",format(round(out2b$coi.b2,2),nsmall=2))
)
#windows(width=9,height=3)
png(file="Figure_4.png",width=1500, height=660,type="quartz",pointsize=25)
par(mfrow=c(1,2),font.axis=6,font.main=8,font.lab=6,mar=c(4,3,3,1),cex.main=1.1,
    mgp=c(1.7,0.35,0),cex.axis=0.9,tcl=-0.15)
plotfun(sce.u1,title="Modest Corrlation rho(0,T)")
legend("topright",bty="n",legend=leg1,text.font=6,text.col="blue",cex=0.85)
plotfun(sce.u2,title="Weak Corrlation rho(0,T)")
legend("topright",bty="n",legend=leg2,text.font=6,text.col="blue",cex=0.85)
dev.off()





