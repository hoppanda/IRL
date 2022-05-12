## Functions -----------------------------------------------------------------
##: (1) To generate simulated data 
simData=function(nbat,ntime,trueSlope=-0.5,varSlope=1,
                 trueInt=100,varInt=1,varRes=2,batname="",
                 seed=12345){
  if(length(ntime)==1) ntime=rep(ntime,nbat)
  if(any(ntime>6)) stop("ntime is too large (>6)")
  tgrid=c(0,6,12,18,24,36)
  
  set.seed(seed)
  slopes=rnorm(nbat,trueSlope,sqrt(varSlope))
  ints=rnorm(nbat,trueInt,sqrt(varInt))
  raw=lapply(1:nbat,function(i){
    y=rnorm(ntime[i],ints[i]+slopes[i]*tgrid[1:ntime[i]],sqrt(varRes))
    data.frame(batch=paste(batname,"Batch",i,sep="-"),time=tgrid[1:ntime[i]],
               res=y)
  })
  do.call("rbind",raw)
}

##: (2) simulation-based calculation of probability mass for multivariate gaussin
pmv=function(mu,S,lower,B=10000,seed=mu[1]*777){
  set.seed(seed)
  temp=rmvnorm(B,mean=mu,sigma=S)
  sum(temp[,1]>=lower[1] & temp[,2]>=lower[2])/B
}
pmv2=function(mu,S,lower,B=10000,seed=mu[1]*777){
  set.seed(seed)
  temp=rmvnorm(B,mean=mu,sigma=S)
  pp1=sum(temp[,1]>=lower[1] & temp[,2]>=lower[2])/B
  pp2=sum(temp[,1]>=lower[1])/B
  pp1/pp2
}
pmv3=function(mu,S,grid,lower2,margin=F){
  if(!margin){
    tt=sapply(grid,function(xx){
      pp1=pmvnorm(lower=c(xx,lower2),upper=c(+Inf,+Inf),mean=mu,sigma=S)
      pp2=1-pnorm(xx,mu[1],sqrt(S[1,1])) 
      pp1/pp2
    })
  }else{
    tt=sapply(grid,function(xx){1-pnorm(xx,mu[1],sqrt(S[1,1]))})
  }
  tt
}
pmv.bayes=function(a,b,sa,sb,se,t=36,lower,grid){
  mu0=a
  muT=a+b*t
  v0=sa^2+se^2
  vT=sa^2+se^2+(sb^2)*(t^2)
  c0T=sa^2
  mu=cbind(mu0,muT)
  S=lapply(seq(mu0),function(i) matrix(c(v0[i],c0T[i],c0T[i],vT[i]),2,2))
  S=simplify2array(S)
  
  allres=matrix(NA,nrow(mu),length(grid))
  for(j in 1:nrow(mu)){
    if(j%%100==0)  print(paste0(j,"/",nrow(mu)))
    set.seed(12345+j)
    temp=pmv3(mu=mu[j,],S=S[,,j],grid=grid,lower2=lower,margin=F) 
    if(any(round(temp,5)==1)|any(temp>1)) temp[min(which(round(temp,5)==1),which(temp>1)):length(temp)]=1
    allres[j,]=temp
  }
  #allres[is.nan(allres)]=1
  fin.res=apply(allres,2,mean,trim=0.1,na.rm=T)
  list(x=grid,y=fin.res)
}

##: (3) to derive values for population model parameters while being conditional on some
setup=function(a=NA,b=NA,b.sd=0.05,a.sd=NA,rho=0.6,p0=NA,p1=NA,lower,sl){
  if(!is.na(a)&!is.na(b)){
    res.sd=sqrt((1-rho)/rho)*a.sd
    m.sl=a+b*sl
    v.0=a.sd^2+res.sd^2
    v.sl=a.sd^2+res.sd^2+b.sd^2*(sl^2)
    ccv=matrix(c(v.0,rep(a.sd^2,2),v.sl),2,2)
    ccr=cov2cor(ccv)
  }
  if(!is.na(p0)&!is.na(p1)){
    if(!is.na(a)){
      a.sd=(lower-a)/qnorm(1-p0)
      res.sd=sqrt((1-rho)/rho)*a.sd
      v.0=a.sd^2+res.sd^2
      v.sl=a.sd^2+res.sd^2+b.sd^2*(sl^2)
      ccv=matrix(c(v.0,rep(a.sd^2,2),v.sl),2,2)
      ccr=cov2cor(ccv)
      b=(lower-sqrt(v.sl)*qnorm(1-p1)-a)/sl
      m.sl=a+b*sl
    }else{
      res.sd=sqrt((1-rho)/rho)*a.sd
      v.0=a.sd^2+res.sd^2
      v.sl=a.sd^2+res.sd^2+b.sd^2*(sl^2)
      ccv=matrix(c(v.0,rep(a.sd^2,2),v.sl),2,2)
      ccr=cov2cor(ccv)
      a=lower-sqrt(v.0)*qnorm(1-p0)
      b=(lower-sqrt(v.sl)*qnorm(1-p1)-a)/sl
      m.sl=a+b*sl
    }
  }
  if(is.na(p0)) p0=1-pnorm(lower,a,sqrt(v.0))
  if(is.na(p1)) p1=1-pnorm(lower,a+b*sl,sqrt(v.sl))
  list(a=a,b=b,a.sd=a.sd,b.sd=b.sd,res.sd=res.sd,m.sl=m.sl,
       v.0=v.0,v.sl=v.sl,ccv=ccv,ccr=ccr,p0=p0,p1=p1)
}

## (4): to generate a curve
curvy=function(sce,eta.grid,choice=1,lower2=95,seed0=sce$a*777,margin=F,spar=0.6){
  if(choice==1){
    set.seed(12345)
    yy=pmv3(mu=c(sce$a,sce$m.sl),S=sce$ccv,grid=eta.grid,lower2=lower2,margin=margin) #B=80000,seed=seed0,
    if(!margin){if(any(round(yy,5)>=1)) yy[min(which(round(yy,5)>=1)):length(yy)]=1}
  }else{
    yy=sapply(eta.grid,function(ee) pmv2(mu=c(sce$a,sce$m.sl),S=sce$ccv,lower=c(ee,lower2),B=50000,seed=ee*777))
    yy[is.nan(yy)]=1
  }
  list(eta=eta.grid,yy=yy) #yy2=predict(smooth.spline(x=eta.grid,y=yy,spar=spar))$y
}

# (5): make the plot with stability data and calculated IRL
plotfun = function(sce.u,title=""){
  dat=simData(nbat=10,ntime=6,trueSlope=sce.u$b,varSlope=sce.u$b.sd^2,
              trueInt=sce.u$a,varInt=sce.u$a.sd^2,varRes=sce.u$res.sd^2,
              batname="BB",seed=sce.u$a*1000)
  spec=c(95,105)
  sl=36
  plotdat=split(dat,f=dat$batch)
  xplot=unique(dat$time)
  yrange=range(dat$res,spec,na.rm=T)
  xrange=range(xplot,sl,na.rm=T)
  colfunc=colorRampPalette(brewer.pal(8,"Dark2"))
  colors=col2rgb(colfunc(length(plotdat)))
  colors=rgb(colors[1,],colors[2,],colors[3,],150,max=255)
  
  plot(0,1,type="n",main=title,xlab="Time in month",ylab="Reportable value",
       xlim=xrange,ylim=c(94,102),yaxt="n",xaxt="n")
  axis(1,at=c(0,6,12,24,36),labels=c(0,6,12,24,36))
  axis(2,las=1)
  for(i in seq(plotdat)){
    ptid1=plotdat[[i]]$time<=24
    ptid2=plotdat[[i]]$time>=24
    points(plotdat[[i]]$res~plotdat[[i]]$time,col=colors[i],pch=16,cex=0.9)
    lines(plotdat[[i]]$res[ptid1]~plotdat[[i]]$time[ptid1],col=colors[i])
    lines(plotdat[[i]]$res[ptid2]~plotdat[[i]]$time[ptid2],col=colors[i],lty=2)
  }
  #abline(v=36,col="gray")
  abline(h=c(95),lty=1,col="coral3")
}
salut=function(eta,pr,op.pr=0.95,PLOT=F){
  tmod=smooth.spline(x=pr,y=eta)
  op.eta=predict(tmod,x=op.pr)$y
  op.eta
}

# (6) derive IRL given example data
output=function(sce.u,xgrid=seq(90,105,length=100),use_prior="nc-0.1"){
  dat=simData(nbat=10,ntime=6,trueSlope=sce.u$b,varSlope=sce.u$b.sd^2,
              trueInt=sce.u$a,varInt=sce.u$a.sd^2,varRes=sce.u$res.sd^2,
              batname="BB",seed=sce.u$a*1000)
  dat=dat[dat$time<=24, ]
  
  #:-> Solution in ADG
  mod.adg1=lm(res~time,dat)
  mod.adg2=lm(res~batch*time,dat)
  slp.se2=summary(mod.adg1)$coef[2,2]^2
  res.var=sigma(mod.adg2)^2
  df1=mod.adg1$df.residual
  df2=mod.adg2$df.residual
  df.adg=((slp.se2+res.var)^2/(slp.se2^2/df1+res.var^2/df2))
  adg.irl=95-coef(mod.adg1)[2]*36+qt(0.95,df.adg)*sqrt(slp.se2*36^2+res.var)
  
  #:-> Solution in CoT
  mod=lmer(res~time+(1+time||batch),dat)
  #summary(mod)
  var.mod=as.data.frame(VarCorr(mod))
  cot.irl=95-(fixef(mod)[2]+qnorm(1-0.8)*var.mod$sdcor[2])*36 
  
  #:-> Solution in CoI Bayes
  if(use_prior=="nc-0.1"){
    fit.b=brm(res~time+(1+time||batch),data=dat,family=gaussian(),
              prior=c(set_prior("normal(100,100)",class="Intercept"),
                      set_prior("normal(0,100)",class="b",coef="time"),
                      set_prior("cauchy(0,0.1)",class="sd")),
              warmup=2000,iter=5000,chains=2,
              control=list(adapt_delta=0.98),
              seed=12345)
  }else if(use_prior=="nc-0.1-info"){
    fit.b=brm(res~time+(1+time||batch),data=dat,family=gaussian(),
              prior=c(set_prior("normal(100,30)",class="Intercept"),
                      set_prior("normal(0,5)",class="b",coef="time"),
                      set_prior("cauchy(0,0.1)",class="sd")),
              warmup=2000,iter=5000,chains=2,
              control=list(adapt_delta=0.98),
              seed=12345)
  }else if(use_prior=="nc-0.01"){
    fit.b=brm(res~time+(1+time||batch),data=dat,family=gaussian(),
              prior=c(set_prior("normal(100,100)",class="Intercept"),
                      set_prior("normal(0,100)",class="b",coef="time"),
                      set_prior("cauchy(0,0.01)",class="sd")),
              warmup=2000,iter=5000,chains=2,
              control=list(adapt_delta=0.98),
              seed=12345)
  }else if(use_prior=="nc-0.01-info"){
    fit.b=brm(res~time+(1+time||batch),data=dat,family=gaussian(),
              prior=c(set_prior("normal(100,30)",class="Intercept"),
                      set_prior("normal(0,5)",class="b",coef="time"),
                      set_prior("cauchy(0,0.01)",class="sd")),
              warmup=2000,iter=5000,chains=2,
              control=list(adapt_delta=0.98),
              seed=12345)
  }
  
  # summary(fit.b)
  # str(fit.b)
  check=as.data.frame(fit.b$fit)
  BayeRes=pmv.bayes(a=check$b_Intercept,
                    b=check$b_time,
                    sa=check$sd_batch__Intercept,
                    sb=check$sd_batch__time,
                    se=check$sigma,
                    t=36,
                    lower=95,
                    grid=xgrid)
  #plot(BayeRes$y~BayeRes$x)
  
  cp.b1=salut(eta=BayeRes$x,pr=BayeRes$y,op.pr=0.95)
  cp.b2.raw=sapply(1:2000,function(i){
    # for(i in 1:2000){
    #if(i%%100 ==0) print(i)
    temp.sce=setup(a=check$b_Intercept[i],b=check$b_time[i],
                   a.sd=check$sd_batch__Intercept[i],b.sd=check$sd_batch__time[i],
                   rho=check$sd_batch__Intercept[i]^2/(check$sd_batch__Intercept[i]^2+check$sigma[i]^2),
                   lower=95,sl=36)
    set.seed(12345+i)
    #xgrid=seq(85,115,length=100)
    temp=pmv3(mu=c(temp.sce$a,temp.sce$m.sl),S=temp.sce$ccv,grid=xgrid,lower2=95,margin=F) 
    if(any(round(temp,5)==1)|any(temp>1)) temp[min(which(round(temp,5)==1),which(temp>1)):length(temp)]=1
    if(temp.sce$p1<0.95){
      pid=min(which(temp>0.95))
      pid= c(pid-1,pid+1)
      if(pid[1]<1) pid[1]=1
      if(pid[1]>length(temp)) pid[2]=length(temp)
      approx(x=temp[pid[1]:pid[2]],y=xgrid[pid[1]:pid[2]],xout=0.95,ties="ordered")$y
    }else{NA}
    # }
  })
  cp.b2.med=median(cp.b2.raw,na.rm=T)
  list(adg.irl=adg.irl, 
       cot.irl=cot.irl,
       coi.b1=cp.b1,
       coi.b2=cp.b2.med,
       mcmc=fit.b)
}

# (7) derive IRL given example data
extMCMC = function(outres){
  extdat = summary(outres$mcmc)
  usedat = rbind(extdat$fixed,
                 extdat$spec_pars,
                 extdat$random$batch)
  usedat
}







