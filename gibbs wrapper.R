run.gibbs=function(dat,ngibbs,a.phi,b.phi,gamma,ncomm,nadapt){
  # rm(list=ls(all=TRUE))
  # library(Rcpp)
  # set.seed(4)
  # 
  # setwd('U:\\presence absence model\\github-presence_absence-SB')
  # source('gibbs functions.R')
  # source('gibbs wrapper.R')
  # sourceCpp('aux1.cpp')
  # 
  # dat=read.csv('fake data y.csv',as.is=T)
  # ngibbs=10000
  # a.phi=0.01
  # b.phi=0.99
  # gamma=0.1
  # ncomm=5
  # nadapt=1000
  
  #convert from a bunch of bernoulli to a single binomial per location
  tmp=aggregate.data(dat)
  y=tmp$dat
  loc.id=tmp$loc.id
  nspp=ncol(y)
  nloc=length(unique(loc.id))
  n=tmp$n
  nmat=matrix(n,nloc,nspp)
  
  #initial values
  theta=matrix(1/ncomm,nloc,ncomm)
  vmat=theta
  vmat[,ncomm]=1
  phi=matrix(0.2,ncomm,nspp,byrow=T)

  #gibbs stuff
  vec.theta=matrix(0,ngibbs,nloc*ncomm)
  vec.phi=matrix(0,ngibbs,ncomm*nspp)
  vec.logl=matrix(NA,ngibbs,1)
  param=list(theta=theta,phi=phi,vmat=vmat)
  
  #for MH algorithm
  jump1=list(vmat=matrix(0.1,nloc,ncomm),
             phi=matrix(0.1,ncomm,nspp))
  accept1=list(vmat=matrix(0,nloc,ncomm),
               phi=matrix(0,ncomm,nspp))
  accept.output=50
  
  options(warn=2)
  count=0
  for (i in 1:ngibbs){
    print(i)
    tmp=update.phi(param=param,jump=jump1$phi,ncomm=ncomm,nspp=nspp,y=y,nmat=nmat,a.phi=a.phi,b.phi=b.phi)
    param$phi=tmp$phi #phi.true #
    accept1$phi=accept1$phi+tmp$accept
    
    tmp=update.theta(param=param,jump=jump1$vmat,nloc=nloc,ncomm=ncomm,y=y,nmat=nmat,gamma=gamma)
    param$theta=tmp$theta
    param$vmat=tmp$v
    accept1$vmat=accept1$vmat+tmp$accept
    
    if (i%%accept.output==0 & i<nadapt){
      k=print.adapt(accept1=accept1,jump1=jump1,accept.output=accept.output)
      accept1=k$accept1
      jump1=k$jump1
    }
    
    #to assess convergence, examine logl
    prob=get.logl(theta=param$theta,phi=param$phi,y=y,nmat=nmat)
    loglikel=sum(prob)+
      sum(dbeta(param$phi,a.phi,b.phi,log=T))+
      sum(dbeta(param$vmat[,-ncomm],1,gamma,log=T))
    
    vec.logl[i]=loglikel
    vec.theta[i,]=param$theta
    vec.phi[i,]=param$phi
  }
  list(logl=vec.logl,theta=vec.theta,phi=vec.phi)
}