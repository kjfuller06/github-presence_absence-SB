run.foldin=function(dat,ngibbs,gamma,phi,nadapt){
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
  ncomm=nrow(phi)
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

  #gibbs stuff
  vec.theta=matrix(0,ngibbs,nloc*ncomm)
  vec.logl=matrix(NA,ngibbs,1)
  param=list(theta=theta,phi=phi,vmat=vmat)
  
  #for MH algorithm
  jump1=list(vmat=matrix(0.1,nloc,ncomm))
  accept1=list(vmat=matrix(0,nloc,ncomm))
  accept.output=50
  
  options(warn=2)
  count=0
  for (i in 1:ngibbs){
    print(i)
    param$phi=phi

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
      sum(dbeta(param$vmat[,-ncomm],1,gamma,log=T))
    
    vec.logl[i]=loglikel
    vec.theta[i,]=param$theta
  }
  list(logl=vec.logl,theta=vec.theta)
}