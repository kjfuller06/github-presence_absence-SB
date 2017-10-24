run.gibbs=function(dat,ngibbs,a.phi,b.phi,gamma,ncomm,nadapt){
  #convert from a bunch of bernoulli to a single binomial per location
  tmp=aggregate.data(dat)
  y=tmp$dat
  loc.id=tmp$loc.id
  nspp=ncol(y)
  nloc=length(unique(loc.id))
  n=tmp$n
  nmat=matrix(n,nloc,nspp)
  
  #get good initial values using k-means clustering
  z=kmeans(y/nmat,ncomm,nstart=20)
  theta=matrix(NA,nloc,ncomm)
  seq1=1:ncomm
  for (i in 1:ncomm){
    cond=z$cluster==i
    theta[cond,i]=0.8
    seq2=seq1[-i]
    theta[cond,seq2]=0.2/9
  }
  vmat=matrix(NA,nloc,ncomm)
  vmat[,1]=theta[,1]
  vmat[,2]=theta[,2]/(1-vmat[,1])
  for (i in 3:(ncomm-1)) vmat[,i]=theta[,i]/apply(1-vmat[,1:(i-1)],1,prod)
  vmat[,ncomm]=1
  
  phi=data.matrix(z$centers)
  cond=phi==0; phi[cond]=0.01
  cond=phi==1; phi[cond]=0.99

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
    param$phi=tmp$phi
    accept1$phi=accept1$phi+tmp$accept
    # param$phi=phi.true
    
    tmp=update.theta(param=param,jump=jump1$vmat,nloc=nloc,ncomm=ncomm,y=y,nmat=nmat,gamma=gamma)
    param$theta=tmp$theta
    param$vmat=tmp$v
    accept1$vmat=accept1$vmat+tmp$accept
    # param$theta=theta.true
    
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