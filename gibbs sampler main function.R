lda.presence.absence=function(dat,id,ncomm,a.phi,b.phi,gamma,ngibbs){
  #useful stuff
  dat1=data.matrix(dat)
  loc.id=id
  nloc=max(loc.id)
  nspp=ncol(dat1)
  nlinhas=nrow(dat1)
  hi=0.999999
  lo=0.000001
  #initial values
  theta=matrix(1/ncomm,nloc,ncomm)
  phi=matrix(0.5,ncomm,nspp)
  z=matrix(sample(1:ncomm,size=nlinhas*nspp,replace=T),nlinhas,nspp)
  #to store outcomes from gibbs sampler
  theta.out=matrix(NA,ngibbs,ncomm*nloc)
  
  phi.out=matrix(NA,ngibbs,ncomm*nspp)
  llk=rep(NA,ngibbs)
  #run gibbs sampler
  options(warn=2)
  for (i in 1:ngibbs){
    print(i)
    #sample z
    rand.u=matrix(runif(nlinhas*nspp),nlinhas,nspp)
    z=samplez(log(theta), log(1-phi), log(phi), dat1, loc.id,rand.u, ncomm, nloc)
    #calculate summaries that depend on z and the data
    tmp=getks(z=z, ncommun=ncomm, dat=dat1)
    nks1=tmp$nks1
    nks0=tmp$nks0
    nlk=getlk(z=z,locid=loc.id, ncommun=ncomm, nloc=nloc)
    #get theta parameters
    theta=get.theta(nlk,gamma,ncomm,nloc)
    theta[theta>hi]=hi; theta[theta<lo]=lo
    #get phi parameters
    phi=matrix(rbeta(nspp*ncomm,nks1+a.phi,nks0+b.phi),ncomm,nspp)
    phi[phi>hi]=hi; phi[phi<lo]=lo
    #calculate probability after integrating out the latent z's
    prob=theta%*%phi
    prob[prob>hi]=hi; prob[prob<lo]=lo
    prob1=prob[loc.id,]
    #calculate logl and store results
    llk[i]=sum(dat1*log(prob1)+(1-dat1)*log(1-prob1))
    theta.out[i,]=theta
    phi.out[i,]=phi
  }
  list(llk=llk,theta=theta.out,phi=phi.out)
}