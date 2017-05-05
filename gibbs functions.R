aggregate.data=function(dat){
  dat1=dat[order(dat$loc.id),]
  nloc=max(dat1$loc.id)
  
  ind=which(colnames(dat1)=='loc.id')
  
  res=matrix(NA,nloc,ncol(dat1)-1)
  n=rep(NA,nloc)
  for (i in 1:nloc){
    cond=dat1$loc.id==i
    dat.tmp=dat1[cond,-ind]
    n[i]=nrow(dat.tmp)
    if (n[i]>1) dat.tmp=apply(dat.tmp,2,sum)
    if (n[i]==1) dat.tmp=as.numeric(dat.tmp)
    res[i,]=dat.tmp
  }
  colnames(res)=colnames(dat1[,-ind])
  loc.id=unique(dat1$loc.id)
  list(dat=res,loc.id=loc.id,n=n)
}
#----------------------------------------------
print.adapt = function(accept1z,jump1z){
  accept1=accept1z; jump1=jump1z; 
  
  for (k in 1:length(accept1)){
    z=accept1[[k]]/accept.output
    print(names(accept1)[k])
    print(mean(z)); print(mean(jump1[[k]]))
  }
  
  for (k in 1:length(jump1)){
    cond=(accept1[[k]]/accept.output)>0.4 & jump1[[k]]<10000
    jump1[[k]][cond] = jump1[[k]][cond]*2       
    cond=(accept1[[k]]/accept.output)<0.2 & jump1[[k]]>0.001
    jump1[[k]][cond] = jump1[[k]][cond]*0.5
    accept1[[k]][]=0
  }
  
  return(list(jump1=jump1,accept1=accept1))
}
#----------------------------
fix.MH=function(lo,hi,old1,new1,jump){
  jold=pnorm(hi,mean=old1,sd=jump)-pnorm(lo,mean=old1,sd=jump)
  jnew=pnorm(hi,mean=new1,sd=jump)-pnorm(lo,mean=new1,sd=jump)
  log(jold)-log(jnew) #add this to pnew
}
#----------------------------------------------------------------------------------------------
tnorm <- function(n,lo,hi,mu,sig){   #generates truncated normal variates based on cumulative normal distribution
  #normal truncated lo and hi
  
  if(length(lo) == 1 & length(mu) > 1)lo <- rep(lo,length(mu))
  if(length(hi) == 1 & length(mu) > 1)hi <- rep(hi,length(mu))
  
  q1 <- pnorm(lo,mu,sig) #cumulative distribution
  q2 <- pnorm(hi,mu,sig) #cumulative distribution
  
  z <- runif(n,q1,q2)
  z <- qnorm(z,mu,sig)
  z[z == -Inf]  <- lo[z == -Inf]
  z[z == Inf]   <- hi[z == Inf]
  z
}
#----------------------------------------------------------------------------------------------
acceptMH <- function(p0,p1,x0,x1,BLOCK){   #accept for M, M-H
  # if BLOCK, then accept as a block,
  # otherwise, accept individually
  
  nz           <- length(x0)  #no. to accept
  if(BLOCK) nz <- 1
  
  a    <- exp(p1 - p0)       #acceptance PR
  z    <- runif(nz,0,1)
  keep <- which(z < a)
  
  if(BLOCK & length(keep) > 0) x0 <- x1
  if(!BLOCK)                   x0[keep] <- x1[keep]           
  accept <- length(keep)        
  
  list(x = x0, accept = accept)
}
#-------------------------------
rmvnorm1=function (n, sigma, pre0.9_9994 = FALSE) 
{
  #   retval <- chol(sigma, pivot = TRUE)
  #   o <- order(attr(retval, "pivot"))
  #   retval <- retval[, o]
  s. <- svd(sigma)
  if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
    warning("sigma is numerically not positive definite")
  }
  R = t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
  retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
  retval
}
#----------------------------
fix.probs=function(probs){
  cond=probs<0.00001
  probs[cond]=0.00001
  cond=probs>0.99999
  probs[cond]=0.99999
  probs
}
#----------------------------
update.phi=function(param,jump){
  phi.orig=phi.old=param$phi
  proposed=matrix(tnorm(nspp*ncomm,lo=0,hi=1,mu=phi.old,sig=jump),ncomm,nspp)
  
  for (i in 1:ncomm){
    phi.new=phi.old
    phi.new[i,]=proposed[i,]
    adj=fix.MH(lo=0,hi=1,phi.old[i,],phi.new[i,],jump[i,])

    prob.old=fix.probs(param$theta%*%phi.old)
    prob.new=fix.probs(param$theta%*%phi.new)
    pold=colSums(dbinom(y,size=nmat,prob=prob.old,log=T))+dbeta(phi.old[i,],a.phi,b.phi,log=T)
    pnew=colSums(dbinom(y,size=nmat,prob=prob.new,log=T))+dbeta(phi.new[i,],a.phi,b.phi,log=T)
    k=acceptMH(pold,pnew+adj,phi.old[i,],phi.new[i,],F)
    phi.old[i,]=k$x
  }
  list(phi=phi.old,accept=phi.orig!=phi.old)
}
#-----------------------------
update.theta=function(param,jump){
  v.orig=v.old=param$vmat
  tmp=tnorm(nloc*(ncomm-1),lo=0,hi=1,mu=v.old[,-ncomm],sig=jump[,-ncomm])
  novos=cbind(matrix(tmp,nloc,ncomm-1),1)
  ajuste=matrix(fix.MH(lo=0,hi=1,v.old,novos,jump),nloc,ncomm)
  
  prior.old=matrix(dbeta(v.old,1,gamma,log=T),nloc,ncomm)
  prior.new=matrix(dbeta(novos,1,gamma,log=T),nloc,ncomm)
  
  for (j in 1:(ncomm-1)){ #last column has to be 1
    v.new=v.old
    v.new[,j]=novos[,j]
    
    theta.old=convertSBtoNormal(vmat=v.old,ncol=ncomm,nrow=nloc,prod=rep(1,nloc))
    theta.new=convertSBtoNormal(vmat=v.new,ncol=ncomm,nrow=nloc,prod=rep(1,nloc))
    
    pold=fix.probs(theta.old%*%param$phi)
    pnew=fix.probs(theta.new%*%param$phi)
    
    p1.old=rowSums(dbinom(y,size=nmat,prob=pold,log=T))
    p1.new=rowSums(dbinom(y,size=nmat,prob=pnew,log=T))
    
    k=acceptMH(p1.old+prior.old[,j],
               p1.new+prior.new[,j]+ajuste[,j],
               v.old[,j],v.new[,j],F)
    v.old[,j]=k$x
  }
  theta=convertSBtoNormal(vmat=v.old,ncol=ncomm,nrow=nloc,prod=rep(1,nloc))
  list(theta=theta,v=v.old,accept=v.old!=v.orig)
}