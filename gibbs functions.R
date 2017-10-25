#----------------------------------------------------------------
get.theta=function(nlk,gamma,ncomm,nloc){
  vmat=matrix(NA,nloc,ncomm)
  for (i in 1:(ncomm-1)){
    if (i==(ncomm-1)) cumsoma=nlk[,ncomm]
    if (i< (ncomm-1)) cumsoma=rowSums(nlk[,(i+1):ncomm])
    vmat[,i]=rbeta(nloc,nlk[,i]+1,cumsoma+gamma)
  }
  vmat[,ncomm]=1
  convertVtoTheta(vmat,rep(1,nloc))
}
#------------------------------------------------------
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
    if (n[i]>1)  dat.tmp=colSums(dat.tmp)
    if (n[i]==1) dat.tmp=as.numeric(dat.tmp)
    res[i,]=dat.tmp
  }
  colnames(res)=colnames(dat1[,-ind])
  loc.id=unique(dat1$loc.id)
  list(dat=res,loc.id=loc.id,n=n)
}
