rm(list=ls(all=TRUE))
library(Rcpp)
library(data.table)
set.seed(4)

setwd('U:\\presence absence model\\github-presence_absence-SB')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

dat1=read.csv('fake data y.csv',as.is=T)
# setwd('U:\\presence absence model\\bbs data\\data florida')
# dat1=read.csv('fifty edited FL.csv',as.is=T)

dat1=dat1[order(dat1$loc.id),]
ind=which(colnames(dat1)%in%c('loc.id'))
y=data.matrix(dat1[,-ind])
loc.id=dat1$loc.id
nspp=ncol(y)
nloc=length(unique(loc.id))
ncomm=5

#hyper-priors
a.phi=0.01
b.phi=0.99
gamma=0.1

#initial values
theta=matrix(1/ncomm,nloc,ncomm)
vmat=theta
vmat[,ncomm]=1
phi=matrix(colMeans(y),ncomm,nspp,byrow=T)
phi[phi==1]=0.99999

#gibbs stuff
ngibbs=10000
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
  tmp=update.phi(param,jump1$phi)
  param$phi=tmp$phi #phi.true #
  accept1$phi=accept1$phi+tmp$accept
  
  tmp=update.theta(param,jump1$vmat)
  param$theta=tmp$theta
  param$vmat=tmp$v
  accept1$vmat=accept1$vmat+tmp$accept
  
  if (i%%accept.output==0 & i<1000){
    k=print.adapt(accept1,jump1)
    accept1=k$accept1
    jump1=k$jump1
  }
  
  #to assess convergence, examine logl
  prob=fix.probs(param$theta%*%param$phi)[loc.id,]
  loglikel=sum(dbinom(y,size=1,prob=prob,log=T))+
           sum(dbeta(param$phi,a.phi,b.phi,log=T))+
           sum(dbeta(param$vmat[,-ncomm],1,gamma,log=T))
  
  vec.logl[i]=loglikel
  vec.theta[i,]=param$theta
  vec.phi[i,]=param$phi
}

seq1=8000:ngibbs
plot(vec.logl[seq1],type='l')

setwd('U:\\presence absence model\\bbs FL results')
theta=matrix(colMeans(vec.theta[seq1,]),nloc,ncomm)
colnames(theta)=paste('comm',1:10,sep='')
write.csv(theta,'bbs FL theta.csv',row.names=F)

phi=matrix(colMeans(vec.phi[seq1,]),ncomm,nspp)
colnames(phi)=colnames(y)
write.csv(phi,'bbs FL phi.csv',row.names=F)
write.csv(vec.logl,'bbs FL logl.csv',row.names=F)