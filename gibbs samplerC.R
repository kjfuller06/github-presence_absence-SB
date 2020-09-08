# rm(list=ls(all=TRUE))
library('Rcpp')
set.seed(4)

#get functions
# setwd('U:\\presence absence model\\github-presence_absence-SB')
source('gibbs functions.R')
sourceCpp('aux1.cpp')

#get data
# setwd('U:\\presence absence model\\simulations')
dat=read.csv('fake data.csv',as.is=T)
tmp=aggregate.data(dat)
y=tmp$dat
n=tmp$n

#useful stuff
ind=which(colnames(dat)=='loc.id')
dat1=data.matrix(dat[,-ind])
loc.id=dat$loc.id
nloc=max(loc.id)
nspp=ncol(dat1)
ncomm=10
nlinhas=nrow(dat1)
hi=0.999999
lo=0.000001

#get good initial values using k-means clustering
# z=kmeans(y/n,ncomm,nstart=20)
# theta=matrix(NA,nloc,ncomm)
# seq1=1:ncomm
# for (i in 1:ncomm){
#   cond=z$cluster==i
#   theta[cond,i]=0.8
#   seq2=seq1[-i]
#   theta[cond,seq2]=0.2/9
# }
# 
# phi=data.matrix(z$centers)
# cond=phi==0; phi[cond]=0.01
# cond=phi==1; phi[cond]=0.99
theta=matrix(1/ncomm,nloc,ncomm)
phi=matrix(0.5,ncomm,nspp)
z=matrix(sample(1:ncomm,size=nlinhas*nspp,replace=T),nlinhas,nspp)

#priors
a.phi=1
b.phi=1
gamma=0.1

#gibbs details
ngibbs=8000
theta.out=matrix(NA,ngibbs,ncomm*nloc)
phi.out=matrix(NA,ngibbs,ncomm*nspp)
llk=rep(NA,ngibbs)
options(warn=2)
for (i in 1:ngibbs){
  # print(i)   
  
  #sample z
  rand.u=matrix(runif(nlinhas*nspp),nlinhas,nspp)
  z=samplez(log(theta), log(1-phi), log(phi), dat1, loc.id,rand.u, ncomm, nloc)

  #calculate summaries
  tmp=getks(z=z, ncommun=ncomm, dat=dat1)
  nks1=tmp$nks1
  nks0=tmp$nks0
  nlk=getlk(z=z,locid=loc.id, ncommun=ncomm, nloc=nloc)
  
  #get parameters  
  theta=get.theta(nlk,gamma,ncomm,nloc) #theta.true#
  theta[theta>hi]=hi; theta[theta<lo]=lo
  phi=matrix(rbeta(nspp*ncomm,nks1+a.phi,nks0+b.phi),ncomm,nspp) #phi.true#
  phi[phi>hi]=hi; phi[phi<lo]=lo
  prob=theta%*%phi
  prob[prob>hi]=hi; prob[prob<lo]=lo
  prob1=prob[loc.id,]

  #store results  
  llk[i]=sum(dat1*log(prob1)+(1-dat1)*log(1-prob1))
  theta.out[i,]=theta
  phi.out[i,]=phi
}

plot(llk,type='l',ylim=range(llk,na.rm=T))
