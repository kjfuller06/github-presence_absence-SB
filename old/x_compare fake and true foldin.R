rm(list=ls(all=TRUE))
set.seed(4)
library(Rcpp)

# setwd('U:\\presence absence model\\bbs data\\data')
# dat=read.csv('foldin data.csv',as.is=T)

setwd('U:\\presence absence model\\github-presence_absence-SB')
source('foldin wrapper.R')
source('gibbs functions.R')
sourceCpp('aux1.cpp')
phi=read.csv('fake phi.csv',as.is=T)
dat=read.csv('fake data y.csv',as.is=T)
nloc=max(dat$loc.id)
ncomm=nrow(phi)

ngibbs=1000
res=run.foldin(dat=dat,ngibbs=ngibbs,gamma=0.1,phi=data.matrix(phi),nadapt=1000)

seq1=900:ngibbs
plot(res$logl[seq1],type='l')

theta=matrix(colMeans(res$theta[seq1,]),nloc,ncomm)
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:10) lines(theta[,i],col=i)