rm(list=ls(all=TRUE))
library(Rcpp)
set.seed(10)

setwd('U:\\presence absence model\\github-presence_absence-SB')
source('gibbs functions.R')
source('gibbs wrapper.R')
sourceCpp('aux1.cpp')

dat=read.csv('fake data y.csv',as.is=T)
ncomm=10
nloc=max(dat$loc.id)
nspp=ncol(dat)-1
ngibbs=1000; nadapt=ngibbs/2
a.phi=1; b.phi=1; gamma=0.1

res=run.gibbs(dat=dat,ngibbs=ngibbs,a.phi=a.phi,b.phi=b.phi,gamma=gamma,ncomm=ncomm,nadapt=nadapt)

