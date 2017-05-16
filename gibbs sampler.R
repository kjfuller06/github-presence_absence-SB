rm(list=ls(all=TRUE))
library(Rcpp)
set.seed(4)

setwd('U:\\presence absence model\\github-presence_absence-SB')
source('gibbs functions.R')
source('gibbs wrapper.R')
sourceCpp('aux1.cpp')

dat=read.csv('fake data y.csv',as.is=T)
nloc=max(dat$loc.id)
nspp=ncol(dat)-1
ncomm=5
res=run.gibbs(dat=dat,ngibbs=1000,a.phi=0.01,b.phi=0.99,gamma=0.1,
              ncomm=ncomm,nadapt=1000)

