rm(list=ls(all=TRUE))
library(Rcpp)
set.seed(4)

setwd('U:\\presence absence model\\github-presence_absence-SB')
source('gibbs functions.R')
source('gibbs wrapper.R')
sourceCpp('aux1.cpp')

dat=read.csv('fake data y.csv',as.is=T)

