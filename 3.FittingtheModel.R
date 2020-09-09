# Fitting the model
library('Rcpp')
set.seed(4)
#get functions
# setwd('U:\\GIT_models\\github-presence_absence-SB')
source('gibbs functions.R')
source('gibbs sampler main function.R')
sourceCpp('aux1.cpp')
#run Gibbs sampler
results=lda.presence.absence(dat=dat,id=loc.id,ncomm=10,
                             a.phi=1,b.phi=1,gamma=0.1,ngibbs=2000)

seq1=1000:2000
par(mfrow=c(2,1),mar=c(3,3,1,1))
plot(results$llk,type='l',xlab='Iterations',ylab='Log-likel.')
plot(x=seq1,y=results$llk[seq1],type='l',xlab='Iterations',ylab='Log-likel.')

#summarize
theta=matrix(colMeans(results$theta[seq1,]),nloc,10)
phi=matrix(colMeans(results$phi[seq1,]),10,nspp)

par(mfrow=c(3,1),mar=c(4,4,4,1))
boxplot(theta,xlab='',ylab=expression(theta[lk]),names=paste('G',1:10,sep=''),
        xlab='Groups')
plot(NA,NA,xlim=c(0,1000),ylim=c(0,1),xlab='Locations',ylab=expression(theta[lk]))
for (i in 1:ncol(theta.true)) lines(theta.true[,i],col='grey',lwd=2)
for (i in 1:10) lines(theta[,i],col=i)
#re-order so that the estimated groups match the true groups
order1=c(4,3,1,2,5)
z=c(0,1)
plot(phi.true,phi[order1,],xlim=z,ylim=z,main=expression(phi[ks]),
     xlab='True',ylab='Estimated')
lines(z,z,col='red',lwd=2)