rm(list=ls(all=TRUE))
set.seed(4)

setwd('U:\\gordon\\presence_absence\\ebird\\data')
dat1=read.csv('edited data.csv',as.is=T)
dat1=dat1[order(dat1$loc.id),]
ind=which(colnames(dat1)%in%c('loc.id'))
y=data.matrix(dat1[,-ind])

loc.id=dat1$loc.id
nspp=ncol(y)
nloc=length(unique(loc.id))
ncommun=5

#generate thetas
base=floor(nloc/(ncommun-2))

x=seq(from=-1,to=1,length.out=base)
y=sqrt(1-(x^2))*0.1
min1=0.0001
y[y<min1]=min1
# plot(x,y)

init=floor(nloc/ncommun)
seq1=c(seq(from=1,to=nloc,by=init),nloc)

theta=matrix(min1,nloc,ncommun)
for (i in 1:ncommun){
  seq2=seq1[i]:(seq1[i]+base-1)
  seq3=seq2[seq2<=nloc]
  theta[seq3,i]=y[1:length(seq3)]
}
theta=theta/matrix(apply(theta,1,sum),nloc,ncommun)
theta.true=theta

plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1))
for (i in 1:ncommun) lines(1:nloc,theta[,i],col=i)

#generate phi's
makesparse=matrix(rbinom(ncommun*nspp,size=1,prob=0.3),ncommun,nspp)
phi.true=phi=matrix(runif(ncommun*nspp,min=0.5,max=1),ncommun,nspp)*makesparse

par(mfrow=c(3,2),mar=rep(1,4))
for (i in 1:ncommun) plot(phi[i,],type='h',ylim=c(0,1))

#generate data
probs=theta%*%phi
hist(probs)
probs1=probs[loc.id,]

#number of observations per location
obs=matrix(rbinom(nrow(probs1)*ncol(probs1),size=1,prob=probs1),nrow(probs1),ncol(probs1))
ind=which(colnames(dat1)%in%c('loc.id'))
colnames(obs)=colnames(dat1)[-ind]
obs1=cbind(obs,loc.id)
colnames(obs1)[ncol(obs1)]='loc.id'

setwd('U:\\gordon\\presence_absence\\new model')
write.csv(obs1,'fake data y.csv',row.names=F)
