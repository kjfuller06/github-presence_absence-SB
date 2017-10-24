rm(list=ls(all=TRUE))
set.seed(5)

nloc=1000
nobs=3
nspp=200
y=matrix(NA,nloc,nspp)

loc.id=rep(1:nloc,each=nobs)
ncommun=10

#generate thetas
base=floor(nloc/(ncommun-4))

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
phi=matrix(rbeta(ncommun*nspp,0.8,0.8),ncommun,nspp)
phi[,1:(ncommun*2)]=cbind(diag(1,ncommun),diag(1,ncommun)) #at least 2 species that are unique to each community
phi.true=phi
apply(phi>0.8,1,sum)
table(apply(phi>0.8,2,sum))

image(phi)

#generate data
probs=theta%*%phi
cond=probs<0.00001; sum(cond)
probs[cond]=0.00001
cond=probs>0.99999; sum(cond)
probs[cond]=0.99999

par(mfrow=c(1,1))
hist(probs)
image(probs)
probs1=probs[loc.id,]

#number of observations per location
obs=matrix(rbinom(nrow(probs1)*ncol(probs1),size=1,prob=probs1),nrow(probs1),ncol(probs1))
colnames(obs)=paste('spp',1:nspp,sep='')
obs1=cbind(obs,loc.id)
colnames(obs1)[ncol(obs1)]='loc.id'

setwd('U:\\presence absence model\\github-presence_absence-SB')
write.csv(obs1,'fake data y.csv',row.names=F)
