# Generating Parameters

rm(list=ls(all=TRUE))
set.seed(4)
nloc=1000 #number of locations
nspp=200 #number of species
ncommun=5 #number of groups
base=floor(nloc/(ncommun-2))
nobs=5 #number of observation per location
#generate thetas
x=seq(from=-1,to=1,length.out=base)
y=sqrt(1-(x^2))*0.1
min1=0.0001
y[y<min1]=min1
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

#visualize thetas
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1),xlab='Locations',ylab=expression(theta[lk]))
for (i in 1:ncommun) lines(1:nloc,theta[,i],col=i)

#generate the id for each location (loc.id)
loc.id=rep(1:nloc,each=nobs) #nobs observations for each location
#generate phi's
phi=matrix(NA,ncommun,nspp)
phi[]=rbeta(ncommun*nspp,0.5,0.5)
phi[,1:ncommun]=diag(1,ncommun) #to ensure that groups have distinct species composition
phi.true=phi