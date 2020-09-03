# Simulating Data

#calculate probabilities after integrating out the z's
probs=theta%*%phi
# image(probs)
#generate data for each location
dat=matrix(NA,nloc*nobs,nspp)
oo=1
for (i in 1:nloc){
  for (j in 1:nobs){
    dat[oo,]=rbinom(nspp,size=1,prob=probs[i,])
    oo=oo+1
  }
}
colnames(dat)=paste('spp',1:nspp,sep='')
setwd('U:\\presence absence model\\appendix tutorial')
dat1=cbind(dat,loc.id)
write.csv(dat1,'fake data.csv',row.names=F)