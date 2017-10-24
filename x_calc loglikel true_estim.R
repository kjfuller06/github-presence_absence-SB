#CALCULATE PROBABILITY FOR TRUE PARAMETERS

#get vmat from theta.true
vmat=matrix(NA,nloc,ncomm-1)
vmat[,1]=theta.true[,1]
vmat[,2]=theta.true[,2]/(1-vmat[,1])
for (i in 3:(ncomm-1)) vmat[,i]=theta.true[,i]/apply(1-vmat[,1:(i-1)],1,prod)

#check to make sure this calculation is correct
theta.tmp=convertSBtoNormal(vmat=cbind(vmat,1),ncol=ncomm,nrow=nloc,prod=rep(1,nloc))
hist(theta.true-theta.tmp)

#get probability
prob=get.logl(theta=theta.true,phi=phi.true,y=y,nmat=nmat)
loglikel=sum(prob)+
  sum(dbeta(phi.true,a.phi,b.phi,log=T))+
  sum(dbeta(vmat,1,gamma,log=T))
loglikel 
#-149,157 + 1125 + 0 = -148,031 (when gamma=1)
#-186,424
#------------------------------------------------------------
#CALCULATE PROBABILITY FOR ESTIMATED PARAMETERS
prob=get.logl(theta=param$theta,phi=param$phi,y=y,nmat=nmat)
loglikel=sum(prob)+
  sum(dbeta(param$phi,a.phi,b.phi,log=T))+
  sum(dbeta(param$vmat[,-ncomm],1,gamma,log=T))
loglikel 
# -184,739 + 4,678 + 0= -180,060 (when gamma=1)
#-219,521

#--------------------------------------------------------------
#Is this an artifact of truncation? No
prob.estim=param$theta%*%param$phi
image(prob.estim,main='estim')
sum(dbinom(y,size=nmat,prob=prob.estim,log=T)) #-126,410

prob.true=theta.true%*%phi.true
image(prob.true,main='true')
sum(dbinom(y,size=nmat,prob=prob.true,log=T)) #-90,828

# filled.contour(abs(prob.estim-prob.true))
image(y/nmat,main='obs')

plot(prob.true,jitter(y/nmat))
table((y/nmat)[50,])

plot(prob.true[50,],jitter(y/nmat)[50,])