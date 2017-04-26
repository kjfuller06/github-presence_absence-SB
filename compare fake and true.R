boxplot(param$theta)

#order=black,red,green,blue,light blue
ind=1:5#c(1,2,5,6,10) #1:5#
theta=param$theta[,ind]
phi=param$phi[ind,]

ind=1:5#c(1,3,4,5,2)
theta=theta[,ind]
plot(NA,NA,ylim=c(0,1),xlim=c(0,nloc))
for (i in 1:ncomm) lines(theta[,i],col=i)

phi=phi[ind,]
rango=range(c(phi.true,phi))
plot(phi.true,phi,xlim=rango,ylim=rango)
lines(rango,rango)

rango=range(c(theta.true,param$theta))
plot(theta.true,param$theta[,ind],xlim=rango,ylim=rango)
lines(rango,rango)
