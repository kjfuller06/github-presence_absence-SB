#order=black,red,green,blue,light blue
ind=1:5#c(1,2,5,6,10) #1:5#
theta=matrix(res$theta[nrow(res$theta),],nloc,ncomm)
theta1=theta[,ind]
phi=matrix(res$phi[nrow(res$phi),],ncomm,nspp)
phi1=phi[ind,]

boxplot(theta1)
plot(NA,NA,ylim=c(0,1),xlim=c(0,nloc))
for (i in 1:ncomm) lines(theta1[,i],col=i)

rango=range(c(phi.true,phi1))
plot(phi.true,phi1,xlim=rango,ylim=rango)
lines(rango,rango)

rango=range(c(theta.true,theta1))
plot(theta.true,theta1[,ind],xlim=rango,ylim=rango)
lines(rango,rango)
