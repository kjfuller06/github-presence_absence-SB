# Generating Parameters

rm(list=ls(all=TRUE))

# Sets the starting number for random number generation so results will be the same every time they are run.
set.seed(4)
nloc=1000 #number of locations
nspp=200 #number of species
ncommun=5 #number of groups
base=floor(nloc/(ncommun-2))
nobs=5 #number of observation per location
#generate thetas-> the proportions of different groups occurring in each location
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
#generate phi's-> the probabilities of species belonging to each group
# Creates a matrix with 200 columns, 5 rows and NAs as values
phi=matrix(NA,ncommun,nspp)
# This generates random numbers from a Beta distribution (rbeta(number of observations, alpha, beta)). alpha and beta are parameters that determine the shape of the distribution. 0.5 for both yields a basically u-shaped density distribution (think normal distribution- this is actually a density distribution showing the distribution of values found in the dataset) with extreme ends as asymptotes- skewing probabilities toward the high and low ends. The outputs are added as the values in the 5 x 200 matrix-> these are the probabilities of a species (column) belonging to a group (row)
phi[]=rbeta(ncommun*nspp,0.5,0.5)
# This selects the first five species (1:5) and makes each sequential group (1:5) have 100% membership probability with the rest have 0%. That way at least one species will be totally different between the groups.
phi[,1:ncommun]=diag(1,ncommun) #to ensure that groups have distinct species composition
phi.true=phi
