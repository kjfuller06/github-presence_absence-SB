# Generating Parameters

rm(list=ls(all=TRUE))

# Sets the starting number for random number generation so results will be the same every time they are run.
set.seed(4)
nloc=1000 #number of locations
nspp=200 #number of species
ncommun=5 #number of groups
# base is 333, why do we want this? floor() calculates the closest possible integer that is lower than the number(s) provided
base=floor(nloc/(ncommun-2))
nobs=5 #number of observation per location
#generate thetas-> the proportions of different groups occurring in each location
# Generates 333 numbers at regular intervals between -1 and 1, why 333?
x=seq(from=-1,to=1,length.out=base)
# Generates a nice, pretty sequence shaped like a rainbow with 1 as the max, zero at the extremes and 333 values. This could help create blending of membership across the locations
y=sqrt(1-(x^2))*0.1
min1=0.0001
# Ensures no values are below 0.0001, which they weren't anyway but that could change with altered numbers of groups, species and locations
y[y<min1]=min1
# init is equal to 200 
init=floor(nloc/ncommun)
# Generates a sequence: 1, 201, 401, 601, 801, 1000
seq1=c(seq(from=1,to=nloc,by=init),nloc)
# Generates a matrix with 5 columns and 1000 rows, with values of 0.0001
theta=matrix(min1,nloc,ncommun)
# Right! So "base" is used to select locations at regular intervals across the 1000 locations. If it were set to 200, there would be no overlap between groups.
# This 1) generates the selection numbers with a bit of overlap 1:333, 201:533, 401:733, 601:933 and 801:1133; 2) truncates the values by the max number of locations (1000) and 3) selects the locations, then assigns the pretty dome shape of values to those locations based on the selected group (group 1 dome is assigned to locations 1:333; group 2 dome is assigned to locations 201:533; group 3 dome is assigned to locations 401:733; etc)
for (i in 1:ncommun){
  seq2=seq1[i]:(seq1[i]+base-1)
  seq3=seq2[seq2<=nloc]
  theta[seq3,i]=y[1:length(seq3)]
}
# This divides the entire theta matrix by the sum of each row- creating a proportion and making the occurrence of different groups behave in relation to the occurrence of other groups at each location, instead of being calculated independently
theta=theta/matrix(apply(theta,1,sum),nloc,ncommun)
theta.true=theta

#visualize thetas
plot(NA,NA,xlim=c(0,nloc),ylim=c(0,1),xlab='Locations',ylab=expression(theta[lk]))
for (i in 1:ncommun) lines(1:nloc,theta[,i],col=i)

#generate the id for each location (loc.id)
# Creates a list of 5000- the locations listed 5 times each
loc.id=rep(1:nloc,each=nobs) #nobs observations for each location

#generate phi's-> the probabilities of species belonging to each group
# Creates a matrix with 200 columns, 5 rows and NAs as values
phi=matrix(NA,ncommun,nspp)
# This generates random numbers from a Beta distribution (rbeta(number of observations, alpha, beta)). alpha and beta are parameters that determine the shape of the distribution. 0.5 for both yields a basically u-shaped density distribution (think normal distribution- this is actually a density distribution showing the distribution of values found in the dataset) with extreme ends as asymptotes- skewing probabilities toward the high and low ends. The outputs are added as the values in the 5 x 200 matrix-> these are the probabilities of a species (column) belonging to a group (row)
phi[]=rbeta(ncommun*nspp,0.5,0.5)
# This selects the first five species (1:5) and makes each sequential group (1:5) have 100% membership probability with the rest have 0%. That way at least one species will be totally different between the groups.
phi[,1:ncommun]=diag(1,ncommun) #to ensure that groups have distinct species composition
phi.true=phi
