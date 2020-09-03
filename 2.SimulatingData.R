# Simulating Data

#calculate probabilities after integrating out the z's
# This multiplies theta by phi
probs=theta%*%phi
# image(probs)
#generate data for each location
# Creates a matrix with 5000 rows and 200 columns, with NA as values
dat=matrix(NA,nloc*nobs,nspp)
oo=1
# This assigns the values in dat. rbinom generates a random 0 or 1 (rbinom(number of observations, number of trials (adds the results together for multiple trials), probablity of success)). Here, it applies the probability of success to all observations at that site using the location and species probability combinations in "probs".

# For each location
for (i in 1:nloc){
  # For each observation within the selected location
  for (j in 1:nobs){
    # Generate a presence or absence using the probability of 1) the group being present and 2) the species belonging to that group
    dat[oo,]=rbinom(nspp,size=1,prob=probs[i,])
    oo=oo+1
  }
}
colnames(dat)=paste('spp',1:nspp,sep='')

# setwd('U:\\presence absence model\\appendix tutorial')

# Generates the final matrix with 201 columns and 5000 rows, with binary values. The columns are the species IDs and one location ID; the rows are unique samplings
dat1=cbind(dat,loc.id)
write.csv(dat1,'fake data.csv',row.names=F)
