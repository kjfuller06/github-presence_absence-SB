strsplit(Sys.getenv("PATH"), split = ";")

# from comments in StackOverflow, it seems the library path I've been using includes a special character ($), which isn't supported by Rcpp. Which is dumb. I have now run RStudio as an administrator to access another library- the one in my C: drive which does not include any "$"'s. I will not install Rcpp into this library and hopefully it will work

install.packages("Rcpp", lib="C:/Program Files/R/R-4.0.2/library")
library(Rcpp)
find.package("Rcpp")
#looks good

evalCpp("2+2")
# YAYAYAYAYAY!!!

# Process did not need to be repeated. Works fine opening RStudio normally now.