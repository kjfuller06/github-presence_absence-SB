x=seq(from=0.0001,to=0.9999,length.out=1000)
y=dbeta(x,0.8,0.8)
plot(x,y,type='l')