library(jpeg)
img = readJPEG("damon.jpg")
par(mar = c(0,0,0,0))
plot(0,0,type = "n")
rasterImage(img, -1,-1,1,1)

ratio=0.4
randomMissfun=function(X,ratio,replace){
  m=nrow(X)
  n=ncol(X)
  A=matrix(runif(m*n,0,1),ncol=n)
  ind=which(A<=ratio)
  vecX=c(X)
  vecX[ind]=replace
  X=matrix(vecX,ncol=n)
  return(X)
}

missImgPlot=randomMissfun(img,ratio,NA)
par(mar = c(0,0,0,0))
plot(0,0,type = "n")
rasterImage(missImgPlot, -1,-1,1,1)

missImg=randomMissfun(img,ratio,1.001)
delta=1
lambda=1
stepSize=0.5

time0=proc.time()
system.time(accProImageRestoration(missImg,lambda,delta,stepSize))
Outcome=accProImageRestoration(missImg,lambda,delta,stepSize)
proc.time()-time0
iteration=Outcome$iteration
reconsImg=Outcome$reconsImg
par(mar = c(0,0,0,0))
plot(0,0,type = "n")
rasterImage(reconsImg, -1,-1,1,1)

