

# set to zero all the negative value of a matrix
pos<-function(A) {
  A[A<0]<-0 
  return(A);
}

# Function that should be minimized
DLQPFunctionSVM<-function(xi,Ai,rho,ui,zi) {
  total<-sum(pos(Ai%*%xi+1))+rho/2*sum(t(xi-zi+ui)%*%(xi-zi+ui))
  return(total)
}

takeStep<-function(xi,zi,ui) {
  # Build Ai matrix multipling x samples for y values (last column).
  Ai<-samples;
  for(i in seq(1,numFeatures -1) ) {
    Ai[,i]<-Ai[,i]*Ai[,numFeatures]
  }
  # calculate the minimum
  res2<-optimx(xi,fn=DLQPFunctionSVM,Ai=Ai,rho=rho,ui=ui,zi=zi, itnmax = 3000)
  #    print(algorithm);
  retBack<-unlist(res2[1,1:numFeatures])
  return(retBack)
}