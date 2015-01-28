


run<-function(epsilon = 0.001, MAX_TIME = 2000, runningPlot = FALSE) {
  logRun<-list();
  history.r_norm<-list();
  history.s_norm<-list();
  history.eps_pri<-list();
  history.eps_dual<-list();
  stopConditionAt<-list();
  XconvRun<-list();
  boydVAL<-list()
  
  timer<<-0;
  x = matrix(1,nrow=numNodes,ncol=numFeatures);
  x_hat = matrix(1,nrow=numNodes,ncol=numFeatures);
  z = array(0,numFeatures);
  u = matrix(0,nrow=numNodes,ncol=numFeatures);
  par(mfrow = c(2,1)) # two rows in the plot
  # loop until the MAX_TIME has been reached
  for(timer in seq(1,MAX_TIME)) {
    
    i<-1;  	
    # collect x vector from each node
    for(i in seq(1,numNodes)) {
      x[i,]<-nodes[[i]]$takeStep(x[i,],z,u[i,]);	
    }
    
    # now calcolate z and u at next step
    i<-1;
    for(i in seq(1,numNodes)) {
      x_hat[i,]<-alpha*x[i,]+(1-alpha)*z
    }  
    zold<-z;
    z<-(( numNodes*rho ) / ( (1/lambda)+numNodes*rho ))*colMeans(x_hat+u) 
    
    zMatrice<- matrix(rep(z,numNodes),nrow=dim(x_hat)[1],byrow=TRUE)
    u<-u+ (x_hat - zMatrice);
    
    # now build some logs and the variables used to check the 
    # stop condition
    logRun[[timer]]<-x;
    history.r_norm[[timer]]<-mtlbNorm(x-zMatrice)
    history.s_norm[[timer]]<-mtlbNorm(-rho*(matrix(rep(z,numNodes),nrow=dim(x_hat)[1],byrow=TRUE) - matrix(rep(zold,numNodes),nrow=dim(x_hat)[1],byrow=TRUE)));
    history.eps_pri[[timer]]<-sqrt(dim(x)[2])*ABSTOL + RELTOL*max(mtlbNorm(x), mtlbNorm(-matrix(rep(z,numNodes),nrow=dim(x_hat)[1],byrow=TRUE)));
    history.eps_dual[[timer]]<-sqrt(dim(x)[2])*ABSTOL + RELTOL*mtlbNorm(rho*u);     
    
    # is the stop condition reached?
    # if so, bye!
    boydVAL[[timer]]<-list();
    boydVAL[[timer]][["got"]] <- 0
    boydVAL[[timer]][["history.eps_pri"]] <- history.eps_pri[[timer]]
    boydVAL[[timer]][["history.r_norm"]] <- history.r_norm[[timer]]
    boydVAL[[timer]][["history.s_norm"]] <- history.s_norm[[timer]]
    boydVAL[[timer]][["history.eps_dual"]] <- history.eps_dual[[timer]]
    
    if ( boydVAL[[timer]][["got"]] == 1 ) boydVAL[[timer]][["got"]] = 1

    # increase the timer    
    lunghezza<-length(logRun)
    if(lunghezza>2) {
      improvement<-mean(colMeans((logRun[[  lunghezza ]]-logRun[[  lunghezza-1 ]])/logRun[[  lunghezza ]]))
      XconvRun[[timer]]<-improvement
    }
    #      SVMimprovement<-mean((colMeans(logRun[[ lunghezza ]])-SVMResult)/SVMResult)
    if(history.r_norm[[timer]] < history.eps_pri[[timer]] & history.s_norm[[timer]] < history.eps_dual[[timer]] ) {
      boydVAL[[timer]][["got"]] <- 1
      if ( boydVAL[[timer]][["got"]] == 1 ) boydVAL[[timer]][["got"]] = 1
      return(list("XconvRun"=XconvRun,"x"=x,"logRun"=logRun,"history.eps_dual"=unlist(history.eps_dual),"history.eps_pri"=unlist(history.eps_pri),"history.s_norm"=unlist(history.s_norm),"history.r_norm"=unlist(history.r_norm),"stopConditionAt"=stopConditionAt))
    }
    
    timer<-timer+1;
  }
  
  # Well, if I'm here is because I reached the MAX_TIME without a converge
  stopConditionAt[["timer"]]<-NA;
  stopConditionAt[["x"]]<-NA;
  return(list("XconvRun"=XconvRun,"x"=x,"logRun"=logRun,"history.eps_dual"=unlist(history.eps_dual),"history.eps_pri"=unlist(history.eps_pri),"history.s_norm"=unlist(history.s_norm),"history.r_norm"=unlist(history.r_norm),"stopConditionAt"=stopConditionAt))    
}
