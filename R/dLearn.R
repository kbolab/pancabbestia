
# ------------------------------------------
# Class dLearn
# ------------------------------------------
dLearn<-function(incomingPathName,lambda=0,rho=0,alpha=0,token = 0, writeTo="screen",basePath = basePath){
  
  # ------------------------------------------
  # Attributes
  # ------------------------------------------
  
  numNodes<-0;
  numTotSamples<-0;
  samples<-list();
  nodes<-list();
  timer<-0;
  MAX_TIME<-0;
  ABSTOL<-0;
  RELTOL<-0;
  lambda<-lambda;
  alpha<-alpha;
  rho<-rho;
  numFeatures<-0;
  image01<-"";
  image02<-"";
  image03<-"";
  image04<-"";
  
  
  # ------------------------------------------
  # Methods
  # ------------------------------------------
  
  # addSamples: add a sample set to node (center)
  addSamples<-function(idCentre,arrSample) {
    nodes[[idCentre]]$addSamples(arrSample);
    numFeatures<<-length(arrSample[1,]);
  }
  
  # addNode: add a new node
  addNode<-function(nodeName) {
    nodes[[length(nodes)+1]]<<-dNode(p_externalName = nodeName,p_externalID = length(nodes)+1, alpha = alpha, rho = rho, lambda = lambda);	
    numNodes<<-length(nodes);
    return(length(nodes));
  }
  
  # getAttributeFromLearner: get an attribute from the class dLearn
  getAttributeFromLearner<-function(whichAttribute) {
    if( whichAttribute == 'nodes' )   return(names(nodes));
  }
  
  # getAttributeFromNode: get an attribute from one of the defined 
  # node. The access to the nodes is mediated by the class dLearn
  # in order to keep the node safe.
  getAttributeFromNode<-function(idNode,whichAttribute) {
    return(nodes[[idNode]]$getAttribute(whichAttribute))
  }
  
  # mtlbNorm: is a particular norm used in the calculus.
  mtlbNorm<-function(matrice) {
    return(max(svd(matrice)$d));
  }
  
  # run: Execute the calculus 
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
#      cat(".")
      boydVAL[[timer]]<-list();
      boydVAL[[timer]][["got"]] <- 0
      if(history.r_norm[[timer]] < history.eps_pri[[timer]] & history.s_norm[[timer]] < history.eps_dual[[timer]] ) {
#        stopConditionAt[["timer"]]<-timer;
#        stopConditionAt[["x"]]<-x;        
#        return(list("x"=x,"logRun"=logRun,"history.eps_dual"=unlist(history.eps_dual),"history.eps_pri"=unlist(history.eps_pri),"history.s_norm"=unlist(history.s_norm),"history.r_norm"=unlist(history.r_norm),"stopConditionAt"=stopConditionAt))
        boydVAL[[timer]][["got"]] <- 1
        
      }

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
#        if( abs(improvement) < abs(epsilon)) { return(list("convRun"=convRun,"x"=x,"logRun"=logRun,"history.eps_dual"=unlist(history.eps_dual),"history.eps_pri"=unlist(history.eps_pri),"history.s_norm"=unlist(history.s_norm),"history.r_norm"=unlist(history.r_norm),"stopConditionAt"=stopConditionAt))}
      }

#      SVMimprovement<-mean((colMeans(logRun[[ lunghezza ]])-SVMResult)/SVMResult)
#      convRun[[timer]]<-SVMimprovement;
#      if( abs(SVMimprovement) < abs(epsilon) ) { return(list("XconvRun"=XconvRun,"x"=x,"logRun"=logRun,"history.eps_dual"=unlist(history.eps_dual),"history.eps_pri"=unlist(history.eps_pri),"history.s_norm"=unlist(history.s_norm),"history.r_norm"=unlist(history.r_norm),"stopConditionAt"=stopConditionAt))}
#      cat(timer,",",SVMimprovement,"LOG",colMeans(logRun[[ lunghezza ]]),"SVM",SVMResult ,"BoydCONV",boydVAL[[timer]]$got,",",boydVAL[[timer]]$history.eps_pri,",",boydVAL[[timer]]$history.r_norm,",",boydVAL[[timer]]$history.s_norm,",",boydVAL[[timer]]$history.eps_dual,"\n")
      cat(timer,",",colMeans(logRun[[ lunghezza ]]) ,"BoydCONV",boydVAL[[timer]]$got,",",boydVAL[[timer]]$history.r_norm,",",boydVAL[[timer]]$history.eps_pri,",",boydVAL[[timer]]$history.s_norm,",",boydVAL[[timer]]$history.eps_dual,"\n")
      # routine for plotting the ongoing graph
      if ( runningPlot == TRUE ) {
        if (timer ==5 ) {  # aggiugner && controllo plot
          
          min.primary <- min(sapply(X = boydVAL, FUN = function(x) return(x$history.eps_pri), simplify = T))
          min.primary <- min.primary / 20
          min.secondary <- min(sapply(X = boydVAL, FUN = function(x) return(x$history.eps_dual), simplify = T))
          min.secondary <- min.secondary / 20
          max.primary <- max(sapply(X = boydVAL, FUN = function(x) return(x$history.r_norm), simplify = T))
          max.primary <- max.primary * 5
          max.secondary <- max(sapply(X = boydVAL, FUN = function(x) return(x$history.s_norm), simplify = T))
          max.secondary <- max.secondary * 5
          if( writeTo == "file") png( image01 )
          plot(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.r_norm), simplify = T), type = "l", col="red", lwd = 2, log="y", ylab = "Primary conv", ylim = c(min.primary, max.primary), xlim=c(0,20), xlab = "Iteration No.")
          lines(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.eps_pri), simplify = T), lty = 4, col="red", lwd = 2)
          if( writeTo == "file") dev.off()
          if( writeTo == "file") png( image02 )
          plot(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.s_norm), simplify = T), type = "l", col="blue", lwd = 2, log="y", ylab = "Secondary conv", ylim = c(min.secondary, max.secondary), xlim=c(0,20), xlab = "Iteration No.")
          lines(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.eps_dual), simplify = T), lty = 4, col="blue", lwd = 2)
          if( writeTo == "file") dev.off()
        }
        if (timer > 5) { # aggiunger && controllo plot
          if( writeTo == "file") png( image01) 
          plot(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.r_norm), simplify = T), type = "l", col="red", lwd = 2, log="y", ylab = "Primary conv", ylim = c(min.primary, max.primary), xlim=c(0,ceiling(timer/20)*20), xlab = "Iteration No.")
          lines(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.eps_pri), simplify = T), lty = 4, col="red", lwd = 2)
          if( writeTo == "file") dev.off()
          if( writeTo == "file") png( image02 )
          plot(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.s_norm), simplify = T), type = "l", col="blue", lwd = 2, log="y", ylab = "Secondary conv", ylim = c(min.secondary, max.secondary), xlim=c(0,ceiling(timer/20)*20), xlab = "Iteration No.")
          lines(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.eps_dual), simplify = T), lty = 4, col="blue", lwd = 2)
          if( writeTo == "file") dev.off()
        }
      }
      if(history.r_norm[[timer]] < history.eps_pri[[timer]] & history.s_norm[[timer]] < history.eps_dual[[timer]] ) {
        if( runningPlot== TRUE  ) {
          if( writeTo == "file") png( image01 )
          plot(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.r_norm), simplify = T), type = "l", col="red", lwd = 2, log="y", ylab = "Primary conv", ylim = c(min.primary, max.primary), xlim=c(0,ceiling(timer/20)*20), xlab = "Iteration No.")
          lines(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.eps_pri), simplify = T), lty = 4, col="red", lwd = 2)
          points(x = timer, y = boydVAL[[timer]]$history.r_norm, pch = 13, cex = 1.5, col = "red", lwd = 2)
          title(main = paste("Algorithm successfully converged after", timer, "iterations!"))
          if( writeTo == "file") dev.off()
          if( writeTo == "file") png( image02 )
          plot(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.s_norm), simplify = T), type = "l", col="blue", lwd = 2, log="y", ylab = "Secondary conv", ylim = c(min.secondary, max.secondary), xlim=c(0,ceiling(timer/20)*20), xlab = "Iteration No.")
          lines(x = c(1:timer), y = sapply(X = boydVAL, FUN = function(x) return(x$history.eps_dual), simplify = T), lty = 4, col="blue", lwd = 2)
          points(x = timer, y = boydVAL[[timer]]$history.s_norm, pch = 13, cex = 1.5, col = "blue", lwd = 2)
          if( writeTo == "file") dev.off()
          Sys.sleep(1)
        }
        return(list("XconvRun"=XconvRun,"x"=x,"logRun"=logRun,"history.eps_dual"=unlist(history.eps_dual),"history.eps_pri"=unlist(history.eps_pri),"history.s_norm"=unlist(history.s_norm),"history.r_norm"=unlist(history.r_norm),"stopConditionAt"=stopConditionAt))
      }

      timer<-timer+1;
    }
    
    # Well, if I'm here is because I reached the MAX_TIME without a converge
    stopConditionAt[["timer"]]<-NA;
    stopConditionAt[["x"]]<-NA;
    return(list("XconvRun"=XconvRun,"x"=x,"logRun"=logRun,"history.eps_dual"=unlist(history.eps_dual),"history.eps_pri"=unlist(history.eps_pri),"history.s_norm"=unlist(history.s_norm),"history.r_norm"=unlist(history.r_norm),"stopConditionAt"=stopConditionAt))    
  }
  
  # createPopulation: Create a pupulation
  # quantiXCentro = number of numSamples cases for each node
  # numFeatures = number of features
  # numCentri = number of nodes
  createPopulation<-function(quantiXCentro,numFeatures,numCentri,meanP = 3,devP = 1.5, meanN=-3, devN=1.5, deltaSDAmongCentroidsAmongCenters = 0) {
#    meanP<-3      # mean of positives
#    devP<-2     # st.dev of positives   (1.5)
#    meanN<--3     # mean of negatives 
#    devN<-2      # st.dev of negatives (1)
    
    Ai<-list()
    arrSD<-c()
    for(i in seq(1:numCentri)) {
      
      uTMP<-list()
      uTMP<-createSubset(
          numSamples = quantiXCentro,
          numFeatures = numFeatures,
          meanP = meanP,
          devP = devP,
          meanN = meanN,
          devN = devN,
          deltaSDAmongCentroidsAmongCenters = deltaSDAmongCentroidsAmongCenters)
      Ai[[i]]<-uTMP$A;      
      arrSD<-c(arrSD,uTMP$sdData)
    }
    return( list( "Ai" = Ai, "mediaSD" = mean(arrSD) ) )  
  }
  
  # createSubset: Create the sub-population of a specific node
  # numSamples = number of clinical cases for each node
  # numFeatures = number of features
  # meanP = mean of positives  
  # devP = st.dev of positives 
  # meanN = mean of negatives  
  # devN = st.dev of negatives 
  createSubset<-function(numSamples,numFeatures,meanP,devP,meanN,devN,deltaSDAmongCentroidsAmongCenters = 0) { 
    
    sbilanciamentoMediaP<-c()
    sbilanciamentoMediaN<-c()
    
    # Positive Samples
    PPoints<-c()

    for( i in seq(1:numFeatures)) {
      sbilanciamentoMediaP <- c(sbilanciamentoMediaP,abs(rnorm( n = 1, mean = 0, sd = deltaSDAmongCentroidsAmongCenters)))
      PPoints<-c( PPoints, rnorm( numSamples, mean = meanP + abs(rnorm( n = 1, mean = 0, sd = deltaSDAmongCentroidsAmongCenters)) , sd = devP )  )
    }
    PMatrix<-matrix(c(PPoints,array(1,numSamples)),ncol=numFeatures+1)
    
    # Negative Samples
    NPoints<-c()
    for( i in seq(1:numFeatures)) {
      sbilanciamentoMediaN <- c(sbilanciamentoMediaN,abs(rnorm( n = 1, mean = 0, sd = deltaSDAmongCentroidsAmongCenters)))
      NPoints<-c( NPoints, rnorm( numSamples, mean = meanN + abs(rnorm( n = 1, mean = 0, sd = deltaSDAmongCentroidsAmongCenters)) , sd = devP )  )
    }
    NMatrix<-matrix(c(NPoints,array(-1,numSamples)),ncol=numFeatures+1)
    
    A<-rbind(NMatrix,PMatrix)
    
    sbilanciamentoMediaP<-mean(sbilanciamentoMediaP)
    sbilanciamentoMediaN<-mean(sbilanciamentoMediaN)

    return( list("A"=A, "sdData"=c(sbilanciamentoMediaP,sbilanciamentoMediaN) ) )
  }
  
  # plotBehaviour: Plot the behaviour of the x vector on time
  # gotFromRun = is the result of the execution of a $run method
  plotBehaviour<-function(res) {
    
    numberOfCenters<-dim(res$logRun[[1]])[1]
    numberOfFeatures<-dim(res$logRun[[1]])[2]
    numberOflog<-length(res$logRun)
    bigArray<-unlist(res$logRun)
    
    DFO<-matrix(0,ncol=numberOflog,nrow=numberOfCenters)
    
    if( writeTo == "file")  png( image03 ) 
    
    for(logNum in seq(1,numberOflog)) {
      for(centNum in seq(1,numberOfCenters)) {
        DFO[centNum,logNum]<-sqrt(res$logRun[[logNum]][centNum,]%*%res$logRun[[logNum]][centNum,])
      }
    }
    plot(0.1,0.1,xlim=c(0,dim(DFO)[2]),ylim=c(min(DFO),max(DFO)),xlab='Step')
    for(i in seq(1,dim(DFO)[1])) {
      points(seq(1,dim(DFO)[2]),DFO[i,],type='l',lty='dotted')
    }
    points(seq(1,dim(DFO)[2]),colMeans(DFO),type='l',col="Red")
    
    if( writeTo == "file") dev.off()
    
    print(c(c(0,dim(DFO)[2]),c(min(DFO),max(DFO)  )))
    return(DFO)
  }
  
  # plotPoints: Plot the point in the space and draws the true SVM
  # gotFromRun = is the result of the execution of a $run method 
  # samplePointsList = is the list with the points
  plotPoints<-function(samplePointsList,gotFromRun,performances=FALSE) {
    if(dim(samplePointsList[[1]])[2]!=3) {
      print("I can plot only a 2 dimensional space, sorry");
      return;
    }
    res<-gotFromRun$x
    logRun<-gotFromRun$logRun
    minVal<-min(unlist(samplePointsList))
    maxVal<-max(unlist(samplePointsList))

    if( writeTo == "file")  png( image04 ) 
    plot(mean(c(minVal,maxVal)),mean(c(minVal,maxVal)),ylim=c(minVal,maxVal),xlim=c(minVal,maxVal),col="Black",pch=3)  
    for(i in seq(1,length(samplePointsList))) {
      for(riga in seq(1,dim(samplePointsList[[i]])[1] )) {

        if( samplePointsList[[i]][riga,3] == 1 ) {
          if( performances$sampleClassification[[i]][riga] ==1 )
            {points(samplePointsList[[i]][riga,1],samplePointsList[[i]][riga,2],col='Red', pch = 1)}
          else
            {points(samplePointsList[[i]][riga,1],samplePointsList[[i]][riga,2],col='Red', pch=16); }
        }
        else {
          if( performances$sampleClassification[[i]][riga] ==1 )
            {points(samplePointsList[[i]][riga,1],samplePointsList[[i]][riga,2],col='Blue', pch = 1)}
          else
            {points(samplePointsList[[i]][riga,1],samplePointsList[[i]][riga,2],col='Blue', pch=16);}
        }
      }
    }
    
    xPoints<-c(minVal,maxVal)

    for(i in seq(1,length(logRun))) {
      vector<-colMeans(logRun[[i]])
      y=-xPoints*(vector[1]/vector[2])-vector[3]/vector[2]
      points(xPoints,y,type='l',lty='dotted',col='Green')
    }  
    vector<-colMeans(res)
    y=-xPoints*(vector[1]/vector[2])-vector[3]/vector[2]
    points(xPoints,y,type='l',col='Black')   
    
    if( writeTo == "file") dev.off()

  }  
  
  # plotAll: Plot the two diagrams
  plotAll<-function(Samples,res,performances=FALSE) {
    
    if(dim(Samples[[1]])[2]==3) {
      par(mfrow=c(2,1))
      plotPoints(samplePointsList = Samples,gotFromRun = res, performances = performances)
    }  else {par(mfrow=c(1,1))}
    
    
#    plotPoints(Samples,res)
    plotBehaviour(res)
  }
  
  # Constructors: it's the method used to initialize
  # attributes which scope is "global" within the object
  constructor<-function() {
    numNodes<<-0;
    numTotSamples<<-0;
    samples<<-list();
    nodes<<-list();
    timer<<-0;	
    if( alpha==0 ) { alpha<<-1.8; } else { alpha<<-alpha }
    if( rho==0 ) { rho<<-1.0; } else { rho<<-rho }
    if( lambda==0 ) { lambda<<-1.0; } else { lambda<<-lambda }
    MAX_TIME<<-100;	
    ABSTOL<<-1e-4;
    RELTOL<<-1e-2;
    numFeatures<<-0;
    writeTo<<-writeTo
    if(token==0)
      {token<<-paste(c(format(Sys.time(),"%Y%m%d"),as.character(runif(1)*10^8),as.character(runif(1)*10^8)),collapse="_"); }
    basePath<<-basePath    
    image01<<-paste(c(basePath,"/","img01_",token,".png"),collapse='')
    image02<<-paste(c(basePath,"/","img02_",token,".png"),collapse='')
    image03<<-paste(c(basePath,"/","img03_",token,".png"),collapse='')
    image04<<-paste(c(basePath,"/","img04_",token,".png"),collapse='')
    
  }
  # invokes functions/methods that must be executed when
  # the object is instantiated   
  constructor();
  return(list(plotAll=plotAll,plotBehaviour=plotBehaviour,plotPoints=plotPoints,createPopulation=createPopulation,run=run,addSamples=addSamples,getAttributeFromNode=getAttributeFromNode,getAttributeFromLearner=getAttributeFromLearner,addNode=addNode,nodes=nodes));
  
}
