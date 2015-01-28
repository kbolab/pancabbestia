
# ------------------------------------------
# Class dNODE
# ------------------------------------------
dNode<-function(p_externalName,p_externalID,algo="SVM", rho = 0, alpha = 0, lambda = 0) {
  # ------------------------------------------
  # Attributes
  # ------------------------------------------
  samples<-list();  
  externalName<-'';
  externalID<-'';
#  lambda<-0;
#  rho<-0;
#  alpha<-0;
  xVar<-0;
  numFeatures<-0;
  algorithm<-algo;
    
# jiojio
  # ------------------------------------------
  # Methods
  # ------------------------------------------
  
  # addSamples: adds an array of samples
  addSamples<-function(arrSample) {
    samples<<-arrSample;
    numFeatures<<-length(arrSample[1,]);
    xVar<<-matrix(1,ncol=numFeatures+1,nrow=1)
  }
  
  # getAttribute: return an attributa value
  getAttribute<-function(whichAttribute) {
    if( whichAttribute == 'externalName' ) return(externalName);
    if( whichAttribute == 'externalID' ) return(externalID);
    if( whichAttribute == 'samples' ) return(samples);
  }
  
  # takeStep: taken x,z,u calculates the new x
  # (n.b: the only variable needed for the computation would be
  # z because u and x could be stored in the node. This flow is done 
  # in order to simplify the code)
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
   
  # Constructors: it's the method used to initialize
  # attributes which scope is "global" within the object
  constructor<-function(algo) {
    externalName<<-p_externalName;
    externalID<<-p_externalID;
    if( lambda == 0 ) {lambda<<-1;} else {lambda<<-lambda}
    if( alpha == 0 ) {alpha<<-1.8;} else {alpha<<-alpha}
    if( rho == 0 ) {rho<<-1;} else {rho<<-rho}
    numFeatures<<-0;
    xVar<<-0;
    algorithm<<-algo;
  }
  
  # invokes functions/methods that must be executed when
  # the object is instantiated 
  constructor(algo);
  return(list(addSamples=addSamples,getAttribute=getAttribute,takeStep=takeStep));
}
