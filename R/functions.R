
# set to zero all the negative value of a matrix
pos<-function(A) {
  A[A<0]<-0 
  return(A);
}

# function to be minimized
DLQPFunctionSVM<-function(xi,Ai,rho,ui,zi) {
  total<-sum(pos(Ai%*%xi+1))+rho/2*sum(t(xi-zi+ui)%*%(xi-zi+ui))
  return(total)
}
DLQPFunctionLogisticRegression<-function(xi,Ai,rho,ui,zi) {
  total<-sum(pos(Ai%*%xi+1))+rho/2*sum(t(xi-zi+ui)%*%(xi-zi+ui))
  return(total)
}

planePointDistance<-function(piano,punto) {
  dist_pP<-sum(head(piano,-1)*punto+piano[length(piano)]) / ( sum(head(piano,-1)^2 ))
  return(dist_pP)
} 
distanzaPianoPunto<-function(piano,punto) {
  dist_pP<-(sum(head(piano,-1)*punto)+piano[length(piano)]) / ( sqrt(sum(head(piano,-1)^2 )))

  #dist_pP<-(punto[2]+punto[1] * piano[1]/piano[2]+piano[3]/piano[2]) / sqrt(1 + (piano[1]/piano[2])^2 )
  return(dist_pP)
} 

writeFilesToFileSystem<-function(listaDati,basePath="/home/kboaria/Desktop/math", baseNameFile="dl_test_rg") {
  
  for( i in seq(1:length(listaDati))) {
    nameFile<-paste(basePath,"/",baseNameFile,"-node",i,".dat",sep="");
    fileConn<-file(nameFile)
    
    # writeLines(toMathSyntax( as.data.frame(Samples[[i]]),"V3",1 ), fileConn)
    writeLines( listaDati[[ i ]], fileConn)
    close(fileConn)
  }
  
}
calcolaPerformances<-function(beta , samples   ) {
  sampleClassification<-c()
  TP<-0; TN<-0; FP<-0; FN<-0;
  totDist<-list()
  sampleClassification<-list()
    
  for(centro in seq(1: length(samples))) {
    totDist[[centro]]<-list()
    sampleClassification[[centro]]<-list()
    for(i in seq(1: dim(samples[[centro]])[1])) {
      distanza <- distanzaPianoPunto( beta , samples[[centro]][i,1:(dim(samples[[centro]])[2]-1)]) 
      totDist[[centro]]<-c(totDist[[centro]],distanza)
      
      if(distanza>0 & samples[[centro]][i,dim(samples[[centro]])[2]] ==-1 ) { sampleClassification[[centro]]<-c(sampleClassification[[centro]],1); TP<-TP+1}
      if(distanza<0 & samples[[centro]][i,dim(samples[[centro]])[2]] ==1 ) { sampleClassification[[centro]]<-c(sampleClassification[[centro]],1);TN<-TN+1}
      if(distanza<0 & samples[[centro]][i,dim(samples[[centro]])[2]] ==-1 ) { sampleClassification[[centro]]<-c(sampleClassification[[centro]],-1);FP<-FP+1}
      if(distanza>0 & samples[[centro]][i,dim(samples[[centro]])[2]] ==1 ) { sampleClassification[[centro]]<-c(sampleClassification[[centro]],-1);FN<-FN+1}
    }
  }

  return(list( "TP"=TP, "TN"=TN, "FP"=FP, "FN"=FN,"sampleClassification"=sampleClassification   ) )
}
experiment<-function(
  meanP = 3,  
  meanN =-3,   
  devN=1.5,  
  devP = 1.5,
  alpha = 1.8,
  rho = 1,
  lambda = 1,
  nodes=3,    
  numFeatures=2,
  samplesPerNode=300,
  deltaSDAmongCentroidsAmongCenters = 0,
  write4Mathematica = FALSE,
  writeToMatlab = FALSE,
  epsilon = .0001,
  plotIt = TRUE,
  returnResult = TRUE,
  runningPlot = TRUE,
  token = 0, 
  basePath = "./outPath",
  writeTo = "screen"
  ) {
  obj<-1
  rm(obj)
  obj<-dLearn( lambda = lambda, rho = rho, alpha = alpha, token = token,writeTo = writeTo,basePath = basePath );

  # create a sample set, just for testing
  uTMP<-list()  
  uTMP<-obj$createPopulation(samplesPerNode,numFeatures,nodes,meanP = meanP,devP = devP, meanN=meanN, devN=devN, deltaSDAmongCentroidsAmongCenters = deltaSDAmongCentroidsAmongCenters );
  Samples<<-uTMP$Ai
  mediaSD<-uTMP$mediaSD

  # ok, now load the created sample set by the creation of the nodes
  for(i in 1:nodes) {
    obj$addNode(as.character(i));
    obj$addSamples(i,Samples[[i]] )
  }
  
  # time to SVM
  #browser();
  tot<-as.data.frame(Samples[[1]], Samples[[2]])
  if (length(Samples)>2)
    for (nicola in 3:length(Samples)) tot<-as.data.frame(rbind(tot, Samples[[nicola]]))
  
  coovariate<-paste(names(tot[,1: (length(tot)-1)  ]),collapse='+')
  dipendente<-names(tot)[length(tot)] 

  
  tot[[dipendente]][tot[[dipendente]]==-1]<-0

  if( write4Mathematica == TRUE ) {
    centro<-list()
    for(i in seq(1:length(Samples))) {
      nomeColonnaOutput<- names(as.data.frame(Samples[[ i ]]))[length(names(as.data.frame(Samples[[ i ]])))]
      centro[[ i ]]<-toMathSyntax( DFIn = as.data.frame(Samples[[ i ]]),column = nomeColonnaOutput, centro= i )
    }
    # calcola i beta
    writeFilesToFileSystem(centro)    
  }  
  if ( writeToMatlab == TRUE ) {
    toMATLABSyntax( Samples )
  }
  
  # and now... Run!
  res<-obj$run(epsilon = epsilon, runningPlot = runningPlot);  

  coefficienti<-colMeans(res$x);
  performances<-calcolaPerformances(  beta = coefficienti, samples = Samples  ) 
  accuracy<-(performances$TP+performances$TN)/(performances$TP+performances$TN+performances$FP+performances$FN);

  if( plotIt == TRUE ) {
    obj$plotAll(Samples,res,performances=performances);
    #obj$plotPoints(Samples,res,performances=performances);
    obs<-seq(max(min(tot$V1),min(tot$V2)),max(max(tot$V1),max(tot$V2)),by=.1)
  }
  

  if( returnResult == TRUE ) {
    dl<-colMeans(res$logRun[[ length(res$logRun) ]])
    return( list("DL"=dl, "log"=res, "mediaSD"=mediaSD,"accuracy"=accuracy) )
  }
}

multiExperiment<-function(
  
  numberOfExperiment = 10000,
  
  fromSD = .1, 
  toSD = 1, 
    
  fromDeltaSDAmongCentroidsAmongCenters = 0,
  toDeltaSDAmongCentroidsAmongCenters  = 0,

  fromNodes=3, 
  toNodes=3, 

  fromNumFeatures=2,
  toNumFeatures=2,

  fromRho = 1,
  toRho = 1,
  fromAlpha = 1.8,
  toAlpha = 1.8,
  fromLambda = 1,
  toLambda = 1,
  
  fromSamplesPerNode=300,
  toSamplesPerNode=300,
  
  plotIt = TRUE,                         # plot It  
  runningPlot = TRUE,
  
  #basePath="/home/kboaria/Desktop/math", # folder 
  basePath="./outPath", # folder 
  baseNameFile="listaResult.cat",        # filename
  buildCSV = FALSE,                       # Do I have to write results on CSV??!?!
  token = 0,                             # token about ID png file name
  writeTo = "screen"
  ) {
  
  nameFile<-paste(basePath,"/",baseNameFile,sep="");
  
  
  stringa<-"i,mediaDistCentrCentroidi,mediaSD,nodes,numFeatures,rho,lambda,alpha,samplesPerNode,accuracy,numberOfIterations,usedTime"
  if( buildCSV == TRUE ) {
    if(!file.exists(nameFile))
      write( stringa, file = nameFile, append=TRUE)
  }
  
  for(i in seq( 1: numberOfExperiment )) {
    # get the values
    devN<-runif(1)*(toSD-fromSD)+fromSD
    devP<-runif(1)*(toSD-fromSD)+fromSD
    
    nodes<-as.integer(  runif(1)*(toNodes-fromNodes+1)+fromNodes )
    numFeatures<-as.integer(  runif(1)*(toNumFeatures-fromNumFeatures+1)+fromNumFeatures )
    samplesPerNode<-as.integer(  runif(1)*(toSamplesPerNode-fromSamplesPerNode+1)+fromSamplesPerNode )
    deltaSDAmongCentroidsAmongCenters<-runif(1)*(toDeltaSDAmongCentroidsAmongCenters-fromDeltaSDAmongCentroidsAmongCenters)+fromDeltaSDAmongCentroidsAmongCenters
    
    rho<-runif(1)*(toRho-fromRho)+fromRho 
    lambda<-runif(1)*(toLambda-fromLambda)+fromLambda 
    alpha<-runif(1)*(toAlpha-fromAlpha)+fromAlpha 
    
    # Experiment!
    usedTime<-system.time(a<-experiment(meanP=-1, meanN=1, devN = devN,devP = devP, nodes = nodes,
                  numFeatures = numFeatures, samplesPerNode = samplesPerNode,
                  rho = rho, lambda = lambda, alpha = alpha,
                  plotIt = plotIt, runningPlot = runningPlot,
                  deltaSDAmongCentroidsAmongCenters = deltaSDAmongCentroidsAmongCenters,
                  token = token, writeTo = writeTo, basePath = basePath
                  ))
    usedTime<-usedTime['elapsed']
    # Compose the results

    numberOfIterations<-length( a$log$logRun )
    errorAmongX<-mean(apply(a$log$logRun[[numberOfIterations ]],2,sd))
    accuracy<-a$accuracy
    
    stringa<-paste(c(
      i,                          # progressivo
      mean(abs(devN),abs(devP)),  # media delle deviazioni standard delle campane dei singoli centri
      a$mediaSD,                  # SD media delle campane rispetto ai centroidi delle 2 classi
      nodes,                      # numero nodi
      numFeatures,                # numero features
      rho,
      lambda,
      alpha,
      samplesPerNode,
      accuracy,
      numberOfIterations,
      usedTime
    ),collapse=',')
    
    # WRITE on Disk!!!! (appending the line....)  
    if( buildCSV == TRUE ) {
      write( stringa, file = nameFile, append=TRUE)
    }
  }
  
}

testDataset<-function( URL = 'http://5.249.147.20:8080/sparql',  arrayOfToken = c(11,12,13) )  {
    ct<-1
    for( token in arrayOfToken) {
      d<-SPARQL(url=URL ,query=paste(c('
        PREFIX db: <http://localhost:8080/resource/>
        PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
        PREFIX owl: <http://www.w3.org/2002/07/owl#>
        PREFIX map: <http://localhost:8080/resource/#>
        PREFIX xsd: <http://www.w3.org/2001/XMLSchema#>
        PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
        PREFIX vocab: <http://localhost:8080/resource/vocab/>
          select ?Query ?URL ?isValid
                where {
                  ?x  <http://localhost:8080/resource/vocab/dataSourceDetails_SPARQLQuery> ?Query.
                  ?x  <http://localhost:8080/resource/vocab/dataSourceDetails_URL> ?URL.
                  ?x  <http://localhost:8080/resource/vocab/dataSourceDetails_isValid> 1.
                  ?x  <http://localhost:8080/resource/vocab/dataSourceDetails_id_dataSourceDetails> 
                ',token,' }'),collapse=''))
      # Retrieve data
      dataSource<-SPARQL(url = d$results$URL,query = d$results$Query);
      # format it as a matrix
      A[[ct]]<-as.matrix(dataSource$res);
      ct<-ct+1
    }
}


#multiExperiment(  numberOfExperiment=1000, fromSD = .2,toSD = .4,fromDeltaSDAmongCentroidsAmongCenters = 0,toDeltaSDAmongCentroidsAmongCenters = .5,
#                  fromNodes = 2,toNodes = 10, fromSamplesPerNode = 50, toSamplesPerNode = 500, plotIt = TRUE, fromNumFeatures = 2, 
#                  toNumFeatures = 2,fromLambda = 1.2, toLambda=1.8, fromAlpha = 1.2, toAlpha = 1.8, fromRho = 1.2, 
#                  toRho = 1.8,buildCSV = FALSE)