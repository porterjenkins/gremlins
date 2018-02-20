#source("prjCtrl.R")

setHRLPrior = function(COV,demoCOV){
  # Create data structure for prior hyperparameters

  # :param COV: (int) Number of covariates in choice model
  # :param demoCOV: (int) Number of demographic covariates in logit for mixture component

  # :return: list of hyperparameters for each prior distribution
  

  prior = list(
  
  slopeBarBar = matrix(0,nrow=COV,1),
  invSlopeBarCov = solve(100*diag(COV)),
  
  
  delta =  matrix(0,nrow=demoCOV,1),
  invDeltaCov = solve(100*diag(demoCOV)),
  
  slopeCovScale = 1,
  slopeCovShape = 1,
  
  sig2Scale = 1,
  sig2Shape = 1,
  sig2GammaScale = 100,
  sig2GammeShape = 5
  
  )
  
  return(prior)
}


setHRLSS = function(){
  
  paramSS = list(
    
    
  )
  
  return(paramSS)
  
}



