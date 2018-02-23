#source("prjCtrl.R")


##### "Roll" Loaded Die #####

loadedDie = function(prob){
  
  u = runif(1)
  state = 1
  cProb = prob[state]
  while(u >= cProb){
    state = state + 1
    cProb = cProb + prob[state]
    
  }
  
  return(state)
}

##### Initialize Prior Data Structure #####
setHRLPrior = function(COV,demoCOV){
  # Create data structure for prior hyperparameters

  # :param COV: (int) Number of covariates in choice model
  # :param demoCOV: (int) Number of demographic covariates in logit for mixture component

  # :return: list of hyperparameters for each prior distribution
  

  prior = list(
  
    slopeBarBar = array(0,dim=c(COV,1)),
    invSlopeBarCov = solve(100*diag(COV)),
    
    
    delta =  array(0,dim=c(demoCOV,1)),
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

##### Initialize SS Parameter Data Structure #####

setHRLSS = function(COV,IND,nS){
  
  paramSS = list(
    
    slope = array(0, dim=c(COV,2,IND)),
    slopeOverSig2 = array(0, dim=c(COV,2,IND)),
    
    slopeBar = array(0,dim=c(COV,2)),
    slopeCov1 = array(0,dim=COV),
    slopeCov2 = array(0,dim=COV),
    sig2 = array(0,dim=c(2,1,nS)),
    nS = array(0,dim=c(nS,2)),
    s = array(0,dim=c(IND,nS)),
    p = array(0,dim=c(IND,nS)),
    N = 0
    
    
  )
  
  return(paramSS)
  
}


##### Initialize Parameter Data Structure #####

setHRLParams = function(COV,IND,nS,useConstraints){
  
  param = list(
    
    slope = array(0, dim=c(COV,1,IND)),
    slopeOverSig2 = array(0, dim=c(COV,1,IND)),
    
    slopeBar = array(0,dim=c(COV,1)),
    slopeCov1 = array(0,dim=COV),
    slopeCov2 = array(0,dim=COV),
    sig2 = array(0,dim=c(2,1,nS)),
    nS = array(0,dim=c(nS,2)),
    s = array(0,dim=c(IND,nS)),
    p = array(0,dim=c(IND,nS)),
    N = 0,
    
    indHitRate = array(0,dim=c(IND,1)),
    indLL =  array(0,dim=c(IND,1)),
    
    TNPARAM = (IND + 2) * COV + 2*nS,
    sSTNREP =  array(0,dim=c(nS,1))

  )
  
  if(useConstraints){
    param[["slopeConstraintFlag"]] = array(0,dim=c(COV,1))
  }
  
  return(param)
  
  
}

##### Generate True Parameter Values #####
genHRLParams = function(prjCtrl,param,dataZ,prior){
  
  slopeTrue = array(0,dim=c(prjCtrl$COV,1,prjCtrl$IND))
  slopeBarTrue = array(0,dim=c(prjCtrl$COV,1))
  slopeCovTrue = diag(prjCtrl$COV)
  
  sig2True = array(0,dim=c(prjCtrl$nS,1))
  sig2True[1,1] = 1
  
  sig2Bump = 20
  for(i in 2:prjCtrl$nS){
    sig2True[i,1] = sig2True[i-1,1] * (1 + sig2Bump*runif(1))
  }
  
  #deltaTrue = t(rmvnorm(n=prjCtrl$nS,mean=rep(0,prjCtrl$demoCOV),sigma=.5*diag(prjCtrl$demoCOV)))
  #deltaTrue = array(0,dim=c(prjCtrl$demoCOV,prjCtrl$nS))
  deltaTrue = array(0,dim=c(prjCtrl$demoCOV,1))
  

  deltaStart = c(-3, .09, 1,-2, .1, -2.5, .5)
  slopeBarStart = c(-1,2,.5,8,-5,2.2,5)
  slopeCovStart = c(9,3,4,5,3,1,18)
  
  slopeBarK = 7
  deltaK = 7
  
  k = 0
  for(j in 1:prjCtrl$COV){
    k = k+1
    if(k>slopeBarK){
      k=1
    }
    slopeBarTrue[j,1] = slopeBarStart[k]
    slopeCovTrue[j,j] = slopeCovStart[k]
  }
  
  k = 0
  for(j in 1:prjCtrl$demoCOV){
    k = k + 1
    if(k>deltaK){
      k=1
    }
    deltaTrue[j,1] = deltaStart[k]
  }
  
  
  
  ## TODO: Consider case where nS > 2. Then segment probability regression shouldb be multinomial. Should delta be a matrix (n_covariates x n_segments)
  #for(j in 1:prjCtrl$nS){
  #  for(i in 1:prjCtrl$demoCOV){
  #    k = k+1
  #    if(k>deltaK){
  #      k=1
  #    }
  #    deltaTrue[i,j] = deltaStart[k]
      
  #  }  
  #}
  
  # Generate true individual state classifications from generated deltas
  
  sTrue = array(0,c(prjCtrl$IND,1))
  pTrue = array(1,c(prjCtrl$nS,1))*(1/prjCtrl$nS)
  nSTrue = array(0,c(prjCtrl$nS,1))
  
  if (prjCtrl$runDelta){
    for(i in 1:prjCtrl$IND){
      
      eta = dataZ[i,] %*% deltaTrue
      sProbs = exp(eta) / (1 + exp(eta))
      pTrue[1] = sProbs
      pTrue[2] = 1 - sProbs
      sTrue[i,1] = loadedDie(pTrue)
      
    }
    
    
  }
  else{
    
    for(i in 1:prjCtrl$IND){
      sTrue[i,1] = loadedDie(pTrue) 
    }
    
  }
  
  for(i in 1:prjCtrl$nS){
    nSTrue[i,1] = length(sTrue[sTrue == i])
  }
  
  print(nSTrue)
  
  
  param$slopeBar = slopeBarTrue
  param$delta = deltaTrue
  param4slopCov = slopeCovTrue
  
  
  
  
  
  return(param)
}





