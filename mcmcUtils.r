library(mvtnorm)

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

##### Slice up function #####

sliceUp = function(xc,mu,sig2,lb){
  u1 = runif(1)
  sqrtPart = sqrt((xc-mu)^2 - 2*sig2*log(u1))
  ub = mu + sqrtPart
  lb = max(mu - sqrtPart, lb)
  
  u2 = runif(1)
  x  = lb + (ub-lb)*u2
  return(x)  
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
  deltaTrue = array(0,dim=c(prjCtrl$demoCOV,prjCtrl$nS))
  #deltaTrue = array(0,dim=c(prjCtrl$demoCOV,1))
  

  deltaStart = c(-3, .09, 1,-2, .1, -2.5, .5,2.1,-1.2,.01,2.2,.75,-.8,-1.9,0,3.6)
  slopeBarStart = c(-1,2,.5,8,-5,2.2,5)
  slopeCovStart = c(9,3,4,5,3,1,18)
  
  slopeBarK = 7
  deltaK = length(deltaStart)
  
  k = 0
  for(j in 1:prjCtrl$COV){
    k = k+1
    if(k>slopeBarK){
      k=1
    }
    slopeBarTrue[j,1] = slopeBarStart[k]
    slopeCovTrue[j,j] = slopeCovStart[k]
  }
  

  for(j in 2:prjCtrl$nS){
    for(i in 1:prjCtrl$demoCOV){

        k = k+1
        if(k>deltaK){
          k=1
        }
        
        deltaTrue[i,j] = deltaStart[k]
      
    }  
  }
  
  
  
  SSE = array(0,dim=c(prjCtrl$COV,prjCtrl$COV,prjCtrl$nS))
  slopeCorrMLE = array(0,dim=c(prjCtrl$COV,prjCtrl$COV,prjCtrl$nS))
  
  
  

  
  # Generate true individual state classifications from generated deltas

  sTrue = array(0,c(prjCtrl$IND,1))
  pTrue = array(1,c(prjCtrl$nS,1))*(1/prjCtrl$nS)
  nSTrue = array(0,c(prjCtrl$nS,1))
  
  if (prjCtrl$runDelta){
    for(i in 1:prjCtrl$IND){
      
      eta = dataZ[i,] %*% deltaTrue
      sProbs = exp(eta) / sum(exp(eta))
      pTrue = t(sProbs)
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
  
  if(prjCtrl$useConstraints){
    atLeastOneNot0 = FALSE
    while(atLeastOneNot0 == FALSE){
      for(j in 1:prjCtrl$COV){
        if(runif(1) < .2){
          if(runif(1) < .5){
            param$slopeConstraintFlag[j,1] = 1
            atLeastOneNot0 = TRUE
          }
          else{
            param$slopeConstraintFlag[j,1] = -1
            atLeastOneNot0 = TRUE
          }
        }
      }
    }
    zer0 = 0
    for(i in 1:prjCtrl$IND){
      #slopTrue: (COV x 1 x IND) - vector of column vectors
      slopeTrue[,1,i] = t(rmvnorm(n=1,mean=slopeBarTrue,sigma=slopeCovTrue))
      for(j in 1:prjCtrl$COV){
        if(param$slopeConstraintFlag[j] != 0){
          if(param$slopeConstraintFlag[j] == 1){
            if(slopeTrue[j,1,i] < 0){
              slopeTrue[j,1,i] = 1
              for(t in 1:25){
                # Slice up
                slopeTrue[j,1,i] = sliceUp(slopeTrue[j,1,i],slopeBarTrue[j],1,zer0)
              }
            }
          }
        }
        else{
          if(slopeTrue[j,1,i] > 0){
            slopeTrue[j,1,i] = -1
            for(t in 1:25){
              # Slice up
              slopeTrue[j,1,i] = sliceUp(slopeTrue[j,1,i],slopeBarTrue[j],1,zer0)
            }
          }
        }
      }
      SSE[,,sTrue[i]] = SSE[,,sTrue[i]] + (slopeTrue[,,i] - slopeBarTrue) %*% t((slopeTrue[,,i] - slopeBarTrue))
    }
  }
  else{
    # No constraints on beta bar
    for(i in 1:prjCtrl$IND){
      slopeTrue[,1,i] = t(rmvnorm(n=1,mean=slopeBarTrue,sigma=slopeCovTrue))
      SSE[,,sTrue[i]] = SSE[,,sTrue[i]] + (slopeTrue[,,i] - slopeBarTrue) %*% t((slopeTrue[,,i] - slopeBarTrue))
    }
  }
  
  for(i in 1:prjCtrl$nS){
    SSE[,,i] = SSE[,,i] / nSTrue[i,1] 
    tMLE = diag(1 / sqrt(diag(SSE[,,i])))
    slopeCorrMLE[,,i] = tMLE%*%SSE[,,i]%*%tMLE
  }
  
  
  param$slope = slopeTrue
  param$slopeBar = slopeBarTrue
  param$slopeCov = slopeCovTrue
  param$sig2 = sig2True
  
  param$delta = deltaTrue
  
  param$nS = nSTrue
  param$p = pTrue
  param$s = sTrue
  

  return(param)
}


genHRLData = function(prjCtrl,param,dataZ){
  print("Generating synthetic data...")
  data = list(
    
    y = array(0,dim=c(prjCtrl$nRep,1,prjCtrl$IND)),
    # X tensor dimensions: (num products x num covariates x num task repitions x n respondents)
    X = array(0,dim=c(prjCtrl$PROD,prjCtrl$COV,prjCtrl$nRep,prjCtrl$IND))
    #Z = array(0,dim=c(prjCtrl$IND,prjCtrl$demoCOV))
    
  )
  
  # Beta covariance for data generation - diagonal covariance only
  cov_X = diag(prjCtrl$COV)
  cov_x_start = c(2,1,3,1,2,1,2)
  cov_x_k = 7
  k = 0
  for(j in 1:prjCtrl$COV){
    k = k + 1
    if(k > cov_x_k){
      k = 1
    }
    cov_X[j,j] = cov_x_start[k]
    
  }
  
  
  
  for(i in 1:prjCtrl$IND){
  
    for(j in 1:prjCtrl$nRep){
      tProb = array(0,c(prjCtrl$PROD,1))
      if(prjCtrl$includeIntercept){
        data$X[,1,j,i] = 1
        data$X[,2:prjCtrl$COV,j,i] = rmvnorm(n=prjCtrl$PROD,mean = rep(0,(prjCtrl$COV-1)),sigma = cov_X[2:prjCtrl$COV,2:prjCtrl$COV])
      }else{
        # No intercept
        data$X[,,j,i] = rmvnorm(n=prjCtrl$PROD,mean = rep(0,prjCtrl$COV),sigma = cov_X)
      }
      
      tProb = (1/param$sig2[param$s[i]])*data$X[,,j,i]%*%param$slope[,1,i]
      maxProb = max(tProb)
      tProb = exp(tProb - maxProb)
      tProb = tProb / sum(tProb)
      data$y[j,1,i] = loadedDie(tProb)
  
    }
  }
  
  if(is.null(dataZ) == FALSE){
    data[['Z']] = dataZ
  }
  
  
  return(data)
}





