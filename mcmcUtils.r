library(mvtnorm)
library(coda)

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
setHRLPrior = function(prjCtrl){
  # Create data structure for prior hyperparameters

  # :param COV: (int) Number of covariates in choice model
  # :param demoCOV: (int) Number of demographic covariates in logit for mixture component

  # :return: list of hyperparameters for each prior distribution
  

  prior = list(
  
    slopeBarBar = array(0,dim=c(prjCtrl$COV,1)),
    invSlopeBarCov = solve(100*diag(prjCtrl$COV)),
    
    
    delta =  array(0,dim=c(prjCtrl$demoCOV,1)),
    invDeltaCov = solve(100*diag(prjCtrl$demoCOV)),
    
    slopeCovScale = 1,
    slopeCovShape = 1,
    
    sig2Scale = 1,
    sig2Shape = 1,
    sig2GammaScale = 100,
    sig2GammaShape = 5
  
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

setHRLParams = function(prjCtrl){
  
  param = list(
    
    slope = array(0, dim=c(prjCtrl$COV,1,prjCtrl$IND)),
    slopeOverSig2 = array(0, dim=c(prjCtrl$COV,1,prjCtrl$IND)),
    
    slopeBar = array(0,dim=c(prjCtrl$COV,1)),
    slopeCov = array(0,dim=c(prjCtrl$COV,prjCtrl$COV)),
    sig2 = array(0,dim=c(2,1,prjCtrl$nS)),
    nS = array(0,dim=c(prjCtrl$nS,2)),
    s = array(0,dim=c(prjCtrl$IND,1)),
    p = array(0,dim=c(prjCtrl$IND,prjCtrl$nS)),
    phi = array(0,dim=c(prjCtrl$IND,prjCtrl$nS)),
    N = 0,
    
    indHitRate = array(0,dim=c(prjCtrl$IND,1)),
    indLL =  array(0,dim=c(prjCtrl$IND,1)),
    
    TNPARAM = (prjCtrl$IND + 2) * prjCtrl$COV + 2*prjCtrl$nS,
    sSTNREP =  array(0,dim=c(prjCtrl$nS,1))

  )
  
  if(prjCtrl$useConstraints){
    param[["slopeConstraintFlag"]] = array(0,dim=c(COV,1))
  }
  if(prjCtrl$runDelta){
    param[["delta"]] = array(0,dim=c(prjCtrl$demoCOV,prjCtrl$nS))
  }

  return(param)
  
}

##### Initialize Parameter Data Draws Structure #####

initParamDraws = function(prjCtrl){
  
  if(prjCtrl$nSample %% prjCtrl$thin == 0){
    n_save_draws = prjCtrl$nSample / prjCtrl$thin
  }else{
    stop('Project control attributes, nSample must be divisible by thin')
  }
  
  paramDraws = list(
    
    slopeBar = array(0,dim = c(prjCtrl$COV,1,n_save_draws)),
    delta = array(0,dim = c(prjCtrl$demoCOV,prjCtrl$nS,n_save_draws)),
    slopeCov = array(0,dim = c(prjCtrl$COV,prjCtrl$COV,n_save_draws)),
  
    slope = array(0,dim = c(prjCtrl$COV,1,prjCtrl$IND,n_save_draws)),
    s = array(0,dim=c(prjCtrl$IND,1,n_save_draws)),
    phi = array(0,dim=c(prjCtrl$IND,prjCtrl$nS,n_save_draws))
  )
  
  
  return(paramDraws)
}

##### Save current parameter draws to parameter draws object  #####

saveCurrentDraw = function(paramDraws,param,n,prjCtrl){
  
  save_idx = (n-prjCtrl$nBurnin)/prjCtrl$thin
  
  paramDraws$slopeBar[,,save_idx] = param$slopeBar
  paramDraws$delta[,,save_idx] = param$delta
  paramDraws$slopeCov[,,save_idx] = param$slopeCov
  
  paramDraws$slope[,,,save_idx] = param$slope
  paramDraws$s[,,save_idx] = param$s
  paramDraws$phi[,,save_idx] = param$phi

  
  return(paramDraws)
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
    X = array(0,dim=c(prjCtrl$PROD,prjCtrl$COV,prjCtrl$nRep,prjCtrl$IND)),
    TNREP = prjCtrl$IND * prjCtrl$nRep
    
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

##### Compute log likelihood of multinomial logit #####

mnlLogLike = function(Y,X,beta,NREP,prjCtrl){
  probs = array(0,dim=c(NREP,prjCtrl$PROD))
  
  for(j in 1:prjCtrl$PROD){
    probs[,j] = t(X[j,,]) %*% beta
  }
  maxProb = array(0,dim=c(NREP,1))
  
  for(n in 1:NREP){
    maxProb[n,1] = max(probs[n,])
  }
  
  for(j in 1:prjCtrl$PROD){
    probs[,j] = probs[,j] - maxProb
  }
  
  log_like = 0
  for(n in 1:NREP){
    for(j in 1:prjCtrl$PROD){
      y_icj_idx = Y[n]
      if(y_icj_idx == j){
        y = 1
      }else{
        y = 0
      }
    
      log_like_icj = y*(probs[n,j] - log(sum(exp(probs[n,]))))
      log_like = log_like + log_like_icj
    }
  }
  
  
  return(log_like)
}

##### Compute log likelihood of multinomial logit for state probabilities #####

mnlLogLikeState = function(stateVector, dataZ, delta, prjCtrl){
  probs = array(0,c(prjCtrl$IND,prjCtrl$nS))
  logLike = 0
  eta = dataZ %*% delta
  
  for(i in 1:prjCtrl$IND){
    max_i = max(eta[i,])
    eta_i = eta[i,] - max_i
    probs[i,] = eta_i - log(sum(exp(eta_i)))
    
    s_i = stateVector[i]
    logLike = logLike + probs[i,s_i] 
    
  }
  
  return(logLike)
}


##### Generate updated beta bar #####
genAggLogit = function(data,prjCtrl,slope,slopeBar,slopeCov){
  jump_size_index = loadedDie(prjCtrl$rwParamProb)
  rw_sig = prjCtrl$rwParamRwSig[jump_size_index]
  p_slope = t(rmvnorm(n=1,mean=slope,sigma = rw_sig*diag(prjCtrl$COV)))
  
  cLL = mnlLogLike(data$Y,data$X,slope,data$TNREP,prjCtrl)
  pLL = mnlLogLike(data$Y,data$X,p_slope,data$TNREP,prjCtrl)
  
  cPrior = -0.5*t(slope-slopeBar)%*%solve(slopeCov)%*%(slope-slopeBar)
  pPrior = -0.5*t(slope-slopeBar)%*%solve(slopeCov)%*%(slope-slopeBar)
  
  lap = pLL + pPrior - cLL - cPrior
  
  a_rate = prjCtrl$aRate
  if(log(runif(1) < lap)){
    slope = p_slope
    a_rate[jump_size_index,1] = a_rate[jump_size_index,1] + 1
    
  }
  
  a_rate[jump_size_index,2] = a_rate[jump_size_index,2] + 1
  return(list(slope,a_rate))
}



##### Get starting values using MLE routine #####
mleHRLStart = function(prjCtrl,data,param,prior){
  
  tSlopeBar = array(0,dim=c(prjCtrl$COV,1))
  ssTSlopeBar = tSlopeBar
  tData = list(
    Y = array(0,dim=c(data$TNREP,1)),
    X = array(0,dim=c(data$TNREP,prjCtrl$COV,prjCtrl$PROD)),
    TNREP = data$TNREP
    
  )
  
  tStart = 1
  tEnd = prjCtrl$nRep
  
  for(i in 1:prjCtrl$IND){
    tData$Y[tStart:tEnd,1] = data$y[,,i]
    for(j in 1:prjCtrl$PROD){
      tData$X[tStart:tEnd,,j] = t(data$X[j,,,i])
    }
    tStart = tEnd + 1
    if(i < prjCtrl$IND){
      tEnd = tStart + prjCtrl$nRep - 1
    }
  }
  
  slopeBarBar = array(0,c(prjCtrl$COV, 1))
  slopeBarCov = 25*diag(prjCtrl$COV)
  
  mle_start_iterations = round(.5*(prjCtrl$startBurnin + prjCtrl$startSample))
  for(n in 1:mle_start_iterations){
    if(prjCtrl$useConstraints){
      #TODO: Implement aggregate logit with constraints
    }else{
        genLogit = genAggLogit(tData,prjCtrl,tSlopeBar,slopeBarBar,slopeBarCov)
        tSlopeBar = genLogit[[1]]
        prjCtrl$aRate = genLogit[[2]]
      
      if(n > prjCtrl$startBurnin){
        ssTSlopeBar = ssTSlopeBar + tSlopeBar
      }
    }
  }
  
  ssTSlopeBar = ssTSlopeBar / prjCtrl$startSample
  
  # Create invidual slope estimates, using aggregate mean
  
  
  
  return(param)
}
##### Get starting values by simulating from priors #####

priorHRLStart = function(prjCtrl,param,prior){
  #[1] "slope"         "slopeOverSig2" "slopeBar"      "slopeCov1"    
  #[5] "slopeCov2"     "sig2"          "nS"            "s"            
  #[9] "p"             "N"             "indHitRate"    "indLL"        
  #[13] "TNPARAM"       "sSTNREP"       "slopeCov"      "delta"
  
  # Simulate beta bar
  param$slopeBar[,1] = rmvnorm(n=1,mean = prior$slopeBarBar,sigma = solve(prior$invSlopeBarCov))
  
  # Simulate Slope Covariance - assume independence, each sig2_j is inv gamma(1,1)
  tSlopeCov = diag(prjCtrl$COV)
  for(p in 1:prjCtrl$COV){
    tSlopeCov[p,p] = tSlopeCov[p,p]*rinvgamma(1,shape=prior$sig2Shape,scale=prior$sig2Scale)
  }
  param$slopeCov = tSlopeCov
  
  
  for(j in 2:prjCtrl$nS){
    # Simulate delta
    param$delta[,j] = rmvnorm(n=1,mean = prior$delta,sigma = solve(prior$invDeltaCov))
    # simulate sig2 (called lambda in paper)
    param$sig2[j] = max(param$sig2[j-1]*1.5,rinvgamma(n=1,shape=prior$sig2GammaShape,rate=prior$sig2GammaScale))
  }
  
  for(i in 1:prjCtrl$IND){
    
    param$slope[,,i] = rmvnorm(n=1,mean=param$slopeBar,sigma=param$slopeCov)
    param$s[i,1] = loadedDie(rep(1,prjCtrl$nS)/prjCtrl$nS)
    param$slopeOverSig2[,,i] = param$slope[,,i] / param$sig2[param$s[i,1]]
    
  }
  
  for(s in 1:prjCtrl$nS){
    param$nS[s,1] = length(param$s[param$s == s])
  }
  
  param$phi = array(1/prjCtrl$nS,dim=c(prjCtrl$IND,prjCtrl$nS))
  
  return(param)
}

##### Generate beta (constrained) from full conditional #####
genHRLSlopeConstraint = function(dataY,dataTNREP,dataX,prjCtrl,sig2,slope,slopeBar,slopeCov,slopeConstraint){
  jump_size_index = loadedDie(prjCtrl$rwParamProb)
  slope = slope*(1/sig2)
  rw_sig = prjCtrl$rwParamRwSig[jump_size_index]
  
  t_lap = 0
  
  p_slope = array(0,dim=c(prjCtrl$COV,1))
  for(j in 1:prjCtrl$COV){
    p_slope[j,1] = rnorm(n=1,mean=slope[j],sd=rw_sig)
    if(slopeConstraint[j,1] != 0){
      if(slopeConstraint[j,1] == 1){
        i = 0
        while(p_slope[j,1]<0 && i<100){
          p_slope[j,1] = rnorm(n=1,mean=slope[j],sd=rw_sig)
          i = i + 1
        }
        if(i == 99){
          p_slope[j,1] = slope[j]
        }else{
          t_lap = t_lap + log(1-pnorm(-slope[j]/rw_sig)) - log(1-pnorm(-p_slope[j,1]/rw_sig))
          }
      }else{
        i = 0
        while(p_slope[j,1]>0 && i<100){
          p_slope[j,1] = rnorm(n=1,mean=slope[j],sd=rw_sig)
          i = i + 1
        }
        if(i == 99){
          p_slope[j,1] = slope[j]
        }else{
          t_lap = t_lap + log(1-pnorm(-slope[j]/rw_sig)) - log(1-pnorm(-p_slope[j,1]/rw_sig))
        }
      }
    }
  }
  
  cLL = mnlLogLike(dataY,dataX,slope,dataTNREP,prjCtrl)
  pLL = mnlLogLike(dataY,dataX,p_slope,dataTNREP,prjCtrl)
  
  cPrior = -0.5*t(slope-slopeBar)%*%solve(slopeCov)%*%(slope-slopeBar)
  pPrior = -0.5*t(slope-slopeBar)%*%solve(slopeCov)%*%(slope-slopeBar)
  
  lap = pLL + pPrior - cLL - cPrior
  
  a_rate = prjCtrl$aRate
  alpha = log(runif(1))
  if(alpha < lap){
    slope = p_slope
    a_rate[jump_size_index,1] = a_rate[jump_size_index,1] + 1
    
  }
  
  slope = slope*sig2
  a_rate[jump_size_index,2] = a_rate[jump_size_index,2] + 1
  
  return(list(slope,a_rate))
}

##### Generate beta (unconstrained) from full conditional #####
genHRLSlope = function(dataY,dataTNREP,dataX,prjCtrl,sig2,slope,slopeBar,slopeCov){
  jump_size_index = loadedDie(prjCtrl$rwParamProb)
  slope = slope*(1/sig2)
  rw_sig = prjCtrl$rwParamRwSig[jump_size_index]
  p_slope = t(rmvnorm(n=1,mean=slope,sigma = rw_sig*diag(prjCtrl$COV)))
  
  cLL = mnlLogLike(dataY,dataX,slope,dataTNREP,prjCtrl)
  pLL = mnlLogLike(dataY,dataX,p_slope,dataTNREP,prjCtrl)
  
  cPrior = -0.5*t(slope-slopeBar)%*%solve(slopeCov)%*%(slope-slopeBar)
  pPrior = -0.5*t(slope-slopeBar)%*%solve(slopeCov)%*%(slope-slopeBar)
  
  lap = pLL + pPrior - cLL - cPrior
  
  a_rate = prjCtrl$aRate
  alpha = log(runif(1))
  if(alpha < lap){
    slope = p_slope
    a_rate[jump_size_index,1] = a_rate[jump_size_index,1] + 1
    
  }
  
  slope = slope*sig2
  a_rate[jump_size_index,2] = a_rate[jump_size_index,2] + 1
  return(list(slope,a_rate))
  
}


genHRLSlopeBar = function(prjCtrl,param,prior){
  sum_slope_i = param$slope[,,1]
  for(i in 1:prjCtrl$IND){
    sum_slope_i = sum_slope_i + param$slope[,,i]
  }
  
  slope_inv_cov = solve(param$slopeCov)
  n = prjCtrl$IND
  
  covariance_mtx = solve(n*slope_inv_cov + prior$invSlopeBarCov)
  mean = slope_inv_cov %*% sum_slope_i + prior$invSlopeBarCov %*% prior$slopeBarBar
  mean = covariance_mtx %*% mean
  
  gen_slope_bar = t(rmvnorm(n=1,mean=mean,sigma=covariance_mtx))
  return(gen_slope_bar)
}


genHRLSlopeCov = function(prjCtrl,param,prior){
  slopeCov = array(0,dim=c(prjCtrl$COV,prjCtrl$COV))
  scale = (param$slope[,,1] - param$slopeBar)**2
  shape = prior$slopeCovShape + .5*prjCtrl$IND
  
  for(i in 2:prjCtrl$IND){
    scale = scale + (param$slope[,,i] - param$slopeBar)**2
  }
  scale = array(prior$slopeCovScale,dim=c(prjCtrl$COV,1)) + .5*scale
  
  for(j in 1:prjCtrl$COV){
    slopeCov[j,j] = 1/rgamma(n=1,shape=shape,scale=1/scale[j])
  }
  return(slopeCov)
}


genHRLDelta = function(stateVector,dataZ,prjCtrl,delta,prior){
  jump_size_index = loadedDie(prjCtrl$rwParamProbDelta)
  rw_sig = prjCtrl$rwParamRwSigDelta[jump_size_index]
  a_rate = prjCtrl$deltaARate
  p_delta = delta
  for(s in 2:prjCtrl$nS){
    p_delta[,s] = rmvnorm(n=1,mean=delta[,s],sigma = rw_sig*diag(prjCtrl$demoCOV))
  
  
    cLL = mnlLogLikeState(stateVector,dataZ,delta,prjCtrl) 
    pLL = mnlLogLikeState(stateVector,dataZ,p_delta,prjCtrl) 
    
  
    cPrior =  -0.5*t(delta[,s]-prior$delta)%*%prior$invDeltaCov%*%(delta[,s]-prior$delta)
    pPrior =  -0.5*t(p_delta[,s]-prior$delta)%*%prior$invDeltaCov%*%(p_delta[,s]-prior$delta)
  
    lap = pLL + pPrior - cLL - cPrior
    
    alpha = log(runif(1))
    if(alpha < lap){
      delta[,s] = p_delta[,s]
      a_rate[jump_size_index,1] = a_rate[jump_size_index,1] + 1
      
    }
    
    a_rate[jump_size_index,2] = a_rate[jump_size_index,2] + 1
  }
  return(list(delta,a_rate))
  
}


genHRLPhi = function(dataZ,delta){
  probs = array(0,c(prjCtrl$IND,prjCtrl$nS))
  eta = dataZ %*% delta
  
  for(i in 1:prjCtrl$IND){
    max_i = max(eta[i,])
    eta_i = eta[i,] - max_i
    probs[i,] = exp(eta_i)/sum(exp(eta_i))
  }
  
  return(probs)
}


genHRLSig2 = function(prjCtrl,data,param){
  sig2 = param$sig2
  for(k in 2:prjCtrl$nS){
    if(k < prjCtrl$nS){
      pSig2 = sig2[k] + ((sig2[k+1]-sig2[k-1])/6)*runif(1)
      while(pSig2 < sig2[k-1] || pSig2 > sig2[k+1]){
        pSig2 = sig2[k] + ((sig2[k+1]-sig2[k-1])/6)*runif(1)
      }
    }else{
      pSig2 = sig2[k] + ((sig2[k]-sig2[1])/(6*prjCtrl$nS))*runif(1)
      while(pSig2 < sig2[k-1]){
        pSig2 = sig2[k] + ((sig2[k]-sig2[1])/(6*prjCtrl$nS))*runif(1)
      }
    }
  
  
    cIndLL = array(0,dim=c(prjCtrl$IND,1))
    pIndLL = array(0,dim=c(prjCtrl$IND,1))
    cIndPrior = array(0,dim=c(prjCtrl$IND,1))
    pIndPrior = array(0,dim=c(prjCtrl$IND,1))
    
    for(i in 1:prjCtrl$IND){
      if(param$s[i,1] == k){
        cSlope = (1/sig2[k])*param$slope[,,i]
        pSlope = (1/pSig2)*param$slope[,,i]
        
        cIndLL[i,1] = mnlLogLike(data$y[,,i],data$X[,,,i],cSlope,prjCtrl$nRep,prjCtrl)
        pIndLL[i,1] = mnlLogLike(data$y[,,i],data$X[,,,i],pSlope,prjCtrl$nRep,prjCtrl)
        
        cIndPrior[i,1] = -0.5*t(param$slope[,,i]-param$slopeBar)%*%solve(param$slopeCov)%*%(param$slope[,,i]-param$slopeBar)
        cIndPrior[i,1] = cIndPrior[i,1] + 2.5*k*(param$nS[k]/prjCtrl$IND)*log(sig2[k]) - sig2[k]/100
        
        pIndPrior[i,1] = -0.5*t(pSig2*pSlope-param$slopeBar)%*%solve(param$slopeCov)%*%(pSig2*pSlope-param$slopeBar)
        pIndPrior[i,1] = pIndPrior[i,1] + 2.5*k*(param$nS[k]/prjCtrl$IND)*log(pSig2) - pSig2/100
        
      }
    }
    cLL = sum(cIndLL)
    pLL = sum(pIndLL)
    cPrior = sum(cIndPrior)
    pPrior = sum(pIndPrior)
    
    lap = pLL + pPrior - cLL - cPrior
    alpha = log(runif(1))
    if(alpha < lap){
      sig2[k] = pSig2
    }
    
  }
  return(sig2)
}


genHRLS = function(prjCtrl,data,param){
  s = array(0,dim=c(prjCtrl$IND,1))
  phi = param$phi
  
  for(i in 1:prjCtrl$IND){
    prob = array(0,prjCtrl$nS)
    diffSlope = param$slope[,,i] - param$slopeBar
    ssSlope = t(diffSlope)%*%solve(param$slopeCov)%*%diffSlope
    
    for(k in 1:prjCtrl$nS){
      cSlope = (1/param$sig2[k])*param$slope[,,i]
      prob[k] = mnlLogLike(data$y[,,i],data$X[,,,i],cSlope,prjCtrl$nRep,prjCtrl)
      prob[k] = prob[k] -.5*ssSlope + log(phi[i,k])
    }
    maxProb = max(prob)
    prob = exp(prob - maxProb)
    prob = prob/sum(prob)
    s[i,1] = loadedDie(prob)
    
  }

  
  return(s)
}


stabilityTest = function(trueParam,paramDraws){
  
  # Cross validate slope bars
  n_slope_bar = length(trueParam$slopeBar)
  slope_bar_in_range = array(0,dim=n_slope_bar)
  for(i in 1:n_slope_bar){
    slope_true_i = trueParam$slopeBar[i]
    slope_bar_i = paramDraws$slopeBar[i,,]
    posterior_range = quantile(slope_bar_i,probs=c(.025,.975))
    if((slope_true_i >= posterior_range[1]) && (slope_true_i <= posterior_range[2])){
      slope_bar_in_range[i] = TRUE
    }
    else{
      slope_bar_in_range[i] = FALSE
    }
  }
  
  print(paste("Percent Slope Bar In Range:",mean(slope_bar_in_range)))
  
  
  
  # cross validate slope i's
  n_slope_bar = dim(trueParam$slope)[1]
  ind = dim(trueParam$slope)[3]
  slope_bar_in_range = array(0,dim=c(n_slope_bar,ind))
  
  for(i in 1:n_slope_bar){
    for(j in 1:ind){
      slope_true_ij = trueParam$slope[i,,j]
      slope_ij = paramDraws$slope[i,,j,]
      posterior_range = quantile(slope_ij,probs=c(.025,.975))
      if((slope_true_ij >= posterior_range[1]) && (slope_true_ij <= posterior_range[2])){
        slope_bar_in_range[i,j] = TRUE
      }
      else{
        slope_bar_in_range[i,j] = FALSE
      }
    
    }  
  }
  
  
  print(paste("Percent Slope i In Range:",mean(slope_bar_in_range)))
  
  # cross validate delta
  n_delta = dim(trueParam$delta)[1]
  n_s = dim(trueParam$delta)[2]
  delta_in_range = array(0,dim=c(n_delta,n_s))
  
  for(i in 1:n_delta){
    for(j in 1:n_s){
      delta_true_ij = trueParam$delta[i,j]
      delta_ij = paramDraws$delta[i,j,]
      posterior_range = quantile(delta_ij,probs=c(.025,.975))
      if((delta_true_ij >= posterior_range[1]) && (delta_true_ij <= posterior_range[2])){
        delta_in_range[i,j] = TRUE
      }
      else{
        delta_in_range[i,j] = FALSE
      }
      
    }  
  }
  
  print(paste("Percent Delta In Range:",mean(delta_in_range)))

  # Cross validate slope covariance - diag only
  n_slope_bar = dim(trueParam$slopeCov)[1]
  slope_cov_in_range = array(0,dim=n_slope_bar)
  for(i in 1:n_slope_bar){
    slope_cov_true_i = trueParam$slopeCov[i,i]
    slope_cov_i = paramDraws$slopeCov[i,i,]
    posterior_range = quantile(slope_cov_i,probs=c(.025,.975))
    #print(c(slope_cov_true_i,posterior_range))
    if((slope_cov_true_i >= posterior_range[1]) && (slope_cov_true_i <= posterior_range[2])){
      slope_cov_in_range[i] = TRUE
    }
    else{
      slope_cov_in_range[i] = FALSE
    }
  }
  
  print(paste("Percent Slope Covariance (diag) In Range:",mean(slope_cov_in_range)))
  
  
}

plotConvergenceDelta = function(paramDraws){
  n_s = dim(paramDraws$delta)[2]
  for(s in 2:n_s){
    delta = mcmc(t(paramDraws$delta[,s,]))
    plot(delta,main=paste("Segment",s))
    
  }
}

plotConvergenceSlopeBar = function(paramDraws){
  slope_bar = mcmc(t(paramDraws$slopeBar[,1,]))
  plot(slope_bar)
  
}

plotConvergenceSlope = function(paramDraws){
  ind = dim(paramDraws$slope)[3]
  for(i in 1:ind){
    slope =  mcmc(t(paramDraws$slope[,1,i,]))
    plot(slope,main=paste("ind",i))
  }
  
}

plotConvergenceSlopeCovDiag = function(paramDraws){
  n_sample = dim(paramDraws$slopeCov)[3]
  n_slope_bar = dim(paramDraws$slopeCov)[2]
  cov_diag = array(0,dim=c(n_sample,n_slope_bar))
  
  for(i in 1:n_sample){
    for(j in 1:n_slope_bar){
      cov_diag[i,j] = paramDraws$slopeCov[j,j,i]
    }
  }
  
  cov_diag = mcmc(cov_diag)
  plot(cov_diag)
  
}