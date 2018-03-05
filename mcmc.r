library(mvtnorm)

source("prjCtrl.r")
source("mcmcUtils.R")

main = function(){
    
  if(prjCtrl$genData){
  
  # Hyper parameters
  prior = setHRLPrior(COV=prjCtrl$COV,demoCOV=prjCtrl$demoCOV)
  if(prjCtrl$saveRData){
    save(prior,file=priorFile)
  }
  
  # Set Parameters
  paramSS = setHRLSS(prjCtrl$COV,prjCtrl$IND,prjCtrl$nS)
  param = setHRLParams(prjCtrl$COV,prjCtrl$IND,prjCtrl$nS,prjCtrl$useConstraints)
  

  if(prjCtrl$useTrueParam){
    # Use true parameter values
    
    if(prjCtrl$runDelta){
      # Generate Z matrix for delta-level regression
      dataZ = rmvnorm(n=prjCtrl$IND,
                      mean = rep(0,prjCtrl$demoCOV),
                      sigma = diag(prjCtrl$demoCOV))
    }
    else{
      dataZ = NULL
    }

    param = genHRLParams(prjCtrl = prjCtrl,
                        param = param,
                        dataZ = dataZ,
                        prior = prior)
    if(prjCtrl$useConstraints){
      paramSS$slopeConstraintFlag = param$slopeConstraintFlag
      write.table(param$slopeConstraintFlag,file=prjCtrl$paramConstraintsFile,col.names = FALSE,row.names=FALSE)
    }
    for(i in 1:prjCtrl$IND){
      param$slopeOverSig2[,,i] = (1/param$sig2[param$s[i]])*param$slope[,,i]
    }
    save(param,file=prjCtrl$trueParamFile)
    
  }
  else{
    # Don't use true parameter values
    load(prjCtrl$trueParamFile)
  }
  
  data = genHRLData(prjCtrl,param,dataZ)
  }  
}




## Execute main MCMC routine
main()

