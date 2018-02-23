library(mvtnorm)

source("prjCtrl.r")
source("mcmcUtils.R")

main = function(){
    
  if(prjCtrl$genData){
    print("Generating synthetic data...")
  }
  
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

    param = genHRLParams(prjCtrl = prjCtrl,
                        param = param,
                        dataZ = dataZ,
                        prior = prior)
    tmp = 0
  }
  else{
    # Don't use true parameter values
    
    
  }
  
  
}




## Execute main MCMC routine
main()

