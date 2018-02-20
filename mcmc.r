source("prjCtrl.r")
source("mcmcUtils.R")

main = function(){
    
  if(genData){
    print("Generating synthetic data...")
  }
  
  # Hyper parameters
  prior = setHRLPrior(COV,demoCOV)
  if(saveRData){
    save(prior,file=priorFile)
  }
  
  # Set Parameters
  paramSS = NULL
  param = setHRLParams()
  

  
}




## Execute main MCMC routine
main()

