library(mvtnorm)

source("prjCtrl.r")
source("mcmcUtils.R")

main = function(){
  
  # INIT core data structures (parameters, paramSS, prior)
  
  # Set Parameters
  # Hyper parameters
  prior = setHRLPrior(COV=prjCtrl$COV,demoCOV=prjCtrl$demoCOV)
  paramSS = setHRLSS(prjCtrl$COV,prjCtrl$IND,prjCtrl$nS)
  param = setHRLParams(prjCtrl$COV,prjCtrl$IND,prjCtrl$nS,prjCtrl$useConstraints)
    
  if(prjCtrl$genData){
  
    if(prjCtrl$saveRData){
      save(prior,file=prjCtrl$priorFile)
    }
  
    
  
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
      if(prjCtrl$saveRData){
        save(param,file=prjCtrl$trueParamFile)
      }
      
    }
    else{
      # Don't use true parameter values
      load(prjCtrl$trueParamFile)
    }
    
    data = genHRLData(prjCtrl,param,dataZ)
    
    if(prjCtrl$saveRData){
      # Save generated data as object
      save(data,file=paste(prjCtrl$projectName,prjCtrl$inputFile,sep=''))
    }else{
      # TODO: re-structure data object as 'flat' file for .csv export
      
    }
    
  }else{
    # read input data from file
  
  }
# Parameter initialization
param_init_start = Sys.time()
print("Getting parameter starting values...")

if(prjCtrl$useConstraints && prjCtrl$genData == FALSE){
  param$slopeConstraintFlag = as.matrix(read.csv2(file = prjCtrl$paramConstraintsFile,header = FALSE))
  paramSS$slopeConstraintsFlag = param$slopeConstraintFlag
}

if(prjCtrl$randomStart){
  if(prjCtrl$debug && prjCtrl$genData){
   
    paramTrue = param 
    
  }
  
  param = mleHRLStart(prjCtrl,data,param,prior)

  
}else{
  # No random parameter start
}


  
  
print("Parameter start complete...")  
print(Sys.time() - param_init_start)
  



# BEGIN MCMC ROUTINE
mcmc_start_time = Sys.time()
print("Starting MCMC...")

  
  



print("MCMC complete")
print(Sys.time() - mcmc_start_time)
}




## Execute main MCMC routine
main()

