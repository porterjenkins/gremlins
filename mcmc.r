library(mvtnorm)
library(invgamma)

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
  print("Initializing parameters at random values...")
  # TODO: Implement John's MLE initialization routine
  #param = mleHRLStart(prjCtrl,data,param,prior)
  param = priorHRLStart(prjCtrl,param,prior)
  
}else{
  # No random parameter start
  if(prjCtrl$genData){
    print("Initializing parameters at true values...")
    paramTrue = param 
  }else{
    
  stop("genData = FALSE, randomStart = FALSE - TODO: write code to read in true values from file")
    
  }
}
  
print("Parameter start complete...")  
print(Sys.time() - param_init_start)
  


# BEGIN MCMC ROUTINE
mcmc_start_time = Sys.time()
print("Starting MCMC...")

for(n in 1:prjCtrl$nBurnin+prjCtrl$nSample){
  if(n %% prjCtrl$printToScreenThin == 0){
    print(paste("iteration:",n))
    #param$indLL = calcIndLogLike()
    #param$cLogLike[1] =
    #param$cLogLike[2] =
    #param$cLogLikeS =
    
    tSlopeBars = array(0,dim=c(prjCtrl$COV,prjCtrl$nS))
    for(i in 1:prjCtrl$IND){
      tSlopeBars[,param$s[i,1]] = tSlopeBars[,param$s[i,1]] + param$slope[,,i]
    }
    
    for(s in 1:prjCtrl$nS){
      tSlopeBars[,s] = tSlopeBars[,s] / param$nS[s,1]
    }
    
    

    
  }
}

  
  



print("MCMC complete")
print(Sys.time() - mcmc_start_time)
}




## Execute main MCMC routine
main()

