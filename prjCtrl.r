# Project control file for MCMC routine
# Create list called prjCtrl that specifies important MCMC meta parameters including: File I/O, Sample Size, Burn-in, model type, etc...

prjCtrl = list(

  projectName = 'tmp_r_debug_',
  
  runDelta = TRUE,
    
    
  nBurnin = 1000,
  nSample = 1000,
  useTrueParam = TRUE,
  randomStart = TRUE,
  genData = TRUE,
  plotDetail = TRUE,
  delayToPrint = FALSE,
  newSeed = TRUE,
  includeIntercept = FALSE,
    
  debug = TRUE,
  detailThin = 10,
  printToScreenThin = 10,
  printFigures = TRUE,
    
  PROD = 5,
  IND = 10,
  COV = 8,
  demoCOV = 10,
  nRep = 10,
  ridgeRegression = .1,
  fullCovar = FALSE,
  saveIndividualDraws = TRUE,
  useConstraints = FALSE,
  useExistingIndSlope = FALSE,
  saveRData = TRUE,
  
  startBurnin = 1000,
  startSample = 1000,
  
  rwParamProb = c(.55,.35,.1),
  rwParamRwSig = c(.1,.65,1.15),
  aRate =  matrix(0,3,2),
  nS = 4, # Number of mixture components
  probKMeansProposal = 0,
  probDecayKMeansProposal = .9,
  probGenS = .5,
  probGenSig2 = .5,

  inputFile = 'data.Rdata'


)

if(prjCtrl$runDelta)
  {
  # Delta I/O
  prjCtrl[["inputDemographicFile"]] = paste(prjCtrl$projectName,"tmp_demo.csv",sep = "")
  prjCtrl[["deltaDetail"]] = paste(prjCtrl$projectName,"deltaDetailFile.txt",sep = "")
  
  # Delta proposal random walk
  prjCtrl[["rwParamProbDelta"]] = c(.25,.5,25)
  prjCtrl[["rwParamRwSigDelta"]] = c(.1,.65,1.15)
  
  # Delta acceptance init
  prjCtrl[["deltaARate"]] = matrix(0,3,2)
  
}

prjCtrl[["indSlopeFile"]] = paste(prjCtrl$projectName,"indSlopeFile.csv",sep = "")
prjCtrl[["detailFile"]] = paste(prjCtrl$projectName,"HRLdetail.txt",sep = "")
prjCtrl[["logLikeFile"]] = paste(prjCtrl$projectName,"HRLLogLike.txt",sep = "")
prjCtrl[["reportFile"]] = paste(prjCtrl$projectName,"HRLreport.txt",sep = "")
prjCtrl[["paramFile"]] = paste(prjCtrl$projectName,"tmp_demo.RData",sep = "") # Save as R Data object
prjCtrl[["paramPostMeanFile"]] = paste(prjCtrl$projectName,"HRLparamPostMean.RData",sep = "") # Save as R Data object
prjCtrl[["paramSSFile"]] = paste(prjCtrl$projectName,"HRLparamSS.RData",sep = "") # Save as R data object
prjCtrl[["trueParamFile"]] = paste(prjCtrl$projectName,"HRLTrueParam.RData",sep = "")
prjCtrl[["trueParamFileTxt"]] = paste(prjCtrl$projectName,"HRLTrueParam.txt",sep = "")
prjCtrl[["priorFile"]] = paste(prjCtrl$projectName,"HRLPrior.RData",sep = "")
prjCtrl[["nParamFile"]] = paste(prjCtrl$projectName,"nParam.txt",sep = "")
prjCtrl[["indResultsFile"]] = paste(prjCtrl$projectName,"IndResults.txt",sep = "")

if(prjCtrl$useConstraints)
{
  prjCtrl[["paramConstraintsFile"]] = paste(prjCtrl$projectName,"ParamConstraints.txt",sep = "")
}





