# Project control file for MCMC routine
# Variables below specify important MCMC meta parameters including: File I/O, Sample Size, Burn-in, model type, etc...


projectName = 'tmp_r_debug_'

runDelta = TRUE
  
  
nBurnin = 1000
nSample = 1000
useTrueParam = TRUE
randomStart = TRUE
genData = TRUE
plotDetail = TRUE
delayToPrint = FALSE
newSeed = TRUE
includeIntercept = FALSE
  
debug = TRUE
detailThin = 10
printToScreenThing = 10
printFigures = TRUE
  
PROD = 5
IND = 100
COV = 8
demoCOV = 10
maxRep = 100
minRep = 100
ridgeRegression = .1
fullCovar = FALSE
saveIndividualDraws = TRUE
useConstraints = FALSE
useExistingIndSlope = FALSE
saveRData = FALSE


inputFile = 'tmp_debug.txt'

if(runDelta)
  {
  # Delta I/O
  inputDemographicFile = paste(projectName,"tmp_demo.csv",sep = "")
  deltaDetail = paste(projectName,"deltaDetailFile.txt",sep = "")
  
  # Delta proposal random walk
  rwParamProbDelta = c(.25,.5,25)
  rwParamRwSigDelta = c(.1,.65,1.15)
  
  # Delta acceptance init
  deltaARate = matrix(0,3,2)
  
}

indSlopeFile = paste(projectName,"indSlopeFile.csv",sep = "")
detailFile = paste(projectName,"HRLdetail.txt",sep = "")
logLikeFile = paste(projectName,"HRLLogLike.txt",sep = "")
reportFile = paste(projectName,"HRLreport.txt",sep = "")
paramFile = paste(projectName,"tmp_demo.RData",sep = "") # Save as R Data object
paramPostMeanFile = paste(projectName,"HRLparamPostMean.RData",sep = "") # Save as R Data object
paramSSFile = paste(projectName,"HRLparamSS.RData",sep = "") # Save as R data object
trueParamFile = paste(projectName,"HRLTrueParam.RData",sep = "")
trueParamFileTxt = paste(projectName,"HRLTrueParam.txt",sep = "")
priorFile = paste(projectName,"HRLPrior.RData",sep = "")
nParamFile = paste(projectName,"nParam.txt",sep = "")
indResultsFile = paste(projectName,"IndResults.txt",sep = "")

if(useConstraints)
{
  paramConstraintsFile = paste(projectName,"ParamConstraints.txt",sep = "")
}

startBurnin = 1000
startSample = 1000

rwParamProb = c(.55,.25,.1)
rwParamRwSig = c(.1,.65,1.15)
aRate =  matrix(0,3,2)
nS = 2 # Number of mixture components
probKMeansProposal = 0
probDecayKMeansProposal = .9
probGenS = .5
probGenSig2 = .5





