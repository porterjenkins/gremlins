
######################################################################
# Create the project control  class
projectControl <- setClass(
  # Set the name for the class
  "Project Control",
  
  # Define the slots
  slots = c(
    
    projectName = 'character',
    inputFile = 'character',
    runDelta = 'logical',
    
    nBurnin = "numeric",
    nSample = "numeric",
    useTrueParam = "logical",
    randomStart = "logical",
    genData = 'logical',
    plotDetail = 'logical',
    delayToPrint = 'logical',
    newSeed = 'logical',
    includeIntercept = 'logical',
    
    debug = 'logical',
    detailThin = 'numeric',
    printToScreenThing = 'numeric',
    printFigures = 'logical',
    
    PROD = 'numeric',
    IND = 'numeric',
    COV = 'numeric',
    demoCOV = 'numeric',
    maxRep = 'numeric',
    minRep = 'numeric',
    ridgeRegression = 'numeric',
    fullCovar = 'logical',
    saveIndividualDraws = 'logical',
    useConstraints = 'logical',
    
    inputDemographicFile = 'character',
    useExistingIndSlope = 'logical',
    indSlopeFile = 'character',
    detailFile = 'character',
    tauFile = 'character',
    loglikeFile = 'character',
    reportFile = 'character',
    paramFile = 'character',
    paramPostMeanFile = 'character',
    paramSSFile = 'character',
    trueParamFileTxt = 'character',
    priorFile = 'character',
    nParamFile = 'character',
    indResultsFile = 'character',
    paramConstraintsFile = 'character',
    deltaDetail = 'character',
    
    startBurnin = 'numeric',
    startSample = 'numeric',
    rwParamProb = 'vector',
    rwParamSig = 'vector',
    rwParamProbDelta = 'vector',
    rwParamSigDelta = 'vector'
    
  
    
  ),
  
  prototype = list(
    projectName = "NULL?]"
  )
  

  
)


setMethod("initialize","ANY",function(.Object,projectName){
  .Object@projectName = projectName
  return(.Object)
  
})



prjCtrl = projectControl()


# Set the default values for the slots. (optional)
tmp=list(
  projectName = 'tmp_r_debug_',
  inputFile = 'tmp_debug.txt',
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
  printToScreenThing = 10,
  printFigures = TRUE,
  
  PROD = 5,
  IND = 100,
  COV = 8,
  demoCOV = 10,
  maxRep = 100,
  minRep = 100,
  ridgeRegression = .1,
  fullCovar = FALSE,
  saveIndividualDraws = TRUE,
  useConstraints = FALSE,
  
  inputDemographicFile = paste(projectName,"tmp_demo")
  
  
  
  
)




