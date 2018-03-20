source("mcmcUtils.R")

load(file='tmp_r_debug_HRLTrueParam.RData')
load(file='tmp_r_debug_paramDraws.RData')

# Perform stability test
stabilityTest(trueParam,paramDraws)

# Get convergence diagnostics
plotConvergenceSlopeCovDiag(paramDraws)
plotConvergenceSlope(paramDraws)
plotConvergenceSlopeBar(paramDraws)
plotConvergenceDelta(paramDraws)

