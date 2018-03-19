source("mcmcUtils.R")

load(file='tmp_r_debug_HRLTrueParam.RData')
load(file='tmp_r_debug_paramDraws.RData')

stabilityTest(trueParam,paramDraws)



