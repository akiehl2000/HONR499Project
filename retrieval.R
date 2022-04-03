model_name <- 'INGARCH_NBinom_WR'

coefs <- readRDS(paste('./Optimize/', model_name, '/coefs.RData', sep = ''))

predictors <- readRDS(paste('./Optimize/', model_name, '/predictors.RData', sep = ''))
test_mses <- readRDS(paste('./Optimize/', model_name, '/test_mses.RData', sep = ''))
tuned <- readRDS(paste('./Optimize/', model_name, '/tuned.RData', sep = ''))
if (grepl('GLARMA', model_name)) {
  modelTuneAIC <- readRDS(paste('./Optimize/', model_name, '/modelTuneAIC.RData', sep = ''))
  varSelectAIC <- readRDS(paste('./Optimize/', model_name, '/varSelectAIC.RData', sep = ''))
} else {
  modelTuneMSE <- readRDS(paste('./Optimize/', model_name, '/modelTuneMSE.RData', sep = ''))
  varSelectMSE <- readRDS(paste('./Optimize/', model_name, '/varSelectMSE.RData', sep = ''))
}
