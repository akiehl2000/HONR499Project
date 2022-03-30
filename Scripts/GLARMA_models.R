# DOCUMENTATION HERE

# SETUP

# install packages
library(tidyverse)
install.packages('glarma', repos = 'http://cran.us.r-project.org')
library(glarma)

# set seed for reproducibility
set.seed(499)

# read in binned data frame
binned_df <- read.csv('../Data/binned_df.csv')

# neuron names and unique trials
pred_names <- c('CA3.1', 'CA3.2', 'CA3.3', 'CA3.4', 'CA3.5', 'CA3.6', 'CA3.7', 
                'CA3.8', 'CA3.9', 'CA3.10', 'CA3.11', 'CA3.12', 'CA3.13', 
                'CA3.14', 'CA3.15', 'CA3.16', 'CA3.17', 'CA3.18', 'CA3.19', 
                'CA3.20', 'CA3.21')
resp_names <- c('CA1.1', 'CA1.2', 'CA1.3', 'CA1.4', 'CA1.5', 'CA1.6', 'CA1.7', 
                'CA1.8', 'CA1.9', 'CA1.10', 'CA1.11', 'CA1.12', 'CA1.13', 
                'CA1.14', 'CA1.15')
trls <- unique(binned_df$trial)

# data selection function
gen_df <- function(train_trials, valid_trials, pred, resp) {
  # train_trails: a vector of trials to be included in the training set
  # valid_trials: a vector of trials to be included in the validation set
  # pred: a vector of predictor variables to be included in the data sets
  # resp: a vector of response variables to be included in the data sets
  
  # generate training set
  train_df <- binned_df %>%
    filter(trial %in% as.vector(train_trials)) %>%
    select(c(resp, pred, trial)) %>%
    data.frame()
  
  # generate validation set
  valid_df <- binned_df %>%
    filter(trial %in% as.vector(valid_trials)) %>%
    select(c(pred, trial)) %>%
    data.frame()
  
  # values to compare predicted values with
  valid_pred <- binned_df %>%
    filter(trial %in% as.vector(valid_trials)) %>%
    select(resp)
  
  # select only non-zero predictor vectors
  zeros <- c()
  zeros <- c(zeros, which(colSums(train_df) == 0))
  zeros <- c(zeros, which(colSums(valid_df) == 0))
  if (length(zeros) > 0) {
    train_df <- train_df[, -zeros]
    valid_df <- valid_df[, -zeros]
    
    print('Removed empty column(s)')
  }
  
  # return data frames
  return(list('train_df' = train_df, 
              'valid_df' = valid_df, 
              'valid_pred' = valid_pred))
}

# cross validation function
cross_valid <- function(pred_list, resp, p, q, dist) {
  # pred_list: a vector of predictor variables to be included in the data sets
  # resp_list: a vector of response variables to be included in the data sets
  # p: the number of past observations to model on
  # q: the number of lags to model the conditional mean on
  # dist: the assumed distribution for the marginal time series process
  
  # define empty vector to store MSE values
  aics <- rep(NA, 4)
  
  for (i in 1:4) {
    # blocking cross-validation method is progressive
    trials <- list((10 * (i - 1) + 1):(10 * i - 2), (10 * i - 1):(10 * i))
    
    # collect data frames
    dfs <- gen_df(trls[trials[[1]]], trls[trials[[2]]], pred_list, resp)
    train_df <- dfs[['train_df']]
    valid_df <- dfs[['valid_df']]
    valid_pred <- unlist(dfs[['valid_pred']])
    
    remove(dfs)
    
    tryCatch(
      {
        # fit GLARMA model
        model <- glarma(train_df[, 1], 
                        as.matrix(train_df[, -1]),
                        phiLags = 1:p,
                        thetaLags = q, 
                        type = dist)
      }, 
      error = function(err) {
        print(paste('Model fitting error on ', resp, 
                    ' (', p, ', ', q, '): ', err, sep = ''))
      }
    )
    
    tryCatch(
      {
        # extract AIC
        aics[i] <- round(model$aic, 4)
      }, 
      error = function(err) {
        print(paste('Model evaluation error on ', resp, 
                    ' (', p, ', ', q, '): ', err, sep = ''))
      }
      
    )
  }
  
  return(mean(aics, na.rm = TRUE))
}

# variable selection function
var_select <- function(select, start, dist, type) {
  # select: the maximum number of predictor variables to model on
  # start: the start time of the model tuning procedure
  # dist: the assumed distribution for the marginal time series process
  # type: 'ir' or 'wr' for inter-region or within-region fitting
  
  working <- TRUE
  
  # define empty list to store best variable subsets
  predictors <<- list()
  var_select_aics <<- list()
  
  for (resp in resp_names) {
    print(paste('Beginning variable selection for:', resp))
    
    working <- TRUE
    
    # collect training data frame
    if (type == 'ir') {
      dfs <- gen_df(trls[1:5], trls[6:10], pred_names, resp)
      train_df <- dfs[['train_df']]
      remove(dfs)
    } else {
      dfs <- gen_df(trls[1:5], trls[6:10], 
                    resp_names[-which(resp_names == resp)], resp)
      train_df <- dfs[['train_df']]
      remove(dfs)
    }
    
    tryCatch(
      {
        # fit full model
        model <- glarma(train_df[, 1],
                        as.matrix(train_df[, -1]), 
                        phiLags = 1:5,
                        thetaLags = 6, 
                        type = dist)
      }, 
      error = function(err) {
        working <- FALSE
        
        print(paste('Model fitting error for', resp, 'full model:', err))
      }
    )
    
    tryCatch(
      {
        # extract estimated coefficients from fitted model
        model_sum <- summary(model)$coefficients1$Estimate
        
        if (dist == 'Poi') {
          est <- model_sum[-length(model_sum)]
        } else{
          est <- model_sum[-c(1:6, length(model_sum) - 1, length(model_sum))]
        }
        
        if (type == 'ir') {
          full_model_select <- data.frame(var = pred_names, coef = est) %>%
            arrange(desc(abs(coef))) %>%
            top_n(select, coef) %>%
            select(var) %>%
            unlist()
        } else {
          full_model_select <- data.frame(var = resp_names[-which(resp_names == resp)], coef = est) %>%
            arrange(desc(abs(coef))) %>%
            top_n(select, coef) %>%
            select(var) %>%
            unlist()
        }
      }, 
      error = function(err) {
        working <- FALSE
        
        print(paste('Model evaluation error for', resp, 'full model:', err))
      }
    )
    
    if (working == TRUE) {
      # define empty vector to store CV AIC results for forward selection
      aics <- rep(NA, select)
      
      # perform forward selection with chosen variables
      for (i in 1:length(full_model_select)) {
        print(paste('selecting top', i, 'variables...'))
        
        preds <- full_model_select[1:i]
        
        aics[i] <- cross_valid(preds, resp, 5, 6, dist)
      }
      
      # store selected best subset 
      forward_select <- data.frame(var_inc = 1:select,
                                   aic = aics)
      predictors[[resp]] <<- full_model_select[1:forward_select$var_inc[which(
        forward_select$aic == min(forward_select$aic))]]
      var_select_aics[[resp]] <<- as.numeric(aics)
    } else {
      print(paste('Could not perform variable selection for', resp))
    }
    
    print(paste('Completed after', difftime(Sys.time(), start, unit = 'mins'), 
                'minutes'))
  }
  
  return(predictors)
}

# parameter tuning function
model_tune <- function(Pmax, Qmax, predictors, start, dist) {
  # Pmax: the maximum number of past observations to model on
  # Qmax: the maximum lagged conditional mean to model on
  # predictors: a list of optimal predictor sets for each response
  # start: the start time of the model tuning procedure
  # dist: the assumed distribution for the marginal time series process
  
  # define empty data frame for storing tuned values
  tuned <<- data.frame(resp = resp_names, 
                       p = rep(NA, length(resp_names)),
                       q = rep(NA, length(resp_names)))
  model_tune_aics <<- list()
  
  # loop over response variables
  for (resp in resp_names) {
    print(paste('Beginning model tuning for:', resp))
    
    # define empty data frame to store MSE results
    tuning <- data.frame(p = rep(1:Pmax, each = Qmax),
                         q = rep(1:Qmax, times = Pmax),
                         aic = rep(NA, (Pmax * Qmax)))
    
    # loop over parameter values of interest
    for (p in 1:Pmax) {
      for (q in 1:Qmax) {
        print(paste('tuning with p=', p, ' & q=', q, '...', sep=''))
        
        if (q > p) {
          # calculate model MSE with given parameters
          if (is.null(predictors[[resp]])) {
            tuning$aic[which(tuning$p == p & tuning$q == q)] <- cross_valid(
              pred_names, resp, p, q, dist)
          } else{
            tuning$aic[which(tuning$p == p & tuning$q == q)] <- cross_valid(
              as.vector(predictors[[resp]]), resp, p, q, dist)
          }
        } 
      }
    }
    
    # save optimal tuned model parameters
    tuning <- tuning %>%
      filter(!is.na(aic))
    tuned_p <- as.numeric(tuning$p[which(tuning$aic == min(tuning$aic, 
                                                           na.rm = TRUE))])
    tuned_q <- as.numeric(tuning$q[which(tuning$aic == min(tuning$aic, 
                                                           na.rm = TRUE))])
    tuned[which(tuned$resp == resp), ] <<- c(resp, tuned_p, tuned_q)
    model_tune_aics[[resp]] <<- as.numeric(tuning$aic)
    
    print(paste('Completed after', difftime(Sys.time(), start, unit = 'mins'), 
                'minutes'))
  }
  
  return(tuned)
}

# optimal model fitting function
fitOpt <- function(predictors, tuned, dist, type) {
  # predictors: a list of optimal predictor sets for each response
  # tuned: a data frame of optimal model parameters for each response
  # dist: the assumed distribution for the marginal time series process
  # type: 'ir' or 'wr' for inter-region or within-region fitting
  
  test_mses <- data.frame(resp = resp_names,
                          mse = rep(NA, length(resp_names)))
  
  for (resp in resp_names) {
    # collect data frames
    dfs <- gen_df(trls[1:40], trls[41:46], predictors[[resp]], resp)
    train_df <- dfs[['train_df']]
    valid_df <- dfs[['valid_df']]
    valid_pred <- dfs[['valid_pred']]
    remove(dfs)
    
    # collect optimized parameter values
    p <- as.numeric(tuned$p[which(tuned$resp == resp)])
    q <- as.numeric(tuned$q[which(tuned$resp == resp)])
    
    print(paste('Fitting optimal model for ', resp, '...', sep = ''))
    
    # fit optimized model
    tryCatch(
      {
        model <- glarma(train_df[, 1],
                        as.matrix(train_df[, -1]),
                        phiLags = 1:p,
                        thetaLags = q, 
                        type = dist)
        
        
        if (dist == 'Poi' & type == 'ir') {
          model %>%
            saveRDS(paste('../Models/GLARMA_Pois_IR/model_', resp, '.rds', 
                          sep = ''))
        } else if (dist == 'Poi' & type == 'wr') {
          model %>%
            saveRDS(paste('../Models/GLARMA_Pois_WR/model_', resp, '.rds', 
                          sep = ''))
        } else if (dist == 'NegBin' & type == 'ir') {
          model %>%
            saveRDS(paste('../Models/GLARMA_NBinom_IR/model_', resp, '.rds', 
                          sep = ''))    
        } else {
          model %>%
            saveRDS(paste('../Models/GLARMA_NBinom_WR/model_', resp, '.rds', 
                          sep = ''))    
        }
      }, 
      error = function(err) {
        print(paste('Model fitting error on optimal ', resp, ': ', err, 
                    sep = ''))
      }
    )
    
    # calculate optimized test MSE
    tryCatch(
      {
        R <- 100
        
        y_pred <- rep(NA, nrow(valid_df))
        
        model <- glarma(train_df[, 1], 
                        as.matrix(train_df[, -1]), 
                        phiLags = 1:p, 
                        thetaLags = q,
                        type = dist)
        
        for (t in 1:nrow(valid_df)) {
          newYs <- rep(NA, R)
          for (r in 1:R) {
            newYs[r] <- forecast(model, n.ahead = 1, newdata = as.matrix(valid_df[t, ]))$Y
          }
          y_pred[t] <- mean(newYs)
          
          print(paste('Forecasting complete for time ', t, '...', sep = ''))
          
          model <- glarma(c(train_df[, 1], valid_pred[t, 1]), 
                          as.matrix(rbind(train_df[, -1], valid_df[t, ])), 
                          phiLags = 1:5, 
                          thetaLags = 6,
                          type = dist)
        }
        
        mse <- mean((y_pred - valid_pred[, 1])^2, na.rm = TRUE)
        test_mses$mse[which(test_mses$resp == resp)] <- mse
      }, 
      error = function(err) {
        print(paste('Model evaluation error on optimal ', resp, ': ', err, 
                    sep = ''))
      }
    )
  }
  
  return(test_mses)
}

# coefficient collection function
getCoefs <- function(predictors, dist, type) {
  # predictors: a list of optimal predictor sets for each response
  # dist: the assumed distribution for the marginal time series process
  # type: 'ir' or 'wr' for inter-region or within-region fitting
  
  working <- TRUE 
  
  # define empty data frame to store coefficients
  coefs <- data.frame(resp = resp_names,
                      phi_1 = rep(NA, 15),
                      phi_2 = rep(NA, 15),
                      phi_3 = rep(NA, 15),
                      phi_4 = rep(NA, 15),
                      phi_5 = rep(NA, 15),
                      theta = rep(NA, 15),
                      CA3.1 = rep(NA, 15),
                      CA3.2 = rep(NA, 15),
                      CA3.3 = rep(NA, 15),
                      CA3.4 = rep(NA, 15),
                      CA3.5 = rep(NA, 15),
                      CA3.6 = rep(NA, 15),
                      CA3.7 = rep(NA, 15),
                      CA3.8 = rep(NA, 15),
                      CA3.9 = rep(NA, 15),
                      CA3.10 = rep(NA, 15),
                      CA3.11 = rep(NA, 15),
                      CA3.12 = rep(NA, 15),
                      CA3.13 = rep(NA, 15),
                      CA3.14 = rep(NA, 15),
                      CA3.15 = rep(NA, 15),
                      CA3.16 = rep(NA, 15),
                      CA3.17 = rep(NA, 15),
                      CA3.18 = rep(NA, 15),
                      CA3.19 = rep(NA, 15),
                      CA3.20 = rep(NA, 15),
                      CA3.21 = rep(NA, 15),
                      trial = rep(NA, 15))
  
  # loop through response variables
  for (resp in resp_names) {
    working <- TRUE
    
    # import trained optimal models
    tryCatch(
      {
        if (dist == 'Poi' & type == 'ir') {
          model <- readRDS(paste('../Models/GLARMA_Pois_IR/model_', resp, 
                                 '.rds', sep = ''))
        } else if (dist == 'Poi' & type =='wr') {
          model <- readRDS(paste('../Models/GLARMA_Pois_WR/model_', resp, 
                                 '.rds', sep = ''))
        } else if (dist == 'NegBin' & type == 'ir') {
          model <- readRDS(paste('../Models/GLARMA_NBinom_IR/model_', resp, 
                                 '.rds', sep = ''))
        } else {
          model <- readRDS(paste('../Models/GLARMA_NBinom_WR/model_', resp, 
                                 '.rds', sep = ''))
        }
      }, 
      error = function(err) {
        print(paste('Could not import trained model for', resp))
      }
    )
    
    # collect model coefficient estimates
    if (working == TRUE & !is.null(predictors[[resp]])) {
      model_sum <- summary(model)
      rowNms <- c(rownames(model_sum$coefficients2),
                  rownames(model_sum$coefficients1))
      est <- c(model_sum$coefficients2$Estimate, 
               model_sum$coefficients1$Estimate)
      
      nms <- rowNms[which(
        !grepl('var', rowNms) & 
          !grepl('Intercept', rowNms) & 
          !grepl('sigma', rowNms) &
          !grepl('trial', rowNms))]
      nms[which(grepl('theta', nms))] <- 'theta'
      nms <- c(nms, predictors[[resp]], 'trial')
      
      coefs[which(coefs$resp == resp), nms] <- est
    }
  }
  
  if (dist == 'Poi' & type == 'ir') {
    coefs %>%
      saveRDS('../Optimize/GLARMA_Pois_IR/coefs.RData')
  } else if (dist == 'Poi' & type == 'wr') {
    coefs %>%
      saveRDS('../Optimize/GLARMA_Pois_WR/coefs.RData')
  } else if (dist == 'NegBin' & type == 'ir') {
    coefs %>%
      saveRDS('../Optimize/GLARMA_NBinom_IR/coefs.RData')
  } else {
    coefs %>%
      saveRDS('../Optimize/GLARMA_NBinom_WR/coefs.RData')
  }
  
  return(coefs)
}

# MODEL FITTING

# fit Poisson inter-region models
print('Fitting Poisson inter-region models...')
tryCatch(
  {
    dist <- 'Poi'
    type <- 'ir'
    
    start <- Sys.time()
    
    # collect optimal predictor subsets
    predictors <- var_select(5, start, dist, type)
    
    # find optimal model parameters
    tuned <- model_tune(5, 5, predictors, start, dist)
    
    # save results
    predictors %>%
      saveRDS('../Optimize/GLARMA_Pois_IR/predictors.RData')
    tuned %>%
      saveRDS('../Optimize/GLARMA_Pois_IR/tuned.RData')
    var_select_aics %>%
      saveRDS('../Optimize/GLARMA_Pois_IR/varSelectAIC.RData')
    model_tune_aics %>%
      saveRDS('../Optimize/GLARMA_Pois_IR/modelTuneAIC.RData')
    
    test_mses <- fitOpt(predictors, tuned, dist, type)
    
    test_mses %>%
      saveRDS('../Optimize/GLARMA_Pois_IR/test_mses.RData')
    
    coefs <- getCoefs(predictors, dist, type)
    
    print(paste('Completed after', 
                difftime(Sys.time(), start, unit = 'mins'), 'minutes'))
  }, error = function(err) {
    print(paste('Error fitting Poisson inter-region models:', err))
  }
)

# fit Poisson within-region models
print('Fitting Poisson within-region models...')
tryCatch(
  {
    type <- 'wr'
    
    start <- Sys.time()
    
    # collect optimal predictor subsets
    predictors <- var_select(5, start, dist, type)
    
    # find optimal model parameters
    tuned <- model_tune(5, 5, predictors, start, dist)
    
    # save results
    predictors %>%
      saveRDS('../Optimize/GLARMA_Pois_WR/predictors.RData')
    tuned %>%
      saveRDS('../Optimize/GLARMA_Pois_WR/tuned.RData')
    var_select_aics %>%
      saveRDS('../Optimize/GLARMA_Pois_WR/varSelectAIC.RData')
    model_tune_aics %>%
      saveRDS('../Optimize/GLARMA_Pois_WR/modelTuneAIC.RData')
    
    test_mses <- fitOpt(predictors, tuned, dist, type)
    
    test_mses %>%
      saveRDS('../Optimize/GLARMA_Pois_WR/test_mses.RData')
    
    coefs <- getCoefs(predictors, dist, type)
    
    print(paste('Completed after', 
                difftime(Sys.time(), start, unit = 'mins'), 'minutes'))
  }, error = function(err) {
    print(paste('Error fitting Poisson within-region models:', err))
  }
)

# # fit Negative Binomial inter-region models
# print('Fitting Negative Binomial inter-region models...')
# tryCatch(
#   {
#     dist <- 'NegBin'
#     type <- 'ir'
#     
#     start <- Sys.time()
#     
#     # collect optimal predictor subsets
#     predictors <- var_select(5, start, dist, type)
#     
#     # find optimal model parameters
#     tuned <- model_tune(5, 5, predictors, start, dist)
#     
#     # save results
#     predictors %>%
#       saveRDS('../Optimize/GLARMA_NBinom_IR/predictors.RData')
#     tuned %>%
#       saveRDS('../Optimize/GLARMA_NBinom_IR/tuned.RData')
#     var_select_aics %>%
#       saveRDS('../Optimize/GLARMA_NBinom_IR/varSelectAIC.RData')
#     model_tune_aics %>%
#       saveRDS('../Optimize/GLARMA_NBinom_IR/modelTuneAIC.RData')
#     
#     test_mses <- fitOpt(predictors, tuned, dist, type)
#     
#     test_mses %>%
#       saveRDS('../Optimize/GLARMA_NBinom_IR/test_mses.RData')
#     
#     coefs <- getCoefs(predictors, dist, type)
#     
#     print(paste('Completed after', 
#                 difftime(Sys.time(), start, unit = 'mins'), 'minutes'))
#   }, error = function(err) {
#     print(paste('Error fitting Negative Binomial inter-region models:', err))
#   }
# )
# 
# # fit Negative Binomial within-region models
# print('Fitting Negative Binomial within-region models...')
# tryCatch(
#   {
#     type <- 'wr'
#     
#     start <- Sys.time()
#     
#     # collect optimal predictor subsets
#     predictors <- var_select(5, start, dist, type)
#     
#     # find optimal model parameters
#     tuned <- model_tune(5, 5, predictors, start, dist)
#     
#     # save results
#     predictors %>%
#       saveRDS('../Optimize/GLARMA_NBinom_WR/predictors.RData')
#     tuned %>%
#       saveRDS('../Optimize/GLARMA_NBinom_WR/tuned.RData')
#     var_select_aics %>%
#       saveRDS('../Optimize/GLARMA_NBinom_WR/varSelectAIC.RData')
#     model_tune_aics %>%
#       saveRDS('../Optimize/GLARMA_NBinom_WR/modelTuneAIC.RData')
#     
#     test_mses <- fitOpt(predictors, tuned, dist, type)
#     
#     test_mses %>%
#       saveRDS('../Optimize/GLARMA_NBinom_WR/test_mses.RData')
#     
#     coefs <- getCoefs(predictors, dist, type)
#     
#     print(paste('Completed after', 
#                 difftime(Sys.time(), start, unit = 'mins'), 'minutes'))
#   }, error = function(err) {
#     print(paste('Error fitting Negative Binomial within-region models:', err))
#   }
# )
