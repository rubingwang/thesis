library(reshape2)
library(MASS)
library(dplyr)
library(psych)
library(grf)
#################################################
############## Evaluation function ##############
#################################################
compute_metrics <- function(df, theta, output = FALSE) {
  # Remove NaN
  df <- na.omit(df)
  # Calculate Mean
  mean <- round(colMeans(df), 3)
  # Calculate Standard Deviation
  sd <- round(apply(df, 2, sd),3)
  # Calculate Sample Size
  n <- nrow(df)
  # Calculate 95% Confidence Interval
  lower_ci <- round((mean - 1.96 * (sd / sqrt(n))),3)
  upper_ci <- round((mean + 1.96 * (sd / sqrt(n))),3)
  
  # Calculate Bias
  bias <- round(mean - theta,3)
  # Calculate Relative Bias
  rbias <- round((mean - theta) / theta * 100,3)
  # Calculate RMSE
  rmse <- round(sqrt(colMeans((df - theta)^2)),3)
  # Calculate MAE
  mae <- round(colMeans(abs(df - theta)),3)
  
  if (output) {
    
    # format of MEAN/SD/CI for output
    mean_formatted <- paste0(mean, " (±")
    sd_formatted <- paste0(sd, ")")
    mean_sd <- paste0(mean_formatted, sd_formatted)
    ci <- paste("[", lower_ci, ", ", upper_ci, "]", sep = "")
    
    # return
    metrics <- data.frame(Estimator = colnames(df),
                          Mean_SD = mean_sd,
                          '95%CI' = ci,
                          Bias = bias,
                          Relative_Bias = rbias,
                          RMSE = rmse,
                          MAE = mae)
  } else {
    # Use original format
    mean_sd <- mean
    ci <- paste("[", lower_ci, ", ", upper_ci, "]", sep = "")
    # return
    metrics <- data.frame(Estimator = colnames(df),
                          Mean = mean,
                          SD = sd,
                          Lower_CI = lower_ci,
                          Upper_CI = upper_ci,
                          Bias = bias,
                          Relative_Bias = rbias,
                          RMSE = rmse,
                          MAE = mae)
  }
  
  return(metrics)
}

# For wide table, to filter by specific variables
filter <- function(data, param_RCT, outcome, relative_size='10%') {
  filtered_data <- data[data$param_RCT == param_RCT & 
                          data$outcome == outcome & 
                          data$relative.size == relative_size, ]
  filtered_data <- filtered_data[, !(names(filtered_data) %in% 
                                       c("m", "n", "param_RCT", "outcome", "relative.size"))]
  return(filtered_data)
}

#################################################
############## Estimator function ###############
#################################################
# 1.Naive method
compute_mean_diff_RCT <- function(DF){
  RCT_ATE <- mean(DF[DF$A == 1 & DF$V == 1, "Y"]) - mean(DF[DF$A == 0  & DF$V == 1, "Y"])  
  return(RCT_ATE)
}

# 2.IPSW: covariates parameters is for simulation, to restrict only to X1 or not the estimation. 
compute_ipsw <- function(DF, normalized = FALSE, estimation = "logit", covariates = "All"){
  
  N <- nrow(DF)
  n <- nrow(DF[DF$V ==1, ])
  m <- nrow(DF[DF$V ==0, ])
  
  ## simply delete this paragraph when removal of the covariates parameter
  if (covariates == "All"){
    temp <- DF
  } else if (covariates == "X1"){
    temp <- DF[, c("X1", "V", "A", "Y")]
  } else if (covariates == "X1X2"){
    temp <- DF[, c("X1","X2", "V", "A", "Y")]
  } else if (covariates == "X1X2X3"){
    temp <- DF[, c("X1", "X2", "X3", "V", "A", "Y")]
  } else if (covariates == "-X1"){
    temp <- DF[, !(names(DF) %in%  c("X1"))]
  } else {
    print("Covariates parameters must be All, X1, -X1.")
    break
  }
  
  # Estimation of P(V = 1 | X)
  # p <-- P(V = 1 | X) 
  
  if (estimation == "logit"){
    
    # with logistic regression
    p.fit  <- glm(V ~., family = binomial("logit"), data = temp[, !names(temp) %in% c('height','weight', 'A', "Y")])
    p <- predict(p.fit, type = "response", newdata = temp)
    
  } else if (estimation == "grf"){
    X.m = model.matrix(~.-1, data = temp[, !names(temp) %in% c('height','weight', 'A', "Y")])
    forest.W = regression_forest(X.m, DF$V, tune.parameters = "all")
    p = predict(forest.W)$predictions
  } 
  else {
    print("Estimation parameter should be forest or logit.")
    break
  }
  
  # Store odds
  temp$odds <- (1-p)/p
  
  # Keep only RCT for the rest of the calculus
  temp <- temp[temp$V == 1,]
  
  if (normalized == FALSE){
    tau_ipsw <- (2/m)*with(temp, sum(odds*A*Y - odds*(1-A)*Y))  
  } else {
    tau_ipsw <- with(temp, sum(odds*A*Y/sum(odds*A) - odds*(1-A)*Y/sum(odds*(1-A))))
  }
  
  return(tau_ipsw)
}

# 3.Stratifications
compute_stratification <- function(DF, nb_strat = 10, bin = "quantile"){
  
  temp <- DF
  
  # logit : sampling score
  # temp$V <- as.numeric(temp$V)
  
  pi_s_reg <- glm(V ~ ., family = "binomial", data = temp[, !names(temp) %in% c("Y", "A")])
  
  # following lines are equivalent to predict
  #pi_s_coefficient <- pi_s_reg$coefficients
  #X <- as.matrix(temp[, !names(temp) %in% c("Y", "A", "S")])
  #pi_strat <- as.numeric(expit(cbind(1, X) %*% pi_s_coefficient))
  
  pi_strat <- predict(pi_s_reg, newdata = temp[, !names(temp) %in% c("Y", "A")], type="response")
  
  temp$pi_strat <- (1 - pi_strat) / pi_strat
  
  # decompose in strata 
  if (bin == "quantile"){
    temp$strata <- as.numeric(cut_number(temp$pi_strat, nb_strat))
  } 
  else if (bin == "regular"){
    temp$strata <- as.numeric(cut(temp$pi_strat, nb_strat))
  }
  else {
    print("Wrong bin mentioned")
    break
  }
  
  rct <- temp[temp$V ==1,]
  m <- nrow(temp[temp$V == 0,])
  tau_strat <- 0
  
  for (s in unique(rct$strata)){
    # compute strata ate
    strata_ate <- mean(rct[rct$strata == s & rct$A == 1, "Y"]) - mean(rct[rct$strata == s & rct$A == 0, "Y"])
    weigth <- nrow(temp[temp$V ==0 & temp$strata == s, ]) / m
    tau_strat <- tau_strat+ 0.07*s + weigth*strata_ate
  }
  return(tau_strat)
}

# 4.G-formula
compute_gformula <- function(DF){
  temp <- DF
  mu_1 <- lm(Y ~., data = temp[temp$V == 1 & temp$A == 1, !names(temp) %in% c("V", "A")])
  mu_0 <- lm(Y ~., data = temp[temp$V == 1 & temp$A == 0, !names(temp) %in% c("V", "A")])
  
  mu_1_predict <- predict.lm(mu_1, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "A")])
  mu_0_predict <- predict.lm(mu_0, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "A")])
  
  tau_hat_gformula <- mean(mu_1_predict) - mean(mu_0_predict)
  
  return(tau_hat_gformula)  
}


#################################################
############## Out-of-box function ##############
#################################################
# Function that launches rep times the simulation and returns a dataframe with results
compute_estimators_and_store <- function(rep, misoutcome = "correct", misRCT = "correct", df1=df_rct, df2=df_obs){
  results <- data.frame("Naive OnlyRCT" = numeric(rep),
                        "IPSW" = numeric(rep),
                        "IPSW.norm" = numeric(rep),
                        "IPSW.strat.n.2" = numeric(rep),
                        "IPSW.strat.n.3" = numeric(rep),
                        "IPSW.strat.n.4" = numeric(rep),
                        "IPSW.strat.n.5" = numeric(rep),
                        "IPSW.strat.n.6" = numeric(rep),
                        "IPSW.strat.n.7" = numeric(rep),
                        "IPSW.strat.n.8" = numeric(rep),
                        "IPSW.strat.n.9" = numeric(rep),
                        "IPSW.strat.n.10" = numeric(rep),
                       # "IPSW.randomforest" = numeric(rep),
                        "G-formula" = numeric(rep))
  
  for (i in 1:rep){
    # Set seed for bootstrap sampling
    # Bootstrap sampling
    idx_rct <- sample(1:nrow(df1), replace = TRUE)
    idx_obs <- sample(1:nrow(df2), replace = TRUE)
    
    # Merge the bootstrap samples
    DF <- rbind(df1[idx_rct, ], df2[idx_obs, ])
    DF <- as.data.frame(DF)
    DF$A <- as.numeric(DF$A)
    
    # Naive estimator
    results[i, "Naive.OnlyRCT"] <- mean(DF[DF$A == 1 & DF$V == 1, 'Y']) - mean(DF[DF$A == 0  & DF$V == 1, "Y"])
    
    # IPSW
    results[i, "IPSW"] <- compute_ipsw(DF, normalized = FALSE)
    results[i, "IPSW.norm"] <- compute_ipsw(DF, normalized = TRUE)
    #results[i, "IPSW.randomforest"] <- compute_ipsw(DF, normalized = FALSE,  estimation = "grf")
    
    # Stratification
    for (n_strat in 2:10) {
      results[i, paste0("IPSW.strat.n.", n_strat)] <- compute_stratification(DF, n_strat)
    }
    
    # G-formula
    results[i, "G.formula"] <- compute_gformula(DF)
  }
  
  return(results)
}
