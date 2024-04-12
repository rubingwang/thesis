library(reshape2)
library(MASS)
library(dplyr)
library(psych)
#################################################
############## Evaluation function ##############
#################################################
compute_metrics <- function(df, theta, output = FALSE) {
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
    p.fit  <- glm(V ~., family = binomial("logit"), data = temp[, !names(temp) %in% c("A", "Y")])
    p <- predict(p.fit, type = "response", newdata = temp)
    
  } 
  else if (estimation == "forest") {
    break
  } 
  else {
    print("Estimation parameter should be forest or logit.")
    break
  }
  
  # Store odds
  temp$odds <- ((1 - p)/p)
  
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
    tau_strat <- tau_strat + weigth*strata_ate
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
############ Data generating function ###########
#################################################

# simulation with continuous outcome
simulate_continuous <- function(n = 1000, m = 49000, 
                                p = 4, mu = rep(1, p), Sigma = diag(p),
                                b0_selection = -2.285, 
                                b_selection = c(-0.53, -0.47, -0.60, -0.55),
                                b0_outcome = - 50, 
                                b_outcome = c(32, 20, 20, 20), 
                                sigma = 1, misRCT = "correct", 
                                misoutcome = "correct", Nested = FALSE) {
  
  # Target population generation 
  covariates <- mvrnorm(n = 50*n, mu, Sigma, tol = 1e-06, empirical = FALSE) # 50*n is roughly the initial population size necessary to have the n
  DF <- as.data.frame(covariates)
  names(DF) <- paste("X", 1:p, sep = "")
  covariates_names <- names(DF)
  
  # RCT probability to sample according to model
  if (misRCT == "correct"){
    etas <- as.vector(covariates %*% b_selection + b0_selection)
  } else if (misRCT == "exponential") {
    # RCT misspecification with exp on all covariates
    etas <- as.vector(exp(covariates) %*% b_selection + b0_selection + 3.58) # 3.58 was found manually to keep same proportion m and n
  } else if (misRCT == "strongbias"){
    b_selection = c(-0.53-1.5, -0.47, -0.60, -0.55)
    etas <- as.vector (covariates %*% b_selection + b0_selection)
  }  else {
    print("Error in RCT specification arguments.")
    break
  }
  
  ps = 1 / (1 + exp(-etas))
  DF$ps <- ps
  
  # from probability to RCT indicator
  RCT_indicator <- rbinom(length(ps), 1, as.vector(ps))
  DF$V <- RCT_indicator 
  
  # random treatment assignement within the RCT
  DF$A <- ifelse(DF$V == 1, rbinom(nrow(DF), 1, 0.5), NA)
  
  # keep only interesting variables
  DF <- DF[, c(covariates_names, "A", "V")]
  
  if (!Nested) {
    
    # drop other data
    DF_rct <- DF[DF$V == 1,] 
    
    # generate new observational data
    covariates_rwe <- mvrnorm(n = m, mu, Sigma, tol = 1e-06, empirical = FALSE) 
    DF_rwe <- as.data.frame(covariates_rwe)
    names(DF_rwe) <- paste("X", 1:p, sep = "")
    DF_rwe$V <- rep(0, m)
    DF_rwe$A <- rep(NA, m)
    
  } else {
    
    #here we need to drop values such that the final data set contains m observational values and n RCT values.
    DF_rct <- DF[DF$V == 1,]
    DF_rwe <- DF[DF$V == 0,]
  }
  
  # stack RCT and RWE
  DF <- rbind(DF_rct, DF_rwe)
  
  # reset row number
  rownames(DF) <- 1:nrow(DF)
  
  # compute Y  
  if (misoutcome == "correct"){
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0_outcome + b_outcome[1]*(DF$A == 1)*DF$X1 + b_outcome[2]*DF$X2 + b_outcome[3]*DF$X3 +
      b_outcome[4]*DF$X4 + error
  } 
  else if (misoutcome == "wrong")
  {
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0_outcome + b_outcome[1]*(DF$A == 1)*DF$X1*DF$X2 +b_outcome[2]*(DF$X2/exp(DF$X3)) + b_outcome[3]*DF$X3 +
      b_outcome[4]*DF$X4 + error
  } 
  else {
    print("Parameters misoutcome is badly specified")
    break
  }
  return(DF)
}


#################################################
############## Out-of-box function ##############
#################################################

# Function that launches rep times the simulation and returns a dataframe with results
compute_estimators_and_store <- function(rep, misoutcome = "correct", misRCT = "correct", N = 50000, m = 49000, n = 1000){
  rct_ate <- c()
  ipsw <- c()
  ipsw_norm <- c()
  strat_2 <- c()
  strat_3 <- c()
  strat_4 <- c()
  strat_5 <- c()
  strat_6 <- c()
  strat_7 <- c()
  strat_8 <- c()
  strat_9 <- c()
  strat_10 <- c()
  gformula <- c()
 
  for (i in 1:rep){
    #set different seeds for each loop [by rubing 8 Mar]
    DF <- simulate_continuous(misoutcome = misoutcome, misRCT = misRCT, m = m, n = n)
    
    # naive estimator
    rct_ate <- c(rct_ate, mean(DF[DF$A == 1 & DF$V == 1, "Y"]) - mean(DF[DF$A == 0  & DF$V == 1, "Y"]))
    
    #ispw
    ipsw  <- c(ipsw, compute_ipsw(DF, normalized = F))
    ipsw_norm <- c(ipsw_norm, compute_ipsw(DF, normalized = T))
    
    #strat
    strat_2 <- c(strat_2, compute_stratification(DF, 2))
    strat_3 <- c(strat_3, compute_stratification(DF, 3))
    strat_4 <- c(strat_4, compute_stratification(DF, 4))
    strat_5 <- c(strat_5, compute_stratification(DF, 5))
    strat_6 <- c(strat_6, compute_stratification(DF, 6))
    strat_7 <- c(strat_7, compute_stratification(DF, 7))
    strat_8 <- c(strat_8, compute_stratification(DF, 8))
    strat_9 <- c(strat_9, compute_stratification(DF, 9))
    strat_10 <- c(strat_10, compute_stratification(DF, 10))
    
    #gformula
    gformula <- c(gformula, compute_gformula(DF))

  }
  
  results <- data.frame("Naive OnlyRCT" = rct_ate,
                        "IPSW" = ipsw,
                        "IPSW norm" = ipsw_norm,
                        "IPSW strat n=2" = strat_2,
                        "IPSW strat n=3" = strat_3,
                        "IPSW strat n=4" = strat_4,
                        "IPSW strat n=5" = strat_5,
                        "IPSW strat n=6" = strat_6,
                        "IPSW strat n=7" = strat_7,
                        "IPSW strat n=8" = strat_8,
                        "IPSW strat n=9" = strat_9,
                        "IPSW strat n=10" = strat_10,
                        "G-formula" = gformula)
}

