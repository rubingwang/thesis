library(reshape2)
library(MASS)
library(dplyr)
library(xgboost)
library(randomForest)

compute_mean_diff_RCT <- function(DF){
  RCT_ATE <- mean(DF[DF$A == 1 & DF$V == 1, "Y"]) - mean(DF[DF$A == 0  & DF$V == 1, "Y"])  
  return(RCT_ATE)
}

# IPSW: covariates parameters is for simulation, to restrict only to X1 or not the estimation. 
compute_ipsw <- function(DF, normalized = FALSE, estimation = "logit", covariates = "All"){
  
  N <- nrow(DF)
  n <- nrow(DF[DF$V ==1, ])
  m <- nrow(DF[DF$V ==0, ])
  
  ## simply delete this paragraph when removal of the covariates parameter
  if (covariates == "All"){
    temp <- DF
  } else if (covariates == "X1"){
    temp <- DF[, c("X1", "V", "A", "Y")]
  } else if (covariates == "-X1"){
    temp <- DF[, !(names(DF) %in%  c("X1"))]
  } else {
    print("Covariates parameters must be All, X1, -X1.")
    break
  }
  
  # Estimation of P(V = 1 | X)
  # p <-- P(V = 1 | X) 
  
  if (estimation == "logit") {
    # with logistic regression
    p.fit  <- glm(V ~., family = binomial("logit"), data = temp[, !names(temp) %in% c("A", "Y")])
    p <- predict(p.fit, type = "response", newdata = temp)
  } 
  
  else if (estimation == "forest") {
    # with Random Forest
    p.fit <- randomForest(factor(V) ~ ., data = temp[, !names(temp) %in% c("A", "Y")])
    p <- predict(p.fit, newdata = temp, type = "response")
  } 
  
  else if (estimation == "boosting") {
    # with XGBOOST
    p.fit <- xgboost(data = as.matrix(temp[, !names(temp) %in% c("A", "Y")]), 
                     label = temp$V, max.depth = 3, eta = 0.3, nrounds = 10, objective = "binary:logistic")
    p <- predict(p.fit, newdata = as.matrix(temp[, !names(temp) %in% c("A", "Y")]))
  } 
  else {
    print("Estimation parameter should be forest, boosting or logit.")
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


# STRATIFICATION
compute_stratification <- function(DF, nb_strat = 10, bin = "quantile", estimation = "logit") {
  
  temp <- DF
  
  # Estimation of P(V = 1 | X)
  if (estimation == "logit") {
    pi_s_reg <- glm(V ~ ., family = "binomial", data = temp[, !names(temp) %in% c("Y", "A")])
    pi_strat <- predict(pi_s_reg, newdata = temp[, !names(temp) %in% c("Y", "A")], type="response")
  } 
  else if (estimation == "forest") {
    library(randomForest)
    pi_s_reg <- randomForest(factor(V) ~ ., data = temp[, !names(temp) %in% c("Y", "A")])
    pi_strat <- predict(pi_s_reg, newdata = temp[, !names(temp) %in% c("Y", "A")], type = "response")
  } 
  else if (estimation == "boosting") {
    library(xgboost)
    pi_s_reg <- xgboost(data = as.matrix(temp[, !names(temp) %in% c("Y", "A")]), 
                        label = temp$V, max.depth = 3, eta = 0.3, nrounds = 10, objective = "binary:logistic")
    pi_strat <- predict(pi_s_reg, newdata = as.matrix(temp[, !names(temp) %in% c("Y", "A")]))
  } 
  else {
    print("Estimation parameter should be logit, forest, or boosting.")
    return(NULL)
  }
  
  temp$pi_strat <- (1 - pi_strat) / pi_strat
  
  # Stratification
  if (bin == "quantile") {
    temp$strata <- as.numeric(cut_number(temp$pi_strat, nb_strat))
  } 
  else if (bin == "regular") {
    temp$strata <- as.numeric(cut(temp$pi_strat, nb_strat))
  }
  else {
    print("Wrong bin mentioned")
    return(NULL)
  }
  
  rct <- temp[temp$V == 1, ]
  m <- nrow(temp[temp$V == 0, ])
  tau_strat <- 0
  
  for (s in unique(rct$strata)) {
    strata_ate <- mean(rct[rct$strata == s & rct$A == 1, "Y"]) - mean(rct[rct$strata == s & rct$A == 0, "Y"])
    weigth <- nrow(temp[temp$V == 0 & temp$strata == s, ]) / m
    tau_strat <- tau_strat + weigth * strata_ate
  }
  
  return(tau_strat)
}

# Simulation setting: transformation to mispecify
transform_X_star <-function(DF){
  DF_transformed <- DF
  DF_transformed$X1 <- exp(DF$X1 / 3)
  DF_transformed$X2 <- DF$X2 / (1+exp(DF$X1)) + 10
  DF_transformed$X3 <- DF$X1 * DF$X3 / 25 + 0.6
  DF_transformed$X4 <- DF$X1 + DF$X4 + 20
  matrix_transformed <- as.matrix(DF_transformed[, c("X1", "X2", "X3", "X4")])
  matrix_transformed <- scale(matrix_transformed, center = c(1,1,1,1))
  return(matrix_transformed)
}

# simulation with continuous outcome
simulate_continuous <- function(n = 1000, m = 49000, p = 4, mu = rep(1, p), Sigma = diag(p), bs = c(-0.5, -0.3, -0.5, -0.4), bs0 = -2.5, beta = c(27.4, 13.7, 13.7, 13.7), b0 = - 100, sigma = 1, misRCT = "correct", misoutcome = "correct", Nested = FALSE) {
  
  # Target population generation 
  covariates <- mvrnorm(n = 50*n, mu, Sigma, tol = 1e-06, empirical = FALSE) # 50*n is roughly the initial population size necessary to have the n
  DF <- as.data.frame(covariates)
  names(DF) <- paste("X", 1:p, sep = "")
  covariates_names <- names(DF)
  
  # RCT probability to sample according to model
  if (misRCT == "correct"){
    etas <- as.vector(covariates %*% bs + bs0)
    
  } else if (misRCT == "exponential") {
    
    # RCT misspecification with exp on all covariates
    etas <- as.vector(exp(covariates) %*% bs + bs0 + 3) # 3 was found manually to keep same proportion m and n
    
  } else if (misRCT == "Partial_X2only"){
    
    # partial misspecification with only X2 affected
    DF_mis <- DF
    DF_mis$X2 <- exp(DF_mis$X2)
    mis_covariatesX2 <- as.matrix(DF_mis[, c("X1", "X2", "X3", "X4")])
    etas <- as.vector(mis_covariatesX2 %*% bs + bs0 + 0.1)
    
  } else if (misRCT == "Partial_X1only"){
    
    # partial misspecification with only X1 affected
    DF_mis <- DF
    DF_mis$X1 <- exp(DF_mis$X1)
    mis_covariatesX1 <- as.matrix(DF_mis[, c("X1", "X2", "X3", "X4")])
    etas <- as.vector(mis_covariatesX1 %*% bs + bs0 + 0.1)
    
  }
  else if(misRCT == "dong"){
    
    miscovariatesDong <- transform_X_star(DF)
    etas <- as.vector (miscovariatesDong %*% bs + bs0)
    
  } else if (misRCT == "strongbias"){
    bs = c(-1.5, -0.3, -0.5, -0.4)
    etas <- as.vector (covariates %*% bs + bs0)
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
    DF$Y = b0 + beta[1]*(DF$A == 1)*DF$X1 + beta[2]*DF$X2 + beta[3]*DF$X3 +
      beta[4]*DF$X4 + error
  } 
  else if (misoutcome == "+a"){
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0 + DF$X1 + beta[2]*DF$X2 + beta[3]*DF$X3 +
      beta[4]*DF$X4 + error + beta[1]*(DF$A == 1)
  }
  else if (misoutcome == "wrong")
  {
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0 + beta[1]*(DF$A == 1)*DF$X1*DF$X2 + beta[2]*DF$X2 + beta[3]*DF$X3 +
      beta[4]*DF$X4 + error
  } 
  else if (misoutcome == "dong") {
    miscovariatesDong <- transform_X_star(DF)
    DF_Xstar <- as.data.frame(miscovariatesDong)
    names(DF_Xstar) <- paste("X", 1:p, sep = "")
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0 + beta[1]*(DF$A == 1)*DF_Xstar$X1 + beta[2]*DF_Xstar$X2 + beta[3]*DF_Xstar$X3 + beta[4]*DF_Xstar$X4 + error
  }
  else {
    print("Parameters misoutcome is badly specified")
    break
  }
  return(DF)
}


# Function that launches rep times the simulation and returns a dataframe with results
compute_estimators_and_store <- function(rep, misoutcome = "correct", misRCT = "correct", N = 50000, m = 49000, n = 1000){
  
  rct_ate <- c()
  ipsw <- c()
  ipsw_norm <- c()
  strat_10 <- c()
  gformula <- c()
  cw <- c()
  aipsw <- c()
  
  for (i in 1:rep){
    #set different seeds for each loop [by rubing 8 Mar]
    DF <- simulate_continuous(misoutcome = misoutcome, misRCT = misRCT, m = m, n = n)
    
    # naive estimator
    rct_ate <- c(rct_ate, mean(DF[DF$A == 1 & DF$V == 1, "Y"]) - mean(DF[DF$A == 0  & DF$V == 1, "Y"]))
    
    #ispw
    ipsw  <- c(ipsw, compute_ipsw(DF, normalized = F))
    ipsw_norm <- c(ipsw_norm, compute_ipsw(DF, normalized = T))
    
    #strat
    strat_10 <- c(strat_10, compute_stratification(DF, 10))
  }
  
  results <- data.frame("RCT" = rct_ate,
                        "IPSW" = ipsw,
                        "IPSW norm" = ipsw_norm,
                        "Stratification n=10" = strat_10,
                        ) 
}

