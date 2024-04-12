# This R code is to determine the specific model parameter of generating function
# Since we need to keep the same m/n proportion to make the estimators comparbale 

# Clear any existing variables
rm(list = ls())
# Set seed for reproducibility
set.seed(411)
# Load implemented functions
source('./function.R')
# Libraries
library(ggplot2) # plots
library(dplyr) # data frame re arrangment
library(table1) # table for baseline
library(xtable) # table for baseline
library(wesanderson) # colors

# Initial setting for stimulation studies
repetitions = 100  # numbers of estimation 
theta =            # true ATE (of the target population)
b_selection = c(-0.53, -0.47, -0.60, -0.55) # beta for selection model
b0_selection = -2.285                     # beta0 for selection model
b_outcome = c(10.5, 20, 20, 20)   # beta for outcome model
b0_outcome  = - 35                     # beta for outcome model

# simulation with continuous outcome
n = 1000
m = 49000
p = 4
mu = rep(1, p)
Sigma = diag(p)
sigma = 1

  
# Target population generation 
covariates <- mvrnorm(n = 50*n, mu, Sigma, tol = 1e-06, empirical = FALSE) # 50*n is roughly the initial population size necessary to have the n
DF <- as.data.frame(covariates)
names(DF) <- paste("X", 1:p, sep = "")
covariates_names <- names(DF)
  
# RCT probability to sample according to model
#  if (misRCT == "correct"){
    etas <- as.vector(covariates %*% b_selection + b0_selection)
#  } else if (misRCT == "exponential") {
    # RCT misspecification with exp on all covariates
    etas <- as.vector(exp(covariates) %*% b_selection + b0_selection + 3.58) # 3 was found manually to keep same proportion m and n
#  } else if (misRCT == "strongbias"){
    b_selection = c(-1.53, -0.47, -0.60, -0.55)
    etas <- as.vector (covariates %*% b_selection + b0_selection)
#  }  
  
  ps = 1 / (1 + exp(-etas))
  DF$ps <- ps
  
  # from probability to RCT indicator
  RCT_indicator <- rbinom(length(ps), 1, as.vector(ps))
  DF$V <- RCT_indicator 
  # random treatment assignement within the RCT
  DF$A <- ifelse(DF$V == 1, rbinom(nrow(DF), 1, 0.5), NA)
  # keep only interesting variables
  DF <- DF[, c(covariates_names, "A", "V")]
    
  # drop other data
  DF_rct <- DF[DF$V == 1,] 
    
  # generate new observational data
  covariates_rwe <- mvrnorm(n = m, mu, Sigma, tol = 1e-06, empirical = FALSE) 
  DF_rwe <- as.data.frame(covariates_rwe)
  names(DF_rwe) <- paste("X", 1:p, sep = "")
  DF_rwe$V <- rep(0, m)
  DF_rwe$A <- rep(NA, m)
  
  # stack RCT and RWE
  DF <- rbind(DF_rct, DF_rwe)
  
  # reset row number
  rownames(DF) <- 1:nrow(DF)
  
  # compute Y  
#  if (misoutcome == "correct"){
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0_outcome + b_outcome[1]*(DF$A == 1)*DF$X1 + b_outcome[2]*DF$X2 + b_outcome[3]*DF$X3 +
      b_outcome[4]*DF$X4 + error
#  } 
#  else if (misoutcome == "wrong")
#  {
    error = rnorm(n = nrow(DF), mean = 0, sd = sigma)
    DF$Y = b0_outcome + b_outcome[1]*(DF$A == 1)*DF$X1*DF$X2 +b_outcome[2]*DF$X2 + b_outcome[3]*DF$X3 +
      b_outcome[4]*DF$X4 + error
#  } 
  
  
########  X1 Shift descriptive ########
# 1.well-specification simulation: weak shift (including all 4 covariate)
one_simulation1 <- simulate_continuous(n = 1000, m = 10000)
one_simulation1$sample <- ifelse(one_simulation1$V == 1, "RCT", "Observational")
baseline_weak <- table1(~ X1 + X2 + X3 + X4 | sample, data = one_simulation1, overall="Total")
# 2.mis-specification simulation: X1 strong shift (including all 4 covariate)
one_shifted_simulation2 <- simulate_continuous(n = 1000, m = 10000, misRCT = "strongbias")
one_shifted_simulation2 $sample <- ifelse(one_shifted_simulation2 $V == 1, "RCT", "Observational")
baseline_x1strong <- table1(~ X1 + X2 + X3 + X4 | sample, data = one_shifted_simulation2, overall="Total")

# Distribution shift plot on X1
x1weak <- one_simulation1 
x1weak$Shift <- rep("X1 shift: weak", nrow(x1weak))
one_shifted_simulation2$Shift <- rep("X1 shift: strong", nrow(one_shifted_simulation2))
shift_comparison1 <- rbind(x1weak, one_shifted_simulation2)

ggplot(shift_comparison1, aes(x = X1, group = sample, fill = sample)) +
  #geom_histogram(binwidth = 0.2, alpha=0.4, position="dodge") + 
  geom_density(alpha=0.4, position="dodge") +
  scale_fill_manual(values=c("darkorchid4", "darkorange1")) +
  theme_classic() +
  theme(legend.title = element_blank(), 
        legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.text = element_text(size=18, face="bold")) +
  ylab("") + 
  theme(axis.text = element_text(vjust = 0.5, hjust=1, size=14, face="bold"), axis.title.x = element_text(size=18, face="bold")) +
  facet_grid(~Shift)  +
  theme(strip.background=element_rect(fill=NA, color=NA),
        strip.text=element_text(size=15, face = "bold")) +
  ggsave("X1strongshift.png", width = 8, height = 8.5, dpi=500)


