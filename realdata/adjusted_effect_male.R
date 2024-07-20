# Clear any existing variables
rm(list = ls())
# Set seed for reproducibility
set.seed(411)
# Load implemented functions
source('~/tmp/thesis/realdata/function_realdata.R')
# Libraries
library(ggplot2) # plots
library(dplyr) # data frame re arrangment
library(table1) # table for baseline
library(xtable) # table for baseline
library(wesanderson) # colors
library(readxl)
library(boot)
rep=10
# Creat file folder
#dir.create("~/tmp/thesis/realdata/realdata_plot")
#dir.create("~/tmp/thesis/realdata/realdata_metrics")

#################################################
############## Analysis in 500 times ############
#################################################

# Read the data
df_obs <- read.csv("/Users/rubing/tmp/thesis/stroke_data.csv")
df_rct <- read_excel("~/tmp/thesis/realdata/RCT_sample.xlsx")
df_rct <- df_rct %>%
  mutate(
    group = ifelse(group == 'treated', '1', '0'),
    sex = ifelse(sex == 'Male', '1', '0'),
    type_stroke = ifelse(type_stroke == 'Hemorrhagic', '1', '0'),
  )
df_obs <- df_obs %>%
  mutate(
    sex = ifelse(sex == 'male', '1', '0'),
    type_stroke = ifelse(type_stroke == 'Hemorrhagic', '1', '0'),
    history_alcoholism = ifelse(history_alcoholism == 'yes', '1', '0')
  )
df_rct <- subset(df_rct, sex == 1)
df_obs <- subset(df_obs, sex == 1) 
df_rct <- df_rct[, c("age", "weight", "height", "BMI", "type_stroke", "group", "history_alcoholism", "d_ankle")]
df_obs <- df_obs[, c("age", "weight", "height", "BMI", "type_stroke",  "history_alcoholism")]


names(df_rct)[names(df_rct) == "group"] <- "A"
names(df_rct)[names(df_rct) == "d_ankle"] <- "Y"
# Add a new column to each data frame to indicate the data source
df_obs$V <- 0
df_rct$V <- 1
df_obs$A <- 99
df_obs$Y <- 99

# set the format
df_rct <- as.data.frame(df_rct)
df_obs <- as.data.frame(df_obs)

######
#i=1
#df1=df_rct
#df2=df_obs
#df <- DF
#################################################
################# new function ##################
#################################################

# Function that launches rep times the simulation and returns a dataframe with results
compute_estimators_sex<- function(rep, df1=df_rct, df2=df_obs){
  results <- data.frame("Naive OnlyRCT" = numeric(rep),
                        "G-formula" = numeric(rep),
                        "IPSW" = numeric(rep),
                         "IPSW.norm" = numeric(rep))
  for (i in 1:rep){
    # Bootstrap sampling
    #idx_rct <- sample(1:nrow(df1), replace = TRUE)
    #idx_obs <- sample(1:nrow(df2), replace = TRUE)
   
    stratum1 <- "type_stroke"
    stratum2 <- "history_alcoholism"
    df1$stratum <- as.factor(interaction(df1[[stratum1]], df1[[stratum2]]))
    df2$stratum <- as.factor(interaction(df2[[stratum1]], df2[[stratum2]]))
    idx_rct <- unlist(tapply(1:nrow(df1), df1$stratum, function(x) sample(x, replace = TRUE)))
    idx_obs <- unlist(tapply(1:nrow(df2), df2$stratum, function(x) sample(x, replace = TRUE)))
    
    # Merge the bootstrap samples
    DF <- rbind(df1[idx_rct, ], df2[idx_obs, ])
    DF <- as.data.frame(DF)
    DF$A <- as.numeric(DF$A)
    
    # IPSW
    results[i, "IPSW"] <- compute_ipsw(DF, normalized = FALSE)
    results[i, "IPSW.norm"] <- compute_ipsw(DF, normalized = TRUE)
    
    # Naive estimator
    results[i, "Naive.OnlyRCT"] <- mean(DF[DF$A == 1 & DF$V == 1, 'Y',]) - mean(DF[DF$A == 0  & DF$V == 1, "Y"])
    # G-formula
    results[i, "G.formula"] <- compute_gformula(DF)
  }
  
  return(results)
}
# G-formula
compute_gformula <- function(DF){
  temp <- DF
  mu_1 <- lm(Y ~., data = temp[temp$V == 1 & temp$A == 1, !names(temp) %in% c("V", "A",'stratum')])
  mu_0 <- lm(Y ~., data = temp[temp$V == 1 & temp$A == 0, !names(temp) %in% c("V", "A",'stratum')])
  
  mu_1_predict <- predict.lm(mu_1, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "A",'stratum')])
  mu_0_predict <- predict.lm(mu_0, newdata = temp[temp$V == 0, !names(temp) %in% c("V", "A",'stratum')])
  
  tau_hat_gformula <- mean(mu_1_predict) - mean(mu_0_predict)
  
  return(tau_hat_gformula)  
}

#################################################
################### 500 times ###################
#################################################

results_list <- list()
seeds <- 1:5
for (seed in seeds) {
  set.seed(seed)  
  results <- compute_estimators_sex(rep=20, df1 = df_rct, df2 = df_obs)
  results_list[[seed]] <- results
}

total_results <- do.call(rbind, results_list)


# plot boxplot
results_melted <- melt(total_results, variable.name = "Estimator", value.name = "Value")
ggplot(results_melted, aes(x = Estimator, y = Value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 4.87371, color = "Naive.OnlyRCT"), 
             size = 0.6, linetype="dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 10)) +
  labs(title = "Performance comparison of the estimators (Sex=male)", x = "Estimator", y = "Estimated PATE")


#################################################
############## Summary Statistics ###############
#################################################
rep=100
estimator_means <- apply(total_results, 2, mean, na.rm = TRUE)
estimator_sds <- apply(total_results, 2, sd, na.rm = TRUE)
rep <- nrow(total_results)

estimator_lower_ci <- estimator_means - 1.96 * (estimator_sds / sqrt(rep))
estimator_upper_ci <- estimator_means + 1.96 * (estimator_sds / sqrt(rep))

estimator_stats <- data.frame(
  Estimator = names(estimator_means),
  Mean = estimator_means,
  SD = estimator_sds,
  CI_lower = round(estimator_lower_ci, 3),
  CI_upper = round(estimator_upper_ci, 3)
)

print(estimator_stats)

write.csv(estimator_stats, file = "~/tmp/thesis/realdata/effect_stats_male.csv", row.names = FALSE)
save(total_results, file = "~/tmp/thesis/realdata/total_results_male.RData")

