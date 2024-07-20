# Clear any existing variables
rm(list = ls())
# Set seed for reproducibility
set.seed(411)
# Load implemented functions
source('~/tmp/thesis/realdata/function_realdata.R')
# Libraries
library(reshape2)
library(ggplot2) # plots
library(dplyr) # data frame re arrangment
library(table1) # table for baseline
library(xtable) # table for baseline
library(wesanderson) # colors
library(readxl)
rep=500
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

df_rct <- df_rct[, c("age", "sex", "weight", "height", "BMI", "type_stroke", "group", "history_alcoholism", "d_ankle")]
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

# total PATE
total_results <- compute_estimators_and_store(rep = 300, df1= df_rct, df2 = df_obs)

# plot boxplot
results_melted <- melt(total_results, variable.name = "Estimator", value.name = "Value")
ggplot(results_melted, aes(x = Estimator, y = Value)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = 7.48468, color = "Naive.OnlyRCT"), 
             size = 0.6, linetype="dashed") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(5, 15)) +
  labs(title = "Performance comparison of the estimators: boxplot based on 500 iterations", x = "Estimator", y = "Estimated PATE")

#################################################
############## Summary Statistics ###############
#################################################
rep=300
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

write.csv(estimator_stats, file = "~/tmp/thesis/realdata/estimator_stats.csv", row.names = FALSE)

