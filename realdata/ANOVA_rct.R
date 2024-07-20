# Clear any existing variables
rm(list = ls())
# Set seed for reproducibility
set.seed(411)
# Load implemented functions
source('~/tmp/thesis/realdata/function_realdata.R')
# Libraries
library(dplyr) # data frame re arrangment


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


model <- lm(Y ~ age + A*sex + BMI + type_stroke , data = df_rct)
summary(model)


initial_model <- lm(Y ~ age + A * sex + BMI + type_stroke, data = df_rct)

# 使用逐步变量选择优化模型
optimized_model <- step(initial_model, direction = "both")

# 打印优化后的模型摘要
summary(optimized_model)




