library(ggplot2)
library(readxl)
library(dplyr)
library(gridExtra)

# Read the data
df_obs <- read.csv("/Users/rubing/tmp/thesis/stroke_data.csv")
df_rct <- read_excel("~/tmp/thesis/realdata/RCT_sample.xlsx")
df_rct <- df_rct %>%
  mutate(history_alcoholism = ifelse(history_alcoholism == 1, 'yes', 'no'))
df_rct <- df_rct %>% select(age, sex, weight, height, BMI, type_stroke, history_alcoholism)

# Add a new column to each data frame to indicate the data source
df_obs$source <- 'Observational'
df_rct$source <- 'RCT'

# Combine the two data frames
combined_df <- rbind(df_obs, df_rct)

# plot densities for each pair of continuous variables in the rct and observational sample
plot_category_distribution <- function(data, category, title, x_label) {
  category_counts <- data %>%
    group_by(source, !!sym(category)) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(percentage = count / sum(count) * 100)
  
  ggplot(category_counts, aes_string(x = category, y = "count", fill = "source")) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5) +
    labs(title = title, x = x_label, y = "Count", fill = "Data Source") +
    theme_minimal() +
    theme(legend.title = element_blank(), 
          legend.position = "bottom", 
          legend.box = "horizontal", 
          legend.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(vjust = 0.5, hjust = 1, size = 12, face = "bold"), 
          axis.title.x = element_text(size = 14, face = "bold"),
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(size = 12, face = "bold"))
}

plot_sex <- plot_category_distribution(combined_df, "sex", "Distribution of Sex between the two samples", "Sex")
plot_typestroke <-plot_category_distribution(combined_df, "type_stroke", "Distribution of Stroke Type between the two samples", "Stroke Type")
plot_alco <-plot_category_distribution(combined_df, "history_alcoholism", "Distribution of History of Alcoholism between the two samples", "History of Alcoholism")
grid.arrange(plot_sex,plot_typestroke,plot_alco,ncol = 3)


# plot densities for each pair of continuous variables in the rct and observational sample
plot_continuous_distribution <- function(data, variable, title, x_label) {
  means <- data %>%
    group_by(source) %>%
    summarize(mean_value = mean(!!sym(variable), na.rm = TRUE))
  
  ggplot(data, aes_string(x = variable, fill = "source", color = "source")) +
    geom_density(alpha = 0.5) +
    geom_vline(data = means, aes(xintercept = mean_value, color = source), 
               linetype = "dashed", size = 0.7) +
    scale_color_manual(values = c("blue", "red")) +
    scale_fill_manual(values = c("blue", "red")) +
    labs(title = title, x = x_label, y = "Density", fill = "Data Source") +
    theme_minimal() +
    theme(legend.title = element_blank(), 
          legend.position = "bottom", 
          legend.box = "horizontal", 
          legend.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(vjust = 0.5, hjust = 1, size = 12, face = "bold"), 
          axis.title.x = element_text(size = 14, face = "bold"),
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(size = 12, face = "bold"))
}

plot_age <- plot_continuous_distribution(combined_df, "age", "Distribution of Age between the two samples (dashlines: two means)", "Age")
plot_weight <- plot_continuous_distribution(combined_df, "weight", "Distribution of Weight between the two samples (dashlines: two means)", "Weight (kg)")
plot_height <- plot_continuous_distribution(combined_df, "height", "Distribution of Height between the two samples (dashlines: two means)", "Height (cm)")
plot_bmi <- plot_continuous_distribution(combined_df, "BMI", "Distribution of BMI between the two samples (dashlines: two means)", "BMI")

# put the four plots in a same plot
shift_plot_continous<-grid.arrange(plot_age, plot_weight, plot_height, plot_bmi, ncol = 2)

#########using statiscs #########
# libaray
library(ggplot2)
library(dplyr)
library(tidyr)

# Kolmogorov-Smirnov
ks_tests <- list(
  age = ks.test(combined_df$age[combined_df$source == "Observational"],
                combined_df$age[combined_df$source == "RCT"]),
  weight = ks.test(combined_df$weight[combined_df$source == "Observational"],
                   combined_df$weight[combined_df$source == "RCT"]),
  height = ks.test(combined_df$height[combined_df$source == "Observational"],
                   combined_df$height[combined_df$source == "RCT"]),
  BMI = ks.test(combined_df$BMI[combined_df$source == "Observational"],
                combined_df$BMI[combined_df$source == "RCT"])
)

ks_tests

# Chi-square
sex_table <- table(combined_df$sex, combined_df$source)
type_stroke_table <- table(combined_df$type_stroke, combined_df$source)
history_alcoholism_table <- table(combined_df$history_alcoholism, combined_df$source)

chisq_tests <- list(
  sex = chisq.test(sex_table),
  type_stroke = chisq.test(type_stroke_table),
  history_alcoholism = chisq.test(history_alcoholism_table)
)

chisq_tests



# 加载必要的库
library(ggplot2)
library(dplyr)
library(gridExtra)

# 函数：绘制连续变量的密度图
plot_density <- function(data, variable, title) {
  ggplot(data, aes_string(x = variable, fill = "source")) +
    geom_density(alpha = 0.5) +
    labs(title = title, x = variable, y = "Density") +
    theme_minimal() +
    theme(legend.title = element_blank(), 
          legend.position = "bottom", 
          legend.box = "horizontal", 
          legend.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(vjust = 0.5, hjust = 1, size = 12, face = "bold"), 
          axis.title.x = element_text(size = 14, face = "bold"),
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(size = 12, face = "bold"))
}

# 函数：绘制分类变量的百分比条形图
plot_category_distribution <- function(data, category, title, x_label) {
  category_counts <- data %>%
    group_by(source, !!sym(category)) %>%
    summarize(count = n()) %>%
    ungroup() %>%
    group_by(source) %>%
    mutate(percentage = count / sum(count) * 100)
  
  ggplot(category_counts, aes_string(x = category, y = "percentage", fill = "source")) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_text(aes(label = paste0(round(percentage, 1), "%")), 
              position = position_dodge(width = 0.9), 
              vjust = -0.5) +
    labs(title = title, x = x_label, y = "Percentage", fill = "Data Source") +
    theme_minimal() +
    theme(legend.title = element_blank(), 
          legend.position = "bottom", 
          legend.box = "horizontal", 
          legend.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(vjust = 0.5, hjust = 1, size = 12, face = "bold"), 
          axis.title.x = element_text(size = 14, face = "bold"),
          strip.background = element_rect(fill = NA, color = NA),
          strip.text = element_text(size = 12, face = "bold"))
}

# 绘制连续变量的密度图
plot_age <- plot_density(combined_df, "age", "Density Plot of Age by Source")
plot_weight <- plot_density(combined_df, "weight", "Density Plot of Weight by Source")
plot_height <- plot_density(combined_df, "height", "Density Plot of Height by Source")
plot_BMI <- plot_density(combined_df, "BMI", "Density Plot of BMI by Source")

# 绘制分类变量的分布图
plot_sex <- plot_category_distribution(combined_df, "sex", "Distribution of Sex between the two samples", "Sex")
plot_typestroke <- plot_category_distribution(combined_df, "type_stroke", "Distribution of Stroke Type between the two samples", "Stroke Type")
plot_alco <- plot_category_distribution(combined_df, "history_alcoholism", "Distribution of History of Alcoholism between the two samples", "History of Alcoholism")

# 安排图表
grid.arrange(plot_age, plot_weight, plot_height, plot_BMI, ncol = 2)
grid.arrange(plot_sex, plot_typestroke, plot_alco, ncol = 3)

