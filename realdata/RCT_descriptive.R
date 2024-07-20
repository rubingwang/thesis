library(readxl)
library(dplyr)
library(effects)
library(effectsize)
library(emmeans)

# Read the RCT sample data
df_rct <- read_excel("~/tmp/thesis/realdata/RCT_sample.xlsx")

# Linear model
lm <- aov(d_ankle ~ group*sex  , data=df_rct)
summary(lm)

# Estimated marginal means
emmeans(lm, ~ group*sex  , data=df_rct)

# Tukey's multiple comparison test
mulcompare <- TukeyHSD(lm)
mulcompare

# Plotting the effects of treatment on ankle dorsiflexion
plot(allEffects(lm), 
     main = "Sex-adjusted treatment effect of TT",
     ylab = "Δ Ankle dorsiflexion (T1-T0)")

# Plot the frequency distribution histogram for age
ggplot(df_rct, aes(x = age)) +
  geom_histogram(binwidth = 5, fill = "skyblue", color = "black") +
  labs(title = "Age Distribution", x = "Age", y = "Frequency") +
  theme_minimal()

# Plot the frequency distribution histogram for weight
ggplot(df_rct, aes(x = weight)) +
  geom_histogram(binwidth = 5, fill = "lightgreen", color = "black") +
  labs(title = "Weight Distribution", x = "Weight (kg)", y = "Frequency") +
  theme_minimal()

# Plot the frequency distribution histogram for height
ggplot(df_rct, aes(x = height)) +
  geom_histogram(binwidth = 5, fill = "salmon", color = "black") +
  labs(title = "Height Distribution", x = "Height (cm)", y = "Frequency") +
  theme_minimal()

# Plot the frequency distribution histogram for BMI
ggplot(df_rct, aes(x = BMI)) +
  geom_histogram(binwidth = 1, fill = "purple", color = "black") +
  labs(title = "BMI Distribution", x = "BMI", y = "Frequency") +
  theme_minimal()

# Plot the bar plot for stroke type distribution
ggplot(df_rct, aes(x = type_stroke)) +
  geom_bar(fill = "orange") +
  labs(title = "Stroke Type Distribution", x = "Stroke Type", y = "Frequency") +
  theme_minimal()

# Plot the bar plot for history of alcoholism distribution
ggplot(df_rct, aes(x = history_alcoholism)) +
  geom_bar(fill = "blue") +
  labs(title = "History of Alcoholism Distribution", x = "History of Alcoholism", y = "Frequency") +
  theme_minimal()

