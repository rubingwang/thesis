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
# Creat file folder
#dir.create("plot")
#dir.create("metrics")

# Initial setting for stimulation studies
repetitions = 100  # numbers of estimation 
theta = 32       # true ATE (of the target population)
# Note: can also change the mis- functions and the constant in log's 'mis-selection'


#################################################
######  Simulation1：systematic analysis ########
#################################################

RCT_param <- c("correct", "strongbias", "exponential")
Outcome_param <- c("correct", "wrong")

total_results <- compute_estimators_and_store(rep = repetitions)
total_results$n = 1000
total_results$m = 49000
total_results$param_RCT = "correct"
total_results$outcome = "correct"

for (m in c(10000)){
  for (rct_param in RCT_param){
    for (outcome_param in Outcome_param){
      results <- compute_estimators_and_store(rep = repetitions, 
                                              n = 1000, m = m, 
                                              misRCT = rct_param, 
                                              misoutcome = outcome_param)
      results$n <- rep(1000, nrow(results))
      results$m <- rep(m, nrow(results))
      results$param_RCT <- rep(rct_param, nrow(results))
      results$outcome <- rep(outcome_param, nrow(results))
      total_results <- rbind(total_results, results)
    }
  }
}
data <- total_results
data$relative.size <- ifelse(data$m == 10000, "10%", "other")

DF <- melt(data , 
           id.vars = c("param_RCT", "outcome", "relative.size"), 
           measure.vars =  c("Naive.OnlyRCT", "IPSW", "IPSW.norm", 
              "IPSW.strat.n.2", "IPSW.strat.n.3", 
              "IPSW.strat.n.4", "IPSW.strat.n.5", 
              "IPSW.strat.n.6", "IPSW.strat.n.7", 
              "IPSW.strat.n.8", "IPSW.strat.n.9", 
              "IPSW.strat.n.10","IPSW.randomforest", "G.formula"))

DF$param_RCT <- ifelse(DF$param_RCT == "correct", 
                       "Well-specified trial selection model", 
                       ifelse(DF$param_RCT == "exponential", 
                              "Mis-specified trial selection model", 
                                     "RCT with X1 strong shift"))

DF$outcome <- ifelse(DF$outcome == "correct", "Well-speciifed outcome model", "Mis-specified outcome model")
DF$outcome <- factor(DF$outcome, levels = c("Well-speciifed outcome model", "Mis-specified outcome model"))

ggplot(data = DF[DF$relative.size == "10%" & 
                   (DF$param_RCT == "Well-specified trial selection model" | 
                      DF$param_RCT == "Mis-specified trial selection model"),], 
       aes(y = variable, x = value)) +
  geom_boxplot(aes(fill=variable)) +
  facet_grid(outcome~param_RCT) +
  theme_bw() +
  geom_vline(aes(xintercept = theta, color = "True ATE"), size = 0.6, linetype = "dashed")+
  xlab("Estimated ATE") +
  ylab("")  +
  theme(legend.title = element_blank(), 
        legend.position="bottom", legend.box = "horizontal") +  # no title in legend
  theme(axis.text = element_text(angle = 45, vjust = 0.5, hjust=1, size=10)) + 
  scale_fill_brewer(palette = "Paired") +
  xlim(0, 50)+
  #annotate("text", x = 0, y = Inf, label = "a)", hjust = 1, vjust = 1, size = 5) +
  ggsave("plot/sim-RCT-outcome-mis.pdf", width = 8.5, height = 10)

# Print performance metrics
# 1.when selection is well-specified, outcome is well-specified
df_wellselection_welloutcome <- filter(data=data, param_RCT = 'correct', outcome = 'correct')
performance_wellselection_welloutcome <- compute_metrics(df_wellselection_welloutcome, theta, output=T)
write.csv(performance_wellselection_welloutcome, file = "metrics/performance_wellselection_welloutcome.csv", row.names = FALSE)

# 2.when selection is well-specified, outcome is mis-specified
df_wellselection_badoutcome <- filter(data=data, param_RCT = 'correct', outcome = 'wrong')
performance_wellselection_badoutcome <- compute_metrics(df_wellselection_badoutcome, theta, output=T)
write.csv(performance_wellselection_badoutcome, file = "metrics/performance_wellselection_badoutcome.csv", row.names = FALSE)

# 3.when selection is mis-specified, outcome is well-specified
df_badselection_welloutcome <- filter(data=data, param_RCT = 'exponential', outcome = 'correct')
performance_badselection_welloutcome <- compute_metrics(df_badselection_welloutcome, theta, output=T)
write.csv(performance_badselection_welloutcome, file = "metrics/performance_badselection_welloutcome.csv", row.names = FALSE)

# 4.when selection is mis-specified, outcome is mis-specified
df_badselection_badoutcome <- filter(data=data, param_RCT = 'exponential', outcome = 'wrong')
performance_badselection_badoutcome <- compute_metrics(df_badselection_badoutcome, theta, output=T)
write.csv(performance_badselection_badoutcome, file = "metrics/performance_badselection_badoutcome.csv", row.names = FALSE)


#################################################
########### Simulation2: X1 Shift ###############
#################################################

# this part focus on the X1 shift, when the outcome model is well-specified
DF2 <- melt(data, id.vars = c("param_RCT", "outcome", "relative.size"),
           measure.vars =  c("Naive.OnlyRCT", "IPSW", "IPSW.norm", 
                             "IPSW.strat.n.2", "IPSW.strat.n.3", 
                             "IPSW.strat.n.4", "IPSW.strat.n.5", 
                             "IPSW.strat.n.6", "IPSW.strat.n.7", 
                             "IPSW.strat.n.8", "IPSW.strat.n.9", 
                             "IPSW.strat.n.10", "G.formula"))

DF2$param_RCT <- ifelse(DF2$param_RCT == "correct", 
                       "X1 shift: weak", 
                       ifelse(DF2$param_RCT == "exponential", 
                              "Mis-specified trial selection model", 
                                     "X1 shift: strong"))

DF2$outcome <- ifelse(DF2$outcome == "correct", "Correct Y", "Mis-specified Y")

ggplot(data = DF2[DF2$relative.size == "10%" & 
                    DF2$outcome == "Correct Y"&
                   (DF2$param_RCT == "X1 shift: weak" | 
                      DF2$param_RCT == "X1 shift: strong"),],
       aes(x = variable, y = value)) +
  geom_boxplot(aes(fill=variable)) +
  facet_wrap(~param_RCT) +
  ##geom_jitter(alpha = 0.2, size = 0.2, width = 0.2)  +
  theme_bw() +
  geom_hline(aes(yintercept = 31.8, color = "True ATE"), size = 0.6, linetype="dashed") +
  xlab("") +
  ylab("Estimated ATE")  +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=14), legend.position="bottom") +
  theme(axis.text = element_text(angle = 45, vjust = 0.5, hjust=1, size=13)) +
  scale_fill_brewer(palette = "Paired") +
  coord_flip() +
  theme(strip.background=element_rect(fill=NA, color=NA),
        strip.text=element_text(size=15, face = "bold")) +
  #scale_y_continuous(limits = c(-3, 45), breaks = seq(0, 40, by = 5))+
  ylim(0, 50)+
  ggsave("plot/sim-X1strong-shift.pdf", width =11.5, height = 12)

# print performance metrics
# 1.when X1 weak shift
df_X1weak <- filter(data=data, param_RCT = 'correct', outcome = 'correct')
performance_X1weak <- compute_metrics(df_X1weak, theta, output=T)
write.csv(performance_X1weak, file = "metrics/performance_X1weak.csv", row.names = FALSE)
# 2.when X1 strong shift
df_X1strong <- filter(data=data, param_RCT = 'strongbias', outcome = 'correct')
performance_X1strong <- compute_metrics(df_X1strong, theta, output=T)
write.csv(performance_X1strong, file = "metrics/performance_X1strong.csv", row.names = FALSE)

# print the covariate descriptive table 
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
  geom_density(alpha=0.5, position="dodge") +
  scale_fill_manual(values=c("darkorchid4", "darkorange1"))   +
  theme_classic() +
  theme(legend.title = element_blank(), 
        legend.position = "bottom", 
        legend.box = "horizontal", 
        legend.text = element_text(size=14, face="bold")) +
  ylab("") + 
  theme(axis.text = element_text(vjust = 0.5, hjust=1, size=13, face="bold"), axis.title.x = element_text(size=18, face="bold")) +
  facet_grid(~Shift)  +
  theme(strip.background=element_rect(fill=NA, color=NA),
        strip.text=element_text(size=14, face = "bold")) +
  ggsave("plot/X1strongshift.pdf", width =8.5, height = 10)


#################################################
######## Simulation4: X1 effect analysis#########
#################################################

rct_ate <- c()
ipsw <- c()
ipsw_norm <- c()
ipsw_x1_only <- c()
ipsw_x1x2_only <- c()
ipsw_x1x2x3_only <- c()
ipsw_wo_x1 <- c()
gformula <- c()
repetitions=100
for (i in 1:repetitions){
  DF <- simulate_continuous(n = 1000, m = 10000)
  
  # naive estimator
  rct_ate <- c(rct_ate, 
               mean(DF[DF$A == 1 & DF$V == 1, "Y"]) - 
                 mean(DF[DF$A == 0  & DF$V == 1, "Y"]))
  
  #ipsw
  ipsw  <- c(ipsw, compute_ipsw(DF, normalized = F))
  ipsw_norm <- c(ipsw_norm, compute_ipsw(DF, normalized = T))
  
  #ipsw with X1 only
  ipsw_x1_only <- c(ipsw_x1_only, compute_ipsw(DF, normalized = FALSE, covariates = "X1"))
  #ipsw with X1 X2 only
  ipsw_x1x2_only <- c(ipsw_x1x2_only, compute_ipsw(DF, normalized = FALSE, covariates = "X1X2"))
  #ipsw with X1 X2 X3 only
  ipsw_x1x2x3_only <- c(ipsw_x1x2x3_only, compute_ipsw(DF, normalized = FALSE, covariates = "X1X2X3"))
  
  #ipsw without X1
  ipsw_wo_x1 <- c(ipsw_wo_x1, compute_ipsw(DF, normalized = FALSE, covariates = "-X1"))
  
  #gformula
  gformula <- c(gformula, compute_gformula(DF))
  
 
}

#ipsw_x1x2_only <- ipsw_x1x2_only[ipsw_x1x2_only < max(ipsw_x1x2_only)]
#ipsw_x1x2x3_only <- ipsw_x1x2x3_only[ipsw_x1x2x3_only < max(ipsw_x1x2x3_only)]
  
results_ipsw <- data.frame("Naive.OnlyRCT" = rct_ate,
                           "IPSW-norm" = ipsw_norm,
                           "IPSW-X1-X2-X3-X4" = ipsw_x1_only,
                           "IPSW-X1" = ipsw,
                           "IPSW-X1-X2" = ipsw_x1x2x3_only,
                           "IPSW-X1-X2-X3" = ipsw_x1x2_only,
                           "IPSW-without-X1" = ipsw_wo_x1,
                           "G.formula" = gformula)

# reorder variables in results_ipsw for a good plot
results_ipsw2 <- results_ipsw[, c("Naive.OnlyRCT","IPSW.without.X1", 
                                  "IPSW.X1", "IPSW.X1.X2", "IPSW.X1.X2.X3", 
                                  "IPSW.X1.X2.X3.X4",
                                  "G.formula")]

ggplot(data = melt(results_ipsw2), aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill=variable)) +
  theme_bw() +
  geom_hline(aes(yintercept = theta, color = "True PATE"), 
             size = 0.6, linetype="dashed") +
  xlab("") +
  ylab("Estimated PATE")  +
  theme(legend.title = element_blank(), legend.text = element_text(size=14)) +  
  theme(axis.text = element_text(angle = 45, vjust = 0.5, hjust=1, size=14)) +          
  coord_flip() +
  ylim(0, 45)+
  ggsave("plot/sim-X1_variation.pdf", width =7.5, height = 8)

# print performance metric
performance_X1variation <- compute_metrics(results_ipsw2, theta=theta, output=T)
write.csv(performance_X1variation, file = "metrics/performance_X1variation.csv", row.names = FALSE)



########################################################################
######## Simulation5: X1 effect analysis for the outcome model#########
########################################################################

rct_ate <- c()
gformula_wo_x1 <- c()
gformula_x1 <- c()
gformula_x1x2 <- c()
gformula_x1x2x3 <- c()
gformula_x1x2x3x4 <- c()
repetitions=100
for (i in 1:repetitions){
  #set.seed(1+i)
  DF <- simulate_continuous(n = 1000, m = 10000)
  # naive estimator
  rct_ate <- c(rct_ate, 
               mean(DF[DF$A == 1 & DF$V == 1, "Y"]) - 
                 mean(DF[DF$A == 0  & DF$V == 1, "Y"]))

  #gformula
  gformula_wo_x1 <- c(gformula_wo_x1, compute_gformula(DF, covariates = "-X1"))
  gformula_x1 <- c(gformula_x1, compute_gformula(DF, covariates = "X1"))
  gformula_x1x2 <- c(gformula_x1x2, compute_gformula(DF, covariates = "X1X2"))
  gformula_x1x2x3 <- c(gformula_x1x2x3, compute_gformula(DF, covariates = "X1X2X3"))
  gformula_x1x2x3x4 <- c(gformula_x1x2x3x4, compute_gformula(DF, covariates = "All"))
  
}

results_g <- data.frame("Naive.OnlyRCT" = rct_ate,
                        "G.formula.without.X1" = gformula_wo_x1,
                        "G.formula.X1" =  gformula_x1,
                        "G.formula.X1.X2" =  gformula_x1x2, 
                        "G.formula.X1.X2.X3" =  gformula_x1x2x3, 
                        "G.formula.X1.X2.X3.X4" =  gformula_x1x2x3x4 )

# reorder variables in results_ipsw for a good plot
results_g2 <- results_g[, c("Naive.OnlyRCT", "G.formula.without.X1",
                            "G.formula.X1","G.formula.X1.X2",
                            "G.formula.X1.X2.X3","G.formula.X1.X2.X3.X4")]

ggplot(data = melt(results_g2), aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill=variable)) +
  theme_bw() +
  geom_hline(aes(yintercept = 32, color = "True PATE"), 
             size = 0.6, linetype="dashed") +
  xlab("") +
  ylab("Estimated PATE")  +
  theme(legend.title = element_blank(), legend.text = element_text(size=14)) +  
  theme(axis.text = element_text(angle = 45, vjust = 0.5, hjust=1, size=14)) +          
  coord_flip() +
  ylim(0, 45)+
  ggsave("plot/sim-X1_variation_gformula.pdf", width =7.5, height = 8)

# print performance metric
performance_X1variation_g <- compute_metrics(results_g2, theta=32, output=T)
write.csv(performance_X1variation_g, file = "metrics/performance_X1variation_gformula.csv", row.names = FALSE)

