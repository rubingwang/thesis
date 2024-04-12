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

# number of repetitions in simulation
repetitions = 100 

########  Simulation：well-specified ######## 
results <- compute_estimators_and_store(rep = repetitions, n = 1000, m = 10000)

ggplot(data = melt(results), aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill=variable)) +
  theme_bw() +
  geom_hline(aes(yintercept = 27.4, color = "True ATE"), 
             size = 0.6, linetype="dashed") +
  xlab("") +
  ylab("Estimated ATE")  +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=14)) +  
  theme(axis.text = element_text(angle = 0, vjust = 0.5, 
                                 hjust=1, size=14)) +          
  scale_fill_brewer(palette = "Paired") +
  coord_flip() +
  ggsave("sim-simple-100.png", width = 7, height = 8, dpi = 400)

# performance_well <- compute_metrics(results,27.4)

########  Simulation：mis-specified ########
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

data <- total_results[2:nrow(total_results),]
data$relative.size <- ifelse(data$m == 10000, "10%", "other")

DF <- melt(data , 
           id.vars = c("param_RCT", "outcome", "relative.size"), 
           measure.vars =  c("Naive.OnlyRCT", "IPSW", "IPSW.norm", 
              "IPSW.strat.n.2", "IPSW.strat.n.3", 
              "IPSW.strat.n.4", "IPSW.strat.n.5", 
              "IPSW.strat.n.6", "IPSW.strat.n.7", 
              "IPSW.strat.n.8", "IPSW.strat.n.9", 
              "IPSW.strat.n.10", "G.formula"))

DF$param_RCT <- ifelse(DF$param_RCT == "correct", 
                       "Well-specified trial selection model", 
                       ifelse(DF$param_RCT == "exponential", 
                              "Mis-specified trial selection model", 
                              "RCT with strong shift"))

DF$outcome <- ifelse(DF$outcome == "correct", "Well-speciifed outcome model", "Mis-specified outcome model")
DF$outcome <- factor(DF$outcome, levels = c("Well-speciifed outcome model", "Mis-specified outcome model"))

ggplot(data = DF[DF$relative.size == "10%" & DF$param_RCT != "RCT with strong shift",], 
       aes(y = variable, x = value)) +
  geom_boxplot(aes(fill=variable)) +
  facet_grid(outcome~param_RCT) +
  theme_bw() +
  geom_vline(aes(xintercept = 27.4, color = "True ATE"), size = 0.6, linetype = "dashed")+
  xlab("Estimated ATE") +
  ylab("")  +
  theme(legend.title = element_blank(), 
        legend.position="bottom", legend.box = "horizontal") +  # no title in legend
  theme(axis.text = element_text(angle = 45, vjust = 0.5, hjust=1, size=10)) + 
  scale_fill_brewer(palette = "Paired") +
  ggsave("sim-RCT-outcome-mis.png", width = 8.5, height = 9.5, dpi = 500)


# When Y is correct, Strong VS Weak Shift
data$relative.size <- ifelse(data$m == 10000, "10%", "other")

DF <- melt(data, id.vars = c("param_RCT", "outcome", "relative.size"),
           measure.vars =  c("Naive.OnlyRCT", "IPSW", "IPSW.norm", 
                             "IPSW.strat.n.2", "IPSW.strat.n.3", 
                             "IPSW.strat.n.4", "IPSW.strat.n.5", 
                             "IPSW.strat.n.6", "IPSW.strat.n.7", 
                             "IPSW.strat.n.8", "IPSW.strat.n.9", 
                             "IPSW.strat.n.10", "G.formula"))

DF$param_RCT <- ifelse(DF$param_RCT == "correct", 
                       "Shift: Weak", 
                       ifelse(DF$param_RCT == "exponential", 
                              "RCT mis-specification", 
                              "Shift: Strong"))
DF$outcome <- ifelse(DF$outcome == "correct", "Correct Y", "Mis-specified Y")

ggplot(data = DF[DF$relative.size == "10%" & 
                   DF$param_RCT !="RCT mis-specification"& 
                   DF$outcome == "Correct Y",], 
       aes(x = variable, y = value)) +
  geom_boxplot(aes(fill=variable)) +
  facet_wrap(~param_RCT) +
  #     #geom_jitter(alpha = 0.2, size = 0.2, width = 0.2)  +
  theme_bw() +
  geom_hline(aes(yintercept = 27.4, color = "True ATE"), size = 0.6, linetype="dashed") +
  xlab("") +
  ylab("Estimated ATE")  +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=9), legend.position="bottom") +
  theme(axis.text = element_text(angle = 0, vjust = 0.5, hjust=1, size=14)) +
  scale_fill_brewer(palette = "Paired") +
  coord_flip() +
  theme(strip.background=element_rect(fill=NA, color=NA),
        strip.text=element_text(size=15, face = "bold")) +
  ggsave("sim-strong-shift.png", width = 9, height = 9.5, dpi=500)


#########################################################
########  Scenario2: distribution shift analysis ########
#########################################################

########  X1 shift: descriptive ########
# 1.well-specification simulation (including all 4 covariate)
one_simulation <- simulate_continuous(n = 1000, m = 10000)
one_simulation$sample <- ifelse(one_simulation$V == 1, "RCT", "Observational")
baseline_weak <- table1(~ X1 + X2 + X3 + X4 | sample, data = one_simulation, overall="Total")
# 2.mis-specification simulation (including all 4 covariate)
one_shifted_simulation <- simulate_continuous(n = 1000, m = 10000, misRCT = "strongbias")
one_shifted_simulation $sample <- ifelse(one_shifted_simulation $V == 1, "RCT", "Observational")
baseline_strong <- table1(~ X1 + X2 + X3 + X4 | sample, data = one_shifted_simulation, overall="Total")

# Distribution shift plot on X1
one_shifted_simulation$Shift <- rep("X1 Shift: Strong", nrow(one_shifted_simulation))
one_simulation$Shift <- rep("X1 Shift: Weak", nrow(one_simulation))
shift_comparison <- rbind(one_simulation, one_shifted_simulation)

ggplot(shift_comparison, aes(x = X1, group = sample, fill = sample)) +
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


####

data$relative.size <- ifelse(data$m == 10000, "10%", "other")

DF <- melt(data, id.vars = c("param_RCT", "outcome", "relative.size"),
           measure.vars = c("RCT", "IPSW", "IPSW.norm",
                            "Stratification.n.10", "G.formula", "AIPSW"))

DF$param_RCT <- ifelse(DF$param_RCT == "correct", 
                       "Shift: Weak", 
                       ifelse(DF$param_RCT == "exponential", 
                              "RCT mis-specification", 
                              "Shift: Strong"))
DF$outcome <- ifelse(DF$outcome == "correct", "Correct Y", "Mis-specified Y")



ggplot(data = DF[DF$relative.size == "10%" & 
                   DF$param_RCT !="RCT mis-specification"& 
                   DF$outcome == "Correct Y",], 
       aes(x = variable, y = value)) +
  geom_boxplot(aes(fill=variable)) +
  facet_wrap(~param_RCT) +
  #     #geom_jitter(alpha = 0.2, size = 0.2, width = 0.2)  +
  theme_bw() +
  geom_hline(aes(yintercept = 27.4, color = "Ture ATE"), size = 0.6, linetype="dashed") +
  xlab("") +
  ylab("Estimated ATE")  +
  theme(legend.title = element_blank(), 
        legend.text = element_text(size=9), legend.position="bottom") +
  theme(axis.text = element_text(angle = 0, vjust = 0.5, hjust=1, size=14)) +
  scale_fill_brewer(palette = "Paired") +
  coord_flip() +
  theme(strip.background=element_rect(fill=NA, color=NA),
        strip.text=element_text(size=15, face = "bold")) +
  ggsave(".sim-strong-shift.png")















# X1 effect


rct_ate <- c()
ipsw <- c()
ipsw_x1_only <- c()
ipsw_wo_x1 <- c()
gformula <- c()

for (i in 1:repetitions){
  DF <- simulate_continuous(n = 1000, m = 10000)
  
  # naive estimator
  rct_ate <- c(rct_ate, 
               mean(DF[DF$A == 1 & DF$V == 1, "Y"]) - 
                 mean(DF[DF$A == 0  & DF$V == 1, "Y"]))
  
  #ipsw
  ipsw  <- c(ipsw, compute_ipsw(DF, normalized = FALSE))
  
  #ipsw with X1 only
  ipsw_x1_only <- c(ipsw_x1_only, compute_ipsw(DF, normalized = FALSE, covariates = "X1"))
  
  #ipsw without X1
  ipsw_wo_x1 <- c(ipsw_wo_x1, compute_ipsw(DF, normalized = FALSE, covariates = "-X1"))
  
  #gformula
  gformula <- c(gformula, compute_gformula(DF))
  
}

results_ipsw <- data.frame("RCT" = rct_ate,
                           "IPSW" = ipsw,
                           "IPSW-X1" = ipsw_x1_only,
                           "IPSW-without-X1" = ipsw_wo_x1,
                           "G.formula" = gformula)




ggplot(data = melt(results_ipsw), aes(x = variable, y = value)) + 
  geom_boxplot(aes(fill=variable)) +
  theme_bw() +
  geom_hline(aes(yintercept = 27.4, color = "Population ATE"), 
             size = 0.6, linetype="dashed") +
  xlab("") +
  ylab("Estimated ATE")  +
  theme(legend.title = element_blank(), legend.text = element_text(size=14)) +  
  theme(axis.text = element_text(angle = 0, vjust = 0.5, hjust=1, size=14)) +          
  coord_flip() +
  ggsave("X1variation.png", dpi=500)


