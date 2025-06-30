rm(list = ls())

# Load necessary libraries
library(summarytools)
library(ggplot2)
library(car)
library(PerformanceAnalytics)
library(dplyr)
library(twang)
library(broom.mixed)
library(mice)
library(miceadds)
library(lme4)
library(lmerTest)

options(scipen=999)
# Load weighted datasets
weighted_data_list <- readRDS("weighted_data_list.rds")
mids_data <- datlist2mids(weighted_data_list)


####### BACKWARD ELIMINATION PROCESS (steps 1-6 below) - BASE CASE ANALYSIS
###### Note, anova() not compatible with MICE imputation so p-values cannot be directly derived for each FE variable
###### Hence, Wald Test used in line with https://stefvanbuuren.name/fimd/sec-stepwise.html (5.4.2)

# Step 1: Fit model with all candidate fixed effects variables
EQindex1 <- with(data = mids_data, expr = lme4::lmer(EQindex ~ Treatment30 
                                                     #Round6+ Stage 
                                                     #Round12+ Pulmonary_comorbidity 
                                                     #Round3+ Durvalumab 
                                                     #Round5+ Chemo_sequence 
                                                     #Round8+ Smoking_status 
                                                     #Round9+ Surgery 
                                                     #Round7+ Adapted 
                                                     #Round14+ TumorLocation 
                                                     #Round15+ TumorCellType 
                                                     #Round4+ Radiation_DoseReceived 
                                                     #Round1+ Gender 
                                                     #Round2+ Age_startRT 
                                                     + WHO_PS 
                                                     + BaselineEQ5D 
                                                     #Round11+ BaselineEQ5D*Treatment30 
                                                     #Round10+ Dysphagia_BL_PHYS_CTCAE
                                                     + Dyspnea_PHYS_BL_CTCAE 
                                                     #Round13+ GTV_volume
                                                     + (1 | PatientID),
                                                     weights = twang_weight, 
                                                     REML = FALSE))
pool(EQindex1)
summary(pool(EQindex1))

# Step 2: Remove one variable at a time (this can be done by placing a # in front of the variable)
EQindex2 <- with(data = mids_data, expr = lme4::lmer(EQindex ~ Treatment30 
                                                     #+ Stage 
                                                     #+ Pulmonary_comorbidity 
                                                     #+ Durvalumab 
                                                     #+ Chemo_sequence 
                                                     #+ Smoking_status 
                                                     #+ Surgery 
                                                     #+ Adapted 
                                                     #+ TumorLocation 
                                                     #+ TumorCellType 
                                                     #+ Radiation_DoseReceived 
                                                     #+ Gender 
                                                     #+ Age_startRT 
                                                     + WHO_PS 
                                                     + BaselineEQ5D 
                                                     #+ BaselineEQ5D*Treatment30 
                                                     #+ Dysphagia_BL_PHYS_CTCAE
                                                     + Dyspnea_PHYS_BL_CTCAE 
                                                     #+ GTV_volume
                                                     + (1 | PatientID),
                                                     weights = twang_weight, 
                                                     REML = FALSE))
pool(EQindex2)
summary(pool(EQindex2))

# Step 3: Run Wald Test to assess significance of removing that variable
D1(EQindex1, EQindex2)
# Step 4: Make a note of the p-value and repeat steps 2 and 3 for each candidate variable (removing one at a time)
## Keep Treatment30 forced in to all models

# Step 5: Identify the least significant variable and remove from EQindex1 and EQindex2 (leave # in front of variable for transparency)
# Step 6: Repeat this process until only significant variables remain (alpha = 0.05) 
#OPTIONAL STEP: Repeat steps 1-6 using "EQ5DVAS" or "GlobalHealthStatus" as outcome variables

## END Selection

### Final Models (base-case) - based on EQindex variable selection
# EQ5D5L
EQindex_Final <- with(data = mids_data, expr = lme4::lmer(EQindex ~ Treatment30 
                                                          + WHO_PS 
                                                          + BaselineEQ5D
                                                          + Dyspnea_PHYS_BL_CTCAE
                                                          + (1 | PatientID),
                                                          weights = twang_weight, 
                                                          REML = FALSE))
pool(EQindex_Final)
summary(pool(EQindex_Final), conf.int = TRUE)

#EQ-VAS
EQVAS_Final <- with(data = mids_data, expr = lme4::lmer(EQ5DVAS ~ Treatment30 
                                                          + WHO_PS 
                                                          + BaselineEQVAS
                                                          + Dyspnea_PHYS_BL_CTCAE
                                                          + (1 | PatientID),
                                                          weights = twang_weight, 
                                                          REML = FALSE))
pool(EQVAS_Final)
summary(pool(EQVAS_Final), conf.int = TRUE)

#EORTC QLQ-C30 Global Health Status
GHS_Final <- with(data = mids_data, expr = lme4::lmer(GlobalHealthStatus ~ Treatment30 
                                                          + WHO_PS 
                                                          + BaselineGHS
                                                          + Dyspnea_PHYS_BL_CTCAE
                                                          + (1 | PatientID),
                                                          weights = twang_weight, 
                                                          REML = FALSE))
pool(GHS_Final)
summary(pool(GHS_Final), conf.int = TRUE)



####Asessing final model
weighted_models <- EQindex_Final$analyses #Adjust EQindex_Final to assess other models
mod1 <- weighted_models[[1]] #Adjust number in brackets to view models for different imputed datasets
# Residuals vs Fitted
plot(mod1, main="Imputation 1: Residuals vs Fitted")
# QQ plot of residuals
qqnorm(resid(mod1)); qqline(resid(mod1), col="red")
# Histogram of residuals
hist(resid(mod1), breaks=30,
     main="Imputation 1: Residual Distribution",
     xlab="Residual")
# Fitted vs Residual scatter
plot(fitted(mod1), resid(mod1),
     main="Imputation 1: Fitted vs Residuals",
     xlab="Fitted values", ylab="Residuals")



### Scenario analysis - proron = <50% fractions as protons
# EQ5D5L
EQindex_SA50 <- with(data = mids_data, expr = lme4::lmer(EQindex ~ Treatment50 
                                                          + WHO_PS 
                                                          + BaselineEQ5D
                                                          + Dyspnea_PHYS_BL_CTCAE
                                                          + (1 | PatientID),
                                                          weights = twang_weight, 
                                                          REML = FALSE))
pool(EQindex_SA50)
summary(pool(EQindex_SA50), conf.int = TRUE)

#EQ-VAS
EQVAS_SA50 <- with(data = mids_data, expr = lme4::lmer(EQ5DVAS ~ Treatment50 
                                                        + WHO_PS 
                                                        + BaselineEQVAS
                                                        + Dyspnea_PHYS_BL_CTCAE
                                                        + (1 | PatientID),
                                                        weights = twang_weight, 
                                                        REML = FALSE))
pool(EQVAS_SA50)
summary(pool(EQVAS_SA50), conf.int = TRUE)

#EORTC QLQ-C30 Global Health Status
GHS_SA50 <- with(data = mids_data, expr = lme4::lmer(GlobalHealthStatus ~ Treatment50 
                                                      + WHO_PS 
                                                      + BaselineGHS
                                                      + Dyspnea_PHYS_BL_CTCAE
                                                      + (1 | PatientID),
                                                      weights = twang_weight, 
                                                      REML = FALSE))
pool(GHS_SA50)
summary(pool(GHS_SA50), conf.int = TRUE)

### Scenario analysis - proton = <80% fractions as protons
# EQ5D5L
EQindex_SA80 <- with(data = mids_data, expr = lme4::lmer(EQindex ~ Treatment80 
                                                         + WHO_PS 
                                                         + BaselineEQ5D
                                                         + Dyspnea_PHYS_BL_CTCAE
                                                         + (1 | PatientID),
                                                         weights = twang_weight, 
                                                         REML = FALSE))
pool(EQindex_SA80)
summary(pool(EQindex_SA80), conf.int = TRUE)

#EQ-VAS
EQVAS_SA80 <- with(data = mids_data, expr = lme4::lmer(EQ5DVAS ~ Treatment80 
                                                       + WHO_PS 
                                                       + BaselineEQVAS
                                                       + Dyspnea_PHYS_BL_CTCAE
                                                       + (1 | PatientID),
                                                       weights = twang_weight, 
                                                       REML = FALSE))
pool(EQVAS_SA80)
summary(pool(EQVAS_SA80), conf.int = TRUE)

#EORTC QLQ-C30 Global Health Status
GHS_SA80 <- with(data = mids_data, expr = lme4::lmer(GlobalHealthStatus ~ Treatment80 
                                                     + WHO_PS 
                                                     + BaselineGHS
                                                     + Dyspnea_PHYS_BL_CTCAE
                                                     + (1 | PatientID),
                                                     weights = twang_weight, 
                                                     REML = FALSE))
pool(GHS_SA80)
summary(pool(GHS_SA80), conf.int = TRUE)


### Scenario analysis - proton = <100% fractions as protons

# EQ5D5L
EQindex_SA100 <- with(data = mids_data, expr = lme4::lmer(EQindex ~ Treatment100 
                                                         + WHO_PS 
                                                         + BaselineEQ5D
                                                         + Dyspnea_PHYS_BL_CTCAE
                                                         + (1 | PatientID),
                                                         weights = twang_weight, 
                                                         REML = FALSE))
pool(EQindex_SA100)
summary(pool(EQindex_SA100), conf.int = TRUE)

#EQ-VAS
EQVAS_SA100 <- with(data = mids_data, expr = lme4::lmer(EQ5DVAS ~ Treatment100 
                                                       + WHO_PS 
                                                       + BaselineEQVAS
                                                       + Dyspnea_PHYS_BL_CTCAE
                                                       + (1 | PatientID),
                                                       weights = twang_weight, 
                                                       REML = FALSE))
pool(EQVAS_SA100)
summary(pool(EQVAS_SA100), conf.int = TRUE)

#EORTC QLQ-C30 Global Health Status
GHS_SA100 <- with(data = mids_data, expr = lme4::lmer(GlobalHealthStatus ~ Treatment100 
                                                     + WHO_PS 
                                                     + BaselineGHS
                                                     + Dyspnea_PHYS_BL_CTCAE
                                                     + (1 | PatientID),
                                                     weights = twang_weight, 
                                                     REML = FALSE))
pool(GHS_SA100)
summary(pool(GHS_SA100), conf.int = TRUE)


### Plotting the trend of scenario analyses
# Extract coefficients for protons
get_coef <- function(pooled_model, pattern) {
  su <- summary(pooled_model)
  idx <- grep(pattern, su$term)
  if(length(idx) != 1) {
    stop("Did not find exactly one term matching ", pattern)
  }
  su$estimate[idx]
}
eq5d_est <- c(
  get_coef(pool(EQindex_Final),  "Treatment30"),
  get_coef(pool(EQindex_SA50),  "Treatment50"),
  get_coef(pool(EQindex_SA80),  "Treatment80"),
  get_coef(pool(EQindex_SA100), "Treatment100")
)
eqvas_est <- c(
  get_coef(pool(EQVAS_Final), "Treatment30"),
  get_coef(pool(EQVAS_SA50), "Treatment50"),
  get_coef(pool(EQVAS_SA80), "Treatment80"),
  get_coef(pool(EQVAS_SA100),"Treatment100")
)
ghs_est <- c(
  get_coef(pool(GHS_Final),   "Treatment30"),
  get_coef(pool(GHS_SA50),   "Treatment50"),
  get_coef(pool(GHS_SA80),   "Treatment80"),
  get_coef(pool(GHS_SA100),  "Treatment100")
)

# Build the 4×3 data.frame
scenario <- factor(c("30%","50%","80%","100%"),
                   levels = c("30%","50%","80%","100%"))
coef_df <- data.frame(
  scenario,
  EQ5D5L = eq5d_est,
  EQVAS  = eqvas_est,
  GHS    = ghs_est
)

#plot line graph
ggplot() +
  # EQ5D5L
  geom_line(data = coef_df,
            aes(x = scenario, y = EQ5D5L, group = 1, color = "EQ5D5L"),
            size = 1) +
  geom_point(data = coef_df,
             aes(x = scenario, y = EQ5D5L, color = "EQ5D5L"),
             size = 2) +
  # EQVAS
  geom_line(data = coef_df,
            aes(x = scenario, y = EQVAS, group = 1, color = "EQVAS"),
            size = 1) +
  geom_point(data = coef_df,
             aes(x = scenario, y = EQVAS, color = "EQVAS"),
             size = 2) +
  # GHS
  geom_line(data = coef_df,
            aes(x = scenario, y = GHS, group = 1, color = "GHS"),
            size = 1) +
  geom_point(data = coef_df,
             aes(x = scenario, y = GHS, color = "GHS"),
             size = 2) +
  # bold horizontal line at y = 0
  geom_hline(yintercept = 0,
             size        = 1.2,      # make it “bold”
             color       = "black",  # or any other colour
             linetype    = "solid") +
  # Legend & labels
  scale_color_discrete(name = "Outcome") +
  labs(
    x = "Proton fractions (%)",
    y = "Estimated Proton coefficient"
  ) +
  theme_minimal()





