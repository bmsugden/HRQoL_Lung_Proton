rm(list = ls())

# Load imputed and matched datasets  
matched_data_list <- readRDS("matched_data_list.rds")

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


####### BACKWARD ELIMINATION PROCESS (steps 1-6 below)
###### Note, anova() not compatible with MICE imputation so p-values cannot be directly derived for each FE variable
###### Hence, Wald Test used in line with https://stefvanbuuren.name/fimd/sec-stepwise.html (5.4.2)

# Step 1: Fit model with all candidate fixed effects variables
EQindex1 <- lapply(matched_data_list, function(df) {
  lmer(
    EQindex ~ Treatment30 
    #Round7+ Stage 
    #Round12+ Pulmonary_comorbidity 
    #Round1+ Durvalumab 
    #Round8+ Chemo_sequence 
    #Round6+ Smoking_status 
    #Round10+ Surgery 
    #Round5+ Adapted 
    #Round14+ TumorLocation 
    #Round15+ TumorCellType 
    #Round4+ Radiation_DoseReceived 
    #Round2+ Gender 
    #Round3+ Age_startRT 
    #Round15+ WHO_PS 
    + BaselineEQ5D 
    #Round11+ BaselineEQ5D*Treatment30 
    #Round9+ Dysphagia_BL_PHYS_CTCAE
    + Dyspnea_PHYS_BL_CTCAE 
    #Round13+ GTV_volume
    + (1 | PatientID),
    data   = df,
    REML   = FALSE
  )
})
# coerce to mira (needed for D1 Wald Test)
class(EQindex1) <- c("mira", "list")
# pool
pool(EQindex1)
# summary
summary(pool(EQindex1), conf.int = TRUE)

# Step 2: Remove one variable at a time (this can be done by placing a # in front of the variable)
EQindex2 <- lapply(matched_data_list, function(df) {
  lmer(
    EQindex ~ Treatment30 
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
    #+ WHO_PS 
    + BaselineEQ5D 
    #+ BaselineEQ5D*Treatment30 
    #+ Dysphagia_BL_PHYS_CTCAE
    + Dyspnea_PHYS_BL_CTCAE 
    #+ GTV_volume
    + (1 | PatientID),
    data   = df,
    REML   = FALSE
  )
})
pool(EQindex2)
class(EQindex2) <- c("mira", "list")
summary(pool(EQindex2), conf.int = TRUE)

# Step 3: Run Wald Test to assess significance of removing that variable
D1(EQindex1, EQindex2)
# Step 4: Make a note of the p-value and repeat steps 2 and 3 for each candidate variable (removing one at a time)
## Keep Treatment30 forced in to all models

# Step 5: Identify the least significant variable and remove from EQindex1 and EQindex2 (leave # in front of variable for transparency)
# Step 6: Repeat this process until only significant variables remain (alpha = 0.05) 
#OPTIONAL STEP: Repeat steps 1-6 using "EQ5DVAS" or "GlobalHealthStatus" as outcome variables

## END Selection

### Final Models (base-case) - based on EQindex variable selection
#EQ5D5L index
EQindex_Final <- lapply(matched_data_list, function(df) {
  lmer(
    EQindex ~ Treatment30 
    + BaselineEQ5D 
    + Dyspnea_PHYS_BL_CTCAE
    + (1 | PatientID),
    data   = df,
    REML   = FALSE
  )
})
pool(EQindex_Final)
class(EQindex_Final) <- c("mira", "list")
summary(pool(EQindex_Final), conf.int = TRUE)
#EQ-VAS
EQVAS_Final <- lapply(matched_data_list, function(df) {
  lmer(
    EQ5DVAS ~ Treatment30 
    + BaselineEQVAS 
    + Dyspnea_PHYS_BL_CTCAE
    + (1 | PatientID),
    data   = df,
    REML   = FALSE
  )
})
pool(EQVAS_Final)
class(EQVAS_Final) <- c("mira", "list")
summary(pool(EQVAS_Final), conf.int = TRUE)
#EORTC QLQ-C3-: Global Health Status
GHS_Final <- lapply(matched_data_list, function(df) {
  lmer(
    GlobalHealthStatus ~ Treatment30 
    + BaselineGHS 
    + Dyspnea_PHYS_BL_CTCAE
    + (1 | PatientID),
    data   = df,
    REML   = FALSE
  )
})
pool(GHS_Final)
class(GHS_Final) <- c("mira", "list")
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