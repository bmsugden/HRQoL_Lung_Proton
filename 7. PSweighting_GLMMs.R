rm(list = ls())

# Load necessary libraries
library(summarytools)
library(ggplot2)
library(car)
library(PerformanceAnalytics)
library(dplyr)
library(survey)
library(twang)
library(broom.mixed)
library(mice)
library(miceadds)
library(mitml)
library(lme4)
library(lmerTest)

options(scipen=999)

# Load weighted datasets
weighted_data_list <- readRDS("weighted_data_list.rds")

mids_data <- datlist2mids(weighted_data_list)

#12April2025 - TESTING (lme4)


# 0) Extract your M completed data‑frames as a LIST
imp_list <- mice::complete(mids_data, action = "all")
# sanity‐check:
str(imp_list)
str(imp_list[[1]])   # must show a data.frame with EQindex and twang_weight columns

# 1) Candidate fixed effects
all_vars <- c(
  "Treatment30", "Stage", "Pulmonary_comorbidity", "Durvalumab",
  "Chemo_sequence", "Smoking_status", "Surgery", "Adapted",
  "TumorLocation", "TumorCellType", "Radiation_DoseReceived",
  "Gender", "Age_startRT", "WHO_PS", "BaselineEQ5D",
  "Dysphagia_BL_PHYS_CTCAE", "Dyspnea_PHYS_BL_CTCAE", "GTV_volume"
)
rand_term <- "(1 | PatientID)"
wt_var    <- "twang_weight"
α_enter   <- 0.05
α_exit    <- 0.05

# 2) Helper to build the right formula
make_formula <- function(vars) {
  rhs <- if (length(vars)) paste(vars, collapse = " + ") else "1"
  as.formula(paste("EQindex ~", rhs, "+", rand_term))
}

# 3) Fit‑one helper: runs lmer() on ONE imputed data.frame
fit_one <- function(df, vars) {
  lmer(
    formula = make_formula(vars),
    data    = df,
    weights = df[[wt_var]]
  )
}

# 4) Turn a vector of vars into a MIRA object
get_mipo <- function(vars) {
  fits <- lapply(imp_list, fit_one, vars = vars)
  as.mira(fits)
}

# 5) Start from the full model
current_vars <- all_vars
mipo_current <- get_mipo(current_vars)

# 6) Bidirectional stepwise
repeat {
  ## 6a) Backward: pooled LRT p for dropping each in‐model term
  rem_p <- if (length(current_vars)) {
    sapply(current_vars, function(v) {
      pool.anova(
        get_mipo(setdiff(current_vars, v)),
        mipo_current
      )[2, "Pr(>Chisq)"]
    })
  } else numeric(0)
  
  ## 6b) Forward: pooled LRT p for adding each out‐of‐model term
  add_cands <- setdiff(all_vars, current_vars)
  add_p <- if (length(add_cands)) {
    sapply(add_cands, function(v) {
      pool.anova(
        mipo_current,
        get_mipo(c(current_vars, v))
      )[2, "Pr(>Chisq)"]
    })
  } else numeric(0)
  
  ## 6c) If any term has drop‑p > α_exit, drop the one with highest p
  if (length(rem_p) && max(rem_p) > α_exit) {
    drop_var    <- names(which.max(rem_p))
    message("Dropping ", drop_var,
            " (p = ", round(rem_p[drop_var], 3), ")")
    current_vars <- setdiff(current_vars, drop_var)
    mipo_current <- get_mipo(current_vars)
    next
  }
  
  ## 6d) Else if any term has add‑p < α_enter, add the one with lowest p
  if (length(add_p) && min(add_p) < α_enter) {
    add_var     <- names(which.min(add_p))
    message("Adding  ", add_var,
            " (p = ", round(add_p[add_var], 3), ")")
    current_vars <- c(current_vars, add_var)
    mipo_current <- get_mipo(current_vars)
    next
  }
  
  ## 6e) Otherwise, stop
  break
}

# 7) Pool & display final model
final_fit <- pool(mipo_current)
summary(final_fit)


















































































# Model selection - stepwise selection using wald test

# Define the number of imputations
m <- length(weighted_data_list)
# Create an empty list to store the models
glm_models <- vector("list", m)
# Loop through each imputed dataset and fit the weighted GLM
for (i in 1:m) {
  # Extract the i-th imputed weighted dataset
  data_i <- weighted_data_list[[i]]
  # Define the survey design using TWANG weights
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  # Fit the weighted GLM model (EQ5D1)
  glm_models[[i]] <- svyglm(
    EQindex ~ Treatment30 + Stage + Pulmonary_comorbidity + Durvalumab +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D1_mira <- as.mira(glm_models)
EQ5D1_pooled <- pool(EQ5D1_mira) # all candidate fixed effects variables
summary(EQ5D1_pooled)

for (i in 1:m) {
  # Extract the i-th imputed weighted dataset
  data_i <- weighted_data_list[[i]]
  # Define the survey design using TWANG weights
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  # Fit the weighted GLM model (EQ5D1)
  glm_models[[i]] <- svyglm(
     ~ Treatment30, 
    design = design_i, family = gaussian()
  )
}
EQ5D1_mira <- as.mira(glm_models)
EQ5D1_pooled <- pool(EQ5D1_mira) # all candidate fixed effects variables
summary(EQ5D1_pooled)

# Fit the weighted GLM model without Stage
glm_models2 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
    design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models2[[i]] <- svyglm(
    EQindex ~ Treatment30 + Pulmonary_comorbidity + Durvalumab +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D2_mira <- as.mira(glm_models2)

# Wald Test - no Stage
D1(EQ5D1_mira, EQ5D2_mira) # => remove stage

# Fit the weighted GLM model without pulmonary comorbidity
glm_models3 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models3[[i]] <- svyglm(
    EQindex ~ Treatment30 + Durvalumab +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D3_mira <- as.mira(glm_models3)

# Wald Test - no Pulmonary comorbidity
D1(EQ5D2_mira, EQ5D3_mira) # => remove Pulmonary comorbidity

# Fit the weighted GLM model without durvalumab
glm_models4 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models4[[i]] <- svyglm(
    EQindex ~ Treatment30 + Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D4_mira <- as.mira(glm_models4)

# Wald Test - no durvalumab
D1(EQ5D3_mira, EQ5D4_mira) # => remove durvalumab

# Fit the weighted GLM model without Chemo sequence
glm_models5 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models5[[i]] <- svyglm(
    EQindex ~ Treatment30 + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D5_mira <- as.mira(glm_models5)

# Wald Test - no chemo sequence
D1(EQ5D4_mira, EQ5D5_mira) # => remove chemo sequence

# Fit the weighted GLM model without smoking status
glm_models6 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models6[[i]] <- svyglm(
    EQindex ~ Treatment30 + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D6_mira <- as.mira(glm_models6)

# Wald Test - no smoking status
D1(EQ5D5_mira, EQ5D6_mira) # => remove smoking status

# Fit the weighted GLM model without surgery
glm_models7 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models7[[i]] <- svyglm(
    EQindex ~ Treatment30 + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D7_mira <- as.mira(glm_models7)

# Wald Test - no surgery
D1(EQ5D6_mira, EQ5D7_mira) # => remove surgery

# Fit the weighted GLM model without Adapted
glm_models8 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models8[[i]] <- svyglm(
    EQindex ~ Treatment30 + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D8_mira <- as.mira(glm_models8)

# Wald Test - no Adapted
D1(EQ5D7_mira, EQ5D8_mira) # => remove Adapted

# Fit the weighted GLM model without tumour location
glm_models9 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models9[[i]] <- svyglm(
    EQindex ~ Treatment30 + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D9_mira <- as.mira(glm_models9)

# Wald Test - no tumour location
D1(EQ5D8_mira, EQ5D9_mira) # => remove tumour location

# Fit the weighted GLM model without tumourcelltype
glm_models10 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models10[[i]] <- svyglm(
    EQindex ~ Treatment30 + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D10_mira <- as.mira(glm_models10)

# Wald Test - no TumorCellType
D1(EQ5D9_mira, EQ5D10_mira) # => remove TumorCellType

# Fit the weighted GLM model without Radiation_DoseReceived
glm_models11 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models11[[i]] <- svyglm(
    EQindex ~ Treatment30 + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D11_mira <- as.mira(glm_models11)

# Wald Test - no Radiation_DoseReceived
D1(EQ5D10_mira, EQ5D11_mira) # => remove Radiation_DoseReceived

# Fit the weighted GLM model without Gender
glm_models12 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models12[[i]] <- svyglm(
    EQindex ~ Treatment30 + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D12_mira <- as.mira(glm_models12)

# Wald Test - no Gender
D1(EQ5D11_mira, EQ5D12_mira) # => remove Gender

# Fit the weighted GLM model without Age_startRT
glm_models13 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models13[[i]] <- svyglm(
    EQindex ~ Treatment30 + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D13_mira <- as.mira(glm_models13)

# Wald Test - no Age_startRT
D1(EQ5D12_mira, EQ5D13_mira) # => remove Age_startRT

# Fit the weighted GLM model without WHO_PS
glm_models14 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models14[[i]] <- svyglm(
    EQindex ~ Treatment30 + BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D14_mira <- as.mira(glm_models14)

# Wald Test - no WHO_PS
D1(EQ5D13_mira, EQ5D14_mira) # => keep WHO_PS (reject EQ5D14_mira)

# Fit the weighted GLM model without BaselineEQ5D
glm_models15 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models15[[i]] <- svyglm(
    EQindex ~ Treatment30 + WHO_PS + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D15_mira <- as.mira(glm_models15)

# Wald Test - no BaselineEQ5D
D1(EQ5D13_mira, EQ5D15_mira) # => keep BaselineEQ5D (reject EQ5D15_mira)

# Fit the weighted GLM model without Dysphagia_BL_PHYS_CTCAE
glm_models16 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models16[[i]] <- svyglm(
    EQindex ~ Treatment30 + WHO_PS + BaselineEQ5D + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D16_mira <- as.mira(glm_models16)

# Wald Test - no Dysphagia_BL_PHYS_CTCAE
D1(EQ5D13_mira, EQ5D16_mira) # => remove Dysphagia_BL_PHYS_CTCAE 


# Fit the weighted GLM model without Dyspnea_PHYS_BL_CTCAE
glm_models17 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models17[[i]] <- svyglm(
    EQindex ~ Treatment30 + WHO_PS + BaselineEQ5D + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQ5D17_mira <- as.mira(glm_models17)

# Wald Test - no GTV_volume
D1(EQ5D16_mira, EQ5D17_mira) # => remove Dyspnea_PHYS_BL_CTCAE 


# Fit the weighted GLM model without GTV_volume
glm_models18 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models18[[i]] <- svyglm(
    EQindex ~ Treatment30 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
EQ5D18_mira <- as.mira(glm_models18)

# Wald Test - no GTV_volume
D1(EQ5D17_mira, EQ5D18_mira) # => remove GTV_volume 

## EQ5D18 = FINAL MODEL
##### END SELECTION

# Final Model
EQ5D18_pooled <- pool(EQ5D18_mira)
summary(EQ5D18_pooled)

###### TESTING the final model - given the order of variable removal was arbitrary
# Fit the weighted GLM model 
glm_modelsTEST <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_modelsTEST[[i]] <- svyglm(
    EQindex ~ Treatment30 + WHO_PS + BaselineEQ5D
    + TumorCellType, # replace this variable and rerun to test order of removal
    design = design_i, family = gaussian()
  )
}
EQ5DTEST_mira <- as.mira(glm_modelsTEST)

# Wald Test 
D1(EQ5DTEST_mira, EQ5D18_mira) # p-value >0.05 suggests additional variable does not significantly improve model fit

# RESULTS OF TESTING: TumorCellType + Radiation_DoseReceived individually both had a significant difference compared to final model
# adding both did not significantly improve vs final model

  



## END Testing


## EQVAS Final Model (selection based on EQ5D5L)
glm_models_EQVAS <- vector("list", m)
for (i in 1:m) {
  # Define the survey design using TWANG weights
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models_EQVAS[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
EQVAS_mira <- as.mira(glm_models_EQVAS)
EQVAS_pooled <- pool(EQVAS_mira)
summary(EQVAS_pooled)


### selection based on EQVAS
#all clinically plausible variables
EQVAS1 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS1[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity + Durvalumab +
  Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
  TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
  BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS1_mira <- as.mira(EQVAS1)
EQVAS1_pooled <- pool(EQVAS1_mira)
summary(EQVAS1_pooled)

#remove stage
EQVAS2 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS2[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Pulmonary_comorbidity + Durvalumab +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS2_mira <- as.mira(EQVAS2)
D1(EQVAS1_mira, EQVAS2_mira) # Wald Test => keep stage => reject EQVAS2_mira

#remove Pulmonary_comorbidity
EQVAS3 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS3[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Durvalumab +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS3_mira <- as.mira(EQVAS3)
D1(EQVAS1_mira, EQVAS3_mira) # Wald Test => keep Pulmonary_comorbidity => reject EQVAS3_mira

#remove Durvalumab
EQVAS4 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS4[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS4_mira <- as.mira(EQVAS4)
D1(EQVAS1_mira, EQVAS4_mira) # Wald Test => remove Durvalumab 

#remove Chemo_sequence
EQVAS5 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS5[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS5_mira <- as.mira(EQVAS5)
D1(EQVAS4_mira, EQVAS5_mira) # Wald Test => keep Chemo_sequence => reject  EQVAS5_mira

#remove Smoking_status
EQVAS6 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS6[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS6_mira <- as.mira(EQVAS6)
D1(EQVAS4_mira, EQVAS6_mira) # Wald Test => remove Smoking_status 

#remove Surgery
EQVAS7 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS7[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS7_mira <- as.mira(EQVAS7)
D1(EQVAS6_mira, EQVAS7_mira) # Wald Test => remove Surgery 

#remove Adapted
EQVAS8 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS8[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS8_mira <- as.mira(EQVAS8)
D1(EQVAS7_mira, EQVAS8_mira) # Wald Test => remove Adapted 

#remove TumorLocation
EQVAS9 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS9[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS9_mira <- as.mira(EQVAS9)
D1(EQVAS8_mira, EQVAS9_mira) # Wald Test => remove TumorLocation

#remove TumorCellType
EQVAS10 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS10[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS10_mira <- as.mira(EQVAS10)
D1(EQVAS9_mira, EQVAS10_mira) # Wald Test => remove TumorCellType

#remove Radiation_DoseReceived
EQVAS11 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS11[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS11_mira <- as.mira(EQVAS11)
D1(EQVAS10_mira, EQVAS11_mira) # Wald Test => keep Radiation_DoseReceived => reject EQVAS11_mira

#remove Gender
EQVAS12 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS12[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Radiation_DoseReceived + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS12_mira <- as.mira(EQVAS12)
D1(EQVAS10_mira, EQVAS12_mira) # Wald Test => remove Gender

#remove Age_startRT
EQVAS13 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS13[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Radiation_DoseReceived + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS13_mira <- as.mira(EQVAS13)
D1(EQVAS12_mira, EQVAS13_mira) # Wald Test => remove Age_startRT

#remove WHO_PS
EQVAS14 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS14[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Radiation_DoseReceived + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS14_mira <- as.mira(EQVAS14)
D1(EQVAS13_mira, EQVAS14_mira) # Wald Test => keep WHO_PS => reject EQVAS14_mira

#remove BaselineEQ5D
EQVAS15 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS15[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Radiation_DoseReceived + WHO_PS + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS15_mira <- as.mira(EQVAS15)
D1(EQVAS13_mira, EQVAS15_mira) # Wald Test => keep BaselineEQ5D => reject EQVAS15_mira

#remove Dysphagia_BL_PHYS_CTCAE
EQVAS16 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS16[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Radiation_DoseReceived + WHO_PS + BaselineEQ5D + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS16_mira <- as.mira(EQVAS16)
D1(EQVAS13_mira, EQVAS16_mira) # Wald Test => remove Dysphagia_BL_PHYS_CTCAE

#remove Dyspnea_PHYS_BL_CTCAE
EQVAS17 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS17[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Radiation_DoseReceived + WHO_PS + BaselineEQ5D + GTV_volume, 
    design = design_i, family = gaussian()
  )
}
EQVAS17_mira <- as.mira(EQVAS17)
D1(EQVAS16_mira, EQVAS17_mira) # Wald Test => remove Dyspnea_PHYS_BL_CTCAE

#remove GTV_volume
EQVAS18 <- vector("list", m) 
for (i in 1:m) {
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS18[[i]] <- svyglm(
    EQ5DVAS ~ Treatment30 + Stage + Pulmonary_comorbidity +
      Chemo_sequence + Radiation_DoseReceived + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
EQVAS18_mira <- as.mira(EQVAS18)
D1(EQVAS17_mira, EQVAS18_mira) # Wald Test => keep GTV_volume => reject EQVAS18_mira

# Final Model = EQVAS17
EQVAS17_pooled <- pool(EQVAS17_mira)
summary(EQVAS17_pooled)




## Global Health Status Final Model
glm_models_GHS <- vector("list", m)
for (i in 1:m) {
  # Define the survey design using TWANG weights
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models_GHS[[i]] <- svyglm(
    GlobalHealthStatus ~ Treatment30 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
GHS_mira <- as.mira(glm_models_GHS)
GHS_pooled <- pool(GHS_mira) 
summary(GHS_pooled)

### Scenario analysis - proron = <50% fractions as protons
# EQ5D5L
glm_models_SA_Tx50 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models_SA_Tx50[[i]] <- svyglm(
    EQindex ~ Treatment50 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
EQ5D_mira_SA_Tx50 <- as.mira(glm_models_SA_Tx50)
EQ5D_SA_Tx50_pooled <- pool(EQ5D_mira_SA_Tx50)
summary(EQ5D_SA_Tx50_pooled)
#EQVAS
EQVAS_SA_Tx50 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS_SA_Tx50[[i]] <- svyglm(
    EQ5DVAS ~ Treatment50 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
EQVAS_SA_Tx50 <- as.mira(EQVAS_SA_Tx50)
EQVAS_SA_Tx50_pooled <- pool(EQVAS_SA_Tx50)
summary(EQVAS_SA_Tx50_pooled)
#GHS
GHS_SA_Tx50 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  GHS_SA_Tx50[[i]] <- svyglm(
    GlobalHealthStatus ~ Treatment50 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
GHS_mira_SA_Tx50 <- as.mira(GHS_SA_Tx50)
GHS_SA_Tx50_pooled <- pool(GHS_mira_SA_Tx50)
summary(GHS_SA_Tx50_pooled)


### Scenario analysis - proron = <80% fractions as protons
# EQ5D5L
glm_models_SA_Tx80 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models_SA_Tx80[[i]] <- svyglm(
    EQindex ~ Treatment80 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
EQ5D_mira_SA_Tx80 <- as.mira(glm_models_SA_Tx80)
EQ5D_SA_Tx80_pooled <- pool(EQ5D_mira_SA_Tx80)
summary(EQ5D_SA_Tx80_pooled)
#EQVAS
EQVAS_SA_Tx80 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS_SA_Tx80[[i]] <- svyglm(
    EQ5DVAS ~ Treatment80 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
EQVAS_SA_Tx80 <- as.mira(EQVAS_SA_Tx80)
EQVAS_SA_Tx80_pooled <- pool(EQVAS_SA_Tx80)
summary(EQVAS_SA_Tx80_pooled)
#GHS
GHS_SA_Tx80 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  GHS_SA_Tx80[[i]] <- svyglm(
    GlobalHealthStatus ~ Treatment80 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
GHS_mira_SA_Tx80 <- as.mira(GHS_SA_Tx80)
GHS_SA_Tx80_pooled <- pool(GHS_mira_SA_Tx80)
summary(GHS_SA_Tx80_pooled)

### Scenario analysis - proron = <100% fractions as protons
# EQ5D5L
glm_models_SA_Tx100 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  glm_models_SA_Tx100[[i]] <- svyglm(
    EQindex ~ Treatment100 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
EQ5D_mira_SA_Tx100 <- as.mira(glm_models_SA_Tx100)
EQ5D_SA_Tx100_pooled <- pool(EQ5D_mira_SA_Tx100)
summary(EQ5D_SA_Tx100_pooled)
#EQVAS
EQVAS_SA_Tx100 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  EQVAS_SA_Tx100[[i]] <- svyglm(
    EQ5DVAS ~ Treatment100 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
EQVAS_SA_Tx100 <- as.mira(EQVAS_SA_Tx100)
EQVAS_SA_Tx100_pooled <- pool(EQVAS_SA_Tx100)
summary(EQVAS_SA_Tx100_pooled)
#GHS
GHS_SA_Tx100 <- vector("list", m)
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  design_i <- svydesign(ids = ~1, weights = ~twang_weight, data = data_i)
  GHS_SA_Tx100[[i]] <- svyglm(
    GlobalHealthStatus ~ Treatment100 + WHO_PS + BaselineEQ5D, 
    design = design_i, family = gaussian()
  )
}
GHS_mira_SA_Tx100 <- as.mira(GHS_SA_Tx100)
GHS_SA_Tx100_pooled <- pool(GHS_mira_SA_Tx100)
summary(GHS_SA_Tx100_pooled)

###################TESTING#######################

rm(list = ls())

# Load libraries
library(WeMix)
library(miceadds)  # For pooling
library(mitml)     # Alternative pooling option if needed

# Load imputed and weighted data
weighted_data_list <- readRDS("weighted_data_list.rds")

m <- length(weighted_data_list) # Number of imputations
mixed_models <- vector("list", m) # Store fitted models

# Loop over each imputed dataset and fit weighted mixed model
for (i in 1:m) {
  data_i <- weighted_data_list[[i]]
  data_i$obs_weight <- 1  # Create an observation-level weight as weighting done on patient-level and WeMix expects length "weights" of 2
  mixed_models[[i]] <- mix(
    formula = EQindex ~ Periode + Treatment30 + Stage + Pulmonary_comorbidity + Durvalumab +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume +
      Treatment30:Periode + (1 | PatientID),
    weights = c("twang_weight", "obs_weight"), 
    data = data_i
  )
}


# Extract coefficient estimates and variance-covariance matrices from each WeMix model
coef_list <- lapply(mixed_models, function(model) model$coef)
vcov_list <- lapply(mixed_models, function(model) model$cov_mat)

# Pool the estimates using miceadds::pool_mi()
WeMix_Pooled <- miceadds:::pool.mi.list(qhat = coef_list, uhat = vcov_list)
summary(WeMix_Pooled)

ls(getNamespace("miceadds"), pattern = "pool")




