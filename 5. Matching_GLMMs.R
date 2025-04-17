rm(list = ls())

# Load imputed datasets from 3. MICE imputation 
matched_data_list <- readRDS("matched_data_list.rds")

library(lme4)
library(lmerTest)
library(summarytools)
library(ggplot2)
library(car)
library(summarytools)
library(PerformanceAnalytics)
library(broom.mixed)
library(mice)
library(miceadds)
library(mitml)
library(dplyr)

options(scipen=999)


m <- length(matched_data_list)

# Create an empty list to store models
glm_models_matched <- vector("list", m)
# Loop through each imputed dataset and fit the GLM
for (i in 1:m) {
  # Extract the i-th imputed matched dataset
  data_i <- matched_data_list[[i]]
  
  # Fit the GLM model on matched data - all candidate covariates
  glm_models_matched[[i]] <- glm(
    EQindex ~ Treatment30 + Stage + Pulmonary_comorbidity + Durvalumab +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D1_matched_mira <- as.mira(glm_models_matched)
EQ5D1_matched_pooled <- pool(EQ5D1_matched_mira)
summary(EQ5D1_matched_pooled)


# Fit matched model without Stage
glm_models_matched2 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched2[[i]] <- glm(
    EQindex ~ Treatment30 + Pulmonary_comorbidity + Durvalumab +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D2_matched_mira <- as.mira(glm_models_matched2) #without stage

#Wald Test (0.05 sig lev)
D1(EQ5D1_matched_mira, EQ5D2_matched_mira) # => remove stage

# Fit matched model without Pulmonary_comorbidity
glm_models_matched3 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched3[[i]] <- glm(
    EQindex ~ Treatment30 + Durvalumab +
      Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D3_matched_mira <- as.mira(glm_models_matched3) #without Pulmonary_comorbidity

#Wald Test (0.05 sig lev)
D1(EQ5D2_matched_mira, EQ5D3_matched_mira) # => remove Pulmonary_comorbidity

# Fit matched model without Durvalumab
glm_models_matched4 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched4[[i]] <- glm(
    EQindex ~ Treatment30 + Chemo_sequence + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D4_matched_mira <- as.mira(glm_models_matched4) #without Durvalumab

#Wald Test (0.05 sig lev)
D1(EQ5D3_matched_mira, EQ5D4_matched_mira) # => remove Durvalumab

# Fit matched model without Chemo_sequence
glm_models_matched5 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched5[[i]] <- glm(
    EQindex ~ Treatment30 + Smoking_status + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D5_matched_mira <- as.mira(glm_models_matched5) #without Chemo_sequence

#Wald Test (0.05 sig lev)
D1(EQ5D4_matched_mira, EQ5D5_matched_mira) # => remove Chemo_sequence

# Fit matched model without Smoking_status
glm_models_matched6 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched6[[i]] <- glm(
    EQindex ~ Treatment30 + Surgery + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D6_matched_mira <- as.mira(glm_models_matched6) #without Smoking_status

#Wald Test (0.05 sig lev)
D1(EQ5D5_matched_mira, EQ5D6_matched_mira) # => remove Smoking_status

# Fit matched model without Surgery
glm_models_matched7 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched7[[i]] <- glm(
    EQindex ~ Treatment30 + Adapted + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D7_matched_mira <- as.mira(glm_models_matched7) #without Surgery

#Wald Test (0.05 sig lev)
D1(EQ5D6_matched_mira, EQ5D7_matched_mira) # => remove Surgery

# Fit matched model without Adapted
glm_models_matched8 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched8[[i]] <- glm(
    EQindex ~ Treatment30 + TumorLocation +
      TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D8_matched_mira <- as.mira(glm_models_matched8) #without Adapted

#Wald Test (0.05 sig lev)
D1(EQ5D7_matched_mira, EQ5D8_matched_mira) # => remove Adapted

# Fit matched model without TumorLocation
glm_models_matched9 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched9[[i]] <- glm(
    EQindex ~ Treatment30 + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D9_matched_mira <- as.mira(glm_models_matched9) #without TumorLocation

#Wald Test (0.05 sig lev)
D1(EQ5D8_matched_mira, EQ5D9_matched_mira) # => remove TumorLocation

# Fit matched model without TumorCellType
glm_models_matched10 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched10[[i]] <- glm(
    EQindex ~ Treatment30 + Radiation_DoseReceived + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D10_matched_mira <- as.mira(glm_models_matched10) #without TumorCellType

#Wald Test (0.05 sig lev)
D1(EQ5D9_matched_mira, EQ5D10_matched_mira) # => remove TumorCellType

# Fit matched model without Radiation_DoseReceived
glm_models_matched11 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched11[[i]] <- glm(
    EQindex ~ Treatment30 + Gender + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D11_matched_mira <- as.mira(glm_models_matched11) #without Radiation_DoseReceived

#Wald Test (0.05 sig lev)
D1(EQ5D10_matched_mira, EQ5D11_matched_mira) # => remove Radiation_DoseReceived

# Fit matched model without Gender
glm_models_matched12 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched12[[i]] <- glm(
    EQindex ~ Treatment30 + Age_startRT + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D12_matched_mira <- as.mira(glm_models_matched12) #without Gender

#Wald Test (0.05 sig lev)
D1(EQ5D11_matched_mira, EQ5D12_matched_mira) # => remove Gender

# Fit matched model without Age_startRT
glm_models_matched13 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched13[[i]] <- glm(
    EQindex ~ Treatment30 + WHO_PS + 
      BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D13_matched_mira <- as.mira(glm_models_matched13) #without Age_startRT

#Wald Test (0.05 sig lev)
D1(EQ5D12_matched_mira, EQ5D13_matched_mira) # => remove Age_startRT

# Fit matched model without WHO_PS
glm_models_matched14 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched14[[i]] <- glm(
    EQindex ~ Treatment30 + BaselineEQ5D + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D14_matched_mira <- as.mira(glm_models_matched14) #without WHO_PS

#Wald Test (0.05 sig lev)
D1(EQ5D13_matched_mira, EQ5D14_matched_mira) # => keep WHO_PS - reject EQ5D14_matched_mira

# Fit matched model without BaselineEQ5D
glm_models_matched15 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched15[[i]] <- glm(
    EQindex ~ Treatment30 + WHO_PS + Dysphagia_BL_PHYS_CTCAE + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D15_matched_mira <- as.mira(glm_models_matched15) #without BaselineEQ5D

#Wald Test (0.05 sig lev)
D1(EQ5D13_matched_mira, EQ5D15_matched_mira) # => keep BaselineEQ5D - reject EQ5D15_matched_mira

# Fit matched model without Dysphagia_BL_PHYS_CTCAE
glm_models_matched16 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched16[[i]] <- glm(
    EQindex ~ Treatment30 + WHO_PS + BaselineEQ5D + Dyspnea_PHYS_BL_CTCAE + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D16_matched_mira <- as.mira(glm_models_matched16) #without Dysphagia_BL_PHYS_CTCAE

#Wald Test (0.05 sig lev)
D1(EQ5D13_matched_mira, EQ5D16_matched_mira) # => remove Dysphagia_BL_PHYS_CTCAE 

# Fit matched model without Dyspnoea_BL_PHYS_CTCAE
glm_models_matched17 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched17[[i]] <- glm(
    EQindex ~ Treatment30 + WHO_PS + BaselineEQ5D + GTV_volume, 
    data = data_i, family = gaussian()
  )
}

EQ5D17_matched_mira <- as.mira(glm_models_matched17) #Keep Dyspnoea_BL_PHYS_CTCAE - reject EQ5D17_matched_mira

#Wald Test (0.05 sig lev)
D1(EQ5D16_matched_mira, EQ5D17_matched_mira) # => remove Dyspnoea_BL_PHYS_CTCAE 


# Fit matched model without GTV_volume
glm_models_matched18 <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matched18[[i]] <- glm(
    EQindex ~ Treatment30 + WHO_PS + BaselineEQ5D + Dyspnea_PHYS_BL_CTCAE, 
    data = data_i, family = gaussian()
  )
}

EQ5D18_matched_mira <- as.mira(glm_models_matched18) #without GTV_volume

#Wald Test (0.05 sig lev)
D1(EQ5D16_matched_mira, EQ5D18_matched_mira) # => remove GTV_volume 

## EQ5D18 = FINAL MODEL
##### END SELECTION

# Final Model
EQ5D18_pooled <- pool(EQ5D18_matched_mira)
summary(EQ5D18_pooled)

###### TESTING the final model - given the order of variable removal was arbitrary

# Fit final model + a removed variable
glm_models_matchedTEST <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_matchedTEST[[i]] <- glm(
    EQindex ~ Treatment30 + WHO_PS + BaselineEQ5D + Dyspnea_PHYS_BL_CTCAE
    + Stage, # replace this variable and rerun to test order of removal
    data = data_i, family = gaussian()
  )
}

EQ5DTEST_matched_mira <- as.mira(glm_models_matchedTEST) 

#Wald Test (0.05 sig lev)
D1(glm_models_matchedTEST, EQ5D18_matched_mira) # p-value >0.05 suggests additional variable does not significantly improve model fit

## END TESTING




## EQVAS Final Model
glm_models_EQVAS <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_EQVAS[[i]] <- glm(
    EQ5DVAS ~ Treatment30 + WHO_PS + BaselineEQ5D + Dyspnea_PHYS_BL_CTCAE, 
    data = data_i, family = gaussian()
  )
}

EQVAS_mira <- as.mira(glm_models_EQVAS)
EQVAS_pooled <- pool(EQVAS_mira)
summary(EQVAS_pooled)

## Global Health Status Final Model
glm_models_GHS <- vector("list", m)
for (i in 1:m) {
  data_i <- matched_data_list[[i]]
  glm_models_GHS[[i]] <- glm(
    GlobalHealthStatus ~ Treatment30 + WHO_PS + BaselineEQ5D + Dyspnea_PHYS_BL_CTCAE, 
    data = data_i, family = gaussian()
  )
}

GHS_mira <- as.mira(glm_models_GHS)
GHS_pooled <- pool(GHS_mira)
summary(GHS_pooled)


