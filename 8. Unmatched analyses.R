rm(list = ls())

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
library(car)

options(scipen=999)

# Load imputed datasets from 3. MICE imputation 
imputed_data <- readRDS("mice_imputed_data.rds")

###############Unmatched analyses###########################
# Univariate models for all considered fixed effects variables 

univmodel_EQ5Dindex_Label <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + (1 | PatientID)))
pool(univmodel_EQ5Dindex_Label)
summary(pool(univmodel_EQ5Dindex_Label)) # Based on indicated treatment (not actually received)

univmodel_EQ5Dindex_Tx30 <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Treatment30 + (1 | PatientID)))
pool(univmodel_EQ5Dindex_Tx30)
summary(pool(univmodel_EQ5Dindex_Tx30))

imputed_long <- complete(imputed_data, action = "long", include = TRUE)
imputed_long$TreatmentSplit <- relevel(factor(imputed_long$TreatmentSplit), ref = "Photons")
imputed_data <- as.mids(imputed_long)
levels( complete(imputed_data, 1)$TreatmentSplit )
univmodel_EQ5Dindex_TxComb <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ TreatmentSplit + (1 | PatientID)))
pool(univmodel_EQ5Dindex_TxComb)
summary(pool(univmodel_EQ5Dindex_TxComb))

univmodel_EQ5Dindex_Stage <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Stage + (1 | PatientID)))
pool(univmodel_EQ5Dindex_Stage)
summary(pool(univmodel_EQ5Dindex_Stage))

univmodel_EQ5Dindex_FracDailFreq <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Fraction_Dailyfrequency + (1 | PatientID)))
pool(univmodel_EQ5Dindex_FracDailFreq)
summary(pool(univmodel_EQ5Dindex_FracDailFreq))

univmodel_EQ5Dindex_PulmComorb <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Pulmonary_comorbidity + (1 | PatientID)))
pool(univmodel_EQ5Dindex_PulmComorb)
summary(pool(univmodel_EQ5Dindex_PulmComorb))

univmodel_EQ5Dindex_Durvalumab <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Durvalumab + (1 | PatientID)))
pool(univmodel_EQ5Dindex_Durvalumab)
summary(pool(univmodel_EQ5Dindex_Durvalumab))

univmodel_EQ5Dindex_ChemoSeq <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Chemo_sequence + (1 | PatientID)))
pool(univmodel_EQ5Dindex_ChemoSeq)
summary(pool(univmodel_EQ5Dindex_ChemoSeq))

univmodel_EQ5Dindex_SmokStat <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Smoking_status + (1 | PatientID)))
pool(univmodel_EQ5Dindex_SmokStat)
summary(pool(univmodel_EQ5Dindex_SmokStat))

univmodel_EQ5Dindex_Surgery <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Surgery + (1 | PatientID)))
pool(univmodel_EQ5Dindex_Surgery)
summary(pool(univmodel_EQ5Dindex_Surgery))

univmodel_EQ5Dindex_Adapted <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Adapted + (1 | PatientID)))
pool(univmodel_EQ5Dindex_Adapted)
summary(pool(univmodel_EQ5Dindex_Adapted))

univmodel_EQ5Dindex_TumLoc <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ TumorLocation + (1 | PatientID)))
pool(univmodel_EQ5Dindex_TumLoc)
summary(pool(univmodel_EQ5Dindex_TumLoc))

univmodel_EQ5Dindex_TumCellType <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ TumorCellType + (1 | PatientID)))
pool(univmodel_EQ5Dindex_TumCellType)
summary(pool(univmodel_EQ5Dindex_TumCellType))

univmodel_EQ5Dindex_DoseReceived <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Radiation_DoseReceived + (1 | PatientID)))
pool(univmodel_EQ5Dindex_DoseReceived)
summary(pool(univmodel_EQ5Dindex_DoseReceived))

univmodel_EQ5Dindex_Gender <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Gender + (1 | PatientID)))
pool(univmodel_EQ5Dindex_Gender)
summary(pool(univmodel_EQ5Dindex_Gender))

univmodel_EQ5Dindex_AgeStartRT <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Age_startRT + (1 | PatientID)))
pool(univmodel_EQ5Dindex_AgeStartRT)
summary(pool(univmodel_EQ5Dindex_AgeStartRT))

univmodel_EQ5Dindex_WHOPS <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ WHO_PS + (1 | PatientID)))
pool(univmodel_EQ5Dindex_WHOPS)
summary(pool(univmodel_EQ5Dindex_WHOPS))

univmodel_EQ5Dindex_BaselineUtility <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ BaselineEQ5D + (1 | PatientID)))
pool(univmodel_EQ5Dindex_BaselineUtility)
summary(pool(univmodel_EQ5Dindex_BaselineUtility))



### EQ5D Utility models (Dutch Tariff)

#Model 1 = all clinically plausible variables    

modelAll_EQ5Dindex <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Treatment30 + Stage + Pulmonary_comorbidity
                                                                 + Durvalumab + Chemo_sequence + Smoking_status + Surgery + Adapted
                                                                 + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                                 + WHO_PS + BaselineEQ5D + GTV_volume + (1 | PatientID)))
pool(modelAll_EQ5Dindex) 
summary(pool(modelAll_EQ5Dindex))


## Multicollinearity Assessment 
## Exploring multicollinearity

lm_EQ5Dindex1a <- with(data = imputed_data, lm(EQindex ~ Label + Stage + Fraction_Dailyfrequency + Pulmonary_comorbidity
                                               + Durvalumab + Chemo_sequence + Smoking_status + Surgery + Adapted
                                               + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                               + WHO_PS + BaselineEQ5D)) #fit linear model for multicollinearity check

lapply(lm_EQ5Dindex1a$analyses, function(model) {
  vif(model)}) # => Fraction_Dailyfrequency and TumorCellType considered high multicollinearity based on GVIF cut off >5 and GVIF^(1/(2*Df) cut off SQRT(10) 

# Explore fraction daily frequency and tumourcelltype
lm_EQ5Dindex1b <- with(data = imputed_data, lm(EQindex ~ Label + Stage + Pulmonary_comorbidity
                                               + Durvalumab + Chemo_sequence + Smoking_status + Surgery + Adapted
                                               + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                               + WHO_PS + BaselineEQ5D)) # Removed Fraction_Dailyfrequency 
lapply(lm_EQ5Dindex1b$analyses, function(model) {
  vif(model)}) # for linear model without Fraction_Dailyfrequency

lm_EQ5Dindex1c <- with(data = imputed_data, lm(EQindex ~ Label + Stage + Fraction_Dailyfrequency + Pulmonary_comorbidity
                                               + Durvalumab + Chemo_sequence + Smoking_status + Surgery + Adapted
                                               + TumorLocation + Radiation_DoseReceived + Gender + Age_startRT
                                               + WHO_PS + BaselineEQ5D)) # Removed tumorcelltype
lapply(lm_EQ5Dindex1c$analyses, function(model) {
  vif(model)}) # for linear model without tumorcelltype

imputed_dataset <- complete(imputed_data, action = 1) #change action number to assess different individual imputed datasets
imputed_dataset_unique <- imputed_dataset[!duplicated(imputed_dataset$PatientID), ]
table(imputed_dataset_unique$Fraction_Dailyfrequency, imputed_dataset_unique$TumorCellType) # => remove fraction daily frequency




##### Unmatched Final Model (Scenario Analysis based on base-case variable selection)

#EQ5D5L
EQindex_Final <- with(
  data = imputed_data,
  expr = lmer(
    EQindex ~ Treatment30
    + WHO_PS
    + BaselineEQ5D
    + Dyspnea_PHYS_BL_CTCAE
    + (1 | PatientID),
    REML = FALSE))
pool(EQindex_Final)
summary(pool(EQindex_Final), conf.int = TRUE)
#EQVAS
EQVAS_Final <- with(
  data = imputed_data,
  expr = lmer(
    EQ5DVAS ~ Treatment30
    + WHO_PS
    + BaselineEQVAS
    + Dyspnea_PHYS_BL_CTCAE
    + (1 | PatientID),
    REML = FALSE))
pool(EQVAS_Final)
summary(pool(EQVAS_Final), conf.int = TRUE)
#EORTC QLQ-C30 Global Health Status
GHS_Final <- with(
  data = imputed_data,
  expr = lmer(
    GlobalHealthStatus ~ Treatment30
    + WHO_PS
    + BaselineGHS
    + Dyspnea_PHYS_BL_CTCAE
    + (1 | PatientID),
    REML = FALSE))
pool(GHS_Final)
summary(pool(GHS_Final), conf.int = TRUE)





#### MODEL ASSESSMENT

plot(EQindex_Final$analyses[[1]]) # From one imputed dataset - Adjust the '1' to view alternatives
qqnorm(residuals(EQindex_Final$analyses[[1]])) # From one imputed dataset - Adjust the '1' to view alternatives
qqline(residuals(EQindex_Final$analyses[[1]])) # From one imputed dataset - Adjust the '1' to view alternatives
hist(residuals(EQindex_Final$analyses[[1]]), breaks = 30, main = "Histogram of Residuals", xlab = "Residuals") # From one imputed dataset - Adjust the '1' to view alternatives

model_residuals <- residuals(EQindex_Final$analyses[[1]])  # From one imputed dataset - Adjust the '1' to view alternatives
hist(model_residuals, main="Residuals Histogram", xlab="Residuals")
plot(fitted(EQindex_Final$analyses[[1]]), model_residuals, main="Fitted vs Residuals",
     xlab="Fitted Values", ylab="Residuals")

### EQ VAS model (based on EQ5D selections)
EQVAS_Final <- with(data = imputed_data, 
                    exp = lme4::lmer(EQindex ~ Label + WHO_PS + BaselineEQVAS + (1 | PatientID))) 
pool(EQVAS_Final)
summary(pool(EQVAS_Final), conf.int = TRUE)


### EORTC QLQ-C30 GHS model (based on EQ5D selections)
GHS_Final <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + WHO_PS + BaselineGHS + (1 | PatientID))) 
pool(GHS_Final)
summary(pool(GHS_Final), conf.int = TRUE)




