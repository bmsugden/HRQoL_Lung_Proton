rm(list = ls())

# Load imputed datasets from 3. MICE imputation 
imputed_data <- readRDS("mice_imputed_data.rds")


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


###############Unmatched analyses###########################
# Univariate models for all considered fixed effects variables 

univmodel_EQ5Dindex_Label <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + (1 | PatientID)))
pool(univmodel_EQ5Dindex_Label)
summary(pool(univmodel_EQ5Dindex_Label))

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

modelAll_EQ5Dindex <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Stage + Fraction_Dailyfrequency + Pulmonary_comorbidity
                                                                 + Durvalumab + Chemo_sequence + Smoking_status + Surgery + Adapted
                                                                 + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                                 + WHO_PS + BaselineEQ5D + (1 | PatientID)))
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

#################################### Backward elimiation
################################### Note, anova() not compatible with MICE imputation so p-values cannot be directly derived for each FE variable
################################## Hence, Wald Test used in line with https://stefvanbuuren.name/fimd/sec-stepwise.html (5.4.2)

#### Backwards elimination: Start 
## 0.05 significance level; Selection approach: Stepwise selection using The Wald Test; force in Label; 
## Starting model based on "multicollinearity assessment" [see above]
## Variables excluded (reason): Fraction_Dailyfrequency (collinearity)
EQ5Dindex_Start <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Stage + Pulmonary_comorbidity
                                                              + Durvalumab + Chemo_sequence + Smoking_status + Surgery + Adapted
                                                              + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                              + WHO_PS + BaselineEQ5D + (1 | PatientID)))

pool(EQ5Dindex_Start)
summary(pool(EQ5Dindex_Start)) 

EQ5Dindex_nostage <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Pulmonary_comorbidity + Durvalumab
                                                                + Chemo_sequence + Smoking_status + Surgery + Adapted
                                                                + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                                + WHO_PS + BaselineEQ5D + (1 | PatientID))) # no Stage
pool(EQ5Dindex_nostage)
summary(pool(EQ5Dindex_nostage)) 

D1(EQ5Dindex_Start, EQ5Dindex_nostage) # Wald Test for StepWise selection => remove stage

EQ5Dindex_noPulmCo <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Durvalumab + Chemo_sequence + Smoking_status + Surgery + Adapted
                                                                 + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                                 + WHO_PS + BaselineEQ5D + (1 | PatientID))) # no Stage/Pulmonary_comorbidity

pool(EQ5Dindex_noPulmCo)
summary(pool(EQ5Dindex_noPulmCo)) 

D1(EQ5Dindex_nostage, EQ5Dindex_noPulmCo) # Wald Test for StepWise selection => remove pulmonary comorbidity

EQ5Dindex_noDurv <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Chemo_sequence + Smoking_status + Surgery + Adapted
                                                               + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                               + WHO_PS + BaselineEQ5D + (1 | PatientID))) # no Stage/Pulmonary_comorbidity/Durvalumab

pool(EQ5Dindex_noDurv)
summary(pool(EQ5Dindex_noDurv)) 

D1(EQ5Dindex_noPulmCo, EQ5Dindex_noDurv) # Wald Test for StepWise selection => remove Durvalumab

EQ5Dindex_ChemSeq <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Smoking_status + Surgery + Adapted
                                                                + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                                + WHO_PS + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence

pool(EQ5Dindex_ChemSeq)
summary(pool(EQ5Dindex_ChemSeq)) 

D1(EQ5Dindex_noDurv, EQ5Dindex_ChemSeq) # Wald Test for StepWise selection => remove Chemo_sequence

EQ5Dindex_SmokStat <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Surgery + Adapted
                                                                 + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                                 + WHO_PS + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status

pool(EQ5Dindex_SmokStat)
summary(pool(EQ5Dindex_SmokStat)) 

D1(EQ5Dindex_ChemSeq, EQ5Dindex_SmokStat) # Wald Test for StepWise selection => remove Smoking_status


EQ5Dindex_Surgery <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Adapted
                                                                + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                                + WHO_PS + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status/Surgery

pool(EQ5Dindex_Surgery)
summary(pool(EQ5Dindex_Surgery)) 

D1(EQ5Dindex_SmokStat, EQ5Dindex_Surgery) # Wald Test for StepWise selection => remove Surgery

EQ5Dindex_Adapted <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + TumorLocation + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                                + WHO_PS + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status/Surgery/Adapted

pool(EQ5Dindex_Adapted)
summary(pool(EQ5Dindex_Adapted)) 

D1(EQ5Dindex_Surgery, EQ5Dindex_Adapted) # Wald Test for StepWise selection => remove Adapted

EQ5Dindex_TumLoc <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + TumorCellType + Radiation_DoseReceived + Gender + Age_startRT
                                                               + WHO_PS + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status/Surgery/TumorLocation

pool(EQ5Dindex_TumLoc)
summary(pool(EQ5Dindex_TumLoc)) 

D1(EQ5Dindex_Adapted, EQ5Dindex_TumLoc) # Wald Test for StepWise selection => remove TumorLocation

EQ5Dindex_TumCellTyp <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Radiation_DoseReceived + Gender + Age_startRT
                                                                   + WHO_PS + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status/Surgery/TumorLocation/TumorCellType

pool(EQ5Dindex_TumCellTyp)
summary(pool(EQ5Dindex_TumCellTyp)) 

D1(EQ5Dindex_TumLoc, EQ5Dindex_TumCellTyp) # Wald Test for StepWise selection => remove TumorCellType

EQ5Dindex_DoseRec <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Gender + Age_startRT
                                                                + WHO_PS + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status/Surgery/TumorLocation/TumorCellType/Radiation_DoseReceived

pool(EQ5Dindex_DoseRec)
summary(pool(EQ5Dindex_DoseRec)) 

D1(EQ5Dindex_TumCellTyp, EQ5Dindex_DoseRec) # Wald Test for StepWise selection => remove Radiation_DoseReceived

EQ5Dindex_Gender <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + Age_startRT
                                                               + WHO_PS + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status/Surgery/TumorLocation/TumorCellType/Radiation_DoseReceived
# .../Gender

pool(EQ5Dindex_Gender)
summary(pool(EQ5Dindex_Gender)) 

D1(EQ5Dindex_DoseRec, EQ5Dindex_Gender) # Wald Test for StepWise selection => remove Gender

EQ5Dindex_Age <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + WHO_PS + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status/Surgery/TumorLocation/TumorCellType/Radiation_DoseReceived
# .../Gender/Age_startRT

pool(EQ5Dindex_Age)
summary(pool(EQ5Dindex_Age)) 

D1(EQ5Dindex_Gender, EQ5Dindex_Age) # Wald Test for StepWise selection => remove Age_startRT

EQ5Dindex_WHOPS <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + BaselineEQ5D + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status/Surgery/TumorLocation/TumorCellType/Radiation_DoseReceived
# .../Gender/Age_startRT/WHO_PS

pool(EQ5Dindex_WHOPS)
summary(pool(EQ5Dindex_WHOPS)) 

D1(EQ5Dindex_Age, EQ5Dindex_WHOPS) # Wald Test for StepWise selection => KEEP WHO_PS

EQ5Dindex_BlineUtil <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + WHO_PS + (1 | PatientID))) 
# no Stage/Pulmonary_comorbidity/Durvalumab/Chemo_sequence/Smoking_status/Surgery/TumorLocation/TumorCellType/Radiation_DoseReceived
# .../Gender/Age_startRT/BaselineEQ5D

pool(EQ5Dindex_BlineUtil)
summary(pool(EQ5Dindex_BlineUtil)) 

D1(EQ5Dindex_Age, EQ5Dindex_BlineUtil) # Wald Test for StepWise selection => KEEP baselineUtil


#### Final Model
EQ5Dindex_FINAL <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + WHO_PS + BaselineEQ5D + (1 | PatientID))) 

pool(EQ5Dindex_FINAL)
summary(pool(EQ5Dindex_FINAL)) 


## CHECKS - add one variable at a time to ensure order did not matter
EQ5Dindex_TEST <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + WHO_PS + BaselineEQ5D 
                                                             + Durvalumab # Adjust this FE variable in addition to final model to assess
                                                             + (1 | PatientID))) 
pool(EQ5Dindex_TEST)
summary(pool(EQ5Dindex_TEST)) 
D1(EQ5Dindex_TEST, EQ5Dindex_FINAL) 

#### Backwards elimination: End

#### FINAL MODEL ASSESSMENT

plot(EQ5Dindex_FINAL$analyses[[1]]) # From one imputed dataset - Adjust the '1' to view alternatives
qqnorm(residuals(EQ5Dindex_FINAL$analyses[[1]])) # From one imputed dataset - Adjust the '1' to view alternatives
qqline(residuals(EQ5Dindex_FINAL$analyses[[1]])) # From one imputed dataset - Adjust the '1' to view alternatives
hist(residuals(EQ5Dindex_FINAL$analyses[[1]]), breaks = 30, main = "Histogram of Residuals", xlab = "Residuals") # From one imputed dataset - Adjust the '1' to view alternatives

model_residuals <- residuals(EQ5Dindex_FINAL$analyses[[1]])  # From one imputed dataset - Adjust the '1' to view alternatives
hist(model_residuals, main="Residuals Histogram", xlab="Residuals")
plot(fitted(EQ5Dindex_FINAL$analyses[[1]]), model_residuals, main="Fitted vs Residuals",
     xlab="Fitted Values", ylab="Residuals")





### EQ VAS model (based on EQ5D selections)
EQ5DVAS <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + WHO_PS + BaselineEQVAS + (1 | PatientID))) 
pool(EQ5DVAS)
summary(pool(EQ5DVAS))

### GHS model (based on EQ5D selections)
GlobalHealthStatus <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Label + WHO_PS + BaselineGHS + (1 | PatientID))) 
pool(GlobalHealthStatus)
summary(pool(GlobalHealthStatus))




## 30% cut of to determine proton eligibility ####### TO BE STRUCTURED BETTER

EQ5D30_Final <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Treatment30 + WHO_PS + BaselineEQ5D + (1 | PatientID)))
pool(EQ5D30_Final)
summary(pool(EQ5D30_Final)) 

EQVAS30_Final <- with(data = imputed_data, exp = lme4::lmer(EQ5DVAS ~ Treatment30 + WHO_PS + BaselineEQVAS + (1 | PatientID)))
pool(EQVAS30_Final)
summary(pool(EQVAS30_Final)) 

EQGHS30_Final <- with(data = imputed_data, exp = lme4::lmer(GlobalHealthStatus ~ Treatment30 + WHO_PS + BaselineGHS + (1 | PatientID)))
pool(EQGHS30_Final)
summary(pool(EQGHS30_Final)) 

#Scenario with 50% cut off
EQ5D40_Final <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Treatment50 + WHO_PS + BaselineEQ5D + (1 | PatientID)))
pool(EQ5D40_Final)
summary(pool(EQ5D40_Final)) 





#with durva forced in (and label forced in - protons determined with 30% cut-off)
EQ5D30_Durv <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Treatment30 + Durvalumab + WHO_PS + BaselineEQ5D + (1 | PatientID)))
pool(EQ5D30_Durv)
summary(pool(EQ5D30_Durv))

#with durva forced in (label not forced in)
EQ5D30_Durv_notreat <- with(data = imputed_data, exp = lme4::lmer(EQindex ~ Durvalumab + WHO_PS + BaselineEQ5D + (1 | PatientID)))
pool(EQ5D30_Durv_notreat)
summary(pool(EQ5D30_Durv_notreat))
