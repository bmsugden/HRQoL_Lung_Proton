rm(list = ls())

# Load imputed data
imputed_data <- readRDS("mice_imputed_data.rds")


# Load required libraries
library(Matching)
library(rgenoud) # used by Matching package for optimisation of matching
library(mice)
library(ggplot2) 
library(dplyr)# For density plots
library(summarytools)
library(car)



# Step 1: Check baseline overlap using ggplot density plots. Matching not recommended if poor overlap. 
## Note// ggplot requires a dataframe (rather than a mids object produced by data imputation). 
## Hence, first create combined dataset with imputation modifier. Densities are plotted for each imputation using facets.

imputed_list <- lapply(1:imputed_data$m, function(i) {
  dat <- complete(imputed_data, i)
  dat$Imputation <- i
  return(dat)
})
combined_data <- bind_rows(imputed_list)

#Summary of baseline imputed data (one row per patient per imputation)
pooled_data <- combined_data %>%
  group_by(PatientID, Imputation) %>%
  summarise(
    # For numeric variables: take the mean across imputations.
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    
    # For character/factor variables: choose the mode (most frequent value)
    across(where(~ is.character(.) || is.factor(.)), ~ {
      tab <- table(.x)
      names(tab)[which.max(tab)]
    }),
    .groups = "drop"
  )

final_pooled_data <- pooled_data %>% #pool across imputations
  group_by(PatientID) %>%
  summarise(
    # Average numeric variables across imputations.
    across(where(is.numeric), ~ mean(.x, na.rm = TRUE)),
    # For character/factor variables: choose the mode across imputations.
    across(where(~ is.character(.) || is.factor(.)), ~ {
      tab <- table(.x)
      names(tab)[which.max(tab)]
    }),
    .groups = "drop"
  )

pooled_data_protons <- final_pooled_data %>%
  filter(grepl("Protons", Treatment30)) #pooled data for proton patients
pooled_data_photons <- final_pooled_data %>%
  filter(grepl("Photons", Treatment30)) # pooled data for photon patients

view(dfSummary(pooled_data_protons)) # pooled proton summary
view(dfSummary(pooled_data_photons)) # pooled photon summary



##Age
ggplot(combined_data, aes(x = Age_startRT, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of Age by Treatment Group by Imputation",
       fill = "Treatment Group")
##Stage
ggplot(combined_data, aes(x = Stage, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of Stage by Treatment Group by Imputation",
       fill = "Treatment Group")
##Gender
ggplot(combined_data, aes(x = Gender, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of Gender by Treatment Group by Imputation",
       fill = "Treatment Group")
##WHO_PS - Grouped variable fo matching
ggplot(combined_data, aes(x = WHO_PS_matching, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of WHO PS by Treatment Group by Imputation",
       fill = "Treatment Group")
##Pulmonary_comorbidity
ggplot(combined_data, aes(x = Pulmonary_comorbidity, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of Pulmonary_comorbidity by Treatment Group by Imputation",
       fill = "Treatment Group")
##Chemo_sequence
ggplot(combined_data, aes(x = Chemo_sequence, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of Chemo_sequence by Treatment Group by Imputation",
       fill = "Treatment Group")
##Smoking_status
ggplot(combined_data, aes(x = Smoking_status, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of Smoking_status by Treatment Group by Imputation",
       fill = "Treatment Group")
##Surgery
ggplot(combined_data, aes(x = Surgery, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of Surgery by Treatment Group by Imputation",
       fill = "Treatment Group")
##TumorLocation
ggplot(combined_data, aes(x = TumorLocation, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of TumorLocation by Treatment Group by Imputation",
       fill = "Treatment Group")
##TumorCellType
ggplot(combined_data, aes(x = TumorCellType, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of TumorCellType by Treatment Group by Imputation",
       fill = "Treatment Group")
##GTV_volume
ggplot(combined_data, aes(x = GTV_volume, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of GTV_volume by Treatment Group by Imputation",
       fill = "Treatment Group")
##Baseline Dysphagia
ggplot(combined_data, aes(x = Dysphagia_BL_PHYS_CTCAE, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of BaselineDysphagia by Treatment Group by Imputation",
       fill = "Treatment Group")
##Baseline Dyspnoea
ggplot(combined_data, aes(x = Dyspnea_PHYS_BL_CTCAE, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  facet_wrap(~ Imputation) +
  labs(title = "Density Plot of BaselineDypnoea by Treatment Group by Imputation",
       fill = "Treatment Group")

##To view for specific dataset - adjust "1" to specific dataset number, adjust "Age_startRT" to variable
ggplot(complete(imputed_data, 1), aes(x = Age_startRT, fill = factor(Treatment30))) +
  geom_density(alpha = 0.5) +
  labs(title = "Density Plot of by specified imputed dataset",
       fill = "Treatment Group")

##MATCHING
m <- imputed_data$m
matched_data_list <- vector("list", m)
model_list <- vector("list", m)  # for storing models later

for(i in 1:m) {
  set.seed(12345 + i)  # Different seed per imputation, but reproducible
  # 1. Extract the full imputed dataset for the i-th imputation (includes all time points)
  full_data_i <- complete(imputed_data, i)
  
  # 2. Collapse the data to one row per PatientID for matching only
  collapsed_data <- full_data_i[!duplicated(full_data_i$PatientID), ]
  
  # Combine Stage III and IV into one category for matching in the collapsed data
  collapsed_data$Stage_combined <- factor(
    ifelse(collapsed_data$Stage %in% c("III", "IV"), "III_IV", as.character(collapsed_data$Stage))
  )
  
  # Define treatment indicator in the collapsed data
  collapsed_data$treatment_indicator <- collapsed_data$Treatment30 == "Protons"
  
  # Define covariates for matching (time invariant)
  covariates <- c("Age_startRT", "Stage_combined", "WHO_PS_matching",
                  "Dyspnea_PHYS_BL_CTCAE", "BaselineEQ5D")
  
  # Propensity score model on collapsed data
  ps_formula <- as.formula(paste("treatment_indicator ~", paste(covariates, collapse = " + ")))
  ps_model <- glm(ps_formula, data = collapsed_data, family = binomial())
  collapsed_data$propensity_score <- predict(ps_model, type = "response")
  
  # Create covariate matrix (including the propensity score) for matching
  X <- model.matrix(~ . - 1, data = collapsed_data[, c(covariates, "propensity_score")])
  
  # Use GenMatch for optimal weights on the collapsed data
  genout <- GenMatch(Tr = collapsed_data$treatment_indicator,
                     X = X,
                     M = 1, 
                     replace = TRUE,
                     pop.size = 300,
                     max.generations = 100,
                     wait.generations = 5)
  
  # Perform matching using the obtained weight matrix on the collapsed data
  match_out <- Match(Tr = collapsed_data$treatment_indicator, 
                     X = X,
                     Weight.matrix = genout)
  
  # Extract the matched indices (from the collapsed data)
  treated_idx <- match_out$index.treated
  control_idx <- match_out$index.control
  matched_indices <- c(treated_idx, control_idx)
  
  # Create a data frame with matched PatientIDs and matching weights
  matched_info <- data.frame(
    PatientID = collapsed_data$PatientID[matched_indices],
    match_weight = match_out$weights[matched_indices]
  )
  
  matched_info <- matched_info[!duplicated(matched_info$PatientID), ] 
  # remove duplicate IDs from matched_info - to ensure extra copies aren't created for merging back with longitudinal data
  
  # 3. Merge the matching information with the full longitudinal dataset
  # This attaches the matching info (and weights) to every time-point record for the matched patients
  full_matched_data <- merge(full_data_i, matched_info, by = "PatientID", all.x = FALSE, all.y = TRUE)
  # Add an imputation indicator
  full_matched_data$Imputation <- i
  # Save the fully re-expanded matched data for this imputation
  matched_data_list[[i]] <- full_matched_data
}
  

### CHECKS
table(duplicated(matched_info$PatientID)) #Total number of uniquepatients matched across inputed datasets
MergedData <- read.csv("MergedData2.csv", stringsAsFactors = TRUE, check.names = FALSE, na.strings = c("", " ","NA", "<NA>", "#NUM!")) #load original data for comparison
sapply(matched_data_list, function(x) length(unique(x$PatientID))) #numbr of matched patients in each inputed dataset
sapply(matched_data_list, nrow) #total number of observations post-matching in each inpiuted dataset
unmatched_counts <- sapply(matched_data_list, function(x) {
  length(setdiff(unique(MergedData$PatientID), unique(x$PatientID)))
})
unmatched_counts # number of patients not matched
matched_counts <- sapply(matched_data_list, function(x) length(unique(x$PatientID)))
total_counts <- matched_counts + unmatched_counts 
total_counts # ensure unmatched + matched equal 242 (i.e., total number of patients)


# Assessing covariate balance
## Covariate balance of continuous and ordinal values with one or more category using Kolmorogov-Smirnov bootstrap p-values
## Covariate balance of variables with dichotomous values - comparison of proportion of patients with each value in each group

## Assessing matching balance using MatchBalance function

balance_results <- MatchBalance(treatment_indicator ~ Age_startRT + Stage_combined + WHO_PS_matching + 
                                Dyspnea_PHYS_BL_CTCAE + BaselineEQ5D,
                                data = collapsed_data,
                                match.out = match_out,
                                nboots = 500)

# Assess individual (matched) inputed datasets
dfSummary(matched_data_list[[3]])  # Number in square brackets corresponds to the inputed dataset




##################################


saveRDS(matched_data_list, file = "matched_data_list.rds")


