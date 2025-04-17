rm(list = ls())

# Load data from CSV file (2. Summaries and Utility calculations)
MergedData <- read.csv("MergedData2.csv", stringsAsFactors = TRUE, check.names = FALSE, na.strings = c("", " ","NA", "<NA>", "#NUM!"))

# library(lme4)
# library(Matrix)
# library(lmerTest)
library(summarytools)
library(ggplot2)
# library(MuMIn)
library(car)
library(summarytools)
library(mice)
library(miceadds)
library(PerformanceAnalytics)
library(dplyr)
library(ggpubr)


str(MergedData)

MergedData$WHO_PS <- as.factor(MergedData$WHO_PS) 
MergedData$Dysphagia_BL_PHYS_CTCAE <- as.factor(MergedData$Dysphagia_BL_PHYS_CTCAE) 
MergedData$Dyspnea_PHYS_BL_CTCAE <- as.factor(MergedData$Dyspnea_PHYS_BL_CTCAE) 
MergedData$Periode <- as.numeric(MergedData$Periode) 
MergedData$Days_StartStop <- as.numeric(MergedData$Days_StartStop) 
MergedData$GTV_volume <- as.numeric(MergedData$GTV_volume) 
dfSummary(MergedData)


sum(rowSums(is.na(MergedData)) > 0) #Number of patients with at least one missing
#PRIOR TO ADDING ADDITIONAL NTCP VARIABLES# 418/763 # => 55 to be used as number of imputed datasets, m, in missing data imputation (i.e., proportion of patients with at least one missing variable)
629/763 # => 82 to be used as number of imputed datasets, m, in missing data imputation (i.e., proportion of observations with at least one missing variable)
str(MergedData)

#Set reference values 
MergedData$Label <- relevel(MergedData$Label, ref = "Photons (PV)")
MergedData$Stage <- relevel(MergedData$Stage, ref = "I/II")
MergedData$Fraction_Dailyfrequency <- relevel(MergedData$Fraction_Dailyfrequency, ref = "Once Daily")
MergedData$Pulmonary_comorbidity <- relevel(MergedData$Pulmonary_comorbidity, ref = "no")
MergedData$Durvalumab <- relevel(MergedData$Durvalumab, ref = "no")
MergedData$Chemo_sequence <- relevel(MergedData$Chemo_sequence, ref = "none")
MergedData$Smoking_status <- relevel(MergedData$Smoking_status, ref = "Never")
MergedData$Surgery <- relevel(MergedData$Surgery, ref = "No")
MergedData$Adapted <- relevel(MergedData$Adapted, ref = "Not adapted")
MergedData$TumorLocation <- relevel(MergedData$TumorLocation, ref = "Other tumours") 
MergedData$TumorCellType <- relevel(MergedData$TumorCellType, ref = "Small Cell") 
MergedData$WHO_PS <- relevel(MergedData$WHO_PS, ref = "0")
MergedData$WHO_PS_matching <- relevel(MergedData$WHO_PS_matching, ref = "0/1") 
MergedData$Dysphagia_BL_PHYS_CTCAE <- relevel(MergedData$Dysphagia_BL_PHYS_CTCAE, ref = "No") 
MergedData$Dyspnea_PHYS_BL_CTCAE <- relevel(MergedData$Dyspnea_PHYS_BL_CTCAE, ref = "No") 


#### DATA IMPUTATION
md.pattern(MergedData, plot = FALSE) #either 0 (missing) or 1 (observed). 
# The first column provides the frequency of each pattern. The last column lists the number of missing entries per pattern. 
# The bottom row provides the number of missing entries per variable, and the total number of missing cells.

#Define imputation methods (2lonly to ensure consistency within-patient)
mice_methods <- make.method(MergedData)
mice_methods["EQindex"] <- "pmm"
mice_methods["EQ5DVAS"] <- "pmm"
mice_methods["GlobalHealthStatus"] <- "pmm"
mice_methods["Radiation_DoseReceived"] <- "2lonly.pmm"
mice_methods["Durvalumab"] <- "2lonly.pmm" 
mice_methods["TumorCellType"] <- "2lonly.pmm"
mice_methods["TumorLocation"] <- "2lonly.pmm"
mice_methods["Surgery"] <- "2lonly.pmm"
mice_methods["Adapted"] <- "2lonly.pmm"
mice_methods["Stage"] <- "2lonly.pmm"
mice_methods["Chemo_sequence"] <- "2lonly.pmm"
mice_methods["GTV_volume"] <- "2lonly.pmm"
mice_methods["Dysphagia_BL_PHYS_CTCAE"] <- "2lonly.pmm"
mice_methods["Dyspnea_PHYS_BL_CTCAE"] <- "2lonly.pmm"
mice_methods["BaselineGHS"] <- "2lonly.pmm"
mice_methods["BaselineEQ5D"] <- "2lonly.pmm"
mice_methods["BaselineEQVAS"] <- "2lonly.pmm"
mice_methods["Fractions_proton"] <- "2lonly.pmm"
mice_methods["Fractions_photon"] <- "2lonly.pmm"
mice_methods["Fractions_total"]  <- "2lonly.pmm"

print(mice_methods)

# Generate initial predictor matrix with maxit = 0 to avoid actual imputation
initial_imputation <- mice(MergedData, method = mice_methods, m = 1, maxit = 0, seed = 123)
predictor_matrix <- initial_imputation$predictorMatrix

print(predictor_matrix)

consistent_vars <- c("Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Treatment30", "Fraction_Dailyfrequency",
                     "Pulmonary_comorbidity", "Smoking_status", "Gender", "Age_startRT", "WHO_PS",  "Adapted", "Stage", 
                     "Chemo_sequence", "Radiation_DoseReceived", "BaselineGHS", "BaselineEQVAS", "BaselineEQ5D", 
                     "Days_StartStop",  "GTV_volume", "Dysphagia_BL_PHYS_CTCAE", "Dyspnea_PHYS_BL_CTCAE")
for (var in consistent_vars) {
predictor_matrix[var, c("Periode", "EQindex", "EQ5DVAS", "GlobalHealthStatus", "Label", "Treatment50")] <- 0  
}
predictor_matrix["Label", ] <- 0  # Exclude 'Label' from predicting any other variable
predictor_matrix["WHO_PS_matching", ] <- 0  # Exclude 'WHO_PS_matching' from predicting any other variable. This is a grouped variable for matching only. 
predictor_matrix["Treatment50", ] <- 0  # Exclude 'Treatment50' from predicting any other variable
predictor_matrix["Treatment80", ] <- 0  # Exclude 'Treatment50' from predicting any other variable
predictor_matrix["Treatment100", ] <- 0  # Exclude 'Treatment50' from predicting any other variable
predictor_matrix[consistent_vars, "PatientID"] <- -2  # to ensure consistent imputation for variables not dependent on time

####### RETROSPECTIVELY ADDING number of fractions as protons and photons 
 ###### MAYBE CHANGE LATER - if so, just add these variables to consistent_vars and remove this
# adding number of fractions as protons and photons (Note - retrospectively adding - imputing two missings)
numfract_consistent_vars <- c("Fractions_proton", "Fractions_photon", "Fractions_total", "TreatmentSplit")
other_vars <- setdiff(names(MergedData), numfract_consistent_vars)
for (var in other_vars) {
  predictor_matrix[var, numfract_consistent_vars] <- 0
} # exclude from predicting other variables 
for (var in numfract_consistent_vars) {
  predictor_matrix[var, "PatientID"] <- -2
} # ensure consistent for a given patient across timepoints
####### END Note for adding 




print(predictor_matrix)


## Perform MICE imputation
imputed_data <- mice(MergedData, method = mice_methods, predictorMatrix = predictor_matrix, m = 82, maxit = 50, seed = 123)
#NOTE// m = number of imputed datasets; maxit = maximum number of iterations; seed = random seed for reproducibility

summary(imputed_data)
plot(imputed_data)
densityplot(imputed_data, main = "Density Plot of Imputed vs. Observed Data") #density plot of observed vs imputed data per variable (numeric variables)

## Assessing individual imputed datasets
dfSummary(MergedData[, c("Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence",
                            "EQindex", "EQ5DVAS", "GlobalHealthStatus", "Radiation_DoseReceived", 
                            "GTV_volume", "Dysphagia_BL_PHYS_CTCAE",
                            "Dyspnea_PHYS_BL_CTCAE")]) #Original data - variables with missing, all observation
imputed_dataset <- complete(imputed_data, action = 1) #change action number to assess different individual imputed datasets
dfSummary(imputed_dataset) 

check_consistency <- function(x) { #CHECK CONSISTENCY OF IMPUTED DATA 
  length(unique(na.omit(x))) == 1
}
consistency_results <- tapply(imputed_dataset$Durvalumab, imputed_dataset$PatientID, check_consistency) #Adjust "Durvalumab" to check other variable consistency
inconsistent_patients <- names(consistency_results[consistency_results == FALSE])
print(inconsistent_patients)


view(dfSummary(imputed_dataset[, c("PatientID", "Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence",
                                   "EQindex", "EQ5DVAS", "GlobalHealthStatus", "Radiation_DoseReceived", 
                                   "GTV_volume", "Dysphagia_BL_PHYS_CTCAE",
                                   "Dyspnea_PHYS_BL_CTCAE")])) # Imputed dataset (observed+imputed), all observations
missing_values <- is.na(MergedData)
imputed_only <- imputed_dataset
imputed_only[!missing_values] <- NA
view(dfSummary(imputed_only)) # imputed values only, all observations

MergedData_unique <- MergedData[!duplicated(MergedData$PatientID), ]
view(dfSummary(MergedData_unique[, c("PatientID", "Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence",
                              "Radiation_DoseReceived", "GTV_volume", "Dysphagia_BL_PHYS_CTCAE",
                              "Dyspnea_PHYS_BL_CTCAE")])) # Observed data, unique patients only

imputed_dataset_unique <- imputed_dataset[!duplicated(imputed_dataset$PatientID), ]
view(dfSummary(imputed_dataset_unique[, c("PatientID", "Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence",
                                       "Radiation_DoseReceived", "GTV_volume", "Dysphagia_BL_PHYS_CTCAE",
                                       "Dyspnea_PHYS_BL_CTCAE")])) #Imputed dataset (observed + imputed), unique patients only
missing_values_unique <- is.na(MergedData_unique)
imputed_only_unique <- imputed_dataset_unique
imputed_only_unique[!missing_values_unique] <- NA
view(dfSummary(imputed_only_unique[, c("PatientID", "Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence",
                                       "Radiation_DoseReceived", "GTV_volume", "Dysphagia_BL_PHYS_CTCAE",
                                       "Dyspnea_PHYS_BL_CTCAE")])) #Imputed values only, unique patients only


# Exploring Durvalumab

Observed_Durvyes <- MergedData_unique[MergedData_unique$Durvalumab == 'yes', ]
view(dfSummary(Observed_Durvyes))
Observed_Durvno <- MergedData_unique[MergedData_unique$Durvalumab == 'no', ]
view(dfSummary(Observed_Durvno))

patients_with_na_durvalumab <- is.na(MergedData_unique$Durvalumab)
imputed_for_na_durvalumab <- imputed_dataset_unique[patients_with_na_durvalumab, ]
view(dfSummary(imputed_for_na_durvalumab[imputed_for_na_durvalumab$Durvalumab == 'yes', ])) # summary of imputed + observed for patients that had imputed Durvalumab yes
view(dfSummary(imputed_for_na_durvalumab[imputed_for_na_durvalumab$Durvalumab == 'no', ])) # summary of imputed + observed for patients that had imputed Durvalumab no


#Continuous variables
imputed_continuous <- imputed_dataset[, c("EQindex", "EQ5DVAS", "GlobalHealthStatus", "Radiation_DoseReceived", 
                                          "GTV_volume", "Fractions_proton", "Fractions_photon", "Fractions_total")]  #numeric variables only
str(imputed_continuous)
chart.Correlation(imputed_continuous, histogram = TRUE, pch = 19)

## Assessing stacked data
long_imputed <- complete(imputed_data, action = 'long') #includes all imputed data stacked
dfSummary(long_imputed)
view(dfSummary(MergedData[, c("Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence",
       "EQindex", "EQ5DVAS", "GlobalHealthStatus", "Radiation_DoseReceived", "GTV_volume", "Dysphagia_BL_PHYS_CTCAE",
       "Dyspnea_PHYS_BL_CTCAE")]))
view(dfSummary(long_imputed[, c("Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence",
                           "EQindex", "EQ5DVAS", "GlobalHealthStatus", "Radiation_DoseReceived", "GTV_volume", "Dysphagia_BL_PHYS_CTCAE",
                           "Dyspnea_PHYS_BL_CTCAE")])) #Includes imputed and original values


#Continuous variables
long_imputed_variables <- imputed_dataset[, c("EQindex", "EQ5DVAS", "GlobalHealthStatus", "Radiation_DoseReceived", 
                                              "GTV_volume")]  #numeric variables only
chart.Correlation(long_imputed_variables, histogram = TRUE, pch = 19)

#factor variables
dfSummary(MergedData[, c("Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence", 
                         "Dysphagia_BL_PHYS_CTCAE", "Dyspnea_PHYS_BL_CTCAE")])
dfSummary(long_imputed[, c("Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence", 
                           "Dysphagia_BL_PHYS_CTCAE", "Dyspnea_PHYS_BL_CTCAE")])
str(long_imputed)
categorical_vars <- c("Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", "Stage", "Chemo_sequence", 
                      "Dysphagia_BL_PHYS_CTCAE", "Dyspnea_PHYS_BL_CTCAE")
for (var in categorical_vars) {
  p <- ggplot() +
    geom_density(data = data.frame(Value = MergedData[[var]], Status = "Original"), 
                 aes(x = Value, fill = Status), alpha = 0.5) +
    geom_density(data = long_imputed, 
                 aes(x = .data[[var]], fill = "Imputed"), alpha = 0.5) +
    labs(title = paste("Density Plot: Original vs. Imputed", var, "Data"),
         x = var,
         y = "Density") +
    scale_fill_manual(values = c("blue", "orange"), 
                      name = "Data Status") +
    theme_minimal()
  
  print(p)  # Print each plot
}

dfSummary(imputed_dataset_unique[, c("PatientID", "Treatment30", "Treatment50", "TreatmentSplit", "Fractions_proton", "Fractions_photon", "Fractions_total",
                                     "Age_startRT", "Gender", "Stage", "TumorCellType", "WHO_PS",  
                                     "TumorLocation", "Pulmonary_comorbidity", "GTV_volume", "Smoking_status", "Chemo_sequence", "Radiation_DoseReceived",
                                     "Days_StartStop", "Fraction_Dailyfrequency", "Surgery", "Adapted", "Durvalumab", "Dysphagia_BL_PHYS_CTCAE",
                                    "Dyspnea_PHYS_BL_CTCAE", "BaselineEQVAS", "BaselineEQ5D", "BaselineGHS", "EQindex", "GlobalHealthStatus", "EQ5DVAS")])

protons_imp_unique <- imputed_dataset_unique %>%
  filter(Treatment30 == "Protons")
dfSummary(protons_imp_unique[, c("PatientID", "Treatment30", "Treatment50", "TreatmentSplit", "Fractions_proton", "Fractions_photon", "Fractions_total",
                                 "Age_startRT", "Gender", "Stage", "TumorCellType", "WHO_PS",  
                                 "TumorLocation", "Pulmonary_comorbidity", "GTV_volume", "Smoking_status", "Chemo_sequence", "Radiation_DoseReceived",
                                 "Days_StartStop", "Fraction_Dailyfrequency", "Surgery", "Adapted", "Durvalumab", "Dysphagia_BL_PHYS_CTCAE",
                                 "Dyspnea_PHYS_BL_CTCAE", "BaselineEQVAS", "BaselineEQ5D", "BaselineGHS", "EQindex", "GlobalHealthStatus", "EQ5DVAS")])
photons_imp_unique <- imputed_dataset_unique %>%
  filter(Treatment30 == "Photons")
dfSummary(photons_imp_unique[, c("PatientID", "Treatment30", "Treatment50", "TreatmentSplit", "Fractions_proton", "Fractions_photon", "Fractions_total",
                                 "Age_startRT", "Gender", "Stage", "TumorCellType", "WHO_PS",  
                                 "TumorLocation", "Pulmonary_comorbidity", "GTV_volume", "Smoking_status", "Chemo_sequence", "Radiation_DoseReceived",
                                 "Days_StartStop", "Fraction_Dailyfrequency", "Surgery", "Adapted", "Durvalumab", "Dysphagia_BL_PHYS_CTCAE",
                                 "Dyspnea_PHYS_BL_CTCAE", "BaselineEQVAS", "BaselineEQ5D", "BaselineGHS", "EQindex", "GlobalHealthStatus", "EQ5DVAS")])


# Exploring Durvalumab (pmm and logreg both result in changed distribution in results)
dfSummary(MergedData$Durvalumab)
dfSummary(long_imputed$Durvalumab)

ggplot() +
  geom_boxplot(data = data.frame(Durvalumab = MergedData$Durvalumab, Status = "Original"),
               aes(x = Status, y = Durvalumab), fill = "blue", alpha = 0.5) +
  geom_boxplot(data = long_imputed, 
               aes(x = "Imputed", y = Durvalumab), fill = "orange", alpha = 0.5) +
  labs(title = "Boxplot: Original vs. Imputed Durvalumab Data",
       x = "Data Status",
       y = "Durvalumab") +
  theme_minimal()

ggplot() +
  geom_density(data = data.frame(Durvalumab = MergedData$Durvalumab, Status = "Original"), 
               aes(x = Durvalumab, fill = Status), alpha = 0.5) +
  geom_density(data = long_imputed, 
               aes(x = Durvalumab, fill = "Imputed"), alpha = 0.5) +
  labs(title = "Density Plot: Original vs. Imputed Durvalumab Data",
       x = "Durvalumab",
       y = "Density") +
  scale_fill_manual(values = c("blue", "orange"), 
                    name = "Data Status") +
  theme_minimal()

ggplot() +
  geom_bar(data = data.frame(Durvalumab = MergedData$Durvalumab, Status = "Original"),
           aes(x = Durvalumab, fill = Status), alpha = 0.5, position = "identity") +
  geom_bar(data = long_imputed,
           aes(x = Durvalumab, fill = "Imputed"), alpha = 0.5, position = "identity") +
  labs(title = "Bar Plot: Original vs. Imputed Durvalumab Data",
       x = "Durvalumab",
       y = "Count") +
  scale_fill_manual(values = c("blue", "orange"), 
                    name = "Data Status") +
  theme_minimal()



### Box plots for imputed data - HRQoL scores over time
imputed_list <- lapply(1:imputed_data$m, function(i) {
  dat <- complete(imputed_data, i)
  dat$Imputation <- i
  return(dat)
})
combined_data <- bind_rows(imputed_list)

give.n <- function(x){    
  return(c(y = mean(x), label = length(x)))
}
palette(rainbow(3, s = 0.5)) #color ramp 
#palette(gray.colors(3)) #color ramp 
combined_data$Periode <- as.factor(combined_data$Periode) 

# Aggregate the data by taking the average across imputations for each subject
aggregated_data <- combined_data %>%
  group_by(PatientID, Periode, Treatment30) %>%  # adjust grouping variables if needed
  summarise(
    avg_EQindex = mean(EQindex, na.rm = TRUE),
    avg_GlobalHealthStatus = mean(GlobalHealthStatus, na.rm = TRUE),
    avg_EQVAS = mean(EQ5DVAS, na.rm = TRUE)
  ) %>%
  ungroup()

# EQ5D5L (using the average EQindex)
bxp_EQindex <- ggplot(aggregated_data, aes(x = Periode, y = avg_EQindex, fill = Treatment30)) +
  scale_fill_manual(values = c("1", "2", "3")) +
  geom_boxplot(alpha = 0.7) + 
  # Remove or update the give.n summary if not needed:
  # stat_summary(fun.data = give.n, geom = "text", position = position_dodge(width = 0.75), size = 2, aes(y = -0.5)) +
  scale_y_continuous(name = "EQ5D5L Utility", limits = c(-0.5, 1)) +
  scale_x_discrete(name = "Time (months)") +
  labs(fill = "") +
  stat_summary(fun = mean, colour = "black", geom = "point", 
               shape = 21, size = 2, show.legend = FALSE, 
               position = position_dodge(width = 0.75))

# Global Health Status
bxp_GHS <- ggplot(aggregated_data, aes(x = Periode, y = avg_GlobalHealthStatus, fill = Treatment30)) +
  scale_fill_manual(values = c("1", "2", "3")) +
  geom_boxplot(alpha = 0.7) + 
  # Optionally add a stat_summary for counts or other annotations if desired:
  scale_y_continuous(name = "Global Health Status Score", limits = c(-0.5, 1)) +
  scale_x_discrete(name = "Time (months)") +
  labs(fill = "") +
  stat_summary(fun = mean, colour = "black", geom = "point", 
               shape = 21, size = 2, show.legend = FALSE, 
               position = position_dodge(width = 0.75))

# EQ VAS 
bxp_EQVAS <- ggplot(aggregated_data, aes(x = Periode, y = avg_EQVAS, fill = Treatment30)) +
  scale_fill_manual(values = c("1", "2", "3")) +
  geom_boxplot(alpha = 0.7) + 
  scale_y_continuous(name = "EQ VAS Utility", limits = c(-0.5, 1)) +
  scale_x_discrete(name = "Time (months)") +
  labs(fill = "") +
  stat_summary(fun = mean, colour = "black", geom = "point", 
               shape = 21, size = 2, show.legend = FALSE, 
               position = position_dodge(width = 0.75))
# Arrange the plots
figure12 <- ggarrange(bxp_EQindex, bxp_EQVAS, bxp_GHS,
                      labels = c("", "", ""),
                      ncol = 1, nrow = 3)


# baseline data summary across imputations

# First, create a list of unique (one row per patient) imputed datasets
imputed_unique_list <- lapply(1:imputed_data$m, function(i) {
  dat <- complete(imputed_data, i)
  dat[!duplicated(dat$PatientID), ]
})

# Subset the unique datasets for protons and photons separately
imputed_unique_list_protons <- lapply(imputed_unique_list, function(data) {
  subset(data, Treatment30 == "Protons")
})
imputed_unique_list_photons <- lapply(imputed_unique_list, function(data) {
  subset(data, Treatment30 == "Photons")
})

# Define the numeric variables you want to summarize (adjust these as needed)
numeric_vars <- c("Age_startRT", "GTV_volume", "Radiation_DoseReceived", 
                  "EQindex", "EQ5DVAS", "GlobalHealthStatus", "Fractions_proton", "Fractions_photon")

# Define a function to get summary statistics for one variable across a list of imputed datasets
get_pooled_numeric <- function(imputed_list, var) {
  # For each imputed dataset, compute mean, sd, median, min, and max
  summaries <- sapply(imputed_list, function(data) {
    c(mean = mean(data[[var]], na.rm = TRUE),
      sd = sd(data[[var]], na.rm = TRUE),
      median = median(data[[var]], na.rm = TRUE),
      min = min(data[[var]], na.rm = TRUE),
      max = max(data[[var]], na.rm = TRUE))
  })
  # Pool by averaging across imputations
  pooled_stats <- rowMeans(summaries)
  return(pooled_stats)
}

# Create pooled summaries for each numeric variable for the overall dataset, protons, and photons
pooled_total <- sapply(numeric_vars, function(var) get_pooled_numeric(imputed_unique_list, var))
pooled_protons <- sapply(numeric_vars, function(var) get_pooled_numeric(imputed_unique_list_protons, var))
pooled_photons <- sapply(numeric_vars, function(var) get_pooled_numeric(imputed_unique_list_photons, var))

## Number of patients with a combination of protons and photons
#Total
count_patients <- function(data) {
  sum(data$Fractions_proton > 0 & data$Fractions_photon > 0, na.rm = TRUE)
}
counts <- sapply(imputed_unique_list, count_patients) # Apply this function to each imputed unique dataset
average_count <- mean(counts) # Compute the average count across all imputations
average_count
#Protons
counts_protons <- sapply(imputed_unique_list_protons, function(data) {
  sum(data$Fractions_proton > 0 & data$Fractions_photon > 0, na.rm = TRUE)
})
average_count_protons <- mean(counts_protons) # Compute the average count across all imputations
average_count_protons
#Photons (no missing for photons so should be the same as pre-imputation)
counts_photons <- sapply(imputed_unique_list_photons, function(data) {
  sum(data$Fractions_proton > 0 & data$Fractions_photon > 0, na.rm = TRUE)
})
average_count_photons <- mean(counts_photons) # Compute the average count across all imputations
average_count_photons

# Print the pooled descriptive statistics
cat("Total (protons + photons):\n")
print(pooled_total)

cat("\nProton patients only:\n")
print(pooled_protons)

cat("\nPhoton patients only:\n")
print(pooled_photons)

categorical_vars <- c("Durvalumab", "TumorCellType", "TumorLocation", "Surgery", "Adapted", 
                      "Stage", "Chemo_sequence", "Smoking_status", "Gender", "Dyspnea_PHYS_BL_CTCAE", 
                      "Dysphagia_BL_PHYS_CTCAE", "WHO_PS", "Pulmonary_comorbidity", "Fraction_Dailyfrequency")
# Function to pool frequency counts for a categorical variable across imputations
get_pooled_categorical <- function(imputed_list, var) {
  # Get all levels present in any of the datasets
  all_levels <- unique(unlist(lapply(imputed_list, function(data) levels(as.factor(data[[var]])))))
  
  # For each level, compute the average frequency across imputations
  pooled_counts <- sapply(all_levels, function(level) {
    mean(sapply(imputed_list, function(data) {
      sum(as.character(data[[var]]) == level, na.rm = TRUE)
    }))
  })
  
  return(pooled_counts)
}

# Assuming imputed_unique_list is already defined as one-row-per-patient datasets
# (see your previous code)
# Create subsets for protons and photons (if not already created)
imputed_unique_list_protons <- lapply(imputed_unique_list, function(data) {
  subset(data, Treatment30 == "Protons")
})
imputed_unique_list_photons <- lapply(imputed_unique_list, function(data) {
  subset(data, Treatment30 == "Photons")
})

# Now loop through each categorical variable and pool the summaries
pooled_categorical_total <- lapply(categorical_vars, function(var) {
  stats <- get_pooled_categorical(imputed_unique_list, var)
  return(stats)
})
names(pooled_categorical_total) <- categorical_vars

pooled_categorical_protons <- lapply(categorical_vars, function(var) {
  stats <- get_pooled_categorical(imputed_unique_list_protons, var)
  return(stats)
})
names(pooled_categorical_protons) <- categorical_vars

pooled_categorical_photons <- lapply(categorical_vars, function(var) {
  stats <- get_pooled_categorical(imputed_unique_list_photons, var)
  return(stats)
})
names(pooled_categorical_photons) <- categorical_vars

# Print the pooled results
cat("Pooled Categorical Variables - Total (Protons + Photons):\n")
print(pooled_categorical_total)

cat("\nPooled Categorical Variables - Proton Patients Only:\n")
print(pooled_categorical_protons)

cat("\nPooled Categorical Variables - Photon Patients Only:\n")
print(pooled_categorical_photons)


### Agregated baseline characteristcs across imputed datasets




# Save imputed data object
saveRDS(imputed_data, file = "mice_imputed_data.rds")

