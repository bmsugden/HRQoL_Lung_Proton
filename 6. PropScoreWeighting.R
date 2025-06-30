rm(list = ls())

# Load imputed data
imputed_data <- readRDS("mice_imputed_data.rds")


# Load required libraries
library(mice)
library(ggplot2) 
library(dplyr)# For density plots
library(summarytools)
library(car)
library(twang)
library(survey) #to use weights in final weighted analysis
library(cobalt) # for love plot to assess covariate balance post-weighting
library(knitr) # for variable influence tables

# Number of imputations
m <- imputed_data$m  
weighted_data_list <- vector("list", m)

# Loop through each imputed dataset
for (i in 1:m) {
  
  # 1. Extract the full imputed dataset for the i-th imputation (includes all time points)
  full_data_i <- complete(imputed_data, i)
  # 2. Collapse the data to one row per patient (select first occurrence)
  collapsed_data <- full_data_i[!duplicated(full_data_i$PatientID), ]
  # Combine Stage III and IV into one category for weighting
  collapsed_data$Stage_combined <- factor(
    ifelse(collapsed_data$Stage %in% c("III", "IV"), "III_IV", as.character(collapsed_data$Stage))
  )
  # Define treatment indicator (note, twang required this to be numeric)
  collapsed_data$treatment_indicator <- as.numeric(collapsed_data$Treatment30 == "Protons")
  # Define covariates for weighting
  covariates <- c("Age_startRT", "Stage_combined", "WHO_PS_matching",
                  "Dyspnea_PHYS_BL_CTCAE", "BaselineEQ5D")
  # 3. Estimate Propensity Scores using TWANG (GBM)
  set.seed(12345) # set seed for reproducibility
  ps_model <- ps(
    formula = as.formula(paste("treatment_indicator ~", paste(covariates, collapse = " + "))),
    data = collapsed_data,
    n.trees = 15000, #GBM builds seq. dec. trees to model PS's. TWANG default 10000 - increase if poor balance
    interaction.depth = 4, #tree depth - controls how many interactions between covariates are selected. Default = 3
    shrinkage = 0.0005, #how much each tree contributes to final model - lower values reduc sensitivity to noise. Prevents overfitting. 
    # If poor balance, reduce shrinkage from 0.01 (default)
    stop.method = c("es.mean", "ks.max"), #stopping rules. effect size mean minimises std. mean diffs across covariates
    # ks.max minimises max kolmogorov-smirnov stat 
    # 2 stopping rules explored. Compare results in balance tables
    estimand = "ATE" 
  )
  # Extract the final stabilised weights
  collapsed_data$twang_weight <- pmin(pmax(ps_model$w$ks.max.ATE, 0.2), 5)  # Final trimmed weights
  # 4. Merge the TWANG weights back into the full longitudinal dataset
  full_weighted_data <- merge(full_data_i, collapsed_data[, c("PatientID", "twang_weight")], 
                              by = "PatientID", all.x = TRUE)
  # Add imputation indicator
  full_weighted_data$Imputation <- i
  # Store in list
  weighted_data_list[[i]] <- full_weighted_data
}
warnings() # warning suggest increasing n.trees. However check this is neccessary by: 
summary(ps_model) # see iterations - n.trees is plenty. Next step: check balance 

#### Assessing balance
# Balance table
## unw = before weighting; es.mean.ATE = effect size mean (after weighting); ks.mean.ATE = kolmogorov-smirnov mean (after weighting)
## 
bal.table(ps_model) #for criteria, see accompanying word file

love.plot(ps_model, threshold = 0.2, stars = "std") #plots std mean diffs for each covariate before and after weighting. 
plot(ps_model) #check convergence
plot(ps_model, plots = 2)  # Check overlap of propensity scores
plot(ps_model, plots = 3)

# Additional check
summary(ps_model$gbm.obj) #shows which covariates were most influential for estimating propensity scores
# Summaries
summary(collapsed_data$twang_weight)
hist(collapsed_data$twang_weight, breaks = 30, main = "Histogram of TWANG Weights")
### Pooled summary of balanced data
# --- 1. Build a list of weighted, one‐row‐per‐patient datasets ----
weighted_unique_list <- lapply(weighted_data_list, function(full_data_i) {
  # collapse to one row per Patient
  unique_pat <- full_data_i[!duplicated(full_data_i$PatientID), ]
  # ensure the weight column is present
  stopifnot("twang_weight" %in% colnames(unique_pat))
  unique_pat
})
# --- 2. Define weighted‐summary function ----
get_pooled_weighted_numeric <- function(data_list, var) {
  # For each imputed dataset, compute weighted stats
  stats_mat <- sapply(data_list, function(df) {
    x <- df[[var]]
    w <- df$twang_weight
    # weighted mean
    m <- weighted.mean(x, w, na.rm = TRUE)
    # weighted variance and sd
    v <- sum(w * (x - m)^2, na.rm = TRUE) / sum(w, na.rm = TRUE)
    s <- sqrt(v)
    # weighted median via Hmisc, or approximate by wtd.quantile
    if (!requireNamespace("Hmisc", quietly = TRUE)) {
      install.packages("Hmisc")
    }
    med <- Hmisc::wtd.quantile(x, weights = w, probs = 0.5, na.rm = TRUE)
    # unweighted min/max (extremes are insensitive to weight)
    mn <- min(x, na.rm = TRUE)
    mx <- max(x, na.rm = TRUE)
    c(mean = m, sd = s, median = med, min = mn, max = mx)
  })
  # Pool by averaging each row (statistic) across imputations
  pooled <- rowMeans(stats_mat)
  return(pooled)
}
# --- 3. Apply to variables ----
numeric_vars <- c("Age_startRT", "GTV_volume", "Radiation_DoseReceived", 
                  "EQindex", "EQ5DVAS", "GlobalHealthStatus", 
                  "Fractions_proton", "Fractions_photon")
# overall pooled weighted summaries
pooled_weighted_total <- sapply(numeric_vars,
                                function(v) get_pooled_weighted_numeric(weighted_unique_list, v))
# separately for protons and photons
weighted_list_protons <- lapply(weighted_unique_list,
                                function(df) subset(df, Treatment30 == "Protons"))
weighted_list_photons <- lapply(weighted_unique_list,
                                function(df) subset(df, Treatment30 == "Photons"))
pooled_weighted_protons <- sapply(numeric_vars,
                                  function(v) get_pooled_weighted_numeric(weighted_list_protons, v))
pooled_weighted_photons <- sapply(numeric_vars,
                                  function(v) get_pooled_weighted_numeric(weighted_list_photons, v))

# --- 4. Present as a table ----
kable(t(pooled_weighted_total),
      caption = "Pooled Weighted Summary Statistics (Overall)",
      digits = 2)

kable(t(pooled_weighted_protons),
      caption = "Pooled Weighted Summary Statistics (Protons)",
      digits = 2)

kable(t(pooled_weighted_photons),
      caption = "Pooled Weighted Summary Statistics (Photons)",
      digits = 2)



# save 
saveRDS(weighted_data_list, file = "weighted_data_list.rds")



